use std::collections::{BTreeMap};
use clap::Parser;
use std::fs::File;
use std::io::Write;
use handlegraph::handle::{Handle, NodeId, Edge};
use handlegraph::handlegraph::*;
use handlegraph::mutablehandlegraph::*;
use handlegraph::pathhandlegraph::{IntoPathIds, GraphPathNames, GraphPathsRef, MutableGraphPaths, GraphPaths};
use handlegraph::hashgraph::HashGraph;
//use handlegraph::pathhandlegraph::PathStep;
use gfa::{gfa::GFA, parser::GFAParser};

#[derive(Parser, Debug)]
#[clap(author, version, about)]
struct Args {
    /// List of GFA file paths.
    #[clap(short, long, value_parser, num_args = 1.., value_delimiter = ' ')]
    gfa_list: Vec<String>,

    /// Output GFA file path
    #[clap(short, long, value_parser)]
    output: String,

    /// Enable debug output
    #[clap(short, long)]
    debug: bool,
}

#[derive(Debug)]
struct RangeInfo {
    start: usize,
    end: usize,
    gfa_id: usize,
    steps: Vec<Handle>,  // Store the path steps for this range
}

// Helper function to check if two ranges are contiguous
fn is_contiguous(r1: &RangeInfo, r2: &RangeInfo) -> bool {
    r1.end == r2.start
}

// Helper function to check if two ranges overlap
fn has_overlap(r1: &RangeInfo, r2: &RangeInfo) -> bool {
    r1.start < r2.end && r2.start < r1.end
}

// Helper function to read GFA file
fn read_gfa(gfa_path: &str, parser: &GFAParser<usize, ()>) -> std::io::Result<GFA<usize, ()>> {
    if gfa_path.ends_with(".gz") {
        let file = std::fs::File::open(gfa_path)?;
        let (mut reader, _format) = niffler::get_reader(Box::new(file)).unwrap();
        
        let mut decompressed = Vec::new();
        reader.read_to_end(&mut decompressed)?;
        
        let temp_path = format!("{}.tmp", gfa_path);
        {
            let mut temp_file = std::fs::File::create(&temp_path)?;
            temp_file.write_all(&decompressed)?;
        }
        
        let result = parser.parse_file(&temp_path).unwrap();
        std::fs::remove_file(&temp_path)?;
        
        Ok(result)
    } else {
        Ok(parser.parse_file(gfa_path).unwrap())
    }
}

fn write_graph_to_gfa(graph: &HashGraph, output_path: &str) -> std::io::Result<()> {
    let mut file = File::create(output_path)?;
    
    // Write GFA version
    writeln!(file, "H\tVN:Z:1.0")?;
    
    // Collect and sort nodes by ID
    let mut nodes: Vec<Handle> = graph.handles().collect();
    nodes.sort_by_key(|handle| handle.id());
    
    // Write sorted nodes (Segments)
    for handle in nodes {
        let sequence = graph.sequence(handle).collect::<Vec<_>>();
        let sequence_str = String::from_utf8(sequence)
            .unwrap_or_else(|_| String::from("N"));
        writeln!(file, "S\t{}\t{}", handle.id(), sequence_str)?;
    }
    
    // Collect and sort edges
    let mut edges: Vec<Edge> = graph.edges().collect();
    edges.sort_by(|a, b| {
        a.0.id().cmp(&b.0.id())
            .then(a.1.id().cmp(&b.1.id()))
    });
    
    // Write sorted edges (Links)
    for edge in edges {
        let from_id = edge.0.id();
        let to_id = edge.1.id();
        let from_orient = if edge.0.is_reverse() { "-" } else { "+" };
        let to_orient = if edge.1.is_reverse() { "-" } else { "+" };
        writeln!(file, "L\t{}\t{}\t{}\t{}\t0M", from_id, from_orient, to_id, to_orient)?;
    }
    
    // Collect and sort paths by name
    let mut paths: Vec<_> = graph.path_ids().collect();
    paths.sort_by_key(|&path_id| {
        graph.get_path_name(path_id)
            .map(|name_iter| name_iter.collect::<Vec<u8>>())
            .unwrap_or_default()
    });
    
    // Write sorted paths
    for path_id in paths {
        if let Some(name_iter) = graph.get_path_name(path_id) {
            let path_name = String::from_utf8(name_iter.collect::<Vec<u8>>())
                .unwrap_or_else(|_| String::from("unknown_path"));
            
            let mut path_elements = Vec::new();
            if let Some(path_ref) = graph.get_path_ref(path_id) {
                for handle in &path_ref.nodes {
                    let orient = if handle.is_reverse() { "-" } else { "+" };
                    path_elements.push(format!("{}{}", handle.id(), orient));
                }
            }
            
            writeln!(file, "P\t{}\t{}\t*", path_name, path_elements.join(","))?;
        }
    }
    
    Ok(())
}
fn main() {
    let args = Args::parse();

    // Create a single combined graph
    let mut combined_graph = HashGraph::new();
    let mut path_key_ranges: BTreeMap<String, Vec<RangeInfo>> = BTreeMap::new();
    let mut id_translations = Vec::new();

    // Process each GFA file
    let parser = GFAParser::new();
    for (gfa_id, gfa_path) in args.gfa_list.iter().enumerate() {
        let gfa = read_gfa(gfa_path, &parser).unwrap();
        let block_graph = HashGraph::from_gfa(&gfa);

        // Record the id translation for this block
        let id_translation = NodeId::from(combined_graph.node_count());
        id_translations.push(id_translation);

        // Add nodes with translated IDs
        for handle in block_graph.handles() {
            let sequence = block_graph.sequence(handle).collect::<Vec<_>>();
            let new_id = id_translation + handle.id().into();
            combined_graph.create_handle(&sequence, new_id);
        }

        // Add edges with translated IDs
        for edge in block_graph.edges() {
            let translated_edge = Edge(
                Handle::pack(id_translation + edge.0.id().into(), edge.0.is_reverse()),
                Handle::pack(id_translation + edge.1.id().into(), edge.1.is_reverse())
            );
            combined_graph.create_edge(translated_edge);
        }
        
        if args.debug {
            eprintln!("GFA file {} ({}) processed: Added {} nodes and {} edges", gfa_id, gfa_path, block_graph.node_count(), block_graph.edge_count());
        }

        // Process paths and collect ranges with their steps
        for path_id in block_graph.path_ids() {
            if let Some(name_iter) = block_graph.get_path_name(path_id) {
                let path_name = String::from_utf8(name_iter.collect::<Vec<u8>>()).unwrap();
                
                if let Some((sample_hap_name, start, end)) = split_path_name(&path_name) {
                    // Get the path steps and translate their IDs
                    let mut translated_steps = Vec::new();
                    if let Some(path_ref) = block_graph.get_path_ref(path_id) {
                        for step in path_ref.nodes.iter() {
                            let translated_id = id_translation + step.id().into();
                            translated_steps.push(Handle::pack(translated_id, step.is_reverse()));
                        }
                    }
                    
                    path_key_ranges.entry(sample_hap_name)
                        .or_default()
                        .push(RangeInfo { 
                            start, 
                            end, 
                            gfa_id,
                            steps: translated_steps,
                        });
                }
            }
        }
    }

    if args.debug {
        eprintln!("Total nodes in combined graph: {}", combined_graph.node_count());
        eprintln!("Total edges in combined graph: {}", combined_graph.edge_count());
        eprintln!("Total path ranges collected: {}", path_key_ranges.len());
    }

    // Sort ranges and create merged paths in the combined graph
    for (path_key, ranges) in path_key_ranges.iter_mut() {
        // Sort ranges by start position
        ranges.sort_by_key(|r| (r.start, r.end));
        
        // Check for overlaps and contiguity
        let mut has_overlaps = false;
        let mut all_contiguous = true;
        
        for window in ranges.windows(2) {
            if has_overlap(&window[0], &window[1]) {
                has_overlaps = true;
            }
            if !is_contiguous(&window[0], &window[1]) {
                all_contiguous = false;
            }
        }
        
        if (has_overlaps || !all_contiguous) && args.debug {
            eprintln!("  Path key '{}' ranges analysis:", path_key);
            
            let mut current_start = ranges[0].start;
            let mut current_end = ranges[0].end;
            let mut current_gfa_ids = vec![ranges[0].gfa_id];
            
            for i in 1..ranges.len() {
                if is_contiguous(&ranges[i-1], &ranges[i]) {
                    // Extend current merged range
                    current_end = ranges[i].end;
                    current_gfa_ids.push(ranges[i].gfa_id);
                } else {
                    // Print current merged range
                    eprintln!("    Merged range: start={}, end={}, gfa_ids={:?}", 
                        current_start, current_end, current_gfa_ids);
                    
                    if !has_overlap(&ranges[i-1], &ranges[i]) {
                        // Calculate and print gap
                        let gap = ranges[i].start - current_end;
                        eprintln!("      Gap to next range: {} positions", gap);
                    } else {
                        // Calculate and print overlap
                        let overlap = current_end - ranges[i].start;
                        eprintln!("      Overlap with next range: {} positions", overlap);
                    }

                    // Start new merged range
                    current_start = ranges[i].start;
                    current_end = ranges[i].end;
                    current_gfa_ids = vec![ranges[i].gfa_id];
                }
            }
            
            // Print final merged range
            eprintln!("    Merged range: start={}, end={}, gfa_ids={:?}", 
                current_start, current_end, current_gfa_ids);
        }

        if all_contiguous && !has_overlaps {
            // Create a single path with the original key
            let path_id = combined_graph.create_path(path_key.as_bytes(), false).unwrap();
            let mut prev_step = None;
            
            // Add all steps from all ranges
            for range in ranges.iter() {
                for step in &range.steps {
                    combined_graph.path_append_step(path_id, *step);
                    
                    if let Some(prev) = prev_step {
                        if !combined_graph.has_edge(prev, *step) {
                            combined_graph.create_edge(Edge(prev, *step));
                        }
                    }
                    prev_step = Some(*step);
                }
            }
        } else {
            // Handle non-contiguous ranges by creating separate paths for each contiguous group
            let mut current_range_idx = 0;
            while current_range_idx < ranges.len() {
                let start_range = &ranges[current_range_idx];
                let mut steps = start_range.steps.clone();
                let mut next_idx = current_range_idx + 1;
                let mut end_range = start_range;
                
                // Merge contiguous ranges
                while next_idx < ranges.len() && is_contiguous(&ranges[next_idx - 1], &ranges[next_idx]) {
                    steps.extend(ranges[next_idx].steps.clone());
                    end_range = &ranges[next_idx];
                    next_idx += 1;
                }
                
                // Create path name with range information
                let path_name = format!("{}:{}-{}", path_key, start_range.start, end_range.end);
                let path_id = combined_graph.create_path(path_name.as_bytes(), false).unwrap();
                
                // Add steps to the path
                let mut prev_step = None;
                for step in steps {
                    combined_graph.path_append_step(path_id, step);
                    
                    if let Some(prev) = prev_step {
                        if !combined_graph.has_edge(prev, step) {
                            combined_graph.create_edge(Edge(prev, step));
                        }
                    }
                    prev_step = Some(step);
                }
                
                current_range_idx = next_idx;
            }
        }
    }

    if args.debug {
        eprintln!("Total paths created: {}", GraphPaths::path_count(&combined_graph));
    }

    // Write the combined graph to GFA file
    match write_graph_to_gfa(&combined_graph, &args.output) {
        Ok(_) => if args.debug {eprintln!("Successfully wrote combined graph to {}", args.output)},
        Err(e) => eprintln!("Error writing GFA file: {}", e),
    }
}

fn split_path_name(path_name: &str) -> Option<(String, usize, usize)> {
    // Find the last ':' to split the range from the key
    if let Some(last_colon) = path_name.rfind(':') {
        let (key, range_str) = path_name.split_at(last_colon);
        // Skip the ':' character
        let range_str = &range_str[1..];
        
        // Find the '-' in the range portion
        if let Some((start_str, end_str)) = range_str.split_once('-') {
            if let (Ok(start), Ok(end)) = (start_str.parse(), end_str.parse()) {
                return Some((key.to_string(), start, end));
            }
        }
    }
    None
}
