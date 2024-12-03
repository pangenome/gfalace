use std::collections::{HashMap, BTreeMap};
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

fn print_graph_details(graph: &HashGraph) {
    println!("\nDetailed Graph Information:");
    
    // Print all node IDs and their sequences
    println!("\nNodes:");
    for handle in graph.handles() {
        let sequence = graph.sequence(handle).collect::<Vec<_>>();
        let sequence_str = String::from_utf8(sequence).unwrap_or_else(|_| "Invalid UTF-8".to_string());
        println!("  Node ID: {}, Sequence: {}", handle.id(), sequence_str);
    }
    
    // Print all edges
    println!("\nEdges:");
    for edge in graph.edges() {
        println!("  Edge: {} ({}) -> {} ({})",
            edge.0.id(),
            if edge.0.is_reverse() { "reverse" } else { "forward" },
            edge.1.id(),
            if edge.1.is_reverse() { "reverse" } else { "forward" }
        );
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
    for (gfa_id, gfa_path) in args.gfa_list.iter().enumerate() {
        let parser = GFAParser::new();
        let gfa: GFA<usize, ()> = parser.parse_file(gfa_path).unwrap();
        let block_graph = HashGraph::from_gfa(&gfa);

        // Record the id translation for this block
        let id_translation = NodeId::from(combined_graph.node_count());
        id_translations.push(id_translation);

        // Add nodes with translated IDs
        for handle in block_graph.handles() {
            let sequence = block_graph.sequence(handle).collect::<Vec<_>>();
            combined_graph.create_handle(&sequence, id_translation + handle.id().into());
        }

        // Add edges with translated IDs
        for edge in block_graph.edges() {
            let translated_edge = Edge(
                Handle::pack(id_translation + edge.0.id().into(), edge.0.is_reverse()),
                Handle::pack(id_translation + edge.1.id().into(), edge.1.is_reverse())
            );
            combined_graph.create_edge(translated_edge);
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
                        .or_insert_with(Vec::new)
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

    // Sort ranges and create merged paths in the combined graph
    for (path_key, ranges) in path_key_ranges.iter_mut() {
        // Sort ranges by start position
        ranges.sort_by_key(|r| (r.start, r.end));
        
        // Create the path in the combined graph
        let path_id = combined_graph.create_path(path_key.as_bytes(), false).unwrap();
        
        // Merge contiguous ranges and add their steps
        let mut current_range_idx = 0;
        while current_range_idx < ranges.len() {
            let mut steps = ranges[current_range_idx].steps.clone();
            let mut next_idx = current_range_idx + 1;
            
            // Keep merging while ranges are contiguous
            while next_idx < ranges.len() && is_contiguous(&ranges[next_idx - 1], &ranges[next_idx]) {
                steps.extend(ranges[next_idx].steps.clone());
                next_idx += 1;
            }
            
            // Add all steps from the merged ranges to the path
            for step in steps {
                combined_graph.path_append_step(path_id, step);
            }
            
            current_range_idx = next_idx;
        }
    }

    println!("Path creation completed");
    println!("Number of paths in combined graph: {}", GraphPaths::path_count(&combined_graph));
    // Print some statistics about the combined graph
    // println!("Combined graph statistics:");
    // println!("  Nodes: {}", combined_graph.node_count());
    // println!("  Edges: {}", combined_graph.edge_count());
    // println!("Number of unique path keys: {}", path_key_ranges.len());

    // Print detailed information about nodes and edges
    //print_graph_details(&combined_graph);

     // Write the combined graph to GFA file
    match write_graph_to_gfa(&combined_graph, &args.output) {
        Ok(_) => println!("Successfully wrote combined graph to {}", args.output),
        Err(e) => eprintln!("Error writing GFA file: {}", e),
    }
}

pub struct PathRange {
    pub path_name: String,
    pub steps: Vec<LacePathStep>,
}

pub struct LacePathStep {
    pub block_id: usize,
    pub node_id: NodeId,
}

pub fn lace_smoothed_blocks(
    smoothed_block_graphs: &[HashGraph],
    path_ranges: &[PathRange],
) -> HashGraph {
    let mut smoothed_graph = HashGraph::new();
    let mut block_id_to_smoothed_id: Vec<HashMap<NodeId, NodeId>> = vec![HashMap::new(); smoothed_block_graphs.len()];

    // Copy nodes and edges from the smoothed blocks to the smoothed graph
    for (block_id, smoothed_block) in smoothed_block_graphs.iter().enumerate() {
        for handle in smoothed_block.handles() {
            let sequence = smoothed_block.sequence(handle).collect::<Vec<_>>();
            let new_node_id = smoothed_graph.create_handle::<usize>(&sequence.into_boxed_slice(), handle.id().into());
            block_id_to_smoothed_id[block_id].insert(handle.id(), new_node_id.into());
        }
        for edge in smoothed_block.edges() {
            smoothed_graph.create_edge(edge);
        }
    }

    // Create paths in the smoothed graph based on the path ranges
    for path_range in path_ranges {
        let path_id = smoothed_graph.create_path(path_range.path_name.as_bytes(), false).unwrap();

        for step in &path_range.steps {
            let smoothed_node_id = block_id_to_smoothed_id[step.block_id][&step.node_id];
            smoothed_graph.path_append_step(path_id, Handle::pack(smoothed_node_id, false));
        }
    }

    // Connect the path steps in the smoothed graph
    let mut edges = Vec::new();
    for path_id in smoothed_graph.path_ids() {
        let mut prev_step = None;
        for step in smoothed_graph.get_path_ref(path_id).unwrap().nodes.iter() {
            if let Some(prev) = prev_step {
                edges.push(Edge(prev, *step));
            }
            prev_step = Some(*step);
        }
    }
    // Create edges collected in previous step
    for edge in edges {
        smoothed_graph.create_edge(edge);
    }
        
    smoothed_graph
}

fn split_path_name(path_name: &str) -> Option<(String, usize, usize)> {
    let parts: Vec<&str> = path_name.split('#').collect();
    if parts.len() == 3 {
        let key_parts = vec![parts[0], parts[1]];
        let chr_range: Vec<&str> = parts[2].split(':').collect();
        if chr_range.len() == 2 {
            let name = chr_range[0];
            let range: Vec<&str> = chr_range[1].split('-').collect();
            if range.len() == 2 {
                let start = range[0].parse().ok()?;
                let end = range[1].parse().ok()?;
                let key = format!("{}#{}", key_parts.join("#"), name);
                return Some((key, start, end));
            }
        }
    }
    None
}
