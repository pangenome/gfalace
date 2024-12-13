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

#[derive(Debug, Clone)]
struct RangeInfo {
    start: usize,
    end: usize,
    gfa_id: usize,
    steps: Vec<Handle>,                    // Path steps for this range
    step_positions: Vec<(usize, usize)>,   // Start and end positions of each step
    step_lengths: Vec<usize>,              // Lengths of each step
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
                    let mut step_positions = Vec::new();
                    let mut step_lengths = Vec::new();
                    let mut cumulative_pos = start;

                    if let Some(path_ref) = block_graph.get_path_ref(path_id) {
                        for step in path_ref.nodes.iter() {
                            let translated_id = id_translation + step.id().into();
                            let translated_step = Handle::pack(translated_id, step.is_reverse());
                            translated_steps.push(translated_step);

                            // Get the sequence length of the node
                            let node_seq = block_graph.sequence(*step).collect::<Vec<_>>();
                            let node_length = node_seq.len();
                            step_lengths.push(node_length);

                            // Record the end position of this step
                            step_positions.push((cumulative_pos, cumulative_pos + node_length));
                            cumulative_pos = cumulative_pos + node_length;
                        }
                    }

                    if !translated_steps.is_empty() {
                        path_key_ranges.entry(sample_hap_name)
                        .or_default()
                        .push(RangeInfo { 
                            start, 
                            end, 
                            gfa_id,
                            steps: translated_steps,
                            step_positions,
                            step_lengths,
                        });
                    } else if args.debug {
                        eprintln!("  Warning: Path '{}' has no steps", path_name);
                    }
                }
            }
        }
    }

    if args.debug {
        eprintln!("Collected {} nodes, {} edges, {} paths, and {} path ranges from all GFA files",
            combined_graph.node_count(), combined_graph.edge_count(), path_key_ranges.len(), path_key_ranges.iter().map(|(_, ranges)| ranges.len()).sum::<usize>());
    }

    let mut next_node_id_value = u64::from(combined_graph.max_node_id()) + 1;

    // Sort ranges and create merged paths in the combined graph
    for (path_key, ranges) in path_key_ranges.iter_mut() {
        // Sort ranges by start position
        ranges.sort_by_key(|r| (r.start, r.end));

        if args.debug {
            eprintln!("Processing path key '{}'", path_key);
            for range in ranges.iter() {
                eprintln!("  Range: start={}, end={}, num.steps={}, gfa_id={}", range.start, range.end, range.steps.len(), range.gfa_id);
            }

            eprintln!("Removing redundant ranges");
        }

        // Remove ranges that are contained within other ranges
        let mut write_idx = 0;
        for read_idx in 1..ranges.len() {
            let (prev_start, prev_end) = (ranges[write_idx].start, ranges[write_idx].end);
            let (curr_start, curr_end) = (ranges[read_idx].start, ranges[read_idx].end);

            if curr_start == prev_start && curr_end == prev_end {
                // Current range is a duplicate of the previous range, skip it
                if args.debug {
                    eprintln!(
                        "  Duplicate range detected: Range [start={}, end={}] is identical to previous range and will be removed.",
                        curr_start, curr_end
                    );
                }
                continue;
            } else if curr_start >= prev_start && curr_end <= prev_end {
                // Current range is contained within previous range, skip it
                continue;
            } else if prev_start >= curr_start && prev_end <= curr_end {
                // Previous range is contained within current range
                if args.debug {
                    eprintln!(
                        "  Redundant range detected: Range [start={}, end={}, len={}] is contained within [start={}, end={} | len={}] and will be removed.",
                        prev_start, prev_end, prev_end - prev_start, curr_start, curr_end, curr_end - curr_start
                    );
                }
                ranges.swap(write_idx, read_idx);
            } else {
                // No containment, advance write_idx and copy current range
                write_idx += 1;
                if write_idx != read_idx {
                    ranges.swap(write_idx, read_idx);
                }
            }
        }
        // Truncate the ranges vector to remove any skipped elements
        ranges.truncate(write_idx + 1);

        if args.debug {
            eprintln!("Path key '{}' without redundancy", path_key);
            for range in ranges.iter() {
                eprintln!("  Range: start={}, end={}, num.steps={}, gfa_id={}", range.start, range.end, range.steps.len(), range.gfa_id);
            }
        }

        if args.debug {
            eprintln!("Trimming overlapping ranges");
        }

        // Process overlaps
        for i in 1..ranges.len() {
            let (left, right) = ranges.split_at_mut(i);
            let r1 = &mut left[left.len()-1];
            let r2 = &mut right[0];

            if has_overlap(r1, r2) {
                // Calculate the overlap region - use max/min to get precise overlap bounds
                let overlap_start = std::cmp::max(r1.start, r2.start);
                let overlap_end = std::cmp::min(r1.end, r2.end);

                if args.debug {
                    let overlap_amount = overlap_end - overlap_start;

                    eprintln!(
                        "  Overlap detected: Range1 [start={}, end={}], Range2 [start={}, end={}], Overlap [start={}, end={}], Overlap size={}",
                        r1.start, r1.end, r2.start, r2.end, overlap_start, overlap_end, overlap_amount
                    );
                }

                // Adjust r2 to remove the overlap
                let mut steps_to_remove = Vec::new();
                let mut step_to_split: Option<usize> = None;
                for (idx, &(step_start, step_end)) in r2.step_positions.iter().enumerate() {
                    if step_end <= overlap_start {
                        // if args.debug {
                        //     eprintln!("    Step {} [start={}, end={}, len={}] before overlap", idx, step_start, step_end, step_end - step_start);
                        // }
                        continue;
                    } else if step_start >= overlap_end {
                        // if args.debug {
                        //     eprintln!("    Step {} [start={}, end={}, len={}] after overlap", idx, step_start, step_end, step_end - step_start);
                        // }
                        break;
                    } else if step_start >= overlap_start && step_end <= overlap_end {
                        // if args.debug {
                        //     eprintln!("    Step {} [start={}, end={}, len={}] fully overlaps", idx, step_start, step_end, step_end - step_start);
                        // }
                        steps_to_remove.push(idx);
                    } else {
                        // if args.debug {
                        //     eprintln!("    Step {} [start={}, end={}, len={}] partially overlaps", idx, step_start, step_end, step_end - step_start);
                        // }
                        if step_to_split.is_some() {
                            panic!("Error: More than one step is partially overlapping, which is not allowed.");
                        }
                        step_to_split = Some(idx);
                    }
                }

                if args.debug && r2.start == 154089 {
                    eprintln!("    Total steps to remove: {}", steps_to_remove.len());
                    eprintln!("    Step to split: {:?}", step_to_split);   
                }

                // Initialize new vectors to store updated steps
                let mut new_steps = Vec::new();
                let mut new_step_positions = Vec::new();
                let mut new_step_lengths = Vec::new();

                // Iterate over the original steps
                for idx in 0..r2.steps.len() {
                    let step_handle = r2.steps[idx];
                    let (step_start, step_end) = r2.step_positions[idx];

                    if steps_to_remove.contains(&idx) {
                        // Skip steps to remove
                        continue;
                    } else if step_to_split == Some(idx) {
                        // Split node for the single partially overlapping step
                        let node_seq = combined_graph.sequence(step_handle).collect::<Vec<_>>();
                        let overlap_within_step_start = std::cmp::max(step_start, overlap_start);
                        let overlap_within_step_end = std::cmp::min(step_end, overlap_end);
                        
                        // Calculate offsets relative to the node sequence
                        let node_len = node_seq.len();
                        // Calculate offsets consistently regardless of strand
                        let overlap_start_offset = (overlap_within_step_start - step_start).min(node_len);
                        let overlap_end_offset = (overlap_within_step_end - step_start).min(node_len);

                        if args.debug {
                            eprintln!("    Splitting step {} [start={}, end={}, len={}] to remove overlap at [start={}, end={}]", 
                                idx, step_start, step_end, step_end - step_start, overlap_within_step_start, overlap_within_step_end);
                            eprintln!("    Overlap offsets: start={}, end={}", overlap_start_offset, overlap_end_offset);
                        }

                        if step_start < overlap_start {
                            if args.debug {
                                eprintln!("    Adding left part of step [start={}, end={}]", step_start, overlap_within_step_start);
                            }
                            assert!(overlap_start_offset > 0);

                            // Keep left part
                            let new_seq = node_seq[0..overlap_start_offset].to_vec();

                            let node_id = NodeId::from(next_node_id_value);
                            next_node_id_value += 1;
                            let new_node = combined_graph.create_handle(&new_seq, node_id);

                            new_steps.push(new_node);
                            new_step_positions.push((step_start, overlap_start));
                            new_step_lengths.push(new_seq.len());
                        } else if step_end > overlap_end {
                            if args.debug {
                                eprintln!("    Adding right part of step [start={}, end={}]", overlap_within_step_end, step_end);
                            }
                            assert!(overlap_end_offset < node_len);

                            // Keep right part
                            let new_seq = node_seq[overlap_end_offset..].to_vec();

                            let node_id = NodeId::from(next_node_id_value);
                            next_node_id_value += 1;
                            let new_node = combined_graph.create_handle(&new_seq, node_id);

                            new_steps.push(new_node);
                            new_step_positions.push((overlap_end, step_end));
                            new_step_lengths.push(new_seq.len());
                        }
                    } else {
                        // Keep steps that are not to be removed or split
                        new_steps.push(step_handle);
                        new_step_positions.push((step_start, step_end));
                        new_step_lengths.push(r2.step_lengths[idx]);
                    }
                }

                // Update r2 with the new steps
                r2.steps = new_steps;
                r2.step_positions = new_step_positions;
                r2.step_lengths = new_step_lengths;

                // Update edges for the modified steps
                for idx in 0..r2.steps.len() {
                    if idx > 0 {
                        let prev_step = r2.steps[idx - 1];
                        if !combined_graph.has_edge(prev_step, r2.steps[idx]) {
                            combined_graph.create_edge(Edge(prev_step, r2.steps[idx]));
                        }
                    }
                }

                // Update r2.start and r2.end based on the new step positions
                if !r2.step_positions.is_empty() {
                    r2.start = r2.step_positions.first().unwrap().0;
                    r2.end = r2.step_positions.last().unwrap().1;
                } else {
                    // If no steps remain, set start and end to overlap_end to effectively remove this range
                    r2.start = overlap_end;
                    r2.end = overlap_end;
                }


                if args.debug {
                    eprintln!("    Updated overlaps: Range2 [start={}, end={}]", r2.start, r2.end);
                }
            }
        }

        if args.debug {
            eprintln!("Path key '{}' without overlaps", path_key);
            for range in ranges.iter() {
                eprintln!("  Range: start={}, end={}, num.steps={}, gfa_id={}", range.start, range.end, range.steps.len(), range.gfa_id);
            }
        }

        // Check for overlaps and contiguity
        let mut all_contiguous = true;
        
        for window in ranges.windows(2) {
            if has_overlap(&window[0], &window[1]) {
                if args.debug {
                    eprintln!("Unresolved overlaps detected between ranges: [start={}, end={}] and [start={}, end={}]", 
                    window[0].start, window[0].end, window[1].start, window[1].end);
                }
                panic!("Unresolved overlaps detected in path key '{}'", path_key);
            }
            if !is_contiguous(&window[0], &window[1]) {
                all_contiguous = false;
            }
        }
                
        if !all_contiguous && args.debug {
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

        if all_contiguous {
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
                println!("Creating path: {} with {} steps", path_name, steps.len());
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
