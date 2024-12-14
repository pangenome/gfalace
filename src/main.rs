use std::{
    collections::BTreeMap,
    fs::File,
    io::{self, Write},
    path::Path,
};
use clap::Parser;
use handlegraph::{
    handle::{Handle, NodeId, Edge},
    handlegraph::*,
    mutablehandlegraph::*,
    pathhandlegraph::{
        IntoPathIds, 
        GraphPathNames, 
        GraphPathsRef,
        MutableGraphPaths,
        GraphPaths
    },
    hashgraph::HashGraph,
};
use gfa::{gfa::GFA, parser::GFAParser};
use bitvec::{bitvec, prelude::BitVec};
use tempfile::NamedTempFile;

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

fn read_gfa(gfa_path: &str, parser: &GFAParser<usize, ()>) -> io::Result<GFA<usize, ()>> {
    if gfa_path.ends_with(".gz") {
        let file = std::fs::File::open(gfa_path).map_err(|e| {
            io::Error::new(
                e.kind(),
                format!("Failed to open gzipped file '{}': {}", gfa_path, e)
            )
        })?;
        
        let (mut reader, _format) = niffler::get_reader(Box::new(file))
            .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;
        
        let mut decompressed = Vec::new();
        reader.read_to_end(&mut decompressed).map_err(|e| {
            io::Error::new(
                e.kind(),
                format!("Failed to decompress file '{}': {}", gfa_path, e)
            )
        })?;
        
        // Create temporary file in the same directory as the input file for better performance
        let parent_dir = Path::new(gfa_path).parent().unwrap_or(Path::new("."));
        let temp_file = NamedTempFile::new_in(parent_dir).map_err(|e| {
            io::Error::new(
                e.kind(),
                format!("Failed to create temporary file: {}", e)
            )
        })?;
        
        // Write decompressed data
        temp_file.as_file().write_all(&decompressed).map_err(|e| {
            io::Error::new(
                e.kind(),
                format!("Failed to write to temporary file: {}", e)
            )
        })?;
        
        // Parse GFA
        parser.parse_file(temp_file.path().to_str().ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "Invalid temporary file path"
            )
        })?).map_err(|e| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("Failed to parse GFA: {}", e)
            )
        })
    } else {
        parser.parse_file(gfa_path).map_err(|e| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("Failed to parse GFA file '{}': {}", gfa_path, e)
            )
        })
    }
}

fn write_graph_to_gfa(graph: &HashGraph, output_path: &str, nodes_to_remove: &BitVec) -> std::io::Result<()> {
    let mut file = File::create(output_path)?;
    
    // Write GFA version
    writeln!(file, "H\tVN:Z:1.0")?;
    
    // Collect nodes, excluding marked nodes
    let mut nodes: Vec<Handle> = graph.handles()
        .filter(|handle| !nodes_to_remove[u64::from(handle.id()) as usize])
        .collect();
    nodes.sort_by_key(|handle| handle.id());
    
    // Create mapping from old IDs to new sequential IDs using a Vec
    // We allocate a vector with capacity for the max node ID
    let max_id = u64::from(nodes.last().unwrap().id()) as usize;
    let mut id_mapping = vec![0; max_id + 1];
    for (new_id, handle) in nodes.iter().enumerate() {
        id_mapping[u64::from(handle.id()) as usize] = new_id + 1; // +1 to start from 1
    }
    
    // Write nodes with new IDs
    for handle in &nodes {
        let sequence = graph.sequence(*handle).collect::<Vec<_>>();
        let sequence_str = String::from_utf8(sequence)
            .unwrap_or_else(|_| String::from("N"));
        let new_id = id_mapping[u64::from(handle.id()) as usize];
        writeln!(file, "S\t{}\t{}", new_id, sequence_str)?;
    }
    
    // Collect and sort edges, excluding those connected to marked nodes
    let mut edges: Vec<Edge> = graph.edges()
        .filter(|edge| {
            !nodes_to_remove[u64::from(edge.0.id()) as usize] && 
            !nodes_to_remove[u64::from(edge.1.id()) as usize]
        })
        .collect();
    edges.sort_by(|a, b| {
        a.0.id().cmp(&b.0.id())
            .then(a.1.id().cmp(&b.1.id()))
    });
    
    // Write edges with new IDs
    for edge in edges {
        let from_id = id_mapping[u64::from(edge.0.id()) as usize];
        let to_id = id_mapping[u64::from(edge.1.id()) as usize];
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

    // Write paths with new IDs
    for path_id in paths {
        if let Some(name_iter) = graph.get_path_name(path_id) {
            let path_name = String::from_utf8(name_iter.collect::<Vec<u8>>())
                .unwrap_or_else(|_| String::from("unknown_path"));
            
            let mut path_elements = Vec::new();
            if let Some(path_ref) = graph.get_path_ref(path_id) {
                for handle in &path_ref.nodes {
                    let new_id = id_mapping[u64::from(handle.id()) as usize];
                    if new_id > 0 { // Only include nodes that weren't removed
                        let orient = if handle.is_reverse() { "-" } else { "+" };
                        path_elements.push(format!("{}{}", new_id, orient));
                    }
                }
            }
            
            if !path_elements.is_empty() {
                writeln!(file, "P\t{}\t{}\t*", path_name, path_elements.join(","))?;
            }
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
                // Skip duplicate range
                // Current range is a duplicate of the previous range, skip it
                if args.debug {
                    eprintln!(
                        "  Duplicate range detected: Range [start={}, end={}] is identical to previous range and will be removed.",
                        curr_start, curr_end
                    );
                }
                continue;
            } else if curr_start >= prev_start && curr_end <= prev_end {
                // Skip range that is fully contained within previous range
                if args.debug {
                    eprintln!(
                        "  Contained range detected: Range [start={}, end={}] is fully contained within previous range [start={}, end={}] and will be removed.",
                        curr_start, curr_end, prev_start, prev_end
                    );
                }
                continue;
            } else if prev_start >= curr_start && prev_end <= curr_end {
                // Previous range is fully contained within current range
                if args.debug {
                    eprintln!(
                        "  Containing range detected: Previous range [start={}, end={}] is fully contained within current range [start={}, end={}] and will be removed.",
                        prev_start, prev_end, curr_start, curr_end
                    );
                }
                ranges.swap(write_idx, read_idx);
            } else if curr_start < prev_end {
                // Handle overlapping ranges - check both previous and next ranges
                let mut should_skip = false;
                
                if read_idx < ranges.len() - 1 {
                    let next_start = ranges[read_idx + 1].start;
                    
                    // Check if current range is significantly overlapped by both neighbors
                    if curr_start > prev_start && next_start < curr_end {
                        let overlap_with_prev = prev_end - curr_start;
                        let overlap_with_next = curr_end - next_start;
                        let range_length = curr_end - curr_start;
                        
                        // Skip if the range is mostly covered by its neighbors
                        if overlap_with_prev + overlap_with_next > range_length {
                            should_skip = true;
                        }
                    }
                }
                
                if !should_skip {
                    if args.debug {
                        eprintln!(
                            "  Overlapping range detected: Range [start={}, end={}] overlaps with previous range [start={}, end={}] and will be kept.",
                            curr_start, curr_end, prev_start, prev_end
                        );
                    }
                    write_idx += 1;
                    if write_idx != read_idx {
                        ranges.swap(write_idx, read_idx);
                    }
                }
            } else {
                // No overlap - keep both ranges
                write_idx += 1;
                if write_idx != read_idx {
                    ranges.swap(write_idx, read_idx);
                }
            }
        }
        ranges.truncate(write_idx + 1);
        

        if args.debug {
            eprintln!("Path key '{}' without redundancy", path_key);
            for range in ranges.iter() {
                eprintln!("  Range: start={}, end={}, num.steps={}, gfa_id={}", range.start, range.end, range.steps.len(), range.gfa_id);
            }
        }

        let mut write_idx = 0;
        for read_idx in 1..ranges.len() {
            let (prev_start, prev_end) = (ranges[write_idx].start, ranges[write_idx].end);
            let (curr_start, curr_end) = (ranges[read_idx].start, ranges[read_idx].end);
        
            if curr_start == prev_start && curr_end == prev_end {
                // Skip duplicate range
                continue;
            } else if curr_start >= prev_start && curr_end <= prev_end {
                // Skip range that is fully contained within previous range
                continue;
            } else if prev_start >= curr_start && prev_end <= curr_end {
                // Previous range is fully contained within current range
                ranges.swap(write_idx, read_idx);
            } else if curr_start < prev_end {
                // Handle overlapping ranges
                if read_idx < ranges.len() - 1 {
                    let next_start = ranges[read_idx + 1].start;
                    
                    // Check if current range is significantly overlapped by both previous and next ranges
                    if curr_start > prev_start && next_start < curr_end {
                        // Calculate overlap percentages
                        let overlap_with_prev = prev_end - curr_start;
                        let overlap_with_next = curr_end - next_start;
                        let range_length = curr_end - curr_start;
                        
                        // If the range is mostly overlapped by both neighbors, skip it
                        if overlap_with_prev + overlap_with_next > range_length {
                            continue;
                        }
                    }
                }
                // Keep current range
                write_idx += 1;
                if write_idx != read_idx {
                    ranges.swap(write_idx, read_idx);
                }
            } else {
                // No overlap - keep both ranges
                write_idx += 1;
                if write_idx != read_idx {
                    ranges.swap(write_idx, read_idx);
                }
            }
        }
        ranges.truncate(write_idx + 1);
        if args.debug {
            eprintln!("Path key '{}' without redundancy 2", path_key);
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
                        // if args.debug && r2.start == 11000 {
                        //     eprintln!("    Step {} [start={}, end={}, len={}] before overlap", idx, step_start, step_end, step_end - step_start);
                        // }
                        continue;
                    } else if step_start >= overlap_end {
                        // if args.debug && r2.start == 11000 {
                        //     eprintln!("    Step {} [start={}, end={}, len={}] after overlap", idx, step_start, step_end, step_end - step_start);
                        // }
                        break;
                    } else if step_start >= overlap_start && step_end <= overlap_end {
                        // if args.debug && r2.start == 11000 {
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

    // After all paths are created but before writing the graph
    if args.debug {
        eprintln!("Total nodes before filtering: {}", combined_graph.node_count());
    }
    
    let nodes_to_skip = mark_nodes_for_removal(&combined_graph);
    
    if args.debug {
        eprintln!("Nodes to be filtered out: {}", nodes_to_skip.count_ones());
    }

    // Write the combined graph to GFA file, skipping unused nodes
    match write_graph_to_gfa(&combined_graph, &args.output, &nodes_to_skip) {
        Ok(_) => if args.debug {eprintln!("Successfully wrote combined graph to {}", args.output)},
        Err(e) => eprintln!("Error writing GFA file: {}", e),
    }
}


fn mark_nodes_for_removal(graph: &HashGraph) -> BitVec {
    // Create a bitvector with all nodes initially marked for removal
    let max_node_id = u64::from(graph.max_node_id());
    let mut nodes_to_remove = bitvec![1; max_node_id as usize + 1];
    
    // Mark nodes used in paths as not to be removed (set bit to 0)
    for path_id in graph.path_ids() {
        if let Some(path_ref) = graph.get_path_ref(path_id) {
            for handle in &path_ref.nodes {
                nodes_to_remove.set(u64::from(handle.id()) as usize, false);
            }
        }
    }
    
    nodes_to_remove
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

#[cfg(test)]
mod tests {
    use super::*;

    // Helper function to create a simple RangeInfo for testing
    fn create_range_info(start: usize, end: usize, gfa_id: usize) -> RangeInfo {
        RangeInfo {
            start,
            end,
            gfa_id,
            steps: vec![],           // Empty steps for testing
            step_positions: vec![],   // Empty positions for testing
            step_lengths: vec![],     // Empty lengths for testing
        }
    }

    #[test]
    fn test_range_containment_removal() {
        // Test cases
        let test_cases = vec![
            // Test case 1: Basic containment
            (
                vec![(10, 50, 0), (20, 30, 1)],  // Input ranges (start, end, gfa_id)
                vec![(10, 50, 0)],               // Expected result
                "Basic containment"
            ),
            
            // Test case 2: No containment
            (
                vec![(10, 20, 0), (30, 40, 1)],
                vec![(10, 20, 0), (30, 40, 1)],
                "No containment"
            ),
            
            // Test case 3: Multiple contained ranges
            (
                vec![(10, 100, 0), (20, 30, 1), (40, 50, 2), (60, 70, 3)],
                vec![(10, 100, 0)],
                "Multiple contained ranges"
            ),
            
            // Test case 4: Identical ranges
            (
                vec![(10, 20, 0), (10, 20, 1)],
                vec![(10, 20, 0)],
                "Identical ranges"
            ),
            
            // Test case 5: Nested containment
            (
                vec![(10, 100, 0), (20, 80, 1), (30, 40, 2)],
                vec![(10, 100, 0)],
                "Nested containment"
            ),
            
            // Test case 6: Partial overlap (should keep both)
            (
                vec![(10, 30, 0), (20, 40, 1)],
                vec![(10, 30, 0), (20, 40, 1)],
                "Partial overlap"
            ),
            
            // Test case 7: Edge cases - touching ranges
            (
                vec![(10, 20, 0), (20, 30, 1)],
                vec![(10, 20, 0), (20, 30, 1)],
                "Touching ranges"
            ),

            // Test case 8: Overlapping ranges from same GFA
            (
                vec![(0, 11742, 0), (9714, 13000, 1), (11000, 19000, 1)],
                vec![(0, 11742, 0), (11000, 19000, 1)],
                "Overlapping ranges from same GFA"
            ),

            // Test case 9: Overlapping ranges with different GFA IDs
            (
                vec![(0, 11742, 0), (9714, 13000, 1), (11000, 19000, 2)],
                vec![(0, 11742, 0), (11000, 19000, 2)],
                "Overlapping ranges"
            ),

            // Test case 10: Overlapping ranges with different GFA IDs 2
            (
                vec![(0, 10, 0), (8, 20, 1), (15, 30, 2)],
                vec![(0, 10, 0), (8, 20, 1), (15, 30, 2)],
                "Overlapping ranges"
            ),

            // Test case 11: Overlapping ranges with different GFA IDs 3
            (
                vec![(8000, 11000, 0), (9694, 12313, 1), (10908, 13908, 2)],
                vec![(8000, 11000, 0), (10908, 13908, 2)],
                "Overlapping ranges"
            ),
        ];

        // Run each test case
        for (case_index, (input_ranges, expected_ranges, case_name)) in test_cases.iter().enumerate() {
            println!("Running test case {}: {}", case_index + 1, case_name);
            
            // Create input ranges
            let mut ranges: Vec<RangeInfo> = input_ranges
                .iter()
                .map(|(start, end, gfa_id)| create_range_info(*start, *end, *gfa_id))
                .collect();
// Sort ranges by start position
ranges.sort_by_key(|r| (r.start, r.end));

let mut write_idx = 0;
for read_idx in 1..ranges.len() {
    let (prev_start, prev_end) = (ranges[write_idx].start, ranges[write_idx].end);
    let (curr_start, curr_end) = (ranges[read_idx].start, ranges[read_idx].end);

    if curr_start == prev_start && curr_end == prev_end {
        // Skip duplicate range
        continue;
    } else if curr_start >= prev_start && curr_end <= prev_end {
        // Skip range that is fully contained within previous range
        continue;
    } else if prev_start >= curr_start && prev_end <= curr_end {
        // Previous range is fully contained within current range
        ranges.swap(write_idx, read_idx);
    } else if curr_start < prev_end {
        // Handle overlapping ranges - check both previous and next ranges
        let mut should_skip = false;
        
        if read_idx < ranges.len() - 1 {
            let next_start = ranges[read_idx + 1].start;
            
            // Check if current range is significantly overlapped by both neighbors
            if curr_start > prev_start && next_start < curr_end {
                let overlap_with_prev = prev_end - curr_start;
                let overlap_with_next = curr_end - next_start;
                let range_length = curr_end - curr_start;
                
                // Skip if the range is mostly covered by its neighbors
                if overlap_with_prev + overlap_with_next > range_length {
                    should_skip = true;
                }
            }
        }
        
        if !should_skip {
            write_idx += 1;
            if write_idx != read_idx {
                ranges.swap(write_idx, read_idx);
            }
        }
    } else {
        // No overlap - keep both ranges
        write_idx += 1;
        if write_idx != read_idx {
            ranges.swap(write_idx, read_idx);
        }
    }
}
ranges.truncate(write_idx + 1);

            // Create expected ranges
            let expected: Vec<RangeInfo> = expected_ranges
                .iter()
                .map(|(start, end, gfa_id)| create_range_info(*start, *end, *gfa_id))
                .collect();

            // Compare results
            assert_eq!(
                ranges.len(),
                expected.len(),
                "Test case '{}': Wrong number of ranges after containment removal",
                case_name
            );

            for (i, (result, expected)) in ranges.iter().zip(expected.iter()).enumerate() {
                assert_eq!(
                    (result.start, result.end, result.gfa_id),
                    (expected.start, expected.end, expected.gfa_id),
                    "Test case '{}': Mismatch at position {}",
                    case_name,
                    i
                );
            }
            
            println!("Test case {} passed: {}", case_index + 1, case_name);
        }
    }
}
