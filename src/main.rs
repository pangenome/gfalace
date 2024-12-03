use std::collections::{HashMap, BTreeMap};
use clap::Parser;
use handlegraph::handle::{Handle, NodeId, Edge};
use handlegraph::handlegraph::*;
use handlegraph::mutablehandlegraph::*;
use handlegraph::pathhandlegraph::{IntoPathIds, GraphPathNames, GraphPathsRef, MutableGraphPaths, PathSteps, PathStep};
use handlegraph::hashgraph::HashGraph;
//use handlegraph::pathhandlegraph::PathStep;
use gfa::{gfa::GFA, parser::GFAParser};

#[derive(Parser, Debug)]
#[clap(author, version, about)]
struct Args {
    /// List of GFA file paths.
    #[clap(short, long, value_parser)]
    gfa_list: Vec<String>,
}

#[derive(Debug)]
struct RangeInfo {
    start: usize,
    end: usize,
    gfa_id: usize,
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

        // Process paths and collect ranges
        for path_id in block_graph.path_ids() {
            if let Some(name_iter) = block_graph.get_path_name(path_id) {
                let path_name = String::from_utf8(name_iter.collect::<Vec<u8>>()).unwrap();
                
                if let Some((sample_hap_name, start, end)) = split_path_name(&path_name) {
                    path_key_ranges.entry(sample_hap_name)
                        .or_insert_with(Vec::new)
                        .push(RangeInfo { start, end, gfa_id });
                }
            }
        }
    }

    // Sort ranges for each key
    for ranges in path_key_ranges.values_mut() {
        ranges.sort_by_key(|r| (r.start, r.end));
    }

    // Print some statistics about the combined graph
    println!("Combined graph statistics:");
    println!("  Nodes: {}", combined_graph.node_count());
    println!("  Edges: {}", combined_graph.edge_count());
    println!("Number of unique path keys: {}", path_key_ranges.len());
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
