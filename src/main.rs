use std::collections::HashMap;
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

fn main() {
    let args = Args::parse();
    let smoothed_block_graphs = args
        .gfa_list
        .iter()
        .map(|gfa_path| {
            let parser = GFAParser::new();
            let gfa: GFA<usize, ()> = parser.parse_file(gfa_path).unwrap();
            HashGraph::from_gfa(&gfa)
        })
        .collect::<Vec<_>>();

    let path_ranges = args
        .gfa_list
        .iter()
        .flat_map(|gfa_path| {
            let parser = GFAParser::new();
            let gfa: GFA<usize, ()> = parser.parse_file(gfa_path).unwrap();
            let graph = HashGraph::from_gfa(&gfa);

            graph.path_ids()
                .filter_map(|path_id| {
                    // Assuming get_path_name returns an Option of an iterator, you can do this
                    let path_name = if let Some(name_iter) = graph.get_path_name(path_id) {
                        name_iter.collect::<Vec<u8>>()
                    } else {
                        // Handle the case where get_path_name returns None, perhaps by returning an error or using a default value
                        Vec::new() // For example, using an empty Vec<u8> as a default
                    };
                    // Then, you can convert the Vec<u8> into a &str
                    let path_name_str = std::str::from_utf8(&path_name).unwrap();
                    if let Some((_genome, _hap, _chr, _start, _end)) = parse_path_name(path_name_str) {
                        let steps: Vec<LacePathStep> = graph
                            .get_path_ref(path_id)?
                            .steps()
                            .map(|step| LacePathStep {
                                block_id: smoothed_block_graphs.len() - 1,
                                node_id: step.handle().id(),
                            })
                            .collect();

                        Some(PathRange {
                            path_name: path_name_str.to_string(),
                            steps,
                        })
                    } else {
                        None
                    }
                })
                .collect::<Vec<_>>()
        })
        //.flatten()
        .collect::<Vec<_>>();

    let _smoothed_graph = lace_smoothed_blocks(&smoothed_block_graphs, &path_ranges);
    // ... (further processing or output of the smoothed graph)
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

fn parse_path_name(path_name: &str) -> Option<(String, usize, String, usize, usize)> {
    let parts: Vec<&str> = path_name.split('#').collect();
    if parts.len() == 3 {
        let genome = parts[0].to_string();
        let hap = parts[1].parse().ok()?;
        let chr_range: Vec<&str> = parts[2].split(':').collect();
        if chr_range.len() == 2 {
            let chr = chr_range[0].to_string();
            let range: Vec<&str> = chr_range[1].split('-').collect();
            if range.len() == 2 {
                let start = range[0].parse().ok()?;
                let end = range[1].parse().ok()?;
                return Some((genome, hap, chr, start, end));
            }
        }
    }
    None
}
