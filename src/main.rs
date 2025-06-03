use bitvec::{bitvec, prelude::BitVec};
use clap::Parser;
use handlegraph::handle::{Handle, NodeId};
use log::{debug, error, info, warn};
use niffler::compression::Format;
use rust_htslib::faidx;
use rustc_hash::{FxHashMap, FxHashSet};
use std::{
    fs::File,
    io::{self, BufRead, BufReader, BufWriter, Read, Seek, SeekFrom, Write},
};
use tempfile::NamedTempFile;

use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use std::num::NonZeroUsize;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};

use gzp::{
    deflate::{Bgzf, Gzip}, // Both Gzip and Bgzf are in deflate module
    par::compress::{ParCompress, ParCompressBuilder},
};
use zstd::stream::Encoder as ZstdEncoder;

// use std::process::Command;

// #[cfg(not(debug_assertions))]
// fn log_memory_usage(stage: &str) {
//     let output = Command::new("ps")
//         .args(&["-o", "rss=", "-p", &std::process::id().to_string()])
//         .output()
//         .expect("Failed to execute ps command");

//     let memory_kb = String::from_utf8_lossy(&output.stdout)
//         .trim()
//         .parse::<u64>()
//         .unwrap_or(0);

//     let memory_mb = memory_kb as f64 / 1024.0;
//     info!("Memory usage at {}: {:.2} MB", stage, memory_mb);
// }

fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&base| match base {
            b'A' | b'a' => b'T',
            b'T' | b't' => b'A',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'U' | b'u' => b'A',
            _ => b'N',
        })
        .collect()
}

#[derive(Parser, Debug)]
#[clap(author, version, about)]
struct Args {
    /// List of GFA file paths to combine (can use wildcards)
    #[clap(short = 'g', long, value_parser, num_args = 1.., value_delimiter = ' ', conflicts_with = "gfa_list")]
    gfa_files: Option<Vec<String>>,

    /// Text file containing list of GFA file paths (one per line)
    #[clap(short = 'l', long, value_parser, conflicts_with = "gfa_files")]
    gfa_list: Option<String>,

    /// Output GFA file path for the combined graph
    #[clap(short, long, value_parser)]
    output: String,

    /// Output compression format: none, gzip (.gz), bgzip (.bgz), or zstd (.zst)
    #[clap(long, value_parser, default_value = "auto")]
    compress: String,

    /// Gap filling mode: 0=none, 1=middle gaps only, 2=all gaps (requires --fasta for end gaps)
    #[clap(long, default_value = "0")]
    fill_gaps: u8,

    /// List of FASTA file paths for gap filling
    #[clap(long, value_parser, num_args = 1.., value_delimiter = ' ', conflicts_with_all = &["fasta_list"])]
    fasta_files: Option<Vec<String>>,

    /// Text file containing list of FASTA file paths (one per line)
    #[clap(long, value_parser, conflicts_with_all = &["fasta_files"])]
    fasta_list: Option<String>,

    /// Temporary directory for temporary files (default: same as input files)
    #[clap(long, value_parser)]
    temp_dir: Option<String>,

    /// Number of threads for parallel processing.
    #[clap(short = 't', long, value_parser, default_value_t = NonZeroUsize::new(4).unwrap())]
    num_threads: NonZeroUsize,

    /// Verbosity level (0 = error, 1 = info, 2 = debug)
    #[clap(short, long, default_value = "0")]
    verbose: u8,
}

fn get_compression_format(compress_arg: &str, output_path: &str) -> Format {
    match compress_arg.to_lowercase().as_str() {
        "none" => Format::No,
        "gzip" | "gz" => Format::Gzip,
        "bgzip" | "bgz" => Format::Bzip,
        "zstd" | "zst" => Format::Zstd,
        "auto" => {
            // Auto-detect based on file extension
            if output_path.ends_with(".gz") {
                Format::Gzip
            } else if output_path.ends_with(".bgz") {
                Format::Bzip
            } else if output_path.ends_with(".zst") {
                Format::Zstd
            } else {
                Format::No
            }
        }
        _ => {
            warn!(
                "Unsupported compression format '{}', using none",
                compress_arg
            );
            Format::No
        }
    }
}

fn main() {
    let args = Args::parse();

    // Initialize logger based on verbosity
    env_logger::Builder::new()
        .filter_level(match args.verbose {
            0 => log::LevelFilter::Error,
            1 => log::LevelFilter::Info,
            _ => log::LevelFilter::Debug,
        })
        .init();

    // Configure thread pool
    ThreadPoolBuilder::new()
        .num_threads(args.num_threads.into())
        .build_global()
        .unwrap();

    // Determine compression format
    let compression_format = get_compression_format(&args.compress, &args.output);

    // Get the list of GFA files
    let gfa_files = match (args.gfa_files, args.gfa_list) {
        (Some(files), None) => files,
        (None, Some(list_file)) => {
            // Read file paths from the list file
            match std::fs::read_to_string(&list_file) {
                Ok(content) => content
                    .lines()
                    .filter(|line| !line.trim().is_empty() && !line.trim().starts_with('#'))
                    .map(|line| line.trim().to_string())
                    .collect(),
                Err(e) => {
                    error!("Failed to read GFA list file '{}': {}", list_file, e);
                    std::process::exit(1);
                }
            }
        }
        (None, None) => {
            error!("Either --gfa-files (-g) or --gfa-list (-l) must be specified");
            std::process::exit(1);
        }
        (Some(_), Some(_)) => {
            error!("Cannot specify both --gfa-files and --gfa-list");
            std::process::exit(1);
        }
    };
    if gfa_files.is_empty() {
        error!("No GFA files specified");
        std::process::exit(1);
    }

    // Build FASTA index if FASTA files are provided
    let fasta_index = if args.fill_gaps > 0 {
        // Get list of FASTA files
        let fasta_files = match (args.fasta_files, args.fasta_list) {
            // Handle --fasta-files option
            (Some(files), None) => files,
            // Handle --fasta-list option
            (None, Some(list_file)) => match std::fs::read_to_string(&list_file) {
                Ok(content) => content
                    .lines()
                    .filter(|line| !line.trim().is_empty() && !line.trim().starts_with('#'))
                    .map(|line| line.trim().to_string())
                    .collect(),
                Err(e) => {
                    error!("Failed to read FASTA list file '{}': {}", list_file, e);
                    std::process::exit(1);
                }
            },
            _ => {
                if args.fill_gaps == 2 {
                    warn!("Gap filling mode 2 requires FASTA files; end gaps will not be filled");
                }
                vec![]
            }
        };

        if fasta_files.is_empty() {
            None
        } else {
            match FastaIndex::build_from_files(&fasta_files) {
                Ok(index) => {
                    info!(
                        "Built FASTA index for {} files with {} sequences",
                        index.fasta_paths.len(),
                        index.path_key_to_fasta.len()
                    );
                    Some(index)
                }
                Err(e) => {
                    error!("Failed to build FASTA index: {}", e);
                    std::process::exit(1);
                }
            }
        }
    } else {
        None
    };

    // log_memory_usage("start");

    // Create a single combined graph without paths and a map of path key to ranges
    info!("Collecting metadata from {} GFA files", gfa_files.len());
    let (combined_graph, mut path_key_ranges) =
        read_gfa_files(&gfa_files, args.temp_dir.as_deref()).unwrap_or_else(|e| {
            error!("Failed to read GFA files: {}", e);
            std::process::exit(1);
        });

    // log_memory_usage("after_reading_files");

    // PASS 1: sort + dedup — purely local, so go parallel
    info!(
        "Sorting and deduplicating {} path ranges from {} path keys",
        path_key_ranges
            .values()
            .map(|ranges| ranges.len())
            .sum::<usize>(),
        path_key_ranges.len()
    );
    path_key_ranges
        .par_iter_mut()
        .for_each(|(_path_key, ranges)| {
            sort_and_filter_ranges(ranges);
        });

    // PASS 2: Process different path keys in parallel with minimal locking
    info!(
        "Trimming overlaps and linking contiguous ranges for {} path ranges from {} path keys",
        path_key_ranges.len(),
        path_key_ranges.len()
    );

    // Wrap graph in Arc<Mutex> for thread-safe access
    let graph_mutex = Arc::new(Mutex::new(combined_graph));

    path_key_ranges
        .par_iter_mut()
        .for_each(|(_path_key, ranges)| {
            trim_range_overlaps(ranges, &graph_mutex);
            link_contiguous_ranges(ranges, &graph_mutex);
        });

    // Unwrap the graph from Arc<Mutex>
    let mut combined_graph = Arc::try_unwrap(graph_mutex)
        .ok()
        .expect("Failed to unwrap graph mutex")
        .into_inner()
        .unwrap();

    info!(
        "Created {} nodes and {} edges",
        combined_graph.node_count,
        combined_graph.edges.len()
    );

    // log_memory_usage("before_writing");

    match write_graph_to_gfa(
        &mut combined_graph,
        &path_key_ranges,
        &args.output,
        compression_format,
        args.fill_gaps,
        &fasta_index,
        args.verbose > 1,
    ) {
        Ok(_) => info!(
            "Successfully wrote the combined graph to {} ({:?} format)",
            args.output, compression_format
        ),
        Err(e) => error!("Error writing the GFA file: {}", e),
    }

    // log_memory_usage("end");
}

// Compact edge representation using bit-packed orientations
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct CompactEdge {
    // Use top 2 bits for orientations, rest for node IDs
    from: u64, // bit 63: orientation, bits 0-62: node ID
    to: u64,   // bit 63: orientation, bits 0-62: node ID
}

impl CompactEdge {
    const ORIENT_MASK: u64 = 1u64 << 63;
    const ID_MASK: u64 = !Self::ORIENT_MASK;

    fn new(from_id: u64, from_rev: bool, to_id: u64, to_rev: bool) -> Self {
        let from = from_id | (if from_rev { Self::ORIENT_MASK } else { 0 });
        let to = to_id | (if to_rev { Self::ORIENT_MASK } else { 0 });
        CompactEdge { from, to }
    }

    fn from_id(&self) -> u64 {
        self.from & Self::ID_MASK
    }
    fn from_rev(&self) -> bool {
        (self.from & Self::ORIENT_MASK) != 0
    }
    fn to_id(&self) -> u64 {
        self.to & Self::ID_MASK
    }
    fn to_rev(&self) -> bool {
        (self.to & Self::ORIENT_MASK) != 0
    }
}

// Sequence storage with memory mapping
struct SequenceStore {
    sequences_file: File,
    offsets: Vec<(u64, u32)>, // (offset, length) for each node
}

impl SequenceStore {
    fn new(temp_dir: Option<&str>) -> io::Result<Self> {
        let temp_file = if let Some(dir) = temp_dir {
            NamedTempFile::new_in(dir)?
        } else {
            NamedTempFile::new()?
        };

        let sequences_file = temp_file.into_file();

        Ok(SequenceStore {
            sequences_file,
            offsets: Vec::new(),
        })
    }

    fn add_sequence(&mut self, seq: &[u8]) -> io::Result<usize> {
        let offset = self.sequences_file.seek(SeekFrom::End(0))?;
        self.sequences_file.write_all(seq)?;

        let idx = self.offsets.len();
        self.offsets.push((offset, seq.len() as u32));
        Ok(idx)
    }

    fn get_sequence(&mut self, idx: usize) -> io::Result<Vec<u8>> {
        if idx >= self.offsets.len() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "Invalid sequence index",
            ));
        }

        let (offset, length) = self.offsets[idx];
        let mut buffer = vec![0u8; length as usize];

        self.sequences_file.seek(SeekFrom::Start(offset))?;
        self.sequences_file.read_exact(&mut buffer)?;

        Ok(buffer)
    }
}

// Simplified graph structure
struct CompactGraph {
    node_count: u64,
    edges: FxHashSet<CompactEdge>,
    sequence_store: SequenceStore,
}

impl CompactGraph {
    fn new(temp_dir: Option<&str>) -> io::Result<Self> {
        Ok(CompactGraph {
            node_count: 0,
            edges: FxHashSet::default(),
            sequence_store: SequenceStore::new(temp_dir)?,
        })
    }

    fn add_node(&mut self, seq: &[u8]) -> io::Result<u64> {
        self.sequence_store.add_sequence(seq)?;
        self.node_count += 1;
        Ok(self.node_count) // Return the actual node ID (1-based as per handlegraph convention)
    }

    fn add_edge(&mut self, from_id: u64, from_rev: bool, to_id: u64, to_rev: bool) {
        self.edges
            .insert(CompactEdge::new(from_id, from_rev, to_id, to_rev));
    }

    fn has_edge(&self, from_id: u64, from_rev: bool, to_id: u64, to_rev: bool) -> bool {
        // Check for the edge in both directions since handlegraph edges are bidirectional
        let edge1 = CompactEdge::new(from_id, from_rev, to_id, to_rev);
        let edge2 = CompactEdge::new(to_id, !to_rev, from_id, !from_rev);
        self.edges.contains(&edge1) || self.edges.contains(&edge2)
    }

    fn get_sequence(&mut self, handle: Handle) -> io::Result<Vec<u8>> {
        let seq = self
            .sequence_store
            .get_sequence((u64::from(handle.id()) - 1) as usize)?;
        Ok(if handle.is_reverse() {
            reverse_complement(&seq)
        } else {
            seq
        })
    }
}

#[derive(Debug, Clone)]
struct RangeInfo {
    start: usize,
    end: usize,
    //gfa_id: usize,          // GFA file ID this range belongs to
    steps: Vec<Handle>, // Path steps for this range
}
impl RangeInfo {
    /// Returns true if this range is immediately followed by another range
    /// with no gap between them
    #[inline]
    fn is_contiguous_with(&self, other: &Self) -> bool {
        self.end == other.start
    }

    /// Returns true if this range overlaps with another range
    /// Two ranges overlap if one starts before the other ends
    #[inline]
    fn overlaps_with(&self, other: &Self) -> bool {
        self.start < other.end && other.start < self.end
    }
}

// Structure to manage multiple FASTA files
struct FastaIndex {
    fasta_paths: Vec<String>,
    path_key_to_fasta: FxHashMap<String, usize>,
}

impl FastaIndex {
    fn new() -> Self {
        FastaIndex {
            fasta_paths: Vec::new(),
            path_key_to_fasta: FxHashMap::default(),
        }
    }

    fn build_from_files(fasta_files: &[String]) -> io::Result<Self> {
        let mut index = FastaIndex::new();

        for (fasta_idx, fasta_path) in fasta_files.iter().enumerate() {
            index.fasta_paths.push(fasta_path.clone());

            // Read the .fai file to get sequence names
            let fai_path = format!("{}.fai", fasta_path);

            // Try to open the .fai file, if it doesn't exist, try to create it
            let fai_content = match std::fs::read_to_string(&fai_path) {
                Ok(content) => content,
                Err(_) => {
                    // Try to create the index using rust-htslib
                    match faidx::Reader::from_path(fasta_path) {
                        Ok(_) => {
                            // Index was created, now read it
                            std::fs::read_to_string(&fai_path)?
                        }
                        Err(e) => {
                            return Err(io::Error::other(
                                format!("Failed to create FASTA index for '{}': {}", fasta_path, e),
                            ));
                        }
                    }
                }
            };

            // Parse the .fai file to get sequence names
            for line in fai_content.lines() {
                if let Some(seq_name) = line.split('\t').next() {
                    if !seq_name.is_empty() {
                        index
                            .path_key_to_fasta
                            .insert(seq_name.to_string(), fasta_idx);
                    }
                }
            }
        }

        Ok(index)
    }

    fn get_fasta_path(&self, path_key: &str) -> Option<&str> {
        self.path_key_to_fasta
            .get(path_key)
            .map(|&idx| self.fasta_paths[idx].as_str())
    }
}

fn read_gfa_files(
    gfa_list: &[String],
    temp_dir: Option<&str>,
) -> io::Result<(CompactGraph, FxHashMap<String, Vec<RangeInfo>>)> {
    // Shared structures protected by Mutex, wrapped in Arc for thread‑safe sharing
    let combined_graph = Arc::new(Mutex::new(CompactGraph::new(temp_dir)?));
    let path_key_ranges: Arc<Mutex<FxHashMap<String, Vec<RangeInfo>>>> =
        Arc::new(Mutex::new(FxHashMap::default()));

    // Cheap thread‑safe counters just for stats
    let num_path_ranges = AtomicUsize::new(0);
    let num_path_range_steps = AtomicUsize::new(0);

    // Process each file in parallel without loading everything in memory
    gfa_list
        .par_iter()
        .enumerate()
        .for_each(|(gfa_id, gfa_path)| {
            if let Ok(reader) = get_gfa_reader(gfa_path) {
                let mut id_translation: FxHashMap<NodeId, NodeId> = FxHashMap::default();
                let mut temp_edges: Vec<CompactEdge> = Vec::new();

                let mut node_count = 0usize;
                let mut edge_count = 0usize;

                for line in reader.lines().flatten() {
                    let line = line.trim();
                    if line.is_empty() || line.starts_with('#') {
                        continue;
                    }

                    let fields: Vec<&str> = line.split('\t').collect();
                    if fields.is_empty() {
                        continue;
                    }

                    match fields[0] {
                        "S" => {
                            // Segment line: S <sid> <seq> [<tag>]*
                            if fields.len() < 3 {
                                warn!("Invalid S line: {}", line);
                                continue;
                            }
                            let node_id: u64 = fields[1].parse().unwrap_or_else(|_| {
                                panic!("Invalid node ID: {}", fields[1]);
                            });
                            let sequence = fields[2].as_bytes();

                            // one node at a time – keeps memory low
                            let new_node_id = {
                                let mut graph = combined_graph.lock().unwrap();
                                graph.add_node(sequence).unwrap()
                            };
                            id_translation.insert(NodeId::from(node_id), NodeId::from(new_node_id));
                            node_count += 1;
                        }
                        "L" => {
                            // Link line: L <sid1> <orient1> <sid2> <orient2> <overlap>
                            if fields.len() < 6 {
                                warn!("Invalid L line: {}", line);
                                continue;
                            }

                            let from_id: u64 = fields[1].parse().unwrap_or_else(|_| {
                                panic!("Invalid from node ID: {}", fields[1]);
                            });
                            let from_rev = fields[2] == "-";
                            let to_id: u64 = fields[3].parse().unwrap_or_else(|_| {
                                panic!("Invalid to node ID: {}", fields[3]);
                            });
                            let to_rev = fields[4] == "-";
                            temp_edges.push(CompactEdge::new(from_id, from_rev, to_id, to_rev));
                            edge_count += 1;
                        }
                        "P" => {
                            // Path line: P <pname> <nodes> <overlaps>
                            if fields.len() < 3 {
                                warn!("Invalid P line: {}", line);
                                continue;
                            }
                            let path_name = fields[1];
                            let nodes_str = fields[2];

                            if let Some((sample_hap_name, start, end)) = split_path_name(path_name)
                            {
                                // Parse path steps
                                let mut translated_steps = Vec::new();
                                for step_str in nodes_str.split(',') {
                                    if step_str.is_empty() {
                                        continue;
                                    }
                                    let (node_str, orient) = if step_str.ends_with('+') {
                                        (&step_str[..step_str.len() - 1], false)
                                    } else if step_str.ends_with('-') {
                                        (&step_str[..step_str.len() - 1], true)
                                    } else {
                                        warn!("Invalid step format: {}", step_str);
                                        continue;
                                    };

                                    let node_id: u64 = node_str.parse().unwrap_or_else(|_| {
                                        panic!("Invalid node ID in path: {}", node_str);
                                    });

                                    // Use the translation map to get the new node ID
                                    if let Some(&translated_id) =
                                        id_translation.get(&NodeId::from(node_id))
                                    {
                                        translated_steps.push(Handle::pack(translated_id, orient));
                                    } else {
                                        warn!(
                                            "Node {} in path {} not found in translation map",
                                            node_id, path_name
                                        );
                                    }
                                }
                                if !translated_steps.is_empty() {
                                    num_path_ranges.fetch_add(1, Ordering::Relaxed);
                                    num_path_range_steps
                                        .fetch_add(translated_steps.len(), Ordering::Relaxed);

                                    let mut map = path_key_ranges.lock().unwrap();
                                    map.entry(sample_hap_name.to_string()).or_default().push(
                                        RangeInfo {
                                            start,
                                            end,
                                            steps: translated_steps,
                                        },
                                    );
                                }
                            }
                        }
                        _ => {
                            // Skip other line types (H, C, etc.)
                        }
                    }
                }

                // Add edges with translated IDs
                for edge in temp_edges {
                    if let (Some(&from_id), Some(&to_id)) = (
                        id_translation.get(&NodeId::from(edge.from_id())),
                        id_translation.get(&NodeId::from(edge.to_id())),
                    ) {
                        let mut graph = combined_graph.lock().unwrap();
                        graph.add_edge(
                            from_id.into(),
                            edge.from_rev(),
                            to_id.into(),
                            edge.to_rev(),
                        );
                    }
                }

                debug!(
                    "GFA file {} ({}) processed: {} nodes, {} edges",
                    gfa_id, gfa_path, node_count, edge_count
                );
            } else {
                error!("Failed to open GFA file '{}'", gfa_path);
                std::process::exit(1);
            }
        });

    // Unwrap the Arc/Mutex wrapper for path_key_ranges now that all threads have finished
    let path_map = Arc::try_unwrap(path_key_ranges)
        .ok()
        .expect("More than one Arc pointer to path map")
        .into_inner()
        .unwrap();

    // Validate all path ranges in parallel using the still-wrapped graph
    info!("Validating path range lengths");
    let validation_errors: Vec<String> = path_map
        .par_iter()
        .flat_map(|(path_key, ranges)| {
            ranges
                .par_iter()
                .filter_map(|range| {
                    let expected_length = range.end - range.start;
                    let mut actual_length = 0;

                    // Compute actual length by summing step lengths
                    for &step_handle in &range.steps {
                        let seq_result = {
                            let mut graph = combined_graph.lock().unwrap();
                            graph.get_sequence(step_handle)
                        };

                        match seq_result {
                            Ok(seq) => actual_length += seq.len(),
                            Err(e) => {
                                return Some(format!(
                                    "Failed to get sequence for node {} in path {}: {}",
                                    step_handle.id(),
                                    path_key,
                                    e
                                ));
                            }
                        }
                    }

                    if expected_length != actual_length {
                        Some(format!(
                            "Path range length mismatch for '{}:{}-{}': \
                             expected length {} but sum of step lengths is {}",
                            path_key, range.start, range.end, expected_length, actual_length
                        ))
                    } else {
                        None
                    }
                })
                .collect::<Vec<String>>()
        })
        .collect();

    // Check if there were any validation errors
    if !validation_errors.is_empty() {
        for error in validation_errors {
            error!("{}", error);
        }
        std::process::exit(1);
    }

    // Now unwrap the graph after validation is complete
    let graph = Arc::try_unwrap(combined_graph)
        .ok()
        .expect("More than one Arc pointer to graph")
        .into_inner()
        .unwrap();

    info!(
        "Collected {} nodes, {} edges, {} path keys, {} path ranges and {} path steps",
        graph.node_count,
        graph.edges.len(),
        path_map.len(),
        num_path_ranges.load(Ordering::Relaxed),
        num_path_range_steps.load(Ordering::Relaxed)
    );

    Ok((graph, path_map))
}

fn get_gfa_reader(gfa_path: &str) -> io::Result<Box<dyn BufRead>> {
    let file = std::fs::File::open(gfa_path).map_err(|e| {
        io::Error::new(
            e.kind(),
            format!("Failed to open file '{}': {}", gfa_path, e),
        )
    })?;

    // Use niffler to handle both compressed and uncompressed files
    let (reader, _format) = niffler::get_reader(Box::new(file)).map_err(|e| {
        io::Error::other(
            format!("Failed to open reader for '{}': {}", gfa_path, e),
        )
    })?;

    // Return a BufReader wrapping the niffler reader
    Ok(Box::new(BufReader::new(reader)))
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

fn sort_and_filter_ranges(ranges: &mut Vec<RangeInfo>) {
    // Sort ranges by start position
    ranges.sort_by_key(|r| (r.start, r.end));

    debug!("  Removing redundant ranges");

    // Remove ranges that are contained within other ranges
    let mut write_idx = 0;
    for read_idx in 1..ranges.len() {
        let (prev_start, prev_end) = (ranges[write_idx].start, ranges[write_idx].end);
        let (curr_start, curr_end) = (ranges[read_idx].start, ranges[read_idx].end);

        if curr_start == prev_start && curr_end == prev_end {
            // Skip duplicate range
            // Current range is a duplicate of the previous range, skip it
            debug!(
                "    Duplicate range detected: Range [start={}, end={}] is identical to previous range and will be removed.",
                curr_start, curr_end
            );

            continue;
        } else if curr_start >= prev_start && curr_end <= prev_end {
            // Skip range that is fully contained within previous range
            debug!(
                "    Contained range detected: Range [start={}, end={}] is fully contained within previous range [start={}, end={}] and will be removed.",
                curr_start, curr_end, prev_start, prev_end
            );

            continue;
        } else if prev_start >= curr_start && prev_end <= curr_end {
            // Previous range is fully contained within current range
            debug!(
                "    Containing range detected: Previous range [start={}, end={}] is fully contained within current range [start={}, end={}] and will be removed.",
                prev_start, prev_end, curr_start, curr_end
            );

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
                debug!(
                    "    Overlapping range detected: Range [start={}, end={}] overlaps with previous range [start={}, end={}] and will be kept.",
                    curr_start, curr_end, prev_start, prev_end
                );
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

    // if debug {
    //     debug!("  Path key '{}' without redundancy", path_key);
    //     for range in ranges.iter() {
    //         //debug!("    Range: start={}, end={}, num.steps={}, gfa_id={}", range.start, range.end, range.steps.len(), range.gfa_id);
    //         debug!("    Range: start={}, end={}, num.steps={}", range.start, range.end, range.steps.len());
    //     }
    // }
}

fn trim_range_overlaps(ranges: &mut [RangeInfo], graph_mutex: &Arc<Mutex<CompactGraph>>) {
    debug!("  Trimming overlapping ranges");

    for i in 1..ranges.len() {
        let (left, right) = ranges.split_at_mut(i);
        let r1 = &mut left[left.len() - 1];
        let r2 = &mut right[0];

        if r1.overlaps_with(r2) {
            // Calculate the overlap region - use max/min to get precise overlap bounds
            let overlap_start = std::cmp::max(r1.start, r2.start);
            let overlap_end = std::cmp::min(r1.end, r2.end);

            debug!(
                "    Overlap detected: Range1 [start={}, end={}], Range2 [start={}, end={}], Overlap [start={}, end={}, size={}]",
                r1.start, r1.end, r2.start, r2.end, overlap_start, overlap_end, overlap_end - overlap_start
            );

            // Adjust r2 to remove the overlap
            let mut steps_to_remove = Vec::new();
            let mut step_to_split: Option<usize> = None;
            let mut cumulative_pos = r2.start;

            // First pass: identify steps to remove/split (read-only, brief lock)
            {
                let mut graph = graph_mutex.lock().unwrap();
                for (idx, &step_handle) in r2.steps.iter().enumerate() {
                    let step_start = cumulative_pos;
                    let node_seq = graph.get_sequence(step_handle).unwrap();
                    let node_length = node_seq.len();
                    cumulative_pos += node_length;
                    let step_end = cumulative_pos;

                    if step_end <= overlap_start {
                        continue;
                    } else if step_start >= overlap_end {
                        break;
                    } else if step_start >= overlap_start && step_end <= overlap_end {
                        steps_to_remove.push(idx);
                    } else {
                        if step_to_split.is_some() {
                            panic!("Error: More than one step is partially overlapping, which is not allowed.");
                        }
                        step_to_split = Some(idx);
                    }
                }
            } // Lock released here

            // Initialize new vectors to store updated steps
            let mut new_steps = Vec::with_capacity(r2.steps.len() / 2);
            let mut range_new_start = None;
            let mut current_pos = None;

            // Reset cumulative position for second pass
            cumulative_pos = r2.start;

            // Second pass: Iterate over the original steps using incrementally computed positions
            for (idx, &step_handle) in r2.steps.iter().enumerate() {
                let step_start = cumulative_pos;

                // Get node sequence with lock
                let node_seq = {
                    let mut graph = graph_mutex.lock().unwrap();
                    graph.get_sequence(step_handle).unwrap()
                }; // Lock released here

                let node_length = node_seq.len();
                cumulative_pos += node_length;
                let step_end = cumulative_pos;

                if steps_to_remove.contains(&idx) {
                    // Skip steps to remove
                    continue;
                } else if step_to_split == Some(idx) {
                    // Split node for the single partially overlapping step
                    let overlap_within_step_start = std::cmp::max(step_start, overlap_start);
                    let overlap_within_step_end = std::cmp::min(step_end, overlap_end);

                    // Calculate offsets relative to the node sequence
                    let node_len = node_seq.len();
                    let overlap_start_offset =
                        (overlap_within_step_start - step_start).min(node_len);
                    let overlap_end_offset = (overlap_within_step_end - step_start).min(node_len);

                    debug!("      Splitting step {} [start={}, end={}, len={}] to remove overlap at [start={}, end={}], Overlap offsets: start={}, end={}",
                        idx, step_start, step_end, step_end - step_start, overlap_within_step_start, overlap_within_step_end, overlap_start_offset, overlap_end_offset);

                    if step_start < overlap_start {
                        debug!(
                            "      Adding left part of step [start={}, end={}]",
                            step_start, overlap_within_step_start
                        );
                        assert!(overlap_start_offset > 0);

                        // Keep left part
                        let new_seq = node_seq[0..overlap_start_offset].to_vec();

                        // Add node with lock
                        let node_id = {
                            let mut graph = graph_mutex.lock().unwrap();
                            graph.add_node(&new_seq).unwrap()
                        }; // Lock released here

                        let new_node = Handle::pack(NodeId::from(node_id), false);

                        new_steps.push(new_node);
                        if range_new_start.is_none() {
                            range_new_start = Some(step_start);
                            current_pos = Some(step_start);
                        }
                    } else if step_end > overlap_end {
                        debug!(
                            "      Adding right part of step [start={}, end={}]",
                            overlap_within_step_end, step_end
                        );
                        assert!(overlap_end_offset < node_len);

                        // Keep right part
                        let new_seq = node_seq[overlap_end_offset..].to_vec();

                        // Add node with lock
                        let node_id = {
                            let mut graph = graph_mutex.lock().unwrap();
                            graph.add_node(&new_seq).unwrap()
                        }; // Lock released here

                        let new_node = Handle::pack(NodeId::from(node_id), false);

                        new_steps.push(new_node);
                        if range_new_start.is_none() {
                            range_new_start = Some(overlap_end);
                            current_pos = Some(overlap_end);
                        }
                        current_pos = Some(current_pos.unwrap() + new_seq.len());
                    }
                } else {
                    // Keep steps that are not to be removed or split
                    new_steps.push(step_handle);
                    if range_new_start.is_none() {
                        range_new_start = Some(step_start);
                        current_pos = Some(step_start);
                    }
                    current_pos = Some(current_pos.unwrap() + node_length);
                }
            }

            // Update r2 with the new steps
            r2.steps = new_steps;

            // Update edges for the modified steps
            for idx in 0..r2.steps.len() {
                if idx > 0 {
                    let prev_step = r2.steps[idx - 1];
                    let curr_step = r2.steps[idx];

                    // Check and add edge with lock
                    let mut graph = graph_mutex.lock().unwrap();
                    if !graph.has_edge(
                        prev_step.id().into(),
                        prev_step.is_reverse(),
                        curr_step.id().into(),
                        curr_step.is_reverse(),
                    ) {
                        debug!(
                            "      Creating edge between steps: {} -> {}",
                            prev_step.id(),
                            curr_step.id()
                        );
                        graph.add_edge(
                            prev_step.id().into(),
                            prev_step.is_reverse(),
                            curr_step.id().into(),
                            curr_step.is_reverse(),
                        );
                    }
                    // Lock released at end of scope
                }
            }

            // Update r2.start and r2.end based on the new step positions
            if !r2.steps.is_empty() {
                r2.start = range_new_start.unwrap();
                r2.end = current_pos.unwrap();
            } else {
                // If no steps remain, set start and end to overlap_end to effectively remove this range
                r2.start = overlap_end;
                r2.end = overlap_end;
            }

            debug!(
                "      Updated overlaps: Range2 [start={}, end={}]",
                r2.start, r2.end
            );
        }
    }
}

fn link_contiguous_ranges(ranges: &[RangeInfo], graph_mutex: &Arc<Mutex<CompactGraph>>) {
    // Trim overlaps
    debug!("  Linking contiguous ranges");

    for i in 1..ranges.len() {
        let r1 = &ranges[i - 1];
        let r2 = &ranges[i];

        // Check if ranges are contiguous
        if r1.is_contiguous_with(r2) {
            // Get last handle from previous range and first handle from current range
            if let (Some(&last_handle), Some(&first_handle)) = (r1.steps.last(), r2.steps.first()) {
                // Lock only for checking and adding edge
                let mut graph = graph_mutex.lock().unwrap();
                if !graph.has_edge(
                    last_handle.id().into(),
                    last_handle.is_reverse(),
                    first_handle.id().into(),
                    first_handle.is_reverse(),
                ) {
                    debug!(
                        "    Creating edge between contiguous ranges at position {}: {} -> {}",
                        r1.end,
                        last_handle.id(),
                        first_handle.id()
                    );
                    graph.add_edge(
                        last_handle.id().into(),
                        last_handle.is_reverse(),
                        first_handle.id().into(),
                        first_handle.is_reverse(),
                    );
                }
                // Lock released at end of scope
            }
        }
    }
}

fn mark_nodes_for_removal(
    node_count: u64,
    path_key_ranges: &FxHashMap<String, Vec<RangeInfo>>,
) -> BitVec {
    // Create a bitvector with all nodes initially marked for removal
    let mut nodes_to_remove = bitvec![1; node_count as usize + 1]; // +1 to account for 0-indexing

    // Mark nodes used in path ranges as not to be removed (set bit to 0)
    for ranges in path_key_ranges.values() {
        for range in ranges {
            for handle in &range.steps {
                nodes_to_remove.set(u64::from(handle.id()) as usize, false);
            }
        }
    }

    nodes_to_remove
}

fn write_graph_to_gfa(
    combined_graph: &mut CompactGraph,
    path_key_ranges: &FxHashMap<String, Vec<RangeInfo>>,
    output_path: &str,
    compression_format: Format,
    fill_gaps: u8,
    fasta_index: &Option<FastaIndex>,
    debug: bool,
) -> std::io::Result<()> {
    info!("Marking unused nodes");
    let nodes_to_remove = mark_nodes_for_removal(combined_graph.node_count, path_key_ranges);
    debug!("Marked {} unused nodes", nodes_to_remove.count_ones() - 1);

    // Create the output file
    let output_file = File::create(output_path)?;

    // Create writer based on compression format
    let writer: Box<dyn Write> = match compression_format {
        Format::Gzip => {
            // Use parallel gzip compression
            let parz: ParCompress<Gzip> = ParCompressBuilder::new()
                .num_threads(rayon::current_num_threads())
                .map_err(|e| {
                    std::io::Error::other(
                        format!("Failed to set threads: {:?}", e),
                    )
                })?
                .compression_level(flate2::Compression::new(6))
                .from_writer(output_file);
            Box::new(parz)
        }
        Format::Bzip => {
            // Use parallel BGZF compression
            let parz: ParCompress<Bgzf> = ParCompressBuilder::new()
                .num_threads(rayon::current_num_threads())
                .map_err(|e| {
                    std::io::Error::other(
                        format!("Failed to set threads: {:?}", e),
                    )
                })?
                .compression_level(flate2::Compression::new(6))
                .from_writer(output_file);
            Box::new(parz)
        }
        Format::Zstd => {
            // Use multi-threaded zstd compression
            let mut encoder = ZstdEncoder::new(output_file, 6)?;
            encoder.multithread(rayon::current_num_threads() as u32)?;
            Box::new(encoder)
        }
        Format::No => {
            // No compression
            Box::new(output_file)
        }
        _ => {
            // Fallback to niffler for other formats
            niffler::get_writer(
                Box::new(output_file),
                compression_format,
                niffler::compression::Level::Six,
            )
            .map_err(std::io::Error::other)?
        }
    };

    let mut file = BufWriter::new(writer);

    // Write GFA version
    writeln!(file, "H\tVN:Z:1.0")?;

    // Write nodes by excluding marked ones and create the id_mapping
    info!("Writing used nodes by compacting their IDs");
    let max_id = combined_graph.node_count as usize;
    let mut id_mapping = vec![0; max_id + 1];
    let mut new_id = 1; // Start from 1

    for node_id in 1..=combined_graph.node_count {
        if !nodes_to_remove[node_id as usize] {
            id_mapping[node_id as usize] = new_id;

            let sequence =
                combined_graph.get_sequence(Handle::pack(NodeId::from(node_id), false))?;
            let sequence_str =
                String::from_utf8(sequence).expect("Node sequence contains invalid UTF-8");
            writeln!(file, "S\t{}\t{}", new_id, sequence_str)?;

            new_id += 1;
        }
    }

    info!("Writing edges connecting used nodes");
    for edge in &combined_graph.edges {
        let from_id = edge.from_id() as usize;
        let to_id = edge.to_id() as usize;

        if !nodes_to_remove[from_id] && !nodes_to_remove[to_id] {
            let from_mapped = id_mapping[from_id];
            let to_mapped = id_mapping[to_id];
            let from_orient = if edge.from_rev() { "-" } else { "+" };
            let to_orient = if edge.to_rev() { "-" } else { "+" };
            writeln!(
                file,
                "L\t{}\t{}\t{}\t{}\t0M",
                from_mapped, from_orient, to_mapped, to_orient
            )?;
        }
    }

    // Write paths by processing ranges directly
    info!("Writing paths by merging contiguous path ranges");
    let mut path_key_vec: Vec<_> = path_key_ranges.keys().collect();
    path_key_vec.par_sort_unstable(); // Sort path keys for consistent output (for path keys, the order of equal elements doesn't matter since they're unique)

    let mut start_gaps = 0;
    let mut middle_gaps = 0;
    let mut end_gaps = 0;

    // Check if a valid FASTA index is provided for end gap filling
    if fill_gaps == 2 && fasta_index.is_none() {
        warn!("Cannot fill end gaps without FASTA files; trailing gaps will be skipped");
    }

    for path_key in path_key_vec {
        let ranges = &path_key_ranges[path_key];
        //if ranges.is_empty() { continue }

        if debug {
            debug!("Processing Path key '{}'", path_key);

            let mut current_start = ranges[0].start;
            let mut current_end = ranges[0].end;

            for i in 1..ranges.len() {
                if ranges[i - 1].is_contiguous_with(&ranges[i]) {
                    // Extend current merged range
                    current_end = ranges[i].end;
                } else {
                    // Print current merged range
                    debug!(
                        "  Merged range: start={}, end={}",
                        current_start, current_end
                    );

                    if !ranges[i - 1].overlaps_with(&ranges[i]) {
                        // Calculate and print gap
                        let gap = ranges[i].start - current_end;
                        debug!("    Gap to next range: {} positions", gap);
                    } else {
                        // Calculate and print overlap (IT SHOULD NOT HAPPEN)
                        let overlap = current_end - ranges[i].start;
                        debug!("    Overlap with next range: {} positions", overlap);
                    }

                    // Start new merged range
                    current_start = ranges[i].start;
                    current_end = ranges[i].end;
                }
            }

            // Print final merged range
            debug!(
                "  Final merged range: start={}, end={}",
                current_start, current_end
            );
        }

        let mut path_elements: Vec<String> = Vec::new();
        let mut path_start = ranges[0].start;
        let mut path_end = ranges[0].start;

        // Handle initial gap if it exists and gap filling is enabled
        if fill_gaps == 2 && ranges[0].start > 0 {
            start_gaps += 1;
            path_elements.push(create_gap_node(
                &mut file,
                (0, ranges[0].start),
                path_key,
                fasta_index,
                None, // No previous element for initial gap
                ranges[0].steps.first(),
                &id_mapping,
                &mut new_id,
            )?);
            path_start = 0;
        }

        // Process subsequent contiguous ranges or add gap nodes
        let mut i = 0;
        while i < ranges.len() {
            // Merge a block of contiguous ranges
            add_range_steps_to_path(&ranges[i], &id_mapping, &mut path_elements);
            let mut last_range_end = ranges[i].end;
            i += 1;

            while i < ranges.len() && ranges[i - 1].is_contiguous_with(&ranges[i]) {
                add_range_steps_to_path(&ranges[i], &id_mapping, &mut path_elements);
                last_range_end = ranges[i].end;
                i += 1;
            }
            path_end = last_range_end;

            // We're now at the end of a contiguous block
            if i < ranges.len() {
                // there is a gap before the next block
                let next_start = ranges[i].start;

                if fill_gaps > 0 {
                    // Bridge the gap
                    middle_gaps += 1;

                    path_elements.push(create_gap_node(
                        &mut file,
                        (last_range_end, next_start),
                        path_key,
                        fasta_index,
                        path_elements.last(),
                        ranges[i].steps.first(),
                        &id_mapping,
                        &mut new_id,
                    )?);
                } else {
                    // Finish current partial path and start a new one
                    if !path_elements.is_empty() {
                        let path_name = format!("{}:{}-{}", path_key, path_start, path_end);
                        writeln!(file, "P\t{}\t{}\t*", path_name, path_elements.join(","))?;
                        path_elements.clear();
                    }
                    path_start = next_start;
                }
            }
        }

        // Handle final gap if it exists and gap filling is enabled
        if fill_gaps == 2 && fasta_index.is_some() {
            // Try to get the FASTA file for this path_key and its total length
            if let Some(fasta_path) = fasta_index.as_ref().unwrap().get_fasta_path(path_key) {
                match faidx::Reader::from_path(fasta_path) {
                    Ok(reader) => {
                        let total_len = reader.fetch_seq_len(path_key) as usize;
                        if path_end < total_len {
                            end_gaps += 1;

                            path_elements.push(create_gap_node(
                                &mut file,
                                (path_end, total_len),
                                path_key,
                                fasta_index,
                                path_elements.last(),
                                None, // No next node for final gap
                                &id_mapping,
                                &mut new_id,
                            )?);
                            path_end = total_len;
                        } else if path_end > total_len {
                            warn!(
                                "Path '{}' extends beyond sequence length ({} > {})",
                                path_key, path_end, total_len
                            );
                        }
                    }
                    Err(e) => {
                        warn!(
                            "Failed to open FASTA file '{}' for path '{}': {}",
                            fasta_path, path_key, e
                        );
                    }
                }
            }
        }

        // Write the remaining path
        if !path_elements.is_empty() {
            let is_full_path = path_start == 0
                && fasta_index
                    .as_ref()
                    .and_then(|idx| idx.get_fasta_path(path_key))
                    .and_then(|fasta_path| faidx::Reader::from_path(fasta_path).ok())
                    .map(|reader| reader.fetch_seq_len(path_key) as usize == path_end)
                    .unwrap_or(true); // Assume path_end matches full length if no FASTA index is available

            let path_name = if is_full_path {
                path_key.to_string()
            } else {
                format!("{}:{}-{}", path_key, path_start, path_end)
            };
            writeln!(file, "P\t{}\t{}\t*", path_name, path_elements.join(","))?;
        }
    }

    if fill_gaps == 2 {
        info!(
            "Filled {} gaps: {} start gaps, {} middle gaps, {} end gaps",
            start_gaps + middle_gaps + end_gaps,
            start_gaps,
            middle_gaps,
            end_gaps
        );
    } else if fill_gaps == 1 {
        info!("Filled {} middle gaps", middle_gaps);
    }

    file.flush()?;
    Ok(())
}

fn create_gap_node<W: Write>(
    file: &mut BufWriter<W>,
    gap_range: (usize, usize),
    path_key: &str,
    fasta_index: &Option<FastaIndex>,
    last_element: Option<&String>,
    next_handle: Option<&Handle>,
    id_mapping: &[usize],
    new_id: &mut usize,
) -> io::Result<String> {
    let (gap_start, gap_end) = gap_range;
    let gap_size = gap_end - gap_start;

    // Get gap sequence either from FASTA or create string of N's
    let gap_sequence = if let Some(index) = fasta_index {
        if let Some(fasta_path) = index.get_fasta_path(path_key) {
            match faidx::Reader::from_path(fasta_path) {
                Ok(reader) => {
                    // Apply the fix for the rust_htslib memory leak bug
                    // https://github.com/rust-bio/rust-htslib/issues/401#issuecomment-1704290171
                    match reader.fetch_seq(path_key, gap_start, gap_end - 1) {
                        Ok(seq) => {
                            let seq_vec = seq.to_vec();
                            // Free the memory allocated by htslib to prevent memory leak
                            unsafe { libc::free(seq.as_ptr() as *mut std::ffi::c_void) };
                            String::from_utf8(seq_vec).unwrap_or_else(|e| {
                                error!("Failed to convert sequence to UTF-8: {}", e);
                                "N".repeat(gap_size)
                            })
                        }
                        Err(e) => {
                            error!("Failed to fetch sequence from '{}': {}", fasta_path, e);
                            "N".repeat(gap_size)
                        }
                    }
                }
                Err(e) => {
                    warn!("Failed to open FASTA file '{}': {}", fasta_path, e);
                    "N".repeat(gap_size)
                }
            }
        } else {
            "N".repeat(gap_size)
        }
    } else {
        "N".repeat(gap_size)
    };

    // Write gap node
    writeln!(file, "S\t{}\t{}", new_id, gap_sequence)?;

    // Add edge from previous node if it exists
    if let Some(last_element) = last_element {
        let last_id = last_element[..last_element.len() - 1]
            .parse::<usize>()
            .unwrap();
        let last_orient = &last_element[last_element.len() - 1..];
        writeln!(file, "L\t{}\t{}\t{}\t+\t0M", last_id, last_orient, new_id)?;
    }

    // Add edge to next node if it exists
    if let Some(handle) = next_handle {
        let next_id = id_mapping[u64::from(handle.id()) as usize];
        let next_orient = if handle.is_reverse() { "-" } else { "+" };
        writeln!(file, "L\t{}\t+\t{}\t{}\t0M", *new_id, next_id, next_orient)?;
    }

    let path_element = format!("{}+", new_id);
    *new_id += 1;

    Ok(path_element)
}

fn add_range_steps_to_path(
    range: &RangeInfo,
    id_mapping: &[usize],
    path_elements: &mut Vec<String>,
) {
    for handle in &range.steps {
        let node_id = id_mapping[u64::from(handle.id()) as usize];
        let orient = if handle.is_reverse() { "-" } else { "+" };
        path_elements.push(format!("{}{}", node_id, orient));
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Helper function to create a simple RangeInfo for testing
    fn create_range_info(start: usize, end: usize) -> RangeInfo {
        RangeInfo {
            start,
            end,
            //gfa_id,
            steps: vec![], // Empty steps for testing
        }
    }

    #[test]
    fn test_range_containment_removal() {
        // Test cases
        let test_cases = vec![
            // Test case 1: Basic containment
            (
                vec![(10, 50), (20, 30)], // Input ranges (start, end)
                vec![(10, 50)],           // Expected result
                "Basic containment",
            ),
            // Test case 2: No containment
            (
                vec![(10, 20), (30, 40)],
                vec![(10, 20), (30, 40)],
                "No containment",
            ),
            // Test case 3: Multiple contained ranges
            (
                vec![(10, 100), (20, 30), (40, 50), (60, 70)],
                vec![(10, 100)],
                "Multiple contained ranges",
            ),
            // Test case 4: Identical ranges
            (vec![(10, 20), (10, 20)], vec![(10, 20)], "Identical ranges"),
            // Test case 5: Nested containment
            (
                vec![(10, 100), (20, 80), (30, 40)],
                vec![(10, 100)],
                "Nested containment",
            ),
            // Test case 6: Partial overlap (should keep both)
            (
                vec![(10, 30), (20, 40)],
                vec![(10, 30), (20, 40)],
                "Partial overlap",
            ),
            // Test case 7: Edge cases - touching ranges
            (
                vec![(10, 20), (20, 30)],
                vec![(10, 20), (20, 30)],
                "Touching ranges",
            ),
            // Test case 8: Overlapping ranges from same GFA
            (
                vec![(0, 11742), (9714, 13000), (11000, 19000)],
                vec![(0, 11742), (11000, 19000)],
                "Overlapping ranges from same GFA",
            ),
            // Test case 9: Overlapping ranges with different GFA IDs
            (
                vec![(0, 11742), (9714, 13000), (11000, 19000)],
                vec![(0, 11742), (11000, 19000)],
                "Overlapping ranges",
            ),
            // Test case 10: Overlapping ranges with different GFA IDs 2
            (
                vec![(0, 10), (8, 20), (15, 30)],
                vec![(0, 10), (8, 20), (15, 30)],
                "Overlapping ranges",
            ),
            // Test case 11: Overlapping ranges with different GFA IDs 3
            (
                vec![(8000, 11000), (9694, 12313), (10908, 13908)],
                vec![(8000, 11000), (10908, 13908)],
                "Overlapping ranges",
            ),
        ];

        // Run each test case
        for (case_index, (input_ranges, expected_ranges, case_name)) in
            test_cases.iter().enumerate()
        {
            println!("Running test case {}: {}", case_index + 1, case_name);

            // Create input ranges
            let mut ranges: Vec<RangeInfo> = input_ranges
                .iter()
                .map(|(start, end)| create_range_info(*start, *end))
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
                .map(|(start, end)| create_range_info(*start, *end))
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
                    (result.start, result.end),
                    (expected.start, expected.end),
                    "Test case '{}': Mismatch at position {}",
                    case_name,
                    i
                );
            }

            println!("Test case {} passed: {}", case_index + 1, case_name);
        }
    }
}
