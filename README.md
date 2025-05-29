# GFALace

GFALace is a Rust tool that combines multiple GFA (Graphical Fragment Assembly) files into a single unified graph. It's designed for working with pangenome graphs that have been split into multiple files, where each file represents a specific pangenomic region.

## Installation

Requires Rust 2021 edition or later. Install using:

```bash
cargo install --git https://github.com/pangenome/gfalace
```

Or build from source:

```bash
git clone https://github.com/pangenome/gfalace
cd gfalace
cargo build --release
```

## Usage

Basic usage:

```bash
gfalace -g *.gfa -o combined.gfa
```

or

```bash
gfalace -g file1.gfa file2.gfa.gz file3.gfa -o combined.gfa
```

of from a file list:

```bash
gfalace -l gfa_list.txt -o combined.gfa
```

You can mix compressed (`.gfa.gz`, `.gfa.bgz`, `.gfa.zst`) and uncompressed (`.gfa`) files in the input.

The input GFA files can be provided in any order. This is because GFALace uses the coordinate information in the path names (`CHROM:START-END`) to determine the correct ordering and relationships between sequences.

## Options

- `-g, --gfa-files`: List of input GFA files (space-separated)
- `-l, --gfa-list`: Text file containing GFA paths (one per line)
- `-o, --output`: Output GFA file path
- `--compress`: Output compression format: `none`, `gzip`, `bgzip`, `zstd`, or `auto` (default: auto-detect from extension)
- `--fill_gaps`: Gap filling mode (0 = `none` [default], 1 = `middle` gaps only, 2 = `all` gaps)
- `--fasta`: FASTA file containing sequences for gap filling
- `--temp-dir`: Directory for temporary files
- `-t, --num-threads`: Number of threads (default: 4)
- `-h, --help`: Show help information
- `-V, --version`: Show version information

## Path Name Format

GFALace expects path names in the format:

```
NAME:START-END
```

Example: `HG002#1#chr20:1000-2000`

The tool uses these coordinates to:
1. Identify which sequences belong together
2. Order the sequences correctly
3. Detect and handle overlaps or gaps

Note: `NAME` can contain ':' characters. When parsing coordinates, GFALace uses the last occurrence of ':' to separate the name from the coordinate range.

## Features

- Combines multiple GFA files while preserving path information
- Parallel processing for improved performance
- Translates node IDs to avoid conflicts
- Creates edges between contiguous path segments
- Handles both contiguous and non-contiguous ranges
- Preserves original sequence and path relationships
- Outputs a standard-compliant GFA 1.0 file

## Post-processing recommendations

After combining the GFA files, the resulting graph will already have compacted node IDs ranging from `1` to the total number of nodes. However, it is strongly recommended to perform post-processing steps using **[ODGI](https://github.com/pangenome/odgi)** to unchop and sort the graph.

```bash
odgi unchop -i combined.gfa -o - -t 16 | \
    odgi sort -i - -o - -p gYs -t 16 | \
    odgi view -i - -g > combined.final.gfa
```

If overlaps were present, and then trimmed during the merging process, it's advisable to run **[GFAffix](https://github.com/marschall-lab/GFAffix)** before the ODGI pipeline to remove redundant nodes introduced by the overlap trimming.

```bash
gfaffix combined.gfa -o combined.fix.gfa &> /dev/null

odgi unchop -i combined.fix.gfa -o - -t 16 | \
    odgi sort -i - -o - -p gYs -t 16 | \
    odgi view -i - -g > combined.final.gfa
```

### Advanced usage

Filling middle gaps with `N`s:

```bash
gfalace -g *.gfa -o combined.gfa --fill_gaps 1
```

Filling all gaps with pangenome sequences:

```bash
gfalace -g *.gfa -o combined.gfa --fill_gaps 2 --fasta pangenome.fasta
```

GFALace provides options to fill gaps between graphs based on the specified gap filling mode:

- **Mode 0 (Default)**: No gap filling. Gaps between segments are left unfilled.
- **Mode 1**: Fills middle gaps between contiguous ranges with `N` characters or with sequences provided with `--fasta`. Useful for connecting segments that are not directly connected but belong to the same path.
- **Mode 2**: Fills all gaps, including start and end gaps, with `N` characters or with sequences provided with `--fasta`. To fill end gaps, the FASTA is required.

## License

MIT License - See LICENSE file for details.
