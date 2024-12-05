# GFALace

GFALace is a Rust tool that combines multiple GFA (Graphical Fragment Assembly) files into a single unified graph. It's designed for working with pangenome graphs that have been split into multiple files, where each file represents a specific genomic region.

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

You can freely mix gzip-compressed (`.gfa.gz`) and uncompressed (`.gfa`) files in the input.

The input GFA files can be provided in any order. This is because GFALace uses the coordinate information embedded in the path names (CHROM:START-END) to determine the correct ordering and relationships between sequences.

Options:
- `-g, --gfa-list`: List of input GFA files (space-separated)
- `-o, --output`: Output GFA file path
- `-d, --debug`: Enable debug output
- `-h, --help`: Show help information
- `-V, --version`: Show version information

## Path Name Format

GFALace expects path names in the format:
```
SAMPLE#HAP#CHROM:START-END
```
Example: `HG002#1#chr20:1000-2000`

The tool uses these coordinates to:
1. Identify which sequences belong together
2. Order the sequences correctly
3. Detect and handle overlaps or gaps

## Features

- Combines multiple GFA files (both .gfa and .gfa.gz) while preserving path information
- Translates node IDs to avoid conflicts
- Creates edges between contiguous path segments
- Handles both contiguous and non-contiguous ranges
- Preserves original sequence and path relationships
- Outputs a standard-compliant GFA 1.0 file

## License

MIT License - See LICENSE file for details.
