# gfalace

`gfalace` is a Rust-based tool designed to process a collection of GFA (Graphical Fragment Assembly) files. It aims to interpret the path names within these files as sub-ranges of a set of genome sequences. Utilizing these positional relationships, GFALace laces together the sub-graphs extracted from the GFA files into a single, comprehensive graph structure. This process facilitates the analysis and visualization of genomic data by providing a unified view of the sequences and their connections.

## Features

- **Parsing GFA Files**: GFALace can parse multiple GFA files, extracting the graph structures and path information contained within.
- **Graph Construction**: Leveraging the `handlegraph` library, it constructs `HashGraph` structures from the parsed GFA data, enabling efficient graph operations.
- **Path Interpretation**: Path names in the GFA files are interpreted as indicating sub-ranges of genome sequences, which GFALace uses to correctly position and lace together sub-graphs.
- **Unified Graph Output**: The tool outputs a single, large graph that represents the laced-together sub-graphs, providing a comprehensive view of the genomic data.

## Usage

To use GFALace, you need to have Rust installed on your system. Once Rust is set up, you can clone the repository and build the project using Cargo, Rust's package manager and build system.

### Basic Command Line Usage

```sh
cargo run --release -- -g <GFA_FILE_PATHS>
```

Here, `<GFA_FILE_PATHS>` should be replaced with the paths to the GFA files you wish to process, separated by spaces.

## Warning

GFALace is currently a work in progress (WIP) and has not been extensively tested. It is being made available publicly as a checkpoint in its development process. Users should be aware that the tool may contain bugs and its output should be verified independently. Contributions, bug reports, and suggestions for improvements are welcome.

## Dependencies

GFALace relies on several Rust crates, including:

- `bstr` for byte string operations.
- `clap` for command-line argument parsing.
- `handlegraph` for graph operations, specifically using the `HashGraph` structure for efficient graph manipulation.
- `gfa` for parsing GFA files.

These dependencies are specified in the `Cargo.toml` file and will be automatically managed by Cargo when building the project.

## Contributing

Contributions to GFALace are welcome. If you have suggestions for improvements, bug fixes, or new features, please feel free to open an issue or submit a pull request.

## License

GFALace is open-source software licensed under the MIT License. See the LICENSE file for more details.
