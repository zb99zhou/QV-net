# VarRange: an aggregated range proof supporting different bases and ranges.

This repository provides an implementation of VarRange, an aggregated range proof that supports different bases and ranges.

An aggregated range proof allows a single prover to prove that multiple committed values lie within specified ranges simultaneously. Given $n$ Pedersen commitments with different commitment keys, VarRange proves that each committed value lies within its respective range.

## Installation

To set up the environment for this project, visit the official Rust website at [https://www.rust-lang.org/tools/install](https://www.rust-lang.org/tools/install) to install Rust.

For Linux and macOS, you can install Rust using the command:

`curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh`

Verify the installation: `rustc --version`

## Build and Test
This project requires the GMP library for high-performance arithmetic operations.

Install GMP
* For Ubuntu/Debian systems:
`sudo apt-get install libgmp-dev`.

* For macOS:
`brew install gmp`

Build the project: `cargo build --release`

Run tests: `cargo test --release`

## File Organization

The VarRange project is structured as follows:

### Root Directory
- **`Cargo.toml`**: Defines the project dependencies, features, and metadata.
- **`README.md`**: Provides an overview of the project, installation instructions, and usage guidelines.
- **`src/`**: Contains the source code for the VarRange implementation.

### `src/` Directory
- **`lib.rs`**: The main entry point for the library. It organizes and exposes the core modules of the project.
- **`proofs/`**: Contains the implementation of the VarRange protocol and related components:
  - **`varrange.rs`**: Implements the core VarRange protocol.
  - **`ipa.rs`**: Implements the inner product argument used in the protocol.
  - **`sigma_pedersen.rs`**: Implements an AND composition for proving knowledge of Pedersen commitment openings.
  - **`transcript.rs`**: Defines the transcript abstraction for managing Fiat-Shamir transformations.
  - **`vec_poly.rs`**: Provides utility functions for working with vector polynomials.
  - **`mod.rs`**: Organizes and re-exports the modules in the `proofs/` directory.