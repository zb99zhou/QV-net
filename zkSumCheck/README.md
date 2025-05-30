# zkSumCheck: a zero-knowledge sum-check protocol

This repository provides an implementation of zkSumCheck, a zero-knowledge sum-check protocol. Except for `zkSumCheck/subroutines/src/poly_iop/zk_sum_check`, other files are directly ported from the [Hyperplonk library](https://github.com/EspressoSystems/hyperplonk), where we delete some unnecessary functions.

## Installation

To set up the environment for this project, visit the [official Rust website](https://www.rust-lang.org/tools/install) to install Rust.

For Linux and macOS, you can install Rust using the command:

`curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh`

Verify the installation: `rustc --version`

## Build and Test

Build the project: `cargo build --release`

Run tests: `cargo test --release`

## File Organization

The zkSumCheck project is structured as follows:

### Root Directory

- **`Cargo.toml`**: Defines the project dependencies, features, and metadata.
- **`README.md`**: Provides an overview of the project, installation instructions, and usage guidelines.
- **`arithmetic/`**: Contains the source code for the efficient polynomial arithmetic primitives implementation.
- **`subroutines/`**: Contains the source code for the zero-knowledge sum-check protocol implementation.
- **`transcript/`**: Contains the source code for the transcript abstraction for managing Fiat-Shamir transformations.
- **`util/`**: Contains utility functions and helper modules used across the project.

### `arithmetic/src/` Directory

- **`lib.rs`**: The main entry point for the library. It organizes and exposes the core modules of the project.
- **`errors.rs`**: Centralized error handling, defines custom error types and conversion logic.
- **`multilinear_polynomial.rs`**: Multilinear polynomial implementation, includes data structures, evaluation, and interpolation algorithms.
- **`virtual_polynomial.rs`**: Virtual polynomial operations, enables efficient polynomial composition, symbolic computation, and specialized operations.

### `subroutines/src/` Directory

- **`lib.rs`**: The main entry point for the library. It organizes and exposes the core modules of the project.

- **`poly_iop/`**: Contains the source code of  polynomial interactive oracle proofs implementation.

  - **`mod.rs`**: Primary module file, declares sub-modules and orchestrates protocol components.

  - **`errors.rs`**: Protocol-specific error handling defines error types.
  - **`structs.rs`**: Core data structures, defines protocol objects.
  - **`prelude.rs`**: Module prelude, aggregates commonly used types/traits for simplified access.

  - **`zk_sum_check/`**: Contains the source code of the zero-knowledge sum-check protocol implementation.
    - **`mod.rs`**: Primary module file, declares sub-modules and orchestrates protocol components.
    - **`prover.rs`**: Prover algorithm implementation, generates ZK proofs for sumcheck protocol.
    - **`verifier.rs`**: Verifier logic implementation, validates proof correctness.
