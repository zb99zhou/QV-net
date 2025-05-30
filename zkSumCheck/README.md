# zkSumCheck: a zero-knowledge sum-check protocol

This repository provides an implementation of zkSumCheck, a zero-knowledge sum-check protocol which is adapted from the [Hyperplonk](https://github.com/EspressoSystems/hyperplonk) library.

Suppose we are given a $\ell$-variate polynomial $g$ defined over a finite field $\mathbb{F}$. A sum-check protocol allows a prover to provide the verifier with the follwing sum:
$$
H:=\sum_{\boldsymbol{b}\in\{0,1\}^\ell}g(\boldsymbol{b}).
$$
And zkSumCheck makes this argumentation process zero-knowledge with the help of random polynomials for mask.

## Installation

To set up the environment for this project, visit the [official Rust website](https://www.rust-lang.org/tools/install) to install Rust.

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

The zkSumCheck project is structured as follows:

### Root Directory

- **`Cargo.toml`**: Defines the project dependencies, features, and metadata.
- **`README.md`**: Provides an overview of the project, installation instructions, and usage guidelines.
- **`arithmetic/`**: Contains the source code for the efficient polynomial arithmetic primitives implementation.
  - **`Cargo.toml`**: Defines the directory dependencies, features, and metadata.
  - **`src/`**: Contains the source code.
- **`subroutines/`**: Contains the source code for the zero-knowledge sum-check protocol implementation.
  - **`Cargo.toml`**: Defines the directory dependencies, features, and metadata.
  - **`src/`**: Contains the source code.
- **`transcript/`**: Contains the source code for the transcript abstraction for managing Fiat-Shamir transformations.
  - **`Cargo.toml`**: Defines the directory dependencies, features, and metadata.
  - **`src/`**: Contains the source code.

### `arithmetic/src/` Directory

- **`lib.rs`**: The main entry point for the library. It organizes and exposes the core modules of the project.
- **`errors.rs`**: Centralized error handling, defines custom error types and conversion logic.
- **`multilinear_polynomial.rs`**: Multilinear polynomial implementation, includes data structures, evaluation, and interpolation algorithms.
- **`virtual_polynomial.rs`**: Virtual polynomial operations, enables efficient polynomial composition, symbolic computation, and specialized operations.

### `subroutines/src/` Directory

- **`lib.rs`**: The main entry point for the library. It organizes and exposes the core modules of the project.

- **`poly_iop/`**: Contains the source code of  polynomial interactive oracle proofs implementation.

  - **`mod.rs`**: Primary module file, declares sub-modules and orchestrates protocol components.

  - **`errors.rs`**: Protocol-specific error handling defines error types (e.g., proof verification failure, invalid parameters).
  - **`structs.rs`**: Core data structures, defines protocol objects (e.g., proof state
  - **`prelude.rs`**: Module prelude, aggregates commonly used types/traits for simplified access.

  - **`zk_sum_check/`**: Contains the source code of the zero-knowledge sum-check protocol implementation.
    - **`mod.rs`**: Primary module file, declares sub-modules and orchestrates protocol components.
    - **`prover.rs`**: Prover algorithm implementation, generates ZK proofs for sumcheck protocol.
    - **`verifier.rs`**: Verifier logic implementation, validates proof correctness and handles challenge-response.

### `transcript/src/` Directory

- **`lib.rs`**: The main entry point for the library. It organizes and exposes the core modules of the project.
- **`errors.rs`**: Centralized error handling, defines custom error types and conversion logic.