# QV-net: Decentralized Self-Tallying Quadratic Voting with Maximal Ballot Secrecy

This repository provides an implementation of QV-net, a decentralized self-tallying quadratic voting with maximal ballot secrecy.

Decentralized e-voting enables secure and transparent elections without relying on trusted authorities, with blockchain emerging as a popular platform. Quadratic voting (QV) is a social choice mechanism where the cost of votes increases quadratically with the number of votes. To cast $n$ votes, a voter needs to spend $n^2$ tokens, which effectively prevents wealthy token holders from dominating the voting outcome while still allowing them to express strong preferences. In Decentralized Autonomous Organizations (DAOs), QV has been adopted by organizations such as MetFi DAO, Karmaverse, Pistachio DAO, Gitcoin, and MoonDAO.

QV-net is the first decentralized quadratic voting scheme, in which voters do not need to trust any external party other than themselves for ballot secrecy. QV-net is self-tallying with maximal ballot secrecy. Self-tallying allows anyone to compute election results once all ballots are cast. Maximal ballot secrecy ensures that what each voter learns from QV-net is nothing more than the tally and their own ballot.

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

The QV-net project is structured as follows:

### Root Directory

- **`Cargo.toml`**: Defines the project dependencies, features, and metadata.
- **`README.md`**: Provides an overview of the project, installation instructions, and usage guidelines.
- **`VarRange/`**: Contains the source code for VarRange, an aggregated range proof that supports different bases and ranges.
- **`zkSumCheck/`**: Contains the source code for zkSumCheck, a zero-knowledge sum-check protocol.
- **`src/`**: Contains the source code for the QV-net implementation.

### `src/` Directory

- **`lib.rs`**: The main entry point for the library. It organizes and exposes the core modules of the project.
- **`transcript.rs`**: Defines the transcript abstraction for managing Fiat-Shamir transformations.
- **`sigma_dl.rs`**: Implements an AND composition for proving knowledge of discrete logarithm.
- **`sigma_dleq.rs`**: Implements an AND composition for proving knowledge of discrete logarithm equality.
- **`zk_eval.rs`**: Implements a zero-knowledge Eval protocol of the discrete-logarithm-based multilinear PCS.
- **`sum_square.rs`**: Implements a succinct ZKAoK for sum-of-square relation.
- **`voting.rs`**: Implements the core QV-net protocol.

