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
  `brew install gmp`. 

  After installation, set the environment variables according to your Mac processor type so that the compiler and linker can find the GMP library files:    
  - For Apple-Silicon macOS:  
    `export LIBRARY_PATH="/opt/homebrew/lib"`

    `export CPATH="/opt/homebrew/include"`

  - For Intel Macs:  
    `export LIBRARY_PATH="/usr/local/lib"`

    `export CPATH="/usr/local/include"`

Build the project: `cargo build --release`

Run tests: `cargo test --release`

## Benchmark Scripts
We provide a set of shell scripts in the  `script` directory to reproduce data reported in the paper. 

Please ensure the `bc` command-line calculator is installed, as some scripts require it for floating-point calculations. 
* On most Linux systems, install it with: 
`sudo apt install bc`. 
* On macOS, install it with: `brew install bc`.

### Usage
1. Navigate to the `script` directory:
   `cd script`

2. Run the corresponding script to generate results for each table or figure:
   - `./real_world_daos.sh` # Table 5.
   - `./shanks.sh` # Table 6.
   - `./test_vary_nc.sh` # Figure 3a - 3d.
   - `./tally_vary_t.sh` # Figure 3f.
   - `./test_vary_nv.sh` # Figure 3g, 3h.
   - `./VarRange_vary_t.sh` # Figure 4.
        
### Notes
All storage overhead is manually calculated based on the storage complexity presented in Table 3, using the field and group element sizes of the Secp256k1 curve (32 bytes for $\mathbb{F}_p$ and 33 bytes for $\mathbb{G}$).

For test cases with a large number of voters, the process may take a considerable amount of time (up to tens of minutes), as all voters are simulated on a single machine. In practice, each voter would execute their voting procedure independently.

For the Shanks' algorithm, we limit the bound to values below $2^{36}$ to avoid excessively long runtimes.


## File Organization

The QV-net project is structured as follows:

### Root Directory

- **`Cargo.toml`**: Defines the project dependencies, features, and metadata.
- **`README.md`**: Provides an overview of the project, installation instructions, and usage guidelines.
- **`VarRange/`**: Contains the source code for VarRange, an aggregated range proof that supports different bases and ranges.
- **`zkSumCheck/`**: Contains the source code for zkSumCheck, a zero-knowledge sum-check protocol.
- **`src/`**: Contains the source code for the QV-net implementation.
- **`script/`**: Contains shell scripts for reproducing the benchmark data reported in the paper. 

### `src/` Directory

- **`lib.rs`**: The main entry point for the library. It organizes and exposes the core modules of the project.
- **`transcript.rs`**: Defines the transcript abstraction for managing Fiat-Shamir transformations.
- **`sigma_dl.rs`**: Implements an AND composition for proving knowledge of discrete logarithm.
- **`sigma_dleq.rs`**: Implements an AND composition for proving knowledge of discrete logarithm equality.
- **`zk_eval.rs`**: Implements a zero-knowledge Eval protocol of the discrete-logarithm-based multilinear PCS.
- **`sum_square.rs`**: Implements a succinct ZKAoK for sum-of-square relation.
- **`voting.rs`**: Implements the core QV-net protocol.

