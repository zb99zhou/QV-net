#![allow(clippy::many_single_char_names, clippy::too_many_arguments)]
#![allow(dead_code)]
#![allow(non_snake_case)]

mod sigma_dl;
mod sigma_dleq;
mod transcript;
mod voting;
mod zk_eval;
mod sum_square;

extern crate serde_derive;
extern crate serde;

extern crate curv;
extern crate generic_array;
extern crate itertools;
extern crate sha2;

#[derive(Copy, PartialEq, Eq, Clone, Debug)]
pub enum Errors {
    SigmaDlProofError,
    SigmaDleqProofError,
    ZkEvalProofError,
    ZkSumSquareArgError,
    VotingError
}
