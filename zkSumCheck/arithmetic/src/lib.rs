// Copyright (c) 2023 Espresso Systems (espressosys.com)
// This file is part of the HyperPlonk library.

// You should have received a copy of the MIT License
// along with the HyperPlonk library. If not, see <https://mit-license.org/>.

mod errors;
mod multilinear_polynomial;
mod virtual_polynomial;

pub use errors::ArithErrors;
pub use multilinear_polynomial::{
    fix_variables,
    DenseMultilinearExtension,
};
pub use virtual_polynomial::{
    build_eq_x_r_vec, VPAuxInfo, VirtualPolynomial,
};
