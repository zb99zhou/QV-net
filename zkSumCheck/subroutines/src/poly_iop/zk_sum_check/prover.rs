use super::ZkSumCheckProver;
use crate::poly_iop::{
    errors::PolyIOPErrors,
    structs::{IOPProverMessage, IOPProverState, RandomMaskPolynomial},
};
use arithmetic::{fix_variables, VirtualPolynomial};
use ark_ff::{batch_inversion, PrimeField};
use ark_poly::DenseMultilinearExtension;
use ark_std::{cfg_into_iter, end_timer, rand::RngCore, start_timer, vec::Vec};
use itertools::max;
use rayon::prelude::{IntoParallelIterator, IntoParallelRefIterator};
use std::sync::Arc;

#[cfg(feature = "parallel")]
use rayon::iter::{IntoParallelRefMutIterator, ParallelIterator};

impl<F: PrimeField> RandomMaskPolynomial<F> {
    pub fn rand<R: RngCore>(
        nv: usize,
        degree: usize,
        rng: &mut R
    ) -> (RandomMaskPolynomial<F>, F) {
        let mut evals = vec![vec![F::zero(); degree+1]; nv];
        for i in 0..nv {
            for j in 1..degree+1 {
                evals[i][j] = F::rand(rng);
            }
        }
        let const_term = F::rand(rng);
        let mask_poly = RandomMaskPolynomial::<F> {
            const_term,
            evaluations: evals
        };
        let mut sum = mask_poly.evaluations.iter().map(|row| row[1]).sum();
        sum *= F::from((1 << (nv-1)) as u64);
        sum += F::from((1 << nv) as u64) * const_term;
        
        (mask_poly, sum)
    }

    pub fn eval(
        &self,
        point: &[F]
    ) -> Result<F, PolyIOPErrors> {
        assert_eq!(point.len(), self.evaluations.len());

        let mut res = F::zero();
        for i in 0..self.evaluations.len() {
            res += interpolate_uni_poly(&self.evaluations[i], point[i])?;
        }

        Ok(res + self.const_term)
    }
}

pub struct ZkSumCheckProverState<F: PrimeField> {
    // sum check prover state
    pub(crate) sum_check_prover_state: IOPProverState<F>,
    // mask polynomial
    pub(crate) mask_poly: RandomMaskPolynomial<F>,
    // `sum_aux[i]` is \sum_{j=i+1, ..., num_variables}g_j(1)
    pub(crate) sum_aux: Vec<F>,
    // current_sum of evaluaions of g_i in random challenge
    pub(crate) current_sum: F
}

impl<F: PrimeField> ZkSumCheckProver<F> for ZkSumCheckProverState<F> {
    type VirtualPolynomial = VirtualPolynomial<F>;
    type ProverMessage = IOPProverMessage<F>;
    type RandomMaskPolynomial = RandomMaskPolynomial<F>;

    /// Initialize the prover state to argue for the sum of the input polynomial
    /// over {0,1}^`num_vars`.
    fn prover_init(polynomial: &Self::VirtualPolynomial, mask_poly: &Self::RandomMaskPolynomial) -> Result<Self, PolyIOPErrors> {
        let start = start_timer!(|| "sum check prover init");
        if polynomial.aux_info.num_variables == 0 {
            return Err(PolyIOPErrors::InvalidParameters(
                "Attempt to prove a constant.".to_string(),
            ));
        } else if mask_poly.evaluations.len() != polynomial.aux_info.num_variables {
            return Err(PolyIOPErrors::InvalidParameters(
                "Number of variables is not match".to_string(),
            ));
        }
        let max_degree = max([polynomial.aux_info.max_degree, mask_poly.evaluations[0].len()-1]).unwrap();

        let sum_check_prover_state = IOPProverState {
            challenges: Vec::with_capacity(polynomial.aux_info.num_variables),
            round: 0,
            poly: polynomial.clone(),
            extrapolation_aux: (1..max_degree)
                .map(|degree| {
                    let points = (0..1 + degree as u64).map(F::from).collect::<Vec<_>>();
                    let weights = barycentric_weights(&points);
                    (points, weights)
                })
                .collect(),
        };

        let mut sum_aux = vec![F::zero(); polynomial.aux_info.num_variables];
        sum_aux[0] = mask_poly.evaluations.iter().map(|row| row[1]).sum();
        for i in 1..sum_aux.len() {
            sum_aux[i] = sum_aux[i-1] - mask_poly.evaluations[i-1][1];
        }
        end_timer!(start);

        Ok(Self {
            sum_check_prover_state,
            mask_poly: mask_poly.clone(),
            sum_aux,
            current_sum: mask_poly.const_term
        })
    }

    /// Receive message from verifier, generate prover message, and proceed to
    /// next round.
    fn prove_round_and_update_state(
        &mut self,
        rho: &F,
        challenge: &Option<F>,
    ) -> Result<Self::ProverMessage, PolyIOPErrors> {
        if self.sum_check_prover_state.round >= self.sum_check_prover_state.poly.aux_info.num_variables {
            return Err(PolyIOPErrors::InvalidProver(
                "Prover is not active".to_string(),
            ));
        }

        // Step 1:
        // fix argument and evaluate f(x) over x_m = r; where r is the challenge
        // for the current round, and m is the round number, indexed from 1
        //
        // i.e.:
        // at round m <= n, for each mle g(x_1, ... x_n) within the flattened_mle
        // which has already been evaluated to
        //
        //    g(r_1, ..., r_{m-1}, x_m ... x_n)
        //
        // eval g over r_m, and mutate g to g(r_1, ... r_m,, x_{m+1}... x_n)
        let mut flattened_ml_extensions: Vec<DenseMultilinearExtension<F>> = self
            .sum_check_prover_state
            .poly
            .flattened_ml_extensions
            .par_iter()
            .map(|x| x.as_ref().clone())
            .collect();

        if let Some(chal) = challenge {
            if self.sum_check_prover_state.round == 0 {
                return Err(PolyIOPErrors::InvalidProver(
                    "first round should be prover first.".to_string(),
                ));
            }
            self.sum_check_prover_state.challenges.push(*chal);

            let r = self.sum_check_prover_state.challenges[self.sum_check_prover_state.round - 1];
            #[cfg(feature = "parallel")]
            flattened_ml_extensions
                .par_iter_mut()
                .for_each(|mle| *mle = fix_variables(mle, &[r]));
            #[cfg(not(feature = "parallel"))]
            flattened_ml_extensions
                .iter_mut()
                .for_each(|mle| *mle = fix_variables(mle, &[r]));
            self.current_sum += interpolate_uni_poly(&self.mask_poly.evaluations[self.sum_check_prover_state.round-1], r)?;
        } else if self.sum_check_prover_state.round > 0 {
            return Err(PolyIOPErrors::InvalidProver(
                "verifier message is empty".to_string(),
            ));
        }
        // end_timer!(fix_argument);

        let mut temp = self.current_sum
            * F::from((1 << (self.sum_check_prover_state.poly.aux_info.num_variables-self.sum_check_prover_state.round-1)) as u64);
        
        if self.sum_check_prover_state.poly.aux_info.num_variables-1 != self.sum_check_prover_state.round {
            temp += self.sum_aux[self.sum_check_prover_state.round+1]
                * F::from((1 << (self.sum_check_prover_state.poly.aux_info.num_variables-self.sum_check_prover_state.round-2)) as u64);
        }

        let mut g_sum: Vec<F> = vec![F::zero(); self.mask_poly.evaluations[0].len()];
        for j in 0..self.mask_poly.evaluations[0].len() {
            g_sum[j] = self.mask_poly.evaluations[self.sum_check_prover_state.round][j]
                * F::from((1 << (self.sum_check_prover_state.poly.aux_info.num_variables-self.sum_check_prover_state.round-1)) as u64)
                + temp;
        }

        self.sum_check_prover_state.round += 1;

        let max_degree = max([self.sum_check_prover_state.poly.aux_info.max_degree, self.mask_poly.evaluations[0].len()-1]).unwrap();

        let products_list = self.sum_check_prover_state.poly.products.clone();
        let mut products_sum = vec![F::zero(); max_degree + 1];

        // Step 2: generate sum for the partial evaluated polynomial:
        // f(r_1, ... r_m,, x_{m+1}... x_n)

        products_list.iter().for_each(|(coefficient, products)| {
            let mut sum = cfg_into_iter!(0..1 << (self.sum_check_prover_state.poly.aux_info.num_variables - self.sum_check_prover_state.round))
                .fold(
                    || {
                        (
                            vec![(F::zero(), F::zero()); products.len()],
                            vec![F::zero(); products.len() + 1],
                        )
                    },
                    |(mut buf, mut acc), b| {
                        buf.iter_mut()
                            .zip(products.iter())
                            .for_each(|((eval, step), f)| {
                                let table = &flattened_ml_extensions[*f];
                                *eval = table[b << 1];
                                *step = table[(b << 1) + 1] - table[b << 1];
                            });
                        acc[0] += buf.iter().map(|(eval, _)| eval).product::<F>();
                        acc[1..].iter_mut().for_each(|acc| {
                            buf.iter_mut().for_each(|(eval, step)| *eval += step as &_);
                            *acc += buf.iter().map(|(eval, _)| eval).product::<F>();
                        });
                        (buf, acc)
                    },
                )
                .map(|(_, partial)| partial)
                .reduce(
                    || vec![F::zero(); products.len() + 1],
                    |mut sum, partial| {
                        sum.iter_mut()
                            .zip(partial.iter())
                            .for_each(|(sum, partial)| *sum += partial);
                        sum
                    },
                );
            sum.iter_mut().for_each(|sum| *sum *= coefficient);
            let extraploation = cfg_into_iter!(0..max_degree - products.len())
                .map(|i| {
                    let (points, weights) = &self.sum_check_prover_state.extrapolation_aux[products.len() - 1];
                    let at = F::from((products.len() + 1 + i) as u64);
                    extrapolate(points, weights, &sum, &at)
                })
                .collect::<Vec<_>>();
            products_sum
                .iter_mut()
                .zip(sum.iter().chain(extraploation.iter()))
                .for_each(|(products_sum, sum)| *products_sum += sum);
        });

        // update prover's state to the partial evaluated polynomial
        self.sum_check_prover_state.poly.flattened_ml_extensions = flattened_ml_extensions
            .par_iter()
            .map(|x| Arc::new(x.clone()))
            .collect();

        assert_eq!(products_sum.len(), g_sum.len());
        for j in 0..products_sum.len() {
            products_sum[j] += *rho * g_sum[j];
        }

        Ok(IOPProverMessage {
            evaluations: products_sum,
        })
    }
}

fn barycentric_weights<F: PrimeField>(points: &[F]) -> Vec<F> {
    let mut weights = points
        .iter()
        .enumerate()
        .map(|(j, point_j)| {
            points
                .iter()
                .enumerate()
                .filter(|&(i, _point_i)| (i != j))
                .map(|(_i, point_i)| *point_j - point_i)
                .reduce(|acc, value| acc * value)
                .unwrap_or_else(F::one)
        })
        .collect::<Vec<_>>();
    batch_inversion(&mut weights);
    weights
}

fn extrapolate<F: PrimeField>(points: &[F], weights: &[F], evals: &[F], at: &F) -> F {
    let (coeffs, sum_inv) = {
        let mut coeffs = points.iter().map(|point| *at - point).collect::<Vec<_>>();
        batch_inversion(&mut coeffs);
        coeffs.iter_mut().zip(weights).for_each(|(coeff, weight)| {
            *coeff *= weight;
        });
        let sum_inv = coeffs.iter().sum::<F>().inverse().unwrap_or_default();
        (coeffs, sum_inv)
    };
    coeffs
        .iter()
        .zip(evals)
        .map(|(coeff, eval)| *coeff * eval)
        .sum::<F>()
        * sum_inv
}

/// Interpolate a uni-variate degree-`p_i.len()-1` polynomial and evaluate this
/// polynomial at `eval_at`:
///
///   \sum_{i=0}^len p_i * (\prod_{j!=i} (eval_at - j)/(i-j) )
///
/// This implementation is linear in number of inputs in terms of field
/// operations. It also has a quadratic term in primitive operations which is
/// negligible compared to field operations.
/// TODO: The quadratic term can be removed by precomputing the lagrange
/// coefficients.
fn interpolate_uni_poly<F: PrimeField>(p_i: &[F], eval_at: F) -> Result<F, PolyIOPErrors> {
    let start = start_timer!(|| "sum check interpolate uni poly opt");

    let len = p_i.len();
    let mut evals = vec![];
    let mut prod = eval_at;
    evals.push(eval_at);

    // `prod = \prod_{j} (eval_at - j)`
    for e in 1..len {
        let tmp = eval_at - F::from(e as u64);
        evals.push(tmp);
        prod *= tmp;
    }
    let mut res = F::zero();
    // we want to compute \prod (j!=i) (i-j) for a given i
    //
    // we start from the last step, which is
    //  denom[len-1] = (len-1) * (len-2) *... * 2 * 1
    // the step before that is
    //  denom[len-2] = (len-2) * (len-3) * ... * 2 * 1 * -1
    // and the step before that is
    //  denom[len-3] = (len-3) * (len-4) * ... * 2 * 1 * -1 * -2
    //
    // i.e., for any i, the one before this will be derived from
    //  denom[i-1] = denom[i] * (len-i) / i
    //
    // that is, we only need to store
    // - the last denom for i = len-1, and
    // - the ratio between current step and fhe last step, which is the product of
    //   (len-i) / i from all previous steps and we store this product as a fraction
    //   number to reduce field divisions.

    // We know
    //  - 2^61 < factorial(20) < 2^62
    //  - 2^122 < factorial(33) < 2^123
    // so we will be able to compute the ratio
    //  - for len <= 20 with i64
    //  - for len <= 33 with i128
    //  - for len >  33 with BigInt
    if p_i.len() <= 20 {
        let last_denominator = F::from(u64_factorial(len - 1));
        let mut ratio_numerator = 1i64;
        let mut ratio_denominator = 1u64;

        for i in (0..len).rev() {
            let ratio_numerator_f = if ratio_numerator < 0 {
                -F::from((-ratio_numerator) as u64)
            } else {
                F::from(ratio_numerator as u64)
            };

            res += p_i[i] * prod * F::from(ratio_denominator)
                / (last_denominator * ratio_numerator_f * evals[i]);

            // compute denom for the next step is current_denom * (len-i)/i
            if i != 0 {
                ratio_numerator *= -(len as i64 - i as i64);
                ratio_denominator *= i as u64;
            }
        }
    } else if p_i.len() <= 33 {
        let last_denominator = F::from(u128_factorial(len - 1));
        let mut ratio_numerator = 1i128;
        let mut ratio_denominator = 1u128;

        for i in (0..len).rev() {
            let ratio_numerator_f = if ratio_numerator < 0 {
                -F::from((-ratio_numerator) as u128)
            } else {
                F::from(ratio_numerator as u128)
            };

            res += p_i[i] * prod * F::from(ratio_denominator)
                / (last_denominator * ratio_numerator_f * evals[i]);

            // compute denom for the next step is current_denom * (len-i)/i
            if i != 0 {
                ratio_numerator *= -(len as i128 - i as i128);
                ratio_denominator *= i as u128;
            }
        }
    } else {
        let mut denom_up = field_factorial::<F>(len - 1);
        let mut denom_down = F::one();

        for i in (0..len).rev() {
            res += p_i[i] * prod * denom_down / (denom_up * evals[i]);

            // compute denom for the next step is current_denom * (len-i)/i
            if i != 0 {
                denom_up *= -F::from((len - i) as u64);
                denom_down *= F::from(i as u64);
            }
        }
    }
    end_timer!(start);
    Ok(res)
}

/// compute the factorial(a) = 1 * 2 * ... * a
#[inline]
fn field_factorial<F: PrimeField>(a: usize) -> F {
    let mut res = F::one();
    for i in 2..=a {
        res *= F::from(i as u64);
    }
    res
}

/// compute the factorial(a) = 1 * 2 * ... * a
#[inline]
fn u128_factorial(a: usize) -> u128 {
    let mut res = 1u128;
    for i in 2..=a {
        res *= i as u128;
    }
    res
}

/// compute the factorial(a) = 1 * 2 * ... * a
#[inline]
fn u64_factorial(a: usize) -> u64 {
    let mut res = 1u64;
    for i in 2..=a {
        res *= i as u64;
    }
    res
}
