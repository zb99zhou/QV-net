#![allow(non_snake_case)]

use std::{marker::PhantomData, sync::Arc};

use arithmetic::{build_eq_x_r_vec, DenseMultilinearExtension, VPAuxInfo, VirtualPolynomial};
use ark_ff::{BigInteger, PrimeField, Zero};
use curv::{arithmetic::{Converter, Modulo}, cryptographic_primitives::hashing::DigestExt, elliptic::curves::{Point, Scalar, Secp256k1}, BigInt};
use ark_secp256k1::Fr;
use merlin::Transcript;
use sha2::{Digest, Sha512};
use subroutines::{IOPProof, PolyIOP, RandomMaskPolynomial, ZkSumCheck};
use crate::{zk_eval::{pad_point_to_power_of_two, pad_scalar_to_power_of_two}, Errors};

use crate::{sigma_dl::generate_random_point, transcript::TranscriptProtocol, zk_eval::ZkEval};

#[derive(Clone)]
pub struct ZkSumSquareArg {
    pub(super) C_g_star: Point<Secp256k1>,
    pub(super) C_R: Point<Secp256k1>,
    pub(super) C_hat: Point<Secp256k1>,
    pub(super) z: Scalar<Secp256k1>,
    pub(super) z_tag: Scalar<Secp256k1>,
    pub(super) f_star_eval: ZkEval,
    pub(super) g_star_eval: ZkEval,
    pub(super) zk_sum_check_proof: IOPProof<Fr>
}

pub fn Scalar2Fr<F: PrimeField>(val: &Scalar<Secp256k1>) -> F {
    if *val == Scalar::<Secp256k1>::zero() {
        return F::zero();
    }
    let bytes = val.to_bigint().to_bytes();
    F::from_be_bytes_mod_order(&bytes)
}

pub fn Fr2Scalar<F: PrimeField>(val: &F) -> Scalar<Secp256k1> {
    let bytes = val.into_bigint().to_bytes_be();
    Scalar::<Secp256k1>::from_bytes(&bytes).unwrap()
}

pub fn poly_eval(
    coeff: &[Scalar<Secp256k1>],
    x_power: &[Scalar<Secp256k1>]
) -> Scalar<Secp256k1> {
    assert_eq!(coeff.len(), x_power.len());

    let out = Scalar::<Secp256k1>::zero();
    let out = coeff.iter().zip(x_power).fold(out, |acc, x| {
        let tmp = x.0 * x.1;
        acc + tmp
    });
    out
}

pub fn build_mle_list<F: PrimeField>(
    multiplicands: Vec<Vec<F>>,
    nv: usize,
    degree: usize,
) -> (Vec<Arc<DenseMultilinearExtension<F>>>, F) {
    assert_eq!(multiplicands.len(), degree);
    assert_eq!(multiplicands[0].len(), 1 << nv);

    let mut sum = F::zero();

    for i in 0..(1 << nv) {
        let mut product = F::one();

        for e in multiplicands.iter() {
            let val = e[i];
            product *= val;
        }
        sum += product;
    }

    let list = multiplicands
        .into_iter()
        .map(|x| Arc::new(DenseMultilinearExtension::from_evaluations_vec(nv, x)))
        .collect();

    (list, sum)
}

impl ZkSumSquareArg {
    pub fn prove(
        transcript: &mut Transcript,
        g_vec: &[Point<Secp256k1>],
        Y_vec: &[Point<Secp256k1>],
        g: &Point<Secp256k1>,
        h: &Point<Secp256k1>,
        B: &Point<Secp256k1>,
        B_hat: &Point<Secp256k1>,
        v_vec: &[Scalar<Secp256k1>],
        x_vec: &[Scalar<Secp256k1>],
        u: &Scalar<Secp256k1>,
        r: &Scalar<Secp256k1>,
        n: usize,
        seed: &BigInt
    ) -> ZkSumSquareArg {
        assert_eq!(g_vec.len(), n);
        assert_eq!(Y_vec.len(), n);
        assert_eq!(v_vec.len(), n);
        assert_eq!(x_vec.len(), n);

        transcript.zk_sum_square_domain_sep(n as u64);
        transcript.append_points_array(b"g_vec", g_vec);
        transcript.append_points_array(b"Y_vec", Y_vec);
        transcript.append_point(b"g", g);
        transcript.append_point(b"h", h);
        transcript.append_point(b"B", B);
        transcript.append_point(b"B_hat", B_hat);

        let mut g_vec = g_vec.to_vec();
        let mut Y_vec = Y_vec.to_vec();
        let mut v_vec = v_vec.to_vec();
        let mut x_vec = x_vec.to_vec();
        pad_point_to_power_of_two(&mut g_vec, seed);
        pad_point_to_power_of_two(&mut Y_vec, &(seed+BigInt::from((g_vec.len()) as u64)));
        pad_scalar_to_power_of_two(&mut v_vec);
        pad_scalar_to_power_of_two(&mut x_vec);
        let mut n = n.next_power_of_two();
        if n == 1 {
            n = 2;
        }
        let seed = &(seed+BigInt::from((2*g_vec.len()) as u64));

        let num_variables = (n as f64).log2() as usize;

        // sample random g_star, a `g_star_const_term` and 6 * \ell `g_star_coeff`
        let g_star_const_term = Scalar::<Secp256k1>::random();
        let g_star_coeff = vec![vec![Scalar::<Secp256k1>::random(); 6]; num_variables];
        // sample `w`, a single scalar, R(X) is (`w`, 0^{n-1})'s MLE 
        let w = Scalar::<Secp256k1>::random();
        
        // evaluate the \ell uni-polys in g_star in 7 points (degree is 6)
        let mut g_star_eval = vec![vec![Scalar::<Secp256k1>::zero(); 7]; num_variables];
        // get [[1^1, 1^2, ..., 1^6],..., [6^1, 6^2, ..., 6^6]]
        // to evaluate the uni-poly, calculate the inner product of
        // x_power's component and non-zero coeff
        let mut x_power = vec![vec![Scalar::<Secp256k1>::zero(); 6]; 6];
        let order = Scalar::<Secp256k1>::group_order();
        for i in 1..7 {
            for j in 1..7 {
                let i_bn = BigInt::from(i as i32);
                let j_bn = BigInt::from(j as i32);
                x_power[i-1][j-1] = Scalar::<Secp256k1>::from_bigint(&BigInt::mod_pow(&i_bn, &j_bn, order));
            }
        }
        for i in 0..g_star_eval.len() {
            for j in 0..6 {
                g_star_eval[i][j+1] = poly_eval(&g_star_coeff[i], &x_power[j]);
            }
        }
        
        // convert the `Scalar<Secp256k1>` to `Fr` and construct the poly
        let mut g_star_eval_fr = vec![vec![Fr::from(0); 7]; num_variables];
        for i in 0..g_star_eval.len() {
            for j in 0..g_star_eval[0].len() {
                g_star_eval_fr[i][j] = Scalar2Fr::<Fr>(&g_star_eval[i][j]);
            }
        }
        let g_star_poly = RandomMaskPolynomial::<Fr>{
            const_term: Scalar2Fr::<Fr>(&g_star_const_term),
            evaluations: g_star_eval_fr
        };
        
        // calculate `u_star` and sample `r_star`
        let mut sum: Fr = g_star_poly.evaluations.iter().map(|row| row[1]).sum();
        sum *= Fr::from((1 << (num_variables-1)) as u64);
        sum += Fr::from((1 << num_variables) as u64) * g_star_poly.const_term;
        let u_star = Fr2Scalar(&sum);
        
        // sample `h_vec` and `Y` to commit `g_star`
        let h_vec = (0..6*num_variables+1)
            .map(|i| {
                let kzen_label_i = BigInt::from(i as u32) + seed;
                let hash_i = Sha512::new().chain_bigint(&kzen_label_i).result_bigint();
                generate_random_point(&Converter::to_bytes(&hash_i))
            })
            .collect::<Vec<Point<Secp256k1>>>();
        let kzen_label = BigInt::from((6*num_variables+1) as i32) + seed;
        let hash = Sha512::new().chain_bigint(&kzen_label).result_bigint();
        let Y = generate_random_point(&Converter::to_bytes(&hash));

        // commitment to `g_star`
        let r_v = Scalar::<Secp256k1>::random();
        let mut v_star: Vec<Scalar<Secp256k1>> = g_star_coeff.into_iter().flatten().collect();
        v_star.insert(0, g_star_const_term);
        let C_g_star = h_vec.iter().zip(v_star.clone()).fold(&Y * r_v.clone(), |acc, x| {
            if x.1 != Scalar::<Secp256k1>::zero(){
                let hi_vi: Point<Secp256k1> = x.0 * x.1;
                acc + &hi_vi
            } else {
                acc
            }
        });

        // commitment to `R`
        let ra_vec: Vec<Scalar<Secp256k1>> = (0..n)
            .map(|_| {
                
                Scalar::<Secp256k1>::random()
            })
            .collect();
        let C_R = Y_vec.iter().zip(ra_vec.clone()).fold(&g_vec[0] * w.clone(), |acc, x| {
            if x.1 != Scalar::<Secp256k1>::zero(){
                let Yi_ri: Point<Secp256k1> = x.0 * x.1;
                acc + &Yi_ri
            } else {
                acc
            }
        });
        
        // commitment to sum
        let r_star = Scalar::<Secp256k1>::random();
        let C_hat = g * u_star.clone() + h * r_star.clone();
        
        // random challenge `rho`
        let rho = transcript.challenge_scalar(b"rho");
        
        // calculate `z` and `z_tag`
        let z = u + rho.clone() * u_star;
        let z_tag = r + rho.clone() * r_star;

        // construct f_star * f_star = f * f + 2 * f * R * Z + R * Z * R * Z
        let mut poly = VirtualPolynomial::new(num_variables);
        let mut sum = Fr::from(0);
        // f
        let mut f_fr = Vec::with_capacity(n);
        for v in v_vec.iter().take(n) {
            f_fr.push(Scalar2Fr::<Fr>(v));
        }
        // R
        let mut R_fr = vec![Fr::from(0); n];
        R_fr[0] = Scalar2Fr(&w);
        // Z
        let mut Z0 = vec![Fr::from(0); n];
        Z0[0] = Fr::from(1);
        let mut Z1 = vec![Fr::from(0); n];
        Z1[n-1] = Fr::from(1);
        // f * f
        let f_f = [f_fr.clone(), f_fr.clone()].to_vec();
        let (product, product_sum) = build_mle_list(f_f, num_variables, 2);
        let coefficient = Fr::from(1);
        poly.add_mle_list(product, coefficient).unwrap();
        sum += product_sum * coefficient;
        // 2 * f * R * Z
        let f_R_Z = [f_fr, R_fr.clone(), Z0.clone(), Z1.clone()].to_vec();
        let (product, product_sum) = build_mle_list(f_R_Z, num_variables, 4);
        let coefficient = Fr::from(2);
        poly.add_mle_list(product, coefficient).unwrap();
        sum += product_sum * coefficient;
        // R * Z * R * Z
        let R_Z_R_Z = [R_fr.clone(), Z0.clone(), Z1.clone(), R_fr, Z0.clone(), Z1.clone()].to_vec();
        let (product, product_sum) = build_mle_list(R_Z_R_Z, num_variables, 6);
        let coefficient = Fr::from(1);
        poly.add_mle_list(product, coefficient).unwrap();
        sum += product_sum * coefficient;

        assert_eq!(Fr2Scalar(&sum), *u);

        // do zk-sum-check
        let mut sum_check_transcript = <PolyIOP<Fr> as ZkSumCheck<Fr>>::init_transcript();
        let proof = <PolyIOP<Fr> as ZkSumCheck<Fr>>::prove(&poly, &g_star_poly, &Scalar2Fr(&rho), &mut sum_check_transcript).unwrap();
        
        // get `challenge_points` in zk-sum-check
        let challenge_points: Vec<Fr> = proof.clone().point;

        // final check by `zkEval`
        // get Z(r)
        let mut Z_poly = VirtualPolynomial::new(num_variables);
        let Z = [Z0, Z1].to_vec();
        let (product, _) = build_mle_list(Z, num_variables, 2);
        let coefficient = Fr::from(1);
        Z_poly.add_mle_list(product, coefficient).unwrap();
        let Z_r: Fr = Z_poly.evaluate(&challenge_points).unwrap();
        let Z_r = Fr2Scalar(&Z_r);

        // get eq(r, x)
        let eq_r = build_eq_x_r_vec(&challenge_points).unwrap();
        let eq_r = (0..n)
            .map(|i| Fr2Scalar(&eq_r[i]))
            .collect::<Vec<Scalar<Secp256k1>>>();

        // get f(X) + R(X) * Z(r)
        let mut f_star_vec = v_vec.to_vec();
        f_star_vec[0] = v_vec[0].clone() + w * Z_r.clone();
        
        // calculate x + Z(r) * r, the exponent of Y in B * C_R ^ Z(r)
        let x_vec = x_vec.to_vec();
        let x_raZ_vec = (0..n)
            .map(|i| x_vec[i].clone() + ra_vec[i].clone() * Z_r.clone())
            .collect::<Vec<Scalar<Secp256k1>>>();
        
        // `zkEval` for `f_star`
        let secret_vec_f_star: Vec<Scalar<Secp256k1>> = f_star_vec.iter().chain(x_raZ_vec.iter()).cloned().collect();
        let public_vec_f_star: Vec<Scalar<Secp256k1>> = eq_r.iter().chain(vec![Scalar::<Secp256k1>::zero(); n].iter()).cloned().collect();
        let base_f_star: Vec<Point<Secp256k1>> = g_vec.iter().chain(Y_vec.iter()).cloned().collect();
        let a_vec = (0..2*n)
            .map(|i| secret_vec_f_star[i].to_bigint())
            .collect::<Vec<BigInt>>();
        let b_vec = (0..2*n)
            .map(|i| public_vec_f_star[i].to_bigint())
            .collect::<Vec<BigInt>>();
        let f_star_eval = ZkEval::eval(transcript, &base_f_star, &(B + &C_R * Z_r), &a_vec, &b_vec, 2*n, &(seed + BigInt::from((6*num_variables+2) as i32)));

        // get `r_star`
        let mut b_vec = Vec::with_capacity(6*num_variables+2);
        b_vec.push(BigInt::from(1));
        for i in 0..challenge_points.len() {
            let temp = Fr2Scalar(&challenge_points[i]).to_bigint();
            for j in 1..7 {
                let temp_power = BigInt::mod_pow(&temp, &BigInt::from(j), order);
                b_vec.push(temp_power);
            }
        }
        b_vec.push(BigInt::zero());
        let mut a_vec = (0..6*num_variables+1)
            .map(|i| v_star[i].to_bigint())
            .collect::<Vec<BigInt>>();
        a_vec.push(r_v.to_bigint());
        let mut base_g_star = h_vec;
        base_g_star.push(Y);
        let g_star_eval = ZkEval::eval(transcript, &base_g_star, &C_g_star, &a_vec, &b_vec, 6*num_variables+2, &(seed + BigInt::from((6*num_variables+3) as i32)));
        
        ZkSumSquareArg {
            C_g_star,
            C_R,
            C_hat,
            z,
            z_tag,
            f_star_eval,
            g_star_eval,
            zk_sum_check_proof: proof
        }
    }

    pub fn verify(
        &self,
        transcript: &mut Transcript,
        g_vec: &[Point<Secp256k1>],
        Y_vec: &[Point<Secp256k1>],
        g: &Point<Secp256k1>,
        h: &Point<Secp256k1>,
        B: &Point<Secp256k1>,
        B_hat: &Point<Secp256k1>,
        n: usize,
        seed: &BigInt
    ) -> Result<(), Errors> {
        assert_eq!(g_vec.len(), n);
        assert_eq!(Y_vec.len(), n);

        transcript.zk_sum_square_domain_sep(n as u64);
        transcript.append_points_array(b"g_vec", g_vec);
        transcript.append_points_array(b"Y_vec", Y_vec);
        transcript.append_point(b"g", g);
        transcript.append_point(b"h", h);
        transcript.append_point(b"B", B);
        transcript.append_point(b"B_hat", B_hat);

        let mut g_vec = g_vec.to_vec();
        let mut Y_vec = Y_vec.to_vec();
        pad_point_to_power_of_two(&mut g_vec, seed);
        pad_point_to_power_of_two(&mut Y_vec, &(seed+BigInt::from((g_vec.len()) as u64)));
        let mut n = n.next_power_of_two();
        if n == 1 {
            n = 2;
        }
        let seed = &(seed+BigInt::from((2*g_vec.len()) as u64));

        let num_variables = (n as f64).log2() as usize;

        // random challenge `rho`
        let rho = transcript.challenge_scalar(b"rho");

        let order = Scalar::<Secp256k1>::group_order();
        
        if B_hat + &self.C_hat * rho.clone() != g * self.z.clone() + h * self.z_tag.clone() {
            return Err(Errors::ZkSumSquareArgError);
        }
        
        let mut sum_check_transcript = <PolyIOP<Fr> as ZkSumCheck<Fr>>::init_transcript();
        let poly_info = VPAuxInfo {
            max_degree: 6,
            num_variables,
            phantom: PhantomData,
        };
        let subclaim = <PolyIOP<Fr> as ZkSumCheck<Fr>>::verify(
            Scalar2Fr(&self.z),
            &self.zk_sum_check_proof,
            &poly_info,
            &mut sum_check_transcript,
            num_variables,
            6
        ).unwrap();
        let challenge_points = subclaim.point;

        // sample `h_vec` and `Y` to commit `g_star`
        let h_vec = (0..6*num_variables+1)
            .map(|i| {
                let kzen_label_i = BigInt::from(i as u32) + seed;
                let hash_i = Sha512::new().chain_bigint(&kzen_label_i).result_bigint();
                generate_random_point(&Converter::to_bytes(&hash_i))
            })
            .collect::<Vec<Point<Secp256k1>>>();
        let kzen_label = BigInt::from((6*num_variables+1) as i32) + seed;
        let hash = Sha512::new().chain_bigint(&kzen_label).result_bigint();
        let Y = generate_random_point(&Converter::to_bytes(&hash));

        // Z
        let mut Z0 = vec![Fr::from(0); n];
        Z0[0] = Fr::from(1);
        let mut Z1 = vec![Fr::from(0); n];
        Z1[n-1] = Fr::from(1);
        // get Z(r)
        let mut Z_poly = VirtualPolynomial::new(num_variables);
        let Z = [Z0, Z1].to_vec();
        let (product, _) = build_mle_list(Z, num_variables, 2);
        let coefficient = Fr::from(1);
        Z_poly.add_mle_list(product, coefficient).unwrap();
        let Z_r: Fr = Z_poly.evaluate(&challenge_points).unwrap();
        let Z_r = Fr2Scalar(&Z_r);

        // get eq(r, x)
        let eq_r = build_eq_x_r_vec(&challenge_points).unwrap();
        let eq_r = (0..n)
            .map(|i| Fr2Scalar(&eq_r[i]))
            .collect::<Vec<Scalar<Secp256k1>>>();

        let public_vec_f_star: Vec<Scalar<Secp256k1>> = eq_r.iter().chain(vec![Scalar::<Secp256k1>::zero(); n].iter()).cloned().collect();
        let base_f_star: Vec<Point<Secp256k1>> = g_vec.iter().chain(Y_vec.iter()).cloned().collect();
        let b_vec = (0..2*n)
            .map(|i| public_vec_f_star[i].to_bigint())
            .collect::<Vec<BigInt>>();

        let f_eval_verifier = self.f_star_eval.verify(
            transcript,
            &base_f_star,
            &(B + &self.C_R * Z_r), 
            &b_vec, 
            2 * n, 
            &(seed + BigInt::from((6*num_variables+2) as i32))
        );

        // get `r_star`
        let mut b_vec = Vec::with_capacity(6*num_variables+2);
        b_vec.push(BigInt::from(1));
        for challenge in challenge_points.iter() {
            let temp = Fr2Scalar(challenge).to_bigint();
            for j in 1..7 {
                let temp_power = BigInt::mod_pow(&temp, &BigInt::from(j), order);
                b_vec.push(temp_power);
            }
        }
        b_vec.push(BigInt::zero());
        let mut base_g_star = h_vec;
        base_g_star.push(Y);
        let g_eval_verifier = self.g_star_eval.verify(
            transcript,
            &base_g_star,
            &self.C_g_star, 
            &b_vec, 
            6*num_variables+2, 
            &(seed + BigInt::from((6*num_variables+3) as i32))
        );

        if f_eval_verifier.is_ok() && g_eval_verifier.is_ok() {
            let v_f = Scalar::<Secp256k1>::from_bigint(&self.f_star_eval.t);
            let v_g = Scalar::<Secp256k1>::from_bigint(&self.g_star_eval.t);
            let eval = v_f.clone() * v_f + rho * v_g;
            if Fr2Scalar(&subclaim.expected_evaluation) == eval {
                Ok(())
            } else {
                Err(Errors::ZkSumSquareArgError)
            }
        } else {
            Err(Errors::ZkSumSquareArgError)
        }
    }
}

mod test {
    use curv::{arithmetic::Converter, cryptographic_primitives::hashing::DigestExt, elliptic::curves::{Point, Scalar, Secp256k1}, BigInt};
    use merlin::Transcript;
    use sha2::{Digest, Sha512};

    use crate::sigma_dl::generate_random_point;

    use super::ZkSumSquareArg;

    pub fn test_helper(seed: &BigInt, n: usize) {
        let g_vec = (0..n)
            .map(|i| {
                let kzen_label_i = BigInt::from(i as u32) + seed;
                let hash_i = Sha512::new().chain_bigint(&kzen_label_i).result_bigint();
                generate_random_point(&Converter::to_bytes(&hash_i))
            })
            .collect::<Vec<Point<Secp256k1>>>();

        let Y_vec = (0..n)
            .map(|i| {
                let kzen_label_j = BigInt::from(n as u32) + BigInt::from(i as u32) + seed;
                let hash_j = Sha512::new().chain_bigint(&kzen_label_j).result_bigint();
                generate_random_point(&Converter::to_bytes(&hash_j))
            })
            .collect::<Vec<Point<Secp256k1>>>();

        let kzen_label = BigInt::from((2*n) as u32) + seed;
        let hash = Sha512::new().chain_bigint(&kzen_label).result_bigint();
        let g = generate_random_point(&Converter::to_bytes(&hash));

        let kzen_label = BigInt::from((2*n+1) as u32) + seed;
        let hash = Sha512::new().chain_bigint(&kzen_label).result_bigint();
        let h = generate_random_point(&Converter::to_bytes(&hash));

        let v_vec: Vec<Scalar<Secp256k1>> = (0..n)
            .map(|_| {
                
                Scalar::<Secp256k1>::random()
            })
            .collect();

        let x_vec: Vec<Scalar<Secp256k1>> = (0..n)
            .map(|_| {
                
                Scalar::<Secp256k1>::random()
            })
            .collect();

        let mut u = Scalar::<Secp256k1>::zero();
        for v in v_vec.iter().take(n) {
            u = u + v * v;
        }

        let r = Scalar::<Secp256k1>::random();

        let B_1 = g_vec.iter().zip(v_vec.clone()).fold(Point::<Secp256k1>::zero(), |acc, x| {
            if x.1 != Scalar::<Secp256k1>::zero(){
                let gi_vi: Point<Secp256k1> = x.0 * x.1;
                acc + &gi_vi
            } else {
                acc
            }
        });

        let B_2 = Y_vec.iter().zip(x_vec.clone()).fold(Point::<Secp256k1>::zero(), |acc, x| {
            if x.1 != Scalar::<Secp256k1>::zero(){
                let gi_vi: Point<Secp256k1> = x.0 * x.1;
                acc + &gi_vi
            } else {
                acc
            }
        });
        let B = B_1 + B_2;

        let B_hat = &g * u.clone() + &h * r.clone();

        let mut verifier = Transcript::new(b"zksstest");
        let proof = ZkSumSquareArg::prove(
            &mut verifier, 
            &g_vec, 
            &Y_vec, 
            &g, 
            &h, 
            &B, 
            &B_hat, 
            &v_vec, 
            &x_vec, 
            &u, 
            &r, 
            n, 
            &(BigInt::from((2*n+2) as u32) + seed)
        );
        let mut verifier = Transcript::new(b"zksstest");
        let result = proof.verify(
            &mut verifier,
             &g_vec, 
             &Y_vec, 
             &g, 
             &h, 
             &B, 
             &B_hat, 
             n, 
             &(BigInt::from((2*n+2) as u32) + seed)
        );
        assert!(result.is_ok());
    }

    #[test]
    pub fn test_zk_ss_4() {
        let KZen: &[u8] = &[75, 90, 101, 110];
        let kzen_label = BigInt::from_bytes(KZen);
        test_helper(&kzen_label, 4);
    }

    #[test]
    pub fn test_zk_ss_8() {
        let KZen: &[u8] = &[75, 90, 101, 110];
        let kzen_label = BigInt::from_bytes(KZen);
        test_helper(&kzen_label, 8);
    }

    #[test]
    pub fn test_zk_ss_16() {
        let KZen: &[u8] = &[75, 90, 101, 110];
        let kzen_label = BigInt::from_bytes(KZen);
        test_helper(&kzen_label, 16);
    }

    #[test]
    pub fn test_zk_ss_32() {
        let KZen: &[u8] = &[75, 90, 101, 110];
        let kzen_label = BigInt::from_bytes(KZen);
        test_helper(&kzen_label, 32);
    }

    #[test]
    pub fn test_zk_ss_64() {
        let KZen: &[u8] = &[75, 90, 101, 110];
        let kzen_label = BigInt::from_bytes(KZen);
        test_helper(&kzen_label, 64);
    }

    #[test]
    pub fn test_zk_ss_3() {
        let KZen: &[u8] = &[75, 90, 101, 110];
        let kzen_label = BigInt::from_bytes(KZen);
        test_helper(&kzen_label, 3);
    }

    #[test]
    pub fn test_zk_ss_7() {
        let KZen: &[u8] = &[75, 90, 101, 110];
        let kzen_label = BigInt::from_bytes(KZen);
        test_helper(&kzen_label, 7);
    }

    #[test]
    pub fn test_zk_ss_15() {
        let KZen: &[u8] = &[75, 90, 101, 110];
        let kzen_label = BigInt::from_bytes(KZen);
        test_helper(&kzen_label, 15);
    }

    #[test]
    pub fn test_zk_ss_31() {
        let KZen: &[u8] = &[75, 90, 101, 110];
        let kzen_label = BigInt::from_bytes(KZen);
        test_helper(&kzen_label, 31);
    }

    #[test]
    pub fn test_zk_ss_55() {
        let KZen: &[u8] = &[75, 90, 101, 110];
        let kzen_label = BigInt::from_bytes(KZen);
        test_helper(&kzen_label, 55);
    }
}