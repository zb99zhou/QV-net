#![allow(non_snake_case)]

use curv::{arithmetic::{Converter, Modulo, Zero}, cryptographic_primitives::hashing::DigestExt, elliptic::curves::{secp256_k1::hash_to_curve::generate_random_point, Point, Scalar, Secp256k1}, BigInt};
use merlin::Transcript;
use sha2::{Digest, Sha512};

use crate::Errors;
use crate::transcript::TranscriptProtocol;

// before invoking zkEval, we should pad the length of vector to power of 2
// pad vectors over the field with 0
pub fn pad_bn_to_power_of_two(bn_vec: &mut Vec<BigInt>) {
    let current_len = bn_vec.len();
    let target_len = current_len.next_power_of_two();
    for _ in current_len..target_len {
        bn_vec.push(BigInt::zero());
    }
}

pub fn pad_scalar_to_power_of_two(scalar_vec: &mut Vec<Scalar<Secp256k1>>) {
    let current_len = scalar_vec.len();
    let mut target_len = current_len.next_power_of_two();
    if target_len == 1 {
        target_len = 2;
    }
    for _ in current_len..target_len {
        scalar_vec.push(Scalar::<Secp256k1>::zero());
    }
}

// pad vectors over the group with random point
pub fn pad_point_to_power_of_two(P_vec: &mut Vec<Point<Secp256k1>>, seed: &BigInt) {
    let current_len = P_vec.len();
    let mut target_len = current_len.next_power_of_two();
    if target_len == 1 {
        target_len = 2;
    }
    for i in current_len..target_len {
        let kzen_label_i = BigInt::from(i as u32) + seed;
        let hash_i = Sha512::new().chain_bigint(&kzen_label_i).result_bigint();
        let r = generate_random_point(&Converter::to_bytes(&hash_i));
        P_vec.push(r);
    }
}

#[derive(Clone)]
pub struct ZkEval {
    pub(super) t: BigInt,
    pub(super) R: Point<Secp256k1>,
    pub(super) s: BigInt,
    pub(super) proof: ZkEvalProof
}

#[derive(Clone)]
pub struct ZkEvalProof {
    pub(super) A_vec: Vec<Point<Secp256k1>>,
    pub(super) B_vec: Vec<Point<Secp256k1>>,
    pub(super) z_vec: Vec<BigInt>
}

pub fn inner_product(a: &[BigInt], b: &[BigInt], modulus: &BigInt) -> BigInt {
    assert_eq!(
        a.len(),
        b.len(),
        "inner_product(a,b): lengths of vectors do not match"
    );
    let out = BigInt::zero();
    let out = a.iter().zip(b).fold(out, |acc, x| {
        let aibi = BigInt::mod_mul(x.0, x.1, modulus);
        BigInt::mod_add(&acc, &aibi, modulus)
    });
    out
}

impl ZkEvalProof {
    pub fn prove(
        transcript: &mut Transcript,
        g_vec: &[Point<Secp256k1>],
        ux: &Point<Secp256k1>,
        b_vec: &[BigInt],
        z_vec: &[BigInt],
        mut A_vec: Vec<Point<Secp256k1>>,
        mut B_vec: Vec<Point<Secp256k1>>
    ) -> ZkEvalProof {
        let n =  g_vec.len();

        // All of the input vectors must have the same length.
        assert_eq!(g_vec.len(), n);
        assert_eq!(b_vec.len(), n);
        assert_eq!(z_vec.len(), n);
        assert!(n.is_power_of_two());

        if n > 4 {
            let order = Scalar::<Secp256k1>::group_order();

            let n = n / 2;
            let (gl_vec, gr_vec) = g_vec.split_at(n);
            let (zl_vec, zr_vec) = z_vec.split_at(n);
            let (bl_vec, br_vec) = b_vec.split_at(n);

            let ux_zl_br = ux * Scalar::<Secp256k1>::from_bigint(&inner_product(zl_vec, br_vec, order));
            let ux_zr_bl = ux * Scalar::<Secp256k1>::from_bigint(&inner_product(zr_vec, bl_vec, order));
            
            let A = gr_vec.iter().zip(zl_vec).fold(ux_zl_br, |acc, x| {
                if x.1 != &BigInt::zero() {
                    let zli = Scalar::<Secp256k1>::from(x.1);
                    let gri_zli: Point<Secp256k1> = x.0 * &zli;
                    acc + &gri_zli
                } else {
                    acc
                }
            });

            let B = gl_vec.iter().zip(zr_vec).fold(ux_zr_bl, |acc, x| {
                if x.1 != &BigInt::zero() {
                    let zri = Scalar::<Secp256k1>::from(x.1);
                    let gli_zri: Point<Secp256k1> = x.0 * &zri;
                    acc + &gli_zri
                } else {
                    acc
                }
            });

            transcript.append_point(b"A", &A);
            transcript.append_point(b"B", &B);

            let ex: Scalar<Secp256k1> = transcript.challenge_scalar(b"ex");
            let ex_bn = ex.to_bigint();

            let z_vec_new = (0..n)
                .map(|i| {
                    let zli = zl_vec[i].clone();
                    let zri_ex = BigInt::mod_mul(&zr_vec[i], &ex_bn, order);
                    BigInt::mod_add(&zli, &zri_ex, order)
                })
                .collect::<Vec<BigInt>>();

            let g_vec_new = (0..n)
                .map(|i| {
                    let gl_ex = &gl_vec[i] * ex.clone();
                    let gr = &gr_vec[i];
                    gl_ex + gr
                })
                .collect::<Vec<Point<Secp256k1>>>();

            let b_vec_new = (0..n)
                .map(|i| {
                    let bli_ex = BigInt::mod_mul(&bl_vec[i], &ex_bn, order);
                    let bri = br_vec[i].clone();
                    BigInt::mod_add(&bli_ex, &bri, order)
                })
                .collect::<Vec<BigInt>>();

            A_vec.push(A);
            B_vec.push(B);
            return ZkEvalProof::prove(transcript, &g_vec_new, ux, &b_vec_new, &z_vec_new, A_vec, B_vec)
        }

        ZkEvalProof {
            A_vec,
            B_vec,
            z_vec: z_vec.to_vec()
        }
    }
    
    pub fn verify(
        &self,
        transcript: &mut Transcript,
        g_vec: &[Point<Secp256k1>],
        ux: &Point<Secp256k1>,
        b_vec: &[BigInt],
        P: &Point<Secp256k1>
    ) -> Result<(), Errors> {
        let n =  g_vec.len();

        // All of the input vectors must have the same length.
        assert_eq!(g_vec.len(), n);
        assert!(n.is_power_of_two());

        let order = Scalar::<Secp256k1>::group_order();

        if n > 4 {
            let n = n / 2;
            let (gl_vec, gr_vec) = g_vec.split_at(n);
            let (bl_vec, br_vec) = b_vec.split_at(n);

            transcript.append_point(b"A", &self.A_vec[0]);
            transcript.append_point(b"B", &self.B_vec[0]);
            let ex: Scalar<Secp256k1> = transcript.challenge_scalar(b"ex");
            let ex_bn = ex.to_bigint();
            let ex_sq_bn = BigInt::mod_mul(&ex_bn, &ex_bn, order);
            let ex_sq_fe = Scalar::<Secp256k1>::from(&ex_sq_bn);

            let g_vec_new = (0..n)
                .map(|i| {
                    let gl_ex = &gl_vec[i] * ex.clone();
                    let gr = &gr_vec[i];
                    gl_ex + gr
                })
                .collect::<Vec<Point<Secp256k1>>>();

            let b_vec_new = (0..n)
                .map(|i| {
                    let bli_ex = BigInt::mod_mul(&bl_vec[i], &ex_bn, order);
                    let bri = br_vec[i].clone();
                    BigInt::mod_add(&bli_ex, &bri, order)
                })
                .collect::<Vec<BigInt>>();

            let Pex = P * ex.clone();
            let Bex_sq = &self.B_vec[0] * ex_sq_fe.clone();
            let P_tag = self.A_vec[0].clone() + Pex + Bex_sq;

            let proof = ZkEvalProof {
                A_vec: self.A_vec[1..].to_vec(),
                B_vec: self.B_vec[1..].to_vec(),
                z_vec: self.z_vec.clone()
            };
            return proof.verify(transcript, &g_vec_new, ux, &b_vec_new, &P_tag)
        }

        let c = inner_product(&self.z_vec, b_vec, order);
        let ux_c = ux * Scalar::<Secp256k1>::from_bigint(&c);

        let g_z_ux_c = g_vec.iter().zip(self.z_vec.clone()).fold(ux_c, |acc, x| {
            if x.1 != BigInt::zero() {
                let zi = Scalar::<Secp256k1>::from(x.1);
                let gi_zi: Point<Secp256k1> = x.0 * &zi;
                acc + &gi_zi
            } else {
                acc
            }
        });

        if P.clone() == g_z_ux_c {
            Ok(())
        } else {
            Err(Errors::ZkEvalProofError)
        }
    }
}

impl ZkEval {
    pub fn eval(
        transcript: &mut Transcript,
        g_vec: &[Point<Secp256k1>],
        C: &Point<Secp256k1>,
        a_vec: &[BigInt],
        b_vec: &[BigInt],
        n: usize,
        seed: &BigInt
    ) -> ZkEval {
        assert_eq!(g_vec.len(), n);
        assert_eq!(a_vec.len(), n);
        assert_eq!(b_vec.len(), n);

        transcript.zk_eval_domain_sep(n as u64);
        transcript.append_point(b"C", C);

        let mut a_vec = a_vec.to_vec();
        let mut b_vec = b_vec.to_vec();
        let mut g_vec = g_vec.to_vec();
        pad_bn_to_power_of_two(&mut a_vec);
        pad_bn_to_power_of_two(&mut b_vec);
        pad_point_to_power_of_two(&mut g_vec, seed);
        let n = n.next_power_of_two();

        let mut r_vec: Vec<Scalar<Secp256k1>> = Vec::with_capacity(n);
        for _j in 0..n {
            r_vec.push(Scalar::<Secp256k1>::random());
        }
        let R = g_vec.iter().zip(r_vec.clone()).fold(Point::<Secp256k1>::zero(), |acc, x| {
            if x.1 != Scalar::<Secp256k1>::zero() {
                let ri = x.1;
                let gi_ri: Point<Secp256k1> = x.0 * &ri;
                acc + &gi_ri
            } else {
                acc
            }
        });

        let r_vec_bn = (0..n)
            .map(|i| r_vec[i].to_bigint())
            .collect::<Vec<BigInt>>();

        let order = Scalar::<Secp256k1>::group_order();

        let s = inner_product(&r_vec_bn, &b_vec, order);

        let challenge_e: Scalar<Secp256k1> = transcript.challenge_scalar(b"e");
        let z_vec = (0..n)
            .map(|i| r_vec[i].clone() + challenge_e.clone() * Scalar::<Secp256k1>::from_bigint(&a_vec[i]))
            .collect::<Vec<Scalar<Secp256k1>>>();
        let z_vec = (0..n)
            .map(|i| z_vec[i].to_bigint())
            .collect::<Vec<BigInt>>();

        let challenge_c: Scalar<Secp256k1> = transcript.challenge_scalar(b"c");

        let kzen_label = &(seed + BigInt::from((g_vec.len()) as u64));
        let hash = Sha512::new().chain_bigint(kzen_label).result_bigint();
        let q = generate_random_point(&Converter::to_bytes(&hash));

        let A_vec: Vec<Point<Secp256k1>> = Vec::with_capacity(n);
        let B_vec: Vec<Point<Secp256k1>> = Vec::with_capacity(n);
        let proof = ZkEvalProof::prove(transcript, &g_vec, &(&q * challenge_c), &b_vec, &z_vec, A_vec, B_vec);

        ZkEval { t: inner_product(&a_vec, &b_vec, order), R, s, proof }
    }

    pub fn verify(
        &self,
        transcript: &mut Transcript,
        g_vec: &[Point<Secp256k1>],
        C: &Point<Secp256k1>,
        b_vec: &[BigInt],
        n: usize,
        seed: &BigInt
    ) -> Result<(), Errors> {
        assert_eq!(g_vec.len(), n);
        assert_eq!(b_vec.len(), n);

        transcript.zk_eval_domain_sep(n as u64);
        transcript.append_point(b"C", C);
        
        let mut b_vec = b_vec.to_vec();
        let mut g_vec = g_vec.to_vec();
        pad_bn_to_power_of_two(&mut b_vec);
        pad_point_to_power_of_two(&mut g_vec, seed);

        let challenge_e: Scalar<Secp256k1> = transcript.challenge_scalar(b"e");
        let challenge_c: Scalar<Secp256k1> = transcript.challenge_scalar(b"c");

        let s_fe = Scalar::<Secp256k1>::from_bigint(&self.s);

        let kzen_label = &(seed + BigInt::from((g_vec.len()) as u64));
        let hash = Sha512::new().chain_bigint(kzen_label).result_bigint();
        let q = generate_random_point(&Converter::to_bytes(&hash));
        let q_c = &q * challenge_c;
        let q_c_s_plus_et = &q_c * (s_fe + challenge_e.clone() * Scalar::<Secp256k1>::from_bigint(&self.t));
        
        self.proof.verify(transcript, &g_vec, &q_c, &b_vec, &(self.R.clone()+C * challenge_e+q_c_s_plus_et))
    }
}

mod test {
    use curv::{arithmetic::{Converter, Zero}, cryptographic_primitives::hashing::DigestExt, elliptic::curves::{Point, Scalar, Secp256k1}, BigInt};
    use merlin::Transcript;
    use sha2::{Digest, Sha512};

    use crate::{sigma_dl::generate_random_point, zk_eval::ZkEval};

    pub fn test_helper(seed: &BigInt, n: usize) {
        let g_vec = (0..n)
            .map(|i| {
                let kzen_label_i = BigInt::from(i as u32) + seed;
                let hash_i = Sha512::new().chain_bigint(&kzen_label_i).result_bigint();
                generate_random_point(&Converter::to_bytes(&hash_i))
            })
            .collect::<Vec<Point<Secp256k1>>>();

        let a_vec: Vec<Scalar<Secp256k1>> = (0..n)
            .map(|_| {
                
                Scalar::<Secp256k1>::random()
            })
            .collect();

        let b_vec: Vec<Scalar<Secp256k1>> = (0..n)
            .map(|_| {
                
                Scalar::<Secp256k1>::random()
            })
            .collect();

        let a_vec = (0..n)
            .map(|i| a_vec[i].to_bigint())
            .collect::<Vec<BigInt>>();
        
        let b_vec = (0..n)
            .map(|i| b_vec[i].to_bigint())
            .collect::<Vec<BigInt>>();

        let C = g_vec.iter().zip(a_vec.clone()).fold(Point::<Secp256k1>::zero(), |acc, x| {
            if x.1 != BigInt::zero() {
                let ai = Scalar::<Secp256k1>::from(x.1);
                let gi_ai: Point<Secp256k1> = x.0 * &ai;
                acc + &gi_ai
            } else {
                acc
            }
        });

        let mut verifier = Transcript::new(b"evaltest");
        let val = ZkEval::eval(&mut verifier, &g_vec, &C, &a_vec, &b_vec, n, &(BigInt::from(n as u32) + seed));
        let mut verifier = Transcript::new(b"evaltest");
        let result = val.verify(&mut verifier, &g_vec, &C, &b_vec, n, &(BigInt::from(n as u32) + seed));
        assert!(result.is_ok());
    }

    #[test]
    pub fn test_zk_eval_2() {
        let KZen: &[u8] = &[75, 90, 101, 110];
        let kzen_label = BigInt::from_bytes(KZen);
        test_helper(&kzen_label, 2);
    }

    #[test]
    pub fn test_zk_eval_4() {
        let KZen: &[u8] = &[75, 90, 101, 110];
        let kzen_label = BigInt::from_bytes(KZen);
        test_helper(&kzen_label, 4);
    }

    #[test]
    pub fn test_zk_eval_8() {
        let KZen: &[u8] = &[75, 90, 101, 110];
        let kzen_label = BigInt::from_bytes(KZen);
        test_helper(&kzen_label, 8);
    }

    #[test]
    pub fn test_zk_eval_16() {
        let KZen: &[u8] = &[75, 90, 101, 110];
        let kzen_label = BigInt::from_bytes(KZen);
        test_helper(&kzen_label, 16);
    }
    
    #[test]
    pub fn test_zk_eval_32() {
        let KZen: &[u8] = &[75, 90, 101, 110];
        let kzen_label = BigInt::from_bytes(KZen);
        test_helper(&kzen_label, 32);
    }

    #[test]
    pub fn test_zk_eval_64() {
        let KZen: &[u8] = &[75, 90, 101, 110];
        let kzen_label = BigInt::from_bytes(KZen);
        test_helper(&kzen_label, 64);
    }
}
