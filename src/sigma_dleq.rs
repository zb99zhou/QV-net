#![allow(non_snake_case)]

/// This file implements an AND composition for proving the discrete logarithm equality.
/// Given h,g,Y,y,B \in \mathbb{G}^n, it allows to prove
/// the knowledge of v,x \in \mathbb{F}^n, such that 
/// for all i\in[n], y_i = h_i^{(x_i)} \land B_i = g_i^{(v_i)} \cdot Y_i^{(x_i)}
/// 
use curv::elliptic::curves::{Point, Scalar, Secp256k1};
use merlin::Transcript;
use crate::{transcript::TranscriptProtocol, Errors::SigmaDleqProofError};
use crate::Errors;

#[derive(Clone)]
pub struct SigmaDleqProof {
    pub(super) R_vec: Vec<Point<Secp256k1>>,
    pub(super) S_vec: Vec<Point<Secp256k1>>,
    pub(super) a_vec: Vec<Scalar<Secp256k1>>,
    pub(super) b_vec: Vec<Scalar<Secp256k1>>
}

impl SigmaDleqProof {
    pub fn prove(
        transcript: &mut Transcript,
        v_vec: &[Scalar<Secp256k1>],
        x_vec: &[Scalar<Secp256k1>],
        y_vec: &[Point<Secp256k1>],
        B_vec: &[Point<Secp256k1>],
        h_vec: &[Point<Secp256k1>],
        g_vec: &[Point<Secp256k1>],
        Y_vec: &[Point<Secp256k1>],
        n: usize
    ) ->SigmaDleqProof {
        assert_eq!(v_vec.len(), n);
        assert_eq!(x_vec.len(), n);
        assert_eq!(h_vec.len(), n);
        assert_eq!(g_vec.len(), n);
        assert_eq!(Y_vec.len(), n);

        transcript.sigma_dleq_domain_sep(n as u64);

        let v_vec = v_vec.to_vec();
        let x_vec = x_vec.to_vec();
        let h_vec = h_vec.to_vec();
        let g_vec = g_vec.to_vec();
        let Y_vec = Y_vec.to_vec();

        transcript.append_points_array(b"h_vec", &h_vec);
        transcript.append_points_array(b"g_vec", &g_vec);
        transcript.append_points_array(b"Y_vec", &Y_vec);
        transcript.append_points_array(b"y_vec", y_vec);
        transcript.append_points_array(b"B_vec", B_vec);

        let mut r_vec: Vec<Scalar<Secp256k1>> = Vec::with_capacity(n);
        let mut s_vec: Vec<Scalar<Secp256k1>> = Vec::with_capacity(n);
        for _j in 0..n {
            r_vec.push(Scalar::<Secp256k1>::random());
            s_vec.push(Scalar::<Secp256k1>::random());
        }

        let R_vec = (0..n)
            .map(|i| &h_vec[i].clone() * r_vec[i].clone())
            .collect::<Vec<Point<Secp256k1>>>();

        let S_vec = (0..n)
            .map(|i| &g_vec[i].clone() * s_vec[i].clone() + &Y_vec[i].clone() * r_vec[i].clone())
            .collect::<Vec<Point<Secp256k1>>>();

        transcript.append_points_array(b"R_vec", &R_vec);
        transcript.append_points_array(b"S_vec", &S_vec);
        let challenge_e: Scalar<Secp256k1> = transcript.challenge_scalar(b"e");

        let a_vec = (0..n)
            .map(|i| r_vec[i].clone() + challenge_e.clone() * x_vec[i].clone())
            .collect::<Vec<Scalar<Secp256k1>>>();
        let b_vec = (0..n)
            .map(|i| s_vec[i].clone() + challenge_e.clone() * v_vec[i].clone())
            .collect::<Vec<Scalar<Secp256k1>>>();
        
        SigmaDleqProof {
            R_vec,
            S_vec,
            a_vec,
            b_vec
        }
    }

    pub fn verify(
        &self,
        transcript: &mut Transcript,
        y_vec: &[Point<Secp256k1>],
        B_vec: &[Point<Secp256k1>],
        h_vec: &[Point<Secp256k1>],
        g_vec: &[Point<Secp256k1>],
        Y_vec: &[Point<Secp256k1>],
        n: usize
    ) -> Result<(), Errors> {
        assert_eq!(h_vec.len(), n);
        assert_eq!(g_vec.len(), n);
        assert_eq!(Y_vec.len(), n);
        assert_eq!(self.R_vec.len(), n);
        assert_eq!(self.S_vec.len(), n);
        assert_eq!(self.a_vec.len(), n);
        assert_eq!(self.b_vec.len(), n);

        transcript.sigma_dleq_domain_sep(n as u64);

        let y_vec = y_vec.to_vec();
        let B_vec = B_vec.to_vec();
        let h_vec = h_vec.to_vec();
        let g_vec = g_vec.to_vec();
        let Y_vec = Y_vec.to_vec();

        transcript.append_points_array(b"h_vec", &h_vec);
        transcript.append_points_array(b"g_vec", &g_vec);
        transcript.append_points_array(b"Y_vec", &Y_vec);
        transcript.append_points_array(b"y_vec", &y_vec);
        transcript.append_points_array(b"B_vec", &B_vec);

        transcript.append_points_array(b"R_vec", &self.R_vec);
        transcript.append_points_array(b"S_vec", &self.S_vec);
        let challenge_e: Scalar<Secp256k1> = transcript.challenge_scalar(b"e");

        let mut flag = true;
        for j in 0..n {
            if &h_vec[j].clone() * self.a_vec[j].clone() != self.R_vec[j].clone() + challenge_e.clone() * y_vec[j].clone() {
                flag = false;
                break;
            }
        }
        for j in 0..n {
            if &g_vec[j].clone() * self.b_vec[j].clone() + &Y_vec[j].clone() * self.a_vec[j].clone() != self.S_vec[j].clone() + challenge_e.clone() * B_vec[j].clone() {
                flag = false;
                break;
            }
        }

        if flag {
            Ok(())
        } else {
            Err(SigmaDleqProofError)
        }
    }
}

mod test {
    use curv::{arithmetic::Converter, cryptographic_primitives::hashing::DigestExt, elliptic::curves::{Point, Scalar, Secp256k1}, BigInt};
    use merlin::Transcript;
    use sha2::{Digest, Sha512};

    use crate::sigma_dl::generate_random_point;

    use super::SigmaDleqProof;

    pub fn test_helper(seed: &BigInt, n: usize) {
        let h_vec = (0..n)
            .map(|i| {
                let kzen_label_i = BigInt::from(i as u32) + seed;
                let hash_i = Sha512::new().chain_bigint(&kzen_label_i).result_bigint();
                generate_random_point(&Converter::to_bytes(&hash_i))
            })
            .collect::<Vec<Point<Secp256k1>>>();

        // can run in parallel to hi:
        let g_vec = (0..n)
            .map(|i| {
                let kzen_label_j = BigInt::from(n as u32) + BigInt::from(i as u32) + seed;
                let hash_j = Sha512::new().chain_bigint(&kzen_label_j).result_bigint();
                generate_random_point(&Converter::to_bytes(&hash_j))
            })
            .collect::<Vec<Point<Secp256k1>>>();

        // can run in parallel to hi:
        let Y_vec = (0..n)
            .map(|i| {
                let kzen_label_j = BigInt::from(n as u32) + BigInt::from(n as u32) + BigInt::from(i as u32) + seed;
                let hash_j = Sha512::new().chain_bigint(&kzen_label_j).result_bigint();
                generate_random_point(&Converter::to_bytes(&hash_j))
            })
            .collect::<Vec<Point<Secp256k1>>>();

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

        let y_vec = (0..n)
            .map(|i| &h_vec[i].clone() * x_vec[i].clone())
            .collect::<Vec<Point<Secp256k1>>>();

        let B_vec = (0..n)
            .map(|i| &g_vec[i].clone() * v_vec[i].clone() + &Y_vec[i].clone() * x_vec[i].clone())
            .collect::<Vec<Point<Secp256k1>>>();

        let mut verifier = Transcript::new(b"sigmatest");
        let proof = SigmaDleqProof::prove(&mut verifier, &v_vec, &x_vec, &y_vec, &B_vec, &h_vec, &g_vec, &Y_vec, n);
        let mut verifier = Transcript::new(b"sigmatest");
        let result = proof.verify(&mut verifier, &y_vec, &B_vec, &h_vec, &g_vec, &Y_vec, n);
        assert!(result.is_ok());
    }


    #[test]
    pub fn test_sigma_dleq_64() {
        let KZen: &[u8] = &[75, 90, 101, 110];
        let kzen_label = BigInt::from_bytes(KZen);
        test_helper(&kzen_label, 64);
    }
}