#![allow(non_snake_case)]

use generic_array::typenum::Unsigned;
use curv::{elliptic::curves::{secp256_k1::Secp256k1, Curve, ECPoint, Point, Scalar}, BigInt};
use curv::arithmetic::traits::*;
use generic_array::GenericArray;
use merlin::Transcript;
use crate::{transcript::TranscriptProtocol, Errors::SigmaDlProofError};
use crate::Errors;

pub fn generate_random_point(bytes: &[u8]) -> Point<Secp256k1> {
    let compressed_point_len =
        <<Secp256k1 as Curve>::Point as ECPoint>::CompressedPointLength::USIZE;
    let truncated = if bytes.len() > compressed_point_len - 1 {
        &bytes[0..compressed_point_len - 1]
    } else {
        bytes
    };
    let mut buffer = GenericArray::<
        u8,
        <<Secp256k1 as Curve>::Point as ECPoint>::CompressedPointLength,
    >::default();
    buffer.as_mut_slice()[0] = 0x2;
    buffer.as_mut_slice()[1..1 + truncated.len()].copy_from_slice(truncated);
    if let Ok(point) = Point::from_bytes(buffer.as_slice()) {
        return point;
    }

    let bn = BigInt::from_bytes(bytes);
    let two = BigInt::from(2);
    let bn_times_two = BigInt::mod_mul(&bn, &two, Scalar::<Secp256k1>::group_order());
    let bytes = BigInt::to_bytes(&bn_times_two);
    generate_random_point(&bytes)
}

#[derive(Clone)]
pub struct SigmaDlProof {
    pub(super) R_vec: Vec<Point<Secp256k1>>,
    pub(super) z_vec: Vec<Scalar<Secp256k1>>,
}

impl SigmaDlProof {
    pub fn prove(
        transcript: &mut Transcript,
        x_vec: &[Scalar<Secp256k1>],
        y_vec: &[Point<Secp256k1>],
        h_vec: &[Point<Secp256k1>],
        n: usize
    ) -> SigmaDlProof {
        assert_eq!(x_vec.len(), n);
        assert_eq!(h_vec.len(), n);

        transcript.sigma_dl_domain_sep(n as u64);
        transcript.append_points_array(b"y_vec", y_vec);

        let x_vec = x_vec.to_vec();
        let h_vec = h_vec.to_vec();

        let mut r_vec: Vec<Scalar<Secp256k1>> = Vec::with_capacity(n);
        for _j in 0..n {
            r_vec.push(Scalar::<Secp256k1>::random());
        }

        let R_vec = (0..n)
            .map(|i| &h_vec[i].clone() * r_vec[i].clone())
            .collect::<Vec<Point<Secp256k1>>>();

        transcript.append_points_array(b"R_vec", &R_vec);
        let challenge_e: Scalar<Secp256k1> = transcript.challenge_scalar(b"e");

        let z_vec = (0..n)
            .map(|i| r_vec[i].clone() + challenge_e.clone() * x_vec[i].clone())
            .collect::<Vec<Scalar<Secp256k1>>>();

        SigmaDlProof {
            R_vec,
            z_vec
        }
    }

    pub fn verify(
        &self,
        transcript: &mut Transcript,
        y_vec: &[Point<Secp256k1>],
        h_vec: &[Point<Secp256k1>],
        n: usize
    ) -> Result<(), Errors> {
        assert_eq!(h_vec.len(), n);
        assert_eq!(self.R_vec.len(), n);
        assert_eq!(self.z_vec.len(), n);

        transcript.sigma_dl_domain_sep(n as u64);
        transcript.append_points_array(b"y_vec", y_vec);

        let y_vec = y_vec.to_vec();
        let h_vec = h_vec.to_vec();

        transcript.append_points_array(b"R_vec", &self.R_vec);
        let challenge_e: Scalar<Secp256k1> = transcript.challenge_scalar(b"e");

        let mut flag = true;
        for j in 0..n {
            if &h_vec[j] * &self.z_vec[j] != &self.R_vec[j] + &challenge_e * &y_vec[j] {
                flag = false;
                break;
            }
        }
        
        if flag {
            Ok(())
        } else {
            Err(SigmaDlProofError)
        }
    } 
}

mod test {
    use curv::{arithmetic::Converter, cryptographic_primitives::hashing::DigestExt, elliptic::curves::{Point, Scalar, Secp256k1}, BigInt};
    use merlin::Transcript;
    use sha2::{Digest, Sha512};

    use crate::sigma_dl::{generate_random_point, SigmaDlProof};

    pub fn test_helper(seed: &BigInt, n: usize) {
        let h_vec = (0..n)
            .map(|i| {
                let kzen_label_i = BigInt::from(i as u32) + seed;
                let hash_i = Sha512::new().chain_bigint(&kzen_label_i).result_bigint();
                generate_random_point(&Converter::to_bytes(&hash_i))
            })
            .collect::<Vec<Point<Secp256k1>>>();

        let x_vec: Vec<Scalar<Secp256k1>> = (0..n)
            .map(|_| {
                let rand = Scalar::<Secp256k1>::random();
                rand
            })
            .collect();

        let y_vec = (0..n)
            .map(|i| &h_vec[i].clone() * x_vec[i].clone())
            .collect::<Vec<Point<Secp256k1>>>();
        
        let mut verifier = Transcript::new(b"sigmatest");
        let proof = SigmaDlProof::prove(&mut verifier, &x_vec, &y_vec, &h_vec, n);
        let mut verifier = Transcript::new(b"sigmatest");
        let result = proof.verify(&mut verifier, &y_vec, &h_vec, n);
        assert!(result.is_ok());
    }

    #[test]
    pub fn test_sigma_dl_4() {
        let KZen: &[u8] = &[75, 90, 101, 110];
        let kzen_label = BigInt::from_bytes(KZen);
        test_helper(&kzen_label, 4);
    }

    #[test]
    pub fn test_sigma_dl_8() {
        let KZen: &[u8] = &[75, 90, 101, 110];
        let kzen_label = BigInt::from_bytes(KZen);
        test_helper(&kzen_label, 8);
    }

    #[test]
    pub fn test_sigma_dl_16() {
        let KZen: &[u8] = &[75, 90, 101, 110];
        let kzen_label = BigInt::from_bytes(KZen);
        test_helper(&kzen_label, 16);
    }
    
    #[test]
    pub fn test_sigma_dl_32() {
        let KZen: &[u8] = &[75, 90, 101, 110];
        let kzen_label = BigInt::from_bytes(KZen);
        test_helper(&kzen_label, 32);
    }

    #[test]
    pub fn test_sigma_dl_64() {
        let KZen: &[u8] = &[75, 90, 101, 110];
        let kzen_label = BigInt::from_bytes(KZen);
        test_helper(&kzen_label, 64);
    }
}