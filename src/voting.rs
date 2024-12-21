#![allow(non_snake_case)]
use crate::Errors::{self, VotingError};
use ark_ff::{One, Zero};
use curv::{arithmetic::{Converter, Roots}, elliptic::curves::{Point, Scalar, Secp256k1}, BigInt};
use merlin::Transcript;
use std::{collections::HashMap, time::Instant};
use VarRange::proofs::varrange::VarRange;

use crate::{sigma_dl::SigmaDlProof, sigma_dleq::SigmaDleqProof, sum_square::ZkSumSquareArg};

#[derive(Clone)]
pub struct Voter {
    VoterID: usize,
    x_vec: Vec<Scalar<Secp256k1>>,
    y_vec: Vec<Point<Secp256k1>>,
    Y_vec: Vec<Point<Secp256k1>>,
    proof_dl: SigmaDlProof
}

pub struct Board {
    pp: PublicParam,
    y_vec: Vec<Vec<Point<Secp256k1>>>,
    ballot_proof: Vec<BallotWithProof>,
    Y_vec_with_global: Vec<Vec<Point<Secp256k1>>>
}

pub struct PublicParam {
    g_vec: Vec<Point<Secp256k1>>,
    h_vec: Vec<Point<Secp256k1>>,
    g: Point<Secp256k1>,
    h: Point<Secp256k1>,
    nv: usize,
    nc: usize
}

#[derive(Clone)]
pub struct BallotWithProof {
    B_vec: Vec<Point<Secp256k1>>,
    B: Point<Secp256k1>,
    proof_dl: SigmaDlProof,
    proof_eq: SigmaDleqProof,
    proof_ss: ZkSumSquareArg,
    proof_ar: VarRange
}

// Shanks' baby-step giant-step algorithm, calculate ind_g{B}
pub fn tally_helper(
    g: &Point<Secp256k1>, 
    B: &Point<Secp256k1>, 
    bound: &BigInt, 
    flag: &mut bool
) -> Scalar<Secp256k1> {
    let n_bn = bound.sqrt() + BigInt::one();
    let n = &Scalar::<Secp256k1>::from_bigint(&n_bn);
    let mut counter = BigInt::one();
    let one = &Scalar::<Secp256k1>::from(1);
    let one_bn = &BigInt::from(1);
    let mut res = &Scalar::<Secp256k1>::zero();

    let mut p = Scalar::<Secp256k1>::from(1);
    let mut q = Scalar::<Secp256k1>::from(0);
    let mut giants: HashMap<Vec<u8>, Scalar<Secp256k1>> = HashMap::new();
    let mut giant: Point<Secp256k1>;
    let mut baby: Point<Secp256k1>;
    
    while counter <= n_bn {
        giant = g * (&p * n);
        giants.insert(giant.x_coord().unwrap().to_bytes(), p.clone());
        p = p + one;
        counter += one_bn;
    }

    counter = BigInt::zero();
    while counter <= n_bn {
        baby = B + g * &q;
        if baby == Point::<Secp256k1>::zero() {
            *flag = true;
            return -q;
        }
        if let Some(giant_p) = giants.get(&baby.x_coord().unwrap().to_bytes()) {
            res = giant_p;
            *flag = true;
            break;
        }
        q = q + one;
        counter += one_bn;
    }

    if g * (res * n - &q) == *B {
        res * n - q
    } else {
        - res * n - q
    }
}

impl PublicParam {
    pub fn new(
        g_vec: &[Point<Secp256k1>],
        h_vec: &[Point<Secp256k1>],
        g: &Point<Secp256k1>,
        h: &Point<Secp256k1>,
        nv: usize,
        nc: usize
    ) -> PublicParam {
        assert_eq!(g_vec.len(), nc);
        assert_eq!(h_vec.len(), nc);
        
        PublicParam {
            g_vec: g_vec.to_vec(),
            h_vec: h_vec.to_vec(),
            g: g.clone(),
            h: h.clone(),
            nv,
            nc
        }
    }
    
}

impl Voter {
    pub fn gen(
        VoterID: usize,
        h_vec: &[Point<Secp256k1>],
        nv: usize,
        nc: usize,
        BulletinBoard: &mut Board
    ) -> Voter {
        assert!(VoterID < nv);
        assert_eq!(h_vec.len(), nc);

        let h_vec = h_vec.to_vec();

        let mut x_vec: Vec<Scalar<Secp256k1>> = Vec::new();
        for _j in 0..nc {
            x_vec.push(Scalar::<Secp256k1>::random());
        }
        
        let start = Instant::now();
        let y_vec = (0..nc)
            .map(|j| &h_vec[j] * &x_vec[j])
            .collect::<Vec<Point<Secp256k1>>>();
        let elapsed = start.elapsed();
        // println!("Time elapsed in generate yij: {:?}", elapsed);
        
        let mut transcript = Transcript::new(b"Proof");
        let start = Instant::now();
        let proof_dl = SigmaDlProof::prove(
            &mut transcript, 
            &x_vec, 
            &y_vec,
            &h_vec, 
            nc
        );
        let elapsed = start.elapsed();
        // println!("Time elapsed in generate pi_i^dl: {:?}", elapsed);

        BulletinBoard.y_vec.push(y_vec.clone());

        Voter {
            VoterID,
            x_vec,
            y_vec,
            Y_vec: vec![],
            proof_dl
        }
    }

    pub fn precompute_for_prove(
        &mut self,
        VoterID: usize,
        BulletinBoard: &Board,
    ) {
        let mut Y_vec: Vec<Point<Secp256k1>> = vec![Point::<Secp256k1>::zero(); BulletinBoard.pp.nc];
        let start = Instant::now();
        for j in 0..BulletinBoard.pp.nc {
            for k in 0..VoterID {
                Y_vec[j] = &Y_vec[j] + &BulletinBoard.y_vec[k][j];
            }
            for k in VoterID + 1..BulletinBoard.pp.nv {
                Y_vec[j] = &Y_vec[j] - &BulletinBoard.y_vec[k][j];
            }
        }
        let elapsed = start.elapsed();
        println!("Time elapsed in generate Yij: {:?}", elapsed);
        self.Y_vec = Y_vec;
    }

    pub fn precompute_for_verify(
        &mut self,
        BulletinBoard: &mut Board,
    ) {
        let mut Y_vec_with_global = Vec::with_capacity(BulletinBoard.pp.nv * BulletinBoard.pp.nc);
        for i in 0..BulletinBoard.pp.nv {
            let mut Yi_vec: Vec<Point<Secp256k1>> = vec![Point::<Secp256k1>::zero(); BulletinBoard.pp.nc];
            for j in 0..BulletinBoard.pp.nc {
                for k in 0..i {
                    Yi_vec[j] = &Yi_vec[j] + &BulletinBoard.y_vec[k][j];
                }
                for k in i + 1..BulletinBoard.pp.nv {
                    Yi_vec[j] = &Yi_vec[j] - &BulletinBoard.y_vec[k][j];
                }
            }
            Y_vec_with_global.push(Yi_vec);
        }
        BulletinBoard.Y_vec_with_global = Y_vec_with_global;
    }

    pub fn vote(
        &self,
        v_vec: &[Scalar<Secp256k1>],
        token: Scalar<Secp256k1>,
        BulletinBoard: &mut Board,
        seed: &BigInt
    ) -> BallotWithProof {
        assert_eq!(v_vec.len(), BulletinBoard.pp.nc);
        assert_eq!(BulletinBoard.y_vec.len(), BulletinBoard.pp.nv);
        assert_eq!(BulletinBoard.y_vec[0].len(), BulletinBoard.pp.nc);
        assert_eq!(BulletinBoard.pp.g_vec.len(), BulletinBoard.pp.nc);
        assert_eq!(BulletinBoard.pp.h_vec.len(), BulletinBoard.pp.nc);

        let start = Instant::now();
        let B_vec = (0..BulletinBoard.pp.nc)
            .map(|j| &BulletinBoard.pp.g_vec[j] * &v_vec[j] + &self.Y_vec[j] * &self.x_vec[j])
            .collect::<Vec<Point<Secp256k1>>>();
        let elapsed = start.elapsed();
        // println!("Time elapsed in generate B_vec: {:?}", elapsed);

        let mut u = Scalar::<Secp256k1>::zero();
        for j in 0..BulletinBoard.pp.nc {
            u = u + &v_vec[j] * &v_vec[j];
        }
        let token_bigint = token.to_bigint();
        assert!(u.to_bigint() <= token_bigint);

        let r = Scalar::<Secp256k1>::random();
        let start = Instant::now();
        let B = &BulletinBoard.pp.g * &u + &BulletinBoard.pp.h * &r;
        let elapsed = start.elapsed();
        // println!("Time elapsed in generate B: {:?}", elapsed);

        let mut transcript = Transcript::new(b"Proof");
        let start = Instant::now();
        let proof_eq = SigmaDleqProof::prove(
            &mut transcript,
            v_vec, 
            &self.x_vec, 
            &self.y_vec,
            &B_vec,
            &BulletinBoard.pp.h_vec, 
            &BulletinBoard.pp.g_vec, 
            &self.Y_vec,
            BulletinBoard.pp.nc
        );
        let elapsed = start.elapsed();
        // println!("Time elapsed in generate pi_i^eq: {:?}", elapsed);

        let mut transcript = Transcript::new(b"Proof");
        let start = Instant::now();
        let proof_ss = ZkSumSquareArg::prove(
            &mut transcript, 
            &BulletinBoard.pp.g_vec, 
            &self.Y_vec,
            &BulletinBoard.pp.g, 
            &BulletinBoard.pp.h, 
            &B_vec.iter().sum(),
            &B, 
            v_vec, 
            &self.x_vec, 
            &u, 
            &r, 
            BulletinBoard.pp.nc, 
            seed
        );
        let elapsed = start.elapsed();
        // println!("Time elapsed in generate pi_i^ss: {:?}", elapsed);

        let mut g_extend = BulletinBoard.pp.g_vec.clone();
        g_extend.push(BulletinBoard.pp.g.clone());
        let mut h_extend = self.Y_vec.clone();
        h_extend.push(BulletinBoard.pp.h.clone());
        let mut v_vec_extend = v_vec.to_vec();
        let mut token_new: Scalar<Secp256k1> = token.clone();
        let mut u_new: Scalar<Secp256k1> = u.clone();
        let mut B_new: Point<Secp256k1> = B.clone();
        let mut r_new: Scalar<Secp256k1> = r.clone();
        if token == Scalar::<Secp256k1>::from(1) {
            token_new = &token * Scalar::<Secp256k1>::from(2);
            u_new = &u * Scalar::<Secp256k1>::from(2);
            B_new = &B * Scalar::<Secp256k1>::from(2);
            r_new = &r * Scalar::<Secp256k1>::from(2);
        }
        v_vec_extend.push(u_new.clone());
        let mut x_vec_extend = self.x_vec.clone();
        x_vec_extend.push(r_new.clone());
        let mut B_vec_extend = B_vec.clone();
        B_vec_extend.push(B_new.clone());
        let s = token_bigint.sqrt();
        for j in 0..v_vec_extend.len()-1 {
            assert!((&v_vec_extend[j] + Scalar::<Secp256k1>::from_bigint(&s)).to_bigint() <= &s + &s);
        }
        assert_eq!(&BulletinBoard.pp.g * u_new.clone() + &BulletinBoard.pp.h * r_new.clone(), B_new);
        assert!(u_new.to_bigint() <= token_new.to_bigint());

        let mut transcript = Transcript::new(b"Proof");
        let start = Instant::now();
        let proof_ar = VarRange::range_prove(
            &mut transcript, 
            &g_extend,
            &h_extend,
            v_vec_extend, 
            x_vec_extend, 
            Scalar::<Secp256k1>::from_bigint(&s), 
            token_new, 
            BulletinBoard.pp.nc, 
            &B_vec_extend, 
            seed
        );
        let elapsed = start.elapsed();
        // println!("Time elapsed in generate pi_i^ar: {:?}", elapsed);

        let ballot_proof = BallotWithProof {
            B_vec,
            B,
            proof_dl: self.proof_dl.clone(),
            proof_eq,
            proof_ss,
            proof_ar
        };
        BulletinBoard.ballot_proof.push(ballot_proof.clone());

        ballot_proof
    }
}

impl Board {
    pub fn setup(
        g_vec: &[Point<Secp256k1>],
        h_vec: &[Point<Secp256k1>],
        g: &Point<Secp256k1>,
        h: &Point<Secp256k1>,
        nv: usize,
        nc: usize
    ) -> Board {
        let y_vec: Vec<Vec<Point<Secp256k1>>> = Vec::with_capacity(nv);
        let pp = PublicParam::new(g_vec, h_vec, g, h, nv, nc);
        let ballot_proof: Vec<BallotWithProof> = Vec::with_capacity(nv);

        Board { pp, y_vec, ballot_proof, Y_vec_with_global: vec![] }
    }

    pub fn precompute_for_verify(&mut self) {
        let mut Y_vec_with_global = Vec::with_capacity(self.pp.nv * self.pp.nc);
        for i in 0..self.pp.nv {
            let mut Yi_vec: Vec<Point<Secp256k1>> = vec![Point::<Secp256k1>::zero(); self.pp.nc];
            for j in 0..self.pp.nc {
                for k in 0..i {
                    Yi_vec[j] = &Yi_vec[j] + &self.y_vec[k][j];
                }
                for k in i + 1..self.pp.nv {
                    Yi_vec[j] = &Yi_vec[j] - &self.y_vec[k][j];
                }
            }
            Y_vec_with_global.push(Yi_vec);
        }
        self.Y_vec_with_global = Y_vec_with_global;
    }

    pub fn batch_verify(
        &self,
        token: Vec<Scalar<Secp256k1>>,
        seed: &BigInt
    ) {
        assert_eq!(self.y_vec.len(), self.pp.nv);
        assert_eq!(self.y_vec[0].len(), self.pp.nc);
        assert_eq!(self.pp.g_vec.len(), self.pp.nc);
        assert_eq!(self.pp.h_vec.len(), self.pp.nc);
        assert_eq!(self.ballot_proof.len(), self.pp.nv);
        assert_eq!(token.len(), self.pp.nv);

        for voteId in 0..self.pp.nv {
            self.verify_single_ballot(&token[voteId], seed, voteId)
        }
    }
    
    pub fn verify_single_ballot(
        &self,
        token: &Scalar<Secp256k1>,
        seed: &BigInt,
        voteId: usize,
    ) {
        // verify the proofs
        let mut transcript = Transcript::new(b"Proof");
        let start = Instant::now();
        let res_dl = self.ballot_proof[voteId].proof_dl.verify(
            &mut transcript,
            &self.y_vec[voteId],
            &self.pp.h_vec,
            self.pp.nc
        );
        let elapsed = start.elapsed();
        println!("Time elapsed in verify pi_i^dl: {:?}", elapsed);
        assert!(res_dl.is_ok());

        let mut transcript = Transcript::new(b"Proof");
        let start = Instant::now();
        let res_dleq = self.ballot_proof[voteId].proof_eq.verify(
            &mut transcript,
            &self.y_vec[voteId],
            &self.ballot_proof[voteId].B_vec,
            &self.pp.h_vec,
            &self.pp.g_vec,
            &self.Y_vec_with_global[voteId],
            self.pp.nc
        );
        let elapsed = start.elapsed();
        println!("Time elapsed in verify pi_i^dleq: {:?}", elapsed);
        assert!(res_dleq.is_ok());

        let mut transcript = Transcript::new(b"Proof");
        let start = Instant::now();
        let res_ss = self.ballot_proof[voteId].proof_ss.verify(
            &mut transcript,
            &self.pp.g_vec,
            &self.Y_vec_with_global[voteId],
            &self.pp.g,
            &self.pp.h,
            &self.ballot_proof[voteId].B_vec.iter().sum(),
            &self.ballot_proof[voteId].B,
            self.pp.nc,
            seed
        );
        let elapsed = start.elapsed();
        println!("Time elapsed in verify pi_i^ss: {:?}", elapsed);
        assert!(res_ss.is_ok());

        let mut gi = self.pp.g_vec.clone();
        gi.push(self.pp.g.clone());
        let mut hi = self.Y_vec_with_global[voteId].clone();
        hi.push(self.pp.h.clone());
        let mut B_vec = self.ballot_proof[voteId].B_vec.clone();
        let mut token_new: Scalar<Secp256k1> = token.clone();
        let mut B_new: Point<Secp256k1> = self.ballot_proof[voteId].B.clone();
        if token == &Scalar::<Secp256k1>::from(1) {
            token_new = token * Scalar::<Secp256k1>::from(2);
            B_new = &self.ballot_proof[voteId].B * Scalar::<Secp256k1>::from(2);
        }
        B_vec.push(B_new);
        let si = BigInt::sqrt(&token.to_bigint());

        let mut transcript = Transcript::new(b"Proof");
        let start = Instant::now();
        let res_ar = self.ballot_proof[voteId].proof_ar.range_verify(
            &mut transcript,
            &gi,
            &hi,
            &B_vec,
            Scalar::<Secp256k1>::from_bigint(&si),
            token_new,
            self.pp.nc,
            seed
        );
        let elapsed = start.elapsed();
        println!("Time elapsed in verify pi_i^ar: {:?}", elapsed);

        assert!(res_ar.is_ok());
    }

    pub fn tally(
        &self,
        bound: &BigInt
    ) -> Result<Vec<Scalar<Secp256k1>>, Errors> {
        assert_eq!(self.y_vec.len(), self.pp.nv);
        assert_eq!(self.y_vec[0].len(), self.pp.nc);
        assert_eq!(self.pp.g_vec.len(), self.pp.nc);
        assert_eq!(self.pp.h_vec.len(), self.pp.nc);
        assert_eq!(self.ballot_proof.len(), self.pp.nv);

        let mut ballots: Vec<Scalar<Secp256k1>> = Vec::with_capacity(self.pp.nc);
        for j in 0..self.pp.nc {
            let mut Bj = Point::<Secp256k1>::zero();
            for i in 0..self.pp.nv {
                Bj = Bj + &self.ballot_proof[i].B_vec[j];
            }
            let mut flag = false;
            let mut Bj_scalar = tally_helper(&self.pp.g_vec[j], &Bj, bound, &mut flag);
            if !flag {
                Bj_scalar = -tally_helper(&(-self.pp.g_vec[j].clone()), &Bj, bound, &mut flag);
            }
            if flag {
                ballots.push(Bj_scalar);
            } else {
                return Err(VotingError)
            }
        }
        Ok(ballots)
    }
}

mod test {
    use std::time::Instant;

    use ark_ff::Zero;
    use curv::{arithmetic::{Converter, Modulo, Roots}, cryptographic_primitives::hashing::DigestExt, elliptic::curves::{Point, Scalar, Secp256k1}, BigInt};
    use sha2::{Digest, Sha512};

    use crate::sigma_dl::generate_random_point;

    use super::{tally_helper, Board, Voter};

    pub fn test_helper_with_verify(seed: &BigInt, nv: usize, nc: usize, token_num: &BigInt) {
        let h_vec = (0..nc)
            .map(|i| {
                let kzen_label_i = BigInt::from(i as u32) + seed;
                let hash_i = Sha512::new().chain_bigint(&kzen_label_i).result_bigint();
                generate_random_point(&Converter::to_bytes(&hash_i))
            })
            .collect::<Vec<Point<Secp256k1>>>();

        // can run in parallel to hi:
        let g_vec = (0..nc)
            .map(|i| {
                let kzen_label_j = BigInt::from(nc as u32) + BigInt::from(i as u32) + seed;
                let hash_j = Sha512::new().chain_bigint(&kzen_label_j).result_bigint();
                generate_random_point(&Converter::to_bytes(&hash_j))
            })
            .collect::<Vec<Point<Secp256k1>>>();

        let kzen_label = BigInt::from(nc as u32) + BigInt::from(nc as u32) + seed;
        let hash = Sha512::new().chain_bigint(&kzen_label).result_bigint();
        let h = generate_random_point(&Converter::to_bytes(&hash));

        // Voting init setup
        let kzen_label =  BigInt::from(1_u32) + BigInt::from(nc as u32) + BigInt::from(nc as u32) + seed;
        let hash = Sha512::new().chain_bigint(&kzen_label).result_bigint();
        let g = generate_random_point(&Converter::to_bytes(&hash));

        let mut Voter_vec: Vec<Voter> = Vec::new();
        let mut board = Board::setup(&g_vec, &h_vec, &g, &h, nv, nc);
        for i in 0..nv {
            let Vi = Voter::gen(i, &h_vec, nv, nc, &mut board);
            Voter_vec.push(Vi);
        }

        // Voting precompute
        for i in 0..nv {
            Voter_vec[i].precompute_for_prove(i, &board);
        }
        board.precompute_for_verify();

        // Voting prove
        let order = Scalar::<Secp256k1>::group_order();
        let mut v_vec: Vec<Vec<Scalar<Secp256k1>>> = Vec::new();
        let mut tokens: Vec<Scalar<Secp256k1>> = Vec::new();
        let bound = BigInt::sqrt(&(token_num / BigInt::from(nc as u64)));
        for i in 0..nv {
            let mut vi_vec: Vec<Scalar<Secp256k1>> = Vec::new();
            let mut square_sum = BigInt::zero();
            for _ in 0..nc {
                // let sign = rng.gen_bool(0.5) as i32;
                // let sign = if sign == 1 { 1 } else { -1 };
                let sign = 1;
                let vij = BigInt::mod_mul(&bound, &(order + BigInt::from(sign)), order);
                vi_vec.push(Scalar::<Secp256k1>::from_bigint(&vij));
                square_sum += &vij * &vij;
            }
            tokens.push(Scalar::<Secp256k1>::from_bigint(&token_num));
            Voter_vec[i].vote(&vi_vec, tokens[i].clone(), &mut board, &(BigInt::from((nc*2+1) as u32) + seed));
            v_vec.push(vi_vec);
        }

        // Voting verify
        board.batch_verify(tokens, &(BigInt::from((nc*2+1) as u32) + seed));
        let start = Instant::now();
        let res = board.tally(&(&bound * BigInt::from(nv as u64))).unwrap();
        let elapsed = start.elapsed();
        println!("Time elapsed in tally: {:?}", elapsed);

        // test correctness
        for j in 0..nc {
            let mut Bj = Scalar::<Secp256k1>::zero();
            for i in 0..nv {
                Bj = Bj + &v_vec[i][j];
            }
            assert_eq!(Bj.to_bigint(), res[j].to_bigint());
        }
    }

    #[test]
    pub fn test_voting_with_verify() {
        let KZen: &[u8] = &[75, 90, 101, 110];
        let kzen_label = BigInt::from_bytes(KZen);
        
        for _i in 0..15 {
           test_helper_with_verify(&kzen_label, 857, 4, &BigInt::from(1_600_000));
        }
        // test_helper_with_verify(&kzen_label, 2800, 15, &BigInt::from(1_200_000));
    }

    #[test]
    pub fn test_shanks() {
        let KZen: &[u8] = &[75, 90, 101, 110];
        let kzen_label = BigInt::from_bytes(KZen);
        let hash = Sha512::new().chain_bigint(&kzen_label).result_bigint();
        let g = generate_random_point(&Converter::to_bytes(&hash));
        
        for exp in 20..=60 {
            let scalar_value = 2_u64.pow(exp);
            let scalar = Scalar::<Secp256k1>::from(scalar_value);

            let point = &g * scalar.clone();
    
            let mut flag = true;
            let start = Instant::now();
            let Bj_scalar = tally_helper(&g, &point, &BigInt::from(scalar_value), &mut flag);
            let elapsed = start.elapsed();
            println!(
                "Time elapsed in Shanks for scalar 2^{}: {:?}",
                exp, elapsed
            );
            // if !flag {
            //     Bj_scalar = tally_helper(&(-g), &point, &BigInt::from((1100000) as u64), &mut flag);
            // }
            assert_eq!(scalar, Bj_scalar);
        }
    }
}