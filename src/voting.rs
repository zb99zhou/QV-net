#![allow(non_snake_case)]

use std::collections::HashMap;

use ark_ff::{One, Zero};
use curv::{arithmetic::{Converter, Roots}, elliptic::curves::{Point, Scalar, Secp256k1}, BigInt};
use merlin::Transcript;
use VarRange::proofs::varrange::VarRange;

use crate::{sigma_dl::SigmaDlProof, sigma_dleq::SigmaDleqProof, sum_square::ZkSumSquareArg};

#[derive(Clone)]
pub struct Voter {
    VoterID: usize,
    xi_vec: Vec<Scalar<Secp256k1>>,
    yi_vec: Vec<Point<Secp256k1>>,
    proofi_dl: SigmaDlProof
}

pub struct Board {
    pp: PublicParam,
    y_vec: Vec<Vec<Point<Secp256k1>>>,
    ballot_proof: Vec<BallotWithProof>
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
    Bi_vec: Vec<Point<Secp256k1>>,
    Bi: Point<Secp256k1>,
    proofi_dl: SigmaDlProof,
    proofi_eq: SigmaDleqProof,
    proofi_ss: ZkSumSquareArg,
    proofi_ar: VarRange
}

// Shanks' baby-step giant-step algorithm, calculate ind_g{B}
pub fn tally_helper(g: &Point<Secp256k1>, B: &Point<Secp256k1>, bound: &BigInt) -> Scalar<Secp256k1> {
    let n_bn = bound.sqrt();
    let n = &Scalar::<Secp256k1>::from_bigint(&n_bn);
    let mut counter = BigInt::one();
    let one = &Scalar::<Secp256k1>::from(1);
    let one_bn = &BigInt::from(1);
    let mut res = Scalar::<Secp256k1>::zero();

    let mut p = Scalar::<Secp256k1>::from(1);
    let mut q = Scalar::<Secp256k1>::from(0);
    let mut giants: HashMap<Vec<u8>, Scalar<Secp256k1>> = HashMap::new();
    let mut giant: Point<Secp256k1>;
    let mut baby: Point<Secp256k1>;
    
    while counter <= n_bn {
        giant = g * (p.clone() * n);
        giants.insert(giant.x_coord().unwrap().to_bytes(), p.clone());
        p = p + one;
        counter += one_bn;
    }

    counter = BigInt::zero();
    while counter <= n_bn {
        baby = B + g * q.clone();
        if let Some(&ref giant_p) = giants.get(&baby.x_coord().unwrap().to_bytes()) {
            res = giant_p.clone();
            break;
        }
        q = q + one;
        counter += one_bn;
    }

    if g * (res.clone() * n - q.clone()) == *B {
        return res * n - q;
    } else {
        return - res * n - q;
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

        let mut xi_vec: Vec<Scalar<Secp256k1>> = Vec::new();
        for _j in 0..nc {
            xi_vec.push(Scalar::<Secp256k1>::random());
        }

        let yi_vec = (0..nc)
            .map(|j| &h_vec[j].clone() * xi_vec[j].clone())
            .collect::<Vec<Point<Secp256k1>>>();
        
        let mut transcript = Transcript::new(b"Proof");
        let proofi_dl = SigmaDlProof::prove(&mut transcript, &xi_vec, &h_vec, nc);

        BulletinBoard.y_vec.push(yi_vec.clone());
        
        Voter {
            VoterID,
            xi_vec,
            yi_vec,
            proofi_dl
        }
    }

    pub fn vote(
        &self,
        vi_vec: &[Scalar<Secp256k1>],
        ti: Scalar<Secp256k1>,
        BulletinBoard: &mut Board,
        seed: &BigInt
    ) -> BallotWithProof {
        assert_eq!(vi_vec.len(), BulletinBoard.pp.nc);
        assert_eq!(BulletinBoard.y_vec.len(), BulletinBoard.pp.nv);
        assert_eq!(BulletinBoard.y_vec[0].len(), BulletinBoard.pp.nc);
        assert_eq!(BulletinBoard.pp.g_vec.len(), BulletinBoard.pp.nc);
        assert_eq!(BulletinBoard.pp.h_vec.len(), BulletinBoard.pp.nc);

        let mut Yi_vec: Vec<Point<Secp256k1>> = vec![Point::<Secp256k1>::zero(); BulletinBoard.pp.nc];
        for j in 0..BulletinBoard.pp.nc {
            for k in 0..self.VoterID {
                Yi_vec[j] = &Yi_vec[j] + &BulletinBoard.y_vec[k][j];
            }
            for k in self.VoterID+1..BulletinBoard.pp.nv {
                Yi_vec[j] = &Yi_vec[j] - &BulletinBoard.y_vec[k][j];
            }
        }

        let Bi_vec = (0..BulletinBoard.pp.nc)
            .map(|j| &BulletinBoard.pp.g_vec[j].clone()*vi_vec[j].clone() + &Yi_vec[j].clone()*self.xi_vec[j].clone())
            .collect::<Vec<Point<Secp256k1>>>();

        let mut ui = Scalar::<Secp256k1>::zero();
        for j in 0..BulletinBoard.pp.nc {
            ui = ui + vi_vec[j].clone() * vi_vec[j].clone();
        }

        let ri = Scalar::<Secp256k1>::random();
        let Bi = &BulletinBoard.pp.g * ui.clone() + &BulletinBoard.pp.h * ri.clone();

        let mut transcript = Transcript::new(b"Proof");
        let proofi_eq = SigmaDleqProof::prove(
            &mut transcript,
            &vi_vec, 
            &self.xi_vec, 
            &BulletinBoard.pp.h_vec, 
            &BulletinBoard.pp.g_vec, 
            &Yi_vec, 
            BulletinBoard.pp.nc
        );

        let mut transcript = Transcript::new(b"Proof");
        let proofi_ss = ZkSumSquareArg::prove(
            &mut transcript, 
            &BulletinBoard.pp.g_vec, 
            &Yi_vec, 
            &BulletinBoard.pp.g, 
            &BulletinBoard.pp.h, 
            &Bi_vec.clone().into_iter().sum(), 
            &Bi, 
            vi_vec, 
            &self.xi_vec, 
            &ui, 
            &ri, 
            BulletinBoard.pp.nc, 
            seed
        );

        let mut gi = BulletinBoard.pp.g_vec.clone();
        gi.push(BulletinBoard.pp.g.clone());
        let mut hi = Yi_vec;
        hi.push(BulletinBoard.pp.h.clone());
        let mut values = vi_vec.to_vec();
        values.push(ui.clone());
        let mut x_vec = self.xi_vec.clone();
        x_vec.push(ri);
        let mut B_vec = Bi_vec.clone();
        B_vec.push(Bi.clone());
        let si = BigInt::sqrt(&ti.to_bigint());
        for j in 0..values.len()-1 {
            assert!(values[j].to_bigint() < si);
        }
        assert!(ui.to_bigint() < ti.to_bigint());

        let mut transcript = Transcript::new(b"Proof");
        let proofi_ar = VarRange::range_prove(
            &mut transcript, 
            &gi,
            &hi,
            values, 
            x_vec, 
            Scalar::<Secp256k1>::from_bigint(&si), 
            ti, 
            BulletinBoard.pp.nc, 
            &B_vec, 
            seed
        );

        let ballot_proof = BallotWithProof {
            Bi_vec,
            Bi,
            proofi_dl: self.proofi_dl.clone(),
            proofi_eq,
            proofi_ss,
            proofi_ar
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

        Board { pp, y_vec, ballot_proof }
    }

    pub fn verify(
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

        for i in 0..self.pp.nv {
            let mut Yi_vec: Vec<Point<Secp256k1>> = vec![Point::<Secp256k1>::zero(); self.pp.nc];
            for j in 0..self.pp.nc {
                for k in 0..i {
                    Yi_vec[j] = &Yi_vec[j] + &self.y_vec[k][j];
                }
                for k in i+1..self.pp.nv {
                    Yi_vec[j] = &Yi_vec[j] - &self.y_vec[k][j];
                }
            }

            // verify the proofs
            let mut transcript = Transcript::new(b"Proof");
            let res_dl = self.ballot_proof[i].proofi_dl.verify(
                &mut transcript, 
                &self.y_vec[i], 
                &self.pp.h_vec, 
                self.pp.nc
            );
            assert!(res_dl.is_ok());

            let mut transcript = Transcript::new(b"Proof");
            let res_dleq = self.ballot_proof[i].proofi_eq.verify(
                &mut transcript, 
                &self.y_vec[i], 
                &self.ballot_proof[i].Bi_vec, 
                &self.pp.h_vec, 
                &self.pp.g_vec, 
                &Yi_vec, 
                self.pp.nc
            );
            assert!(res_dleq.is_ok());

            let mut transcript = Transcript::new(b"Proof");
            let res_ss = self.ballot_proof[i].proofi_ss.verify(
                &mut transcript,
                &self.pp.g_vec, 
                &Yi_vec,
                &self.pp.g, 
                &self.pp.h, 
                &self.ballot_proof[i].Bi_vec.clone().into_iter().sum(), 
                &self.ballot_proof[i].Bi, 
                self.pp.nc, 
                seed
            );
            assert!(res_ss.is_ok());

            let mut gi = self.pp.g_vec.clone();
            gi.push(self.pp.g.clone());
            let mut hi = Yi_vec;
            hi.push(self.pp.h.clone());
            let mut B_vec = self.ballot_proof[i].Bi_vec.clone();
            B_vec.push(self.ballot_proof[i].Bi.clone());
            let si = BigInt::sqrt(&token[i].clone().to_bigint());

            let mut transcript = Transcript::new(b"Proof");
            let res_ar = self.ballot_proof[i].proofi_ar.range_verify(
                &mut transcript, 
                &gi, 
                &hi,
                &B_vec, 
                Scalar::<Secp256k1>::from_bigint(&si), 
                token[i].clone(), 
                self.pp.nc,
                seed
            );
            assert!(res_ar.is_ok());
        }
    }

    pub fn tally(
        &self
    ) -> Vec<Scalar<Secp256k1>> {
        assert_eq!(self.y_vec.len(), self.pp.nv);
        assert_eq!(self.y_vec[0].len(), self.pp.nc);
        assert_eq!(self.pp.g_vec.len(), self.pp.nc);
        assert_eq!(self.pp.h_vec.len(), self.pp.nc);
        assert_eq!(self.ballot_proof.len(), self.pp.nv);

        let mut ballots: Vec<Scalar<Secp256k1>> = Vec::with_capacity(self.pp.nc);
        for j in 0..self.pp.nc {
            let mut Bj = Point::<Secp256k1>::zero();
            for i in 0..self.pp.nv {
                Bj = Bj + self.ballot_proof[i].Bi_vec[j].clone();
            }
            ballots.push(tally_helper(&self.pp.g_vec[j], &Bj, &BigInt::from((1_u64 << 16) as u64)));
        }
        ballots
    }
}

mod test {
    use curv::{arithmetic::Converter, cryptographic_primitives::hashing::DigestExt, elliptic::curves::{Point, Scalar, Secp256k1}, BigInt};
    use sha2::{Digest, Sha512};
    use rand::Rng;

    use crate::{sigma_dl::generate_random_point, voting::tally_helper};

    use super::{Board, Voter};

    pub fn test_helper_with_verify(seed: &BigInt, nv: usize, nc: usize) {
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

        let kzen_label =  BigInt::from(1 as u32) + BigInt::from(nc as u32) + BigInt::from(nc as u32) + seed;
        let hash = Sha512::new().chain_bigint(&kzen_label).result_bigint();
        let g = generate_random_point(&Converter::to_bytes(&hash));

        let mut Voter_vec: Vec<Voter> = Vec::new();
        let mut y_vec: Vec<Vec<Point<Secp256k1>>> = Vec::new();
        let mut board = Board::setup(&g_vec, &h_vec, &g, &h, nv, nc);
        for i in 0..nv {
            let Vi = Voter::gen(i, &h_vec, nv, nc, &mut board);
            y_vec.push(Vi.yi_vec.clone());
            Voter_vec.push(Vi.clone());
        }
        
        let mut rng = rand::thread_rng();
        let start: u64 = 1;
        let end: u64 = 10;
        let mut v_vec: Vec<Vec<Scalar<Secp256k1>>> = Vec::new();
        let mut tokens: Vec<Scalar<Secp256k1>> = Vec::new();
        for i in 0..nv {
            let mut vi_vec: Vec<Scalar<Secp256k1>> = Vec::new();
            let mut square_sum = 0;
            for _ in 0..nc {
                let vij = rng.gen_range(start..=end);
                vi_vec.push(Scalar::<Secp256k1>::from(vij));
                square_sum += vij * vij;
            }
            square_sum += rng.gen_range(start..=end);
            tokens.push(Scalar::<Secp256k1>::from(square_sum));
            Voter_vec[i].vote(&vi_vec, tokens[i].clone(), &mut board, &(BigInt::from((nc*2+1) as u32) + seed));
            v_vec.push(vi_vec);
        }
        
        board.verify(tokens, &(BigInt::from((nc*2+1) as u32) + seed));
        let res = board.tally();

        // test correctness
        for j in 0..nc {
            let mut Bj = Scalar::<Secp256k1>::zero();
            for i in 0..nv {
                Bj = Bj + v_vec[i][j].clone();
            }
            assert_eq!(Bj, res[j]);
        }
    }

    #[test]
    pub fn test_voting_with_verify() {
        let KZen: &[u8] = &[75, 90, 101, 110];
        let kzen_label = BigInt::from_bytes(KZen);
        test_helper_with_verify(&kzen_label, 3, 10);
    }

    #[test]
    pub fn test_shanks() {
        let KZen: &[u8] = &[75, 90, 101, 110];
        let kzen_label = BigInt::from_bytes(KZen);
        let hash = Sha512::new().chain_bigint(&kzen_label).result_bigint();
        let g = generate_random_point(&Converter::to_bytes(&hash));
        let scalar = Scalar::<Secp256k1>::from(65564);
        let point = &g * scalar.clone();
        assert_eq!(scalar, tally_helper(&g, &point, &BigInt::from((1_u64 << 32) as u64)));
    }
}