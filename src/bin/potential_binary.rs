//! Checks the access lemma for sam-model's binary multi-write pointer.
//!
//! Claim: with Phi = sum over internal nodes v of log2 s(v) (s = aliases below v),
//! every dereference satisfies  reads <= 2 log2 f + 1 + Phi_before - Phi_after,
//! for every access order (the bound is deterministic and holds adaptively).
use rand::{rngs::StdRng, Rng, SeedableRng};
use sam_model::pointer::{MultiWriteCell, MultiWritePointer};
use sam_model::{AccessPolicy, Address, DryRunSam, SingleAccessMachine};
use std::collections::HashMap;

fn shape(sam: &DryRunSam<MultiWriteCell<u64>>) -> (f64, HashMap<u64, usize>) {
    let snap = sam.snapshot();
    let mut parent = HashMap::new();
    let mut has_child: HashMap<u64, bool> = HashMap::new();
    let mut root = 0;
    for b in &snap.blocks {
        match &b.value {
            MultiWriteCell::Root(_) => root = b.identifier,
            MultiWriteCell::Inner { parent: Address::Oblivious(p), .. } => {
                parent.insert(b.identifier, *p);
                has_child.insert(*p, true);
            }
            _ => panic!(),
        }
    }
    let mut s: HashMap<u64, usize> = HashMap::new();
    let mut depth = HashMap::new();
    for (&a, _) in parent.iter().filter(|(a, _)| !has_child.contains_key(a)) {
        let (mut d, mut c) = (0, a);
        while c != root {
            c = parent[&c];
            d += 1;
            *s.entry(c).or_default() += 1;
        }
        depth.insert(a, d);
    }
    (s.values().map(|&x| (x as f64).log2()).sum(), depth)
}

fn main() {
    let steps = 3000;
    println!("f order | mean_reads max_reads | max_slack(reads - 2log2f - 1 - dPhi) | Phi0-(f-1) | T/k  bound/k");
    let mut violations = 0;
    for f in [2usize, 3, 5, 16, 100, 952, 3000] {
        for order in ["random", "cyclic", "reverse", "zipf", "bursts", "deepest"] {
            let mut sam: DryRunSam<MultiWriteCell<u64>> = DryRunSam::new(AccessPolicy::MULTI_WRITE);
            let mut first = MultiWritePointer::new(&mut sam, 7u64).unwrap();
            let mut ptrs = vec![];
            if f > 1 {
                let copies = first.copy_many(&mut sam, f - 1).unwrap();
                ptrs.push(first);
                ptrs.extend(copies);
            } else {
                ptrs.push(first);
            }
            let (mut phi, mut depth) = shape(&sam);
            let phi0 = phi;
            let lg = (f as f64).log2();
            let mut rng = StdRng::seed_from_u64(3);
            let (mut tot, mut mx, mut worst) = (0u64, 0u64, f64::MIN);
            for step in 0..steps {
                let i = match order {
                    "random" => rng.gen_range(0..f),
                    "cyclic" => step % f,
                    "reverse" => f - 1 - step % f,
                    "zipf" => { let u: f64 = rng.gen(); ((f as f64).powf(u) as usize - 1).min(f - 1) }
                    "bursts" => (step / 20) % f,
                    _ => (0..f).max_by_key(|&j| { let Some(Address::Oblivious(h)) = ptrs[j].head() else { unreachable!() }; depth[&h] }).unwrap(),
                };
                let r0 = sam.stats().operations.reads;
                ptrs[i].get(&mut sam).unwrap();
                let reads = sam.stats().operations.reads - r0;
                let (nphi, nd) = shape(&sam);
                let slack = reads as f64 - (2.0 * lg + 1.0) - (phi - nphi);
                if slack > 1e-9 { violations += 1; }
                worst = worst.max(slack);
                tot += reads; mx = mx.max(reads);
                phi = nphi; depth = nd;
            }
            let bound = steps as f64 * (2.0 * lg + 1.0) + phi0 - phi;
            println!("{f} {order} | {:.2} {mx} | {:+.3} | {:.1} | {:.2} {:.2}", tot as f64 / steps as f64, worst, phi0 - (f as f64 - 1.0), tot as f64 / steps as f64, bound / steps as f64);
        }
    }
    println!("violations of the per-step access lemma: {violations}");
}
