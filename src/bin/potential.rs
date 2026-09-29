//! Numerical test of amortized bounds for r-ary alias walks.
//!
//! After every `get` we rebuild the alias tree from the dry-run SAM's live
//! cells and evaluate S = sum over internal nodes v of log2 s(v), where s(v)
//! is the number of aliases below v. For a scale c, the amortized cost of
//! step i is  reads_i + c (S_i - S_{i-1}).  A valid potential needs
//! max_i amortized_i <= A(f) with A growing like log f.
use rand::{rngs::StdRng, Rng, SeedableRng};
use sam_model::pointer::{RaryCell, RaryPointer};
use sam_model::{AccessPolicy, Address, DryRunSam, SingleAccessMachine};
use std::collections::HashMap;

struct Shape {
    /// sum_v log2 s(v) over internal nodes (nodes with children)
    s_log: f64,
    /// depth (edges to root) of each alias leaf, by address
    depth: HashMap<u64, usize>,
    internal: usize,
}

fn shape(sam: &DryRunSam<RaryCell<u64>>) -> Shape {
    let snap = sam.snapshot();
    let mut parent: HashMap<u64, u64> = HashMap::new();
    let mut children: HashMap<u64, usize> = HashMap::new();
    let mut root = 0;
    for b in &snap.blocks {
        match &b.value {
            RaryCell::Root(_) => root = b.identifier,
            RaryCell::Node { parent: Address::Oblivious(p), .. } => {
                parent.insert(b.identifier, *p);
                *children.entry(*p).or_default() += 1;
            }
            _ => panic!("unexpected cell"),
        }
    }
    let mut s: HashMap<u64, usize> = HashMap::new();
    let mut depth = HashMap::new();
    for (&a, _) in parent.iter().filter(|(a, _)| !children.contains_key(a)) {
        // alias leaf: walk to root
        let mut d = 0;
        let mut cur = a;
        while cur != root {
            cur = parent[&cur];
            d += 1;
            *s.entry(cur).or_default() += 1;
        }
        depth.insert(a, d);
    }
    let s_log = s.values().map(|&x| (x as f64).log2()).sum();
    Shape { s_log, depth, internal: s.len() }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let steps: usize = args.get(1).and_then(|x| x.parse().ok()).unwrap_or(4000);
    let cs = [0.0, 0.5, 1.0, 2.0, 3.0, 4.0, 6.0];
    print!("b f order steps mean_reads max_reads range_S");
    for c in cs { print!(" maxam_c{c}"); }
    println!();
    for (b, f) in [(4usize, 64usize), (4, 300), (8, 300), (16, 200), (16, 952), (30, 100), (30, 952), (30, 3000)] {
        for order in ["random", "cyclic", "deepest", "zipf"] {
            let mut sam: DryRunSam<RaryCell<u64>> = DryRunSam::new(AccessPolicy::MULTI_WRITE);
            let mut ptrs = RaryPointer::install(&mut sam, 7u64, f, b).unwrap();
            let mut rng = StdRng::seed_from_u64(11);
            let mut sh = shape(&sam);
            let (mut smin, mut smax) = (sh.s_log, sh.s_log);
            let mut maxam = vec![f64::MIN; cs.len()];
            let (mut sum, mut mx) = (0u64, 0u64);
            for step in 0..steps {
                let i = match order {
                    "random" => rng.gen_range(0..f),
                    "cyclic" => step % f,
                    "zipf" => { let u: f64 = rng.gen(); ((f as f64).powf(u) as usize - 1).min(f - 1) }
                    _ => {
                        // adversary: an alias of maximum current depth
                        let best = ptrs.iter().enumerate().map(|(j, p)| {
                            let Some(Address::Oblivious(h)) = p.head() else { unreachable!() };
                            (sh.depth[&h], j)
                        }).max().unwrap();
                        best.1
                    }
                };
                let r0 = sam.stats().operations.reads;
                ptrs[i].get(&mut sam).unwrap();
                let reads = sam.stats().operations.reads - r0;
                let new = shape(&sam);
                for (k, c) in cs.iter().enumerate() {
                    let am = reads as f64 + c * (new.s_log - sh.s_log);
                    maxam[k] = maxam[k].max(am);
                }
                sum += reads; mx = mx.max(reads);
                smin = smin.min(new.s_log); smax = smax.max(new.s_log);
                sh = new;
            }
            print!("{b} {f} {order} {steps} {:.2} {mx} {:.1}", sum as f64 / steps as f64, smax - smin);
            for m in maxam { print!(" {m:.1}"); }
            println!("  internal={}", sh.internal);
        }
    }
}
