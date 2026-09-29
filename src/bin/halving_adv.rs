//! Greedy adaptive adversary against r-ary halving: at each step, try the
//! deepest alias plus 40 random ones on a clone and pick the one with the
//! largest amortized cost reads + beta*dPhi (beta = 2/log2 b), or the largest
//! reads (mode 1). Reports the worst amortized step and the mean reads.
use rand::{rngs::StdRng, Rng, SeedableRng};
use efficient_agro::aliasmodel::AliasTree;
fn main() {
    for &(b, f) in &[(4usize, 952usize), (8, 952), (30, 952), (30, 3000), (4, 3000)] {
        for mode in 0..2 {
            let beta = 2.0 / (b as f64).log2();
            let mut h = AliasTree::install(f, b, false, 0);
            let mut rng = StdRng::seed_from_u64(5);
            let (mut worst, mut tot, mut maxr) = (f64::NEG_INFINITY, 0.0, 0u64);
            let steps = 1500;
            for _ in 0..steps {
                let d = h.depths();
                let mut cands: Vec<usize> = (0..40).map(|_| rng.gen_range(0..f)).collect();
                cands.push((0..f).max_by_key(|&j| d[j]).unwrap());
                let p0 = h.potential();
                let mut best = (f64::NEG_INFINITY, 0usize);
                for &c in &cands {
                    let mut t = h.clone();
                    let r = t.access_halving(c) as f64;
                    let score = if mode == 0 { r + beta * (t.potential() - p0) } else { r };
                    if score > best.0 { best = (score, c); }
                }
                let r = h.access_halving(best.1);
                let a = r as f64 + beta * (h.potential() - p0);
                worst = worst.max(a); tot += r as f64; maxr = maxr.max(r);
            }
            println!("b={b:>2} f={f:>5} adversary={} | worst amortized (beta=2/log2 b) {worst:.2} | 2 log_b f = {:.2} | mean reads {:.2} max {maxr} | 2log2 f+1 = {:.1}",
                if mode == 0 { "max-amortized" } else { "max-reads    " }, 2.0 * (f as f64).ln() / (b as f64).ln(), tot / steps as f64, 2.0 * (f as f64).log2() + 1.0);
        }
    }
}
