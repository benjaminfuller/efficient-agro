//! Empirical amortized cost of r-ary halving under the potential beta*Phi,
//! Phi = sum log2 size: A = max_t reads_t + beta*(Phi_{t+1} - Phi_t).
use rand::{rngs::StdRng, Rng, SeedableRng};
use efficient_agro::aliasmodel::AliasTree;
fn main() {
    let steps = 4000;
    println!("  b     f | log_b f | A(beta = c/log2 b) for c = 0.5 1 2 4 | A(beta=1)");
    for &b in &[4usize, 8, 30] {
        for &f in &[100usize, 952, 3000] {
            let cs = [0.5, 1.0, 2.0, 4.0];
            let mut worst = vec![f64::NEG_INFINITY; 5];
            for sched in 0..6 {
                let mut h = AliasTree::install(f, b, false, 0);
                let mut rng = StdRng::seed_from_u64(99 + sched as u64);
                let mut burst = 0;
                for t in 0..steps {
                    let i = match sched {
                        0 => rng.gen_range(0..f), 1 => t % f, 2 => (f - 1) - t % f,
                        3 => { if t % 20 == 0 { burst = rng.gen_range(0..f); } burst }
                        4 => (t * 7919) % f,
                        _ => { let d = h.depths(); (0..f).max_by_key(|&j| d[j]).unwrap() }
                    };
                    let p0 = h.potential(); let r = h.access_halving(i) as f64; let dp = h.potential() - p0;
                    for (k, &c) in cs.iter().enumerate() { worst[k] = worst[k].max(r + c / (b as f64).log2() * dp); }
                    worst[4] = worst[4].max(r + dp);
                }
            }
            println!("{b:>3} {f:>5} | {:>7.2} | {:>6.2} {:>6.2} {:>6.2} {:>6.2} | {:>6.2}", (f as f64).ln() / (b as f64).ln(), worst[0], worst[1], worst[2], worst[3], worst[4]);
        }
    }
}
