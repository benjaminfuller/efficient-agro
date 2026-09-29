//! Oblivious adversary vs randomized layout: for a schedule fixed in advance,
//! the distribution (over layout randomness) of the total reads of k
//! consecutive dereferences, after a warm-up of 2f accesses.
use rand::{rngs::StdRng, Rng, SeedableRng};
use efficient_agro::aliasmodel::AliasTree;

fn schedule(name: &str, f: usize, b: usize, len: usize) -> Vec<usize> {
    let mut rng = StdRng::seed_from_u64(424242);
    match name {
        "cyclic" => (0..len).map(|i| i % f).collect(),
        "stride" => (0..len).map(|i| (i * 7919) % f).collect(),
        "uniform" => (0..len).map(|_| rng.gen_range(0..f)).collect(),
        "zipf" => (0..len).map(|_| { let u: f64 = rng.gen(); ((f as f64).powf(u) as usize - 1).min(f - 1) }).collect(),
        // b hot aliases hammered, then one full scan, repeated
        "hotscan" => { let mut v = Vec::new(); while v.len() < len { for r in 0..4 * f { v.push(r % b.min(f)); } for i in 0..f { v.push(i); } } v.truncate(len); v }
        // each alias repeated 20 times in turn
        "bursts" => (0..len).map(|i| (i / 20) % f).collect(),
        _ => unreachable!(),
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let seeds: u64 = args.get(1).and_then(|x| x.parse().ok()).unwrap_or(100);
    println!("b f schedule k | det_mean | rnd_mean rnd_sd rnd_max (max-mean)/sd sd/sqrt(k) | single_max_det single_max_rnd");
    for (b, f) in [(4usize, 300usize), (16, 952), (30, 952), (30, 3000)] {
        for name in ["cyclic", "stride", "uniform", "zipf", "hotscan", "bursts"] {
            for k in [1000usize, 10000] {
                let warm = 2 * f;
                let sch = schedule(name, f, b, warm + k);
                let run = |randomized: bool, seed: u64| {
                    let mut t = AliasTree::install(f, b, randomized, seed);
                    let mut tot = 0u64; let mut mx = 0u64;
                    for (i, &a) in sch.iter().enumerate() {
                        let r = t.access(a);
                        if i >= warm { tot += r; mx = mx.max(r); }
                    }
                    (tot, mx)
                };
                let (det, det_mx) = run(false, 0);
                let res: Vec<(u64, u64)> = (0..seeds).map(|s| run(true, 1000 + s)).collect();
                let n = res.len() as f64;
                let mean = res.iter().map(|r| r.0 as f64).sum::<f64>() / n;
                let sd = (res.iter().map(|r| (r.0 as f64 - mean).powi(2)).sum::<f64>() / (n - 1.0)).sqrt();
                let mx = res.iter().map(|r| r.0).max().unwrap() as f64;
                let smx = res.iter().map(|r| r.1).max().unwrap();
                println!("{b} {f} {name} {k} | {:.3} | {:.3} {:.4} {:.3} {:.2} {:.3} | {det_mx} {smx}",
                    det as f64 / k as f64, mean / k as f64, sd / k as f64, mx / k as f64,
                    if sd > 0.0 { (mx - mean) / sd } else { 0.0 }, sd / (k as f64).sqrt());
            }
        }
    }
}
