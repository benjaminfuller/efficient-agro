//! Checks r-ary path halving: access lemma reads <= 2 log2 f + 1 + Phi - Phi',
//! depth growth <= 1 for every other alias, flat regime; compares mean reads
//! with the current rule.
use rand::{rngs::StdRng, Rng, SeedableRng};
use efficient_agro::aliasmodel::AliasTree;
fn zipf(rng: &mut StdRng, n: usize) -> usize {
    let h: f64 = (1..=n).map(|i| 1.0 / i as f64).sum();
    let mut u = rng.gen::<f64>() * h;
    for i in 1..=n { u -= 1.0 / i as f64; if u <= 0.0 { return i - 1; } }
    n - 1
}
fn main() {
    let steps = 3000;
    let mut total = 0u64; let mut worst_slack = f64::INFINITY; let mut max_growth = 0i64;
    println!("{:>3} {:>5} {:>9} | {:>8} {:>8} | {:>8} {:>8} | {:>6}", "b", "f", "schedule", "halv mean", "halv max", "curr mean", "curr max", "2lg f+1");
    for &b in &[4usize, 8, 30] {
        for &f in &[16usize, 100, 952, 3000] {
            for sched in ["uniform", "cyclic", "zipf", "deepest", "bursty"] {
                let mut h = AliasTree::install(f, b, false, 0);
                let mut c = AliasTree::install(f, b, false, 0);
                let mut rng = StdRng::seed_from_u64(f as u64 * 31 + b as u64);
                let (mut hs, mut hm, mut cs, mut cm) = (0u64, 0u64, 0u64, 0u64);
                let bound = 2.0 * (f as f64).log2() + 1.0;
                let mut burst = 0usize;
                for t in 0..steps {
                    let i = match sched {
                        "uniform" => rng.gen_range(0..f),
                        "cyclic" => t % f,
                        "zipf" => zipf(&mut rng, f),
                        "bursty" => { if t % 20 == 0 { burst = rng.gen_range(0..f); } burst }
                        _ => { let d = h.depths(); (0..f).max_by_key(|&j| d[j]).unwrap() }
                    };
                    let before = h.depths(); let phi0 = h.potential();
                    let r = h.access_halving(i);
                    let phi1 = h.potential(); let after = h.depths();
                    for j in 0..f { if j != i { max_growth = max_growth.max(after[j] as i64 - before[j] as i64); } else { assert!(after[j] <= before[j].max(1)); } }
                    let slack = bound + phi0 - phi1 - r as f64;
                    worst_slack = worst_slack.min(slack);
                    assert!(slack > -1e-9, "access lemma violated b={b} f={f} {sched} t={t}: reads {r}, bound {bound}, dphi {}", phi1 - phi0);
                    if f <= b { assert_eq!(r, 2); }
                    hs += r; hm = hm.max(r);
                    let ci = if sched == "deepest" { let d = c.depths(); (0..f).max_by_key(|&j| d[j]).unwrap() } else { i };
                    let rc = c.access(ci); cs += rc; cm = cm.max(rc);
                    total += 1;
                }
                println!("{b:>3} {f:>5} {sched:>9} | {:>8.2} {hm:>8} | {:>8.2} {cm:>8} | {bound:>6.1}", hs as f64 / steps as f64, cs as f64 / steps as f64);
            }
        }
    }
    println!("checked {total} dereferences: access lemma holds (min slack {worst_slack:.3}), max depth growth of another alias = {max_growth}");
}
