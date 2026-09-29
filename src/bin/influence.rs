//! Influence of one step's layout randomness on the expected future cost.
//!
//! Doob-martingale concentration needs, for every history and step t,
//!   |E[T | rho_1..rho_t] - E[T | rho_1..rho_{t-1}, rho_t']| <= c.
//! We estimate it: run to step t, branch into two copies that draw different
//! randomness at step t, continue both with the *same* per-step random
//! streams (common random numbers) for L steps, and average the cost
//! difference over M futures.
use rand::{rngs::StdRng, Rng, SeedableRng};
use efficient_agro::aliasmodel::AliasTree;

fn mix(a: u64, b: u64, c: u64) -> u64 {
    a.wrapping_mul(0x9E3779B97F4A7C15) ^ b.wrapping_mul(0xC2B2AE3D27D4EB4F) ^ c.wrapping_mul(0x165667B19E3779F9)
}

fn main() {
    println!("b f schedule L | mean|I| max|I| | mean|D| max|D| p99|D| | coalesced_frac");
    let histories = 20u64;
    let futures = 40u64;
    for (b, f) in [(4usize, 300usize), (16, 952), (30, 952)] {
        for name in ["uniform", "cyclic", "bursts", "hotscan"] {
            let len = 12 * f + 4000;
            let mut rng = StdRng::seed_from_u64(99);
            let sch: Vec<usize> = match name {
                "uniform" => (0..len).map(|_| rng.gen_range(0..f)).collect(),
                "cyclic" => (0..len).map(|i| i % f).collect(),
                "bursts" => (0..len).map(|i| (i / 20) % f).collect(),
                _ => { let mut v = Vec::new(); while v.len() < len { for r in 0..4 * f { v.push(r % b.min(f)); } for i in 0..f { v.push(i); } } v.truncate(len); v }
            };
            for lookahead in [200usize, 2000] {
                let mut infl = Vec::new();
                let mut diffs = Vec::new();
                let mut coalesced = 0usize;
                for h in 0..histories {
                    let t = 2 * f + (h as usize * 997) % (8 * f);
                    let mut base = AliasTree::install(f, b, true, mix(h, 0, 1));
                    for (i, &a) in sch[..t].iter().enumerate() {
                        base.reseed(mix(h, i as u64, 2));
                        base.access(a);
                    }
                    let mut sum = 0.0;
                    for m in 0..futures {
                        let mut x = base.clone();
                        let mut y = base.clone();
                        x.reseed(mix(h, m, 3));
                        y.reseed(mix(h, m, 4));
                        let (mut cx, mut cy) = (x.access(sch[t]), y.access(sch[t]));
                        for (j, &a) in sch[t + 1..t + 1 + lookahead].iter().enumerate() {
                            let s = mix(h * 1000 + m, j as u64, 5);
                            x.reseed(s);
                            y.reseed(s);
                            cx += x.access(a);
                            cy += y.access(a);
                        }
                        let d = cx as f64 - cy as f64;
                        diffs.push(d.abs());
                        sum += d;
                        if x.depths() == y.depths() {
                            coalesced += 1;
                        }
                    }
                    infl.push((sum / futures as f64).abs());
                }
                diffs.sort_by(|a, b| a.partial_cmp(b).unwrap());
                let p99 = diffs[(diffs.len() * 99) / 100];
                let n = infl.len() as f64;
                println!(
                    "{b} {f} {name} {lookahead} | {:.2} {:.2} | {:.2} {:.0} {:.0} | {:.2}",
                    infl.iter().sum::<f64>() / n,
                    infl.iter().cloned().fold(0.0, f64::max),
                    diffs.iter().sum::<f64>() / diffs.len() as f64,
                    diffs.last().unwrap(),
                    p99,
                    coalesced as f64 / diffs.len() as f64
                );
            }
        }
    }
}
