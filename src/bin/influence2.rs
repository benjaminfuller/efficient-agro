//! Signed influence estimates with standard errors: is E[D] distinguishable
//! from zero, and does it grow with the look-ahead L?
use rand::{rngs::StdRng, Rng, SeedableRng};
use efficient_agro::aliasmodel::AliasTree;
fn mix(a: u64, b: u64, c: u64) -> u64 {
    a.wrapping_mul(0x9E3779B97F4A7C15) ^ b.wrapping_mul(0xC2B2AE3D27D4EB4F) ^ c.wrapping_mul(0x165667B19E3779F9)
}
fn main() {
    let (b, f) = (30usize, 952usize);
    let futures = 400u64;
    for name in ["uniform", "cyclic", "hotscan"] {
        let len = 12 * f + 8000;
        let mut rng = StdRng::seed_from_u64(99);
        let sch: Vec<usize> = match name {
            "uniform" => (0..len).map(|_| rng.gen_range(0..f)).collect(),
            "cyclic" => (0..len).map(|i| i % f).collect(),
            _ => { let mut v = Vec::new(); while v.len() < len { for r in 0..4 * f { v.push(r % b); } for i in 0..f { v.push(i); } } v.truncate(len); v }
        };
        for h in 0..6u64 {
            let t = 2 * f + (h as usize * 1597) % (8 * f);
            let mut base = AliasTree::install(f, b, true, mix(h, 0, 1));
            for (i, &a) in sch[..t].iter().enumerate() { base.reseed(mix(h, i as u64, 2)); base.access(a); }
            let checkpoints = [1usize, 10, 100, 1000, 4000];
            let mut d = vec![Vec::new(); checkpoints.len()];
            for m in 0..futures {
                let (mut x, mut y) = (base.clone(), base.clone());
                x.reseed(mix(h, m, 3)); y.reseed(mix(h, m, 4));
                let mut diff = x.access(sch[t]) as i64 - y.access(sch[t]) as i64;
                let mut ci = 0;
                for j in 0..*checkpoints.last().unwrap() {
                    if j + 1 == checkpoints[ci] { d[ci].push(diff as f64); ci += 1; if ci == checkpoints.len() { break; } }
                    let s = mix(h * 7919 + m, j as u64, 5);
                    x.reseed(s); y.reseed(s);
                    let a = sch[t + 1 + j];
                    diff += x.access(a) as i64 - y.access(a) as i64;
                }
                if ci < checkpoints.len() { d[ci].push(diff as f64); }
            }
            print!("{name} hist{h}:");
            for (k, v) in d.iter().enumerate() {
                let n = v.len() as f64;
                let mean = v.iter().sum::<f64>() / n;
                let sd = (v.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / (n - 1.0)).sqrt();
                print!("  L={}: {:+.2}±{:.2}", checkpoints[k], mean, sd / n.sqrt());
            }
            println!();
        }
    }
}
