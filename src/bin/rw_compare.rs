//! Reads + writes per dereference (the quantity OSAM+ leaks) for the current
//! r-ary rebuild versus a shape-preserving rebuild of the same installed tree.
use rand::{rngs::StdRng, Rng, SeedableRng};
use efficient_agro::aliasmodel::AliasTree;
/// Evenly split b-ary tree as built by RaryPointer::install: returns, for each
/// leaf, (depth, sum of the fanouts of the nodes on its path).
fn static_costs(f: usize, b: usize) -> Vec<(u64, u64)> {
    fn rec(count: usize, b: usize, depth: u64, fan: u64, out: &mut Vec<(u64, u64)>) {
        if count == 1 { out.push((depth, fan)); return; }
        if count <= b { for _ in 0..count { out.push((depth + 1, fan + count as u64)); } return; }
        let (base, rem) = (count / b, count % b);
        for i in 0..b { rec(base + usize::from(i < rem), b, depth + 1, fan + b as u64, out); }
    }
    let mut out = Vec::new(); rec(f, b, 0, 0, &mut out); out
}
fn main() {
    println!("{:>4} {:>5} | {:>22} | {:>28} | {:>10}", "b", "f", "static R+W (min/max)", "current rule R+W mean/max(1 dereference)", "binary 3R-1 cap");
    for &f in &[8usize, 30, 64, 256, 1024, 4096] {
        for &b in &[2usize, 3, 4, 8, 16, 30] {
            let st = static_costs(f, b);
            let c: Vec<u64> = st.iter().map(|&(h, s)| (h + 1) + s + 1).collect();
            let (smin, smax) = (*c.iter().min().unwrap(), *c.iter().max().unwrap());
            let mut t = AliasTree::install(f, b, false, 0);
            let mut rng = StdRng::seed_from_u64(7);
            let (mut tot, mut mx, n) = (0u64, 0u64, 20 * f.max(50));
            for step in 0..n {
                let i = if step % 2 == 0 { rng.gen_range(0..f) } else { step % f };
                let r = t.access(i); let w = t.last_writes();
                tot += r + w; mx = mx.max(r + w);
            }
            let bin = 6.0 * (f as f64).log2() + 2.0;
            println!("{b:>4} {f:>5} | {smin:>10} / {smax:<9} | {:>14.1} / {mx:<11} | {bin:>10.1}", tot as f64 / n as f64);
        }
    }
}
