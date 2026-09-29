//! Numerical check of the structural lemmas:
//!  L1: with f <= b aliases every dereference costs exactly 2 reads, forever.
//!  L2: one dereference increases the depth of any alias by at most 1, and
//!      does not increase the depth of the dereferenced alias.
use rand::{rngs::StdRng, Rng, SeedableRng};
use efficient_agro::aliasmodel::AliasTree;
fn main() {
    let mut checked = 0u64;
    let mut worst_inc = 0i64;
    let mut worst_self = i64::MIN;
    for randomized in [false, true] {
        for (b, f) in [(2usize, 2usize), (2, 9), (4, 4), (4, 5), (4, 64), (4, 300), (6, 100), (16, 16), (16, 17), (16, 952), (30, 30), (30, 31), (30, 952), (30, 3000)] {
            for order in 0..4 {
                let mut t = AliasTree::install(f, b, randomized, 5);
                let mut rng = StdRng::seed_from_u64(f as u64);
                let steps = if f > 1000 { 3000 } else { 6000 };
                let mut before = t.depths();
                for s in 0..steps {
                    let a = match order { 0 => rng.gen_range(0..f), 1 => s % f, 2 => (s / 20) % f, _ => { // deepest
                        (0..f).max_by_key(|&j| before[j]).unwrap() } };
                    let r = t.access(a);
                    let after = t.depths();
                    if f <= b { assert_eq!(r, if f == 1 { 1 } else { 2 }, "L1 b={b} f={f}"); }
                    for j in 0..f {
                        let inc = after[j] as i64 - before[j] as i64;
                        if j == a { worst_self = worst_self.max(inc); assert!(inc <= 0, "L2 self"); }
                        else { worst_inc = worst_inc.max(inc); assert!(inc <= 1, "L2 b={b} f={f} inc={inc}"); }
                    }
                    checked += 1;
                    before = after;
                }
            }
        }
    }
    println!("checked {checked} dereferences: max depth increase of others = {worst_inc}, of the accessed alias = {worst_self}");
}
