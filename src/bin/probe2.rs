//! Access-order dependence of r-ary alias walks: random vs cyclic vs
//! "least recently used first" (cyclic over a shuffled order) vs a skewed
//! (Zipf-like) workload.
use rand::{rngs::StdRng, Rng, SeedableRng};
fn main() {
    let k = 30000;
    println!("b f order mean max win10 win100 win1000");
    for (b, f) in [(30usize, 100usize), (30, 952), (30, 3000), (16, 952)] {
        let mut rng = StdRng::seed_from_u64(7);
        let random: Vec<usize> = (0..k).map(|_| rng.gen_range(0..f)).collect();
        let cyclic: Vec<usize> = (0..k).map(|i| i % f).collect();
        let rev: Vec<usize> = (0..k).map(|i| f - 1 - (i % f)).collect();
        let zipf: Vec<usize> = (0..k).map(|_| { let u: f64 = rng.gen(); ((f as f64).powf(u) as usize - 1).min(f - 1) }).collect();
        for (name, order) in [("random", &random), ("cyclic", &cyclic), ("reverse", &rev), ("zipf", &zipf)] {
            let t = efficient_agro::probe::alias_trace_order(b, f, order);
            let t = &t[f.min(t.len()/2)..]; // drop warm-up
            let mean = t.iter().sum::<u64>() as f64 / t.len() as f64;
            let w = |k: usize| t.windows(k).map(|w| w.iter().sum::<u64>()).max().unwrap() as f64 / k as f64;
            println!("{b} {f} {name} {mean:.2} {} {:.2} {:.2} {:.2}", t.iter().max().unwrap(), w(10), w(100), w(1000));
        }
    }
}
