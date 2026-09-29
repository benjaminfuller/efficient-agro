//! Per-get vs windowed (whole-trace) read counts for r-ary aliases.
fn main() {
    let trials = 20000;
    println!("b f mean max_single max_window_k10/10 max_window_k100/100 max_window_k1000/1000 var");
    for b in [16usize, 30] {
        for f in [2usize, 16, 31, 100, 952, 3000] {
            let t = efficient_agro::probe::alias_trace(b, f, trials, 1);
            let mean = t.iter().sum::<u64>() as f64 / t.len() as f64;
            let var = t.iter().map(|&x| (x as f64 - mean).powi(2)).sum::<f64>() / t.len() as f64;
            let win = |k: usize| {
                t.windows(k).map(|w| w.iter().sum::<u64>()).max().unwrap() as f64 / k as f64
            };
            println!("{b} {f} {mean:.2} {} {:.2} {:.2} {:.2} {var:.2}", t.iter().max().unwrap(), win(10), win(100), win(1000));
        }
    }
}
