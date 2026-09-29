//! Measures r-ary alias-walk length (reads per get) under random gets.
use rand::{rngs::StdRng, Rng, SeedableRng};
use sam_model::pointer::{RaryCell, RaryPointer};
use sam_model::{AccessPolicy, DryRunSam, SingleAccessMachine};

/// Returns per-get read counts for `trials` uniformly random gets on one
/// object with `f` aliases (after one warm-up pass over every alias).
pub fn alias_trace(b: usize, f: usize, trials: usize, seed: u64) -> Vec<u64> {
    let mut sam: DryRunSam<RaryCell<u64>> = DryRunSam::new(AccessPolicy::MULTI_WRITE);
    let mut ptrs = RaryPointer::install(&mut sam, 7u64, f, b).unwrap();
    for p in ptrs.iter_mut() {
        p.get(&mut sam).unwrap();
    }
    let mut rng = StdRng::seed_from_u64(seed);
    (0..trials)
        .map(|_| {
            let i = rng.gen_range(0..f);
            let r0 = sam.stats().operations.reads;
            assert_eq!(ptrs[i].get(&mut sam).unwrap(), 7);
            sam.stats().operations.reads - r0
        })
        .collect()
}

/// Per-get read counts for a deterministic access order `order` (alias ids).
pub fn alias_trace_order(b: usize, f: usize, order: &[usize]) -> Vec<u64> {
    let mut sam: DryRunSam<RaryCell<u64>> = DryRunSam::new(AccessPolicy::MULTI_WRITE);
    let mut ptrs = RaryPointer::install(&mut sam, 7u64, f, b).unwrap();
    order
        .iter()
        .map(|&i| {
            let r0 = sam.stats().operations.reads;
            assert_eq!(ptrs[i].get(&mut sam).unwrap(), 7);
            sam.stats().operations.reads - r0
        })
        .collect()
}
