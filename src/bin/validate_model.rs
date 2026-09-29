//! Checks the standalone alias-tree model against sam-model's RaryPointer,
//! read count by read count, on several access orders.
use rand::{rngs::StdRng, Rng, SeedableRng};
use sam_model::pointer::{RaryCell, RaryPointer};
use sam_model::{AccessPolicy, DryRunSam, SingleAccessMachine};
use efficient_agro::aliasmodel::AliasTree;
fn main() {
    let mut total = 0;
    for (b, f) in [(2usize, 7usize), (2, 50), (4, 64), (4, 300), (6, 77), (8, 300), (16, 952), (30, 31), (30, 952), (30, 3000)] {
        for order in 0..3 {
            let mut sam: DryRunSam<RaryCell<u64>> = DryRunSam::new(AccessPolicy::MULTI_WRITE);
            let mut ptrs = RaryPointer::install(&mut sam, 7u64, f, b).unwrap();
            let mut model = AliasTree::install(f, b, false, 0);
            let mut rng = StdRng::seed_from_u64(b as u64 * 1000 + f as u64);
            for step in 0..5000 {
                let i = match order { 0 => rng.gen_range(0..f), 1 => step % f, _ => (step * 7919) % f };
                let r0 = sam.stats().operations.reads;
                let w0 = sam.stats().operations.writes;
                ptrs[i].get(&mut sam).unwrap();
                let lib = sam.stats().operations.reads - r0;
                let libw = sam.stats().operations.writes - w0;
                let m = model.access(i);
                assert_eq!(lib, m, "b={b} f={f} order={order} step={step}: library {lib} reads, model {m}");
                assert_eq!(libw, model.last_writes(), "b={b} f={f} order={order} step={step}: library {libw} writes, model {}", model.last_writes());
                total += 1;
            }
        }
    }
    println!("model matches library (reads and writes) on {total} dereferences");
}
