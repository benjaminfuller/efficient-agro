//! Randomized correctness and trace-length tests on small synthetic data.

use rand::{rngs::StdRng, Rng, SeedableRng};
use sam_model::{AccessPolicy, DryRunSam};
use efficient_agro::ctx::{dummy_reads, reads, Cell, Ctx};
use efficient_agro::data::{oracle, random_queries, Dataset, Point};
use efficient_agro::layered::{Kind, Layered, LayeredConfig, Tails};
use efficient_agro::sparq::Sparq;

fn synthetic(d: usize, n: usize, domain: i64, seed: u64) -> Dataset {
    let mut rng = StdRng::seed_from_u64(seed);
    let mut seen = std::collections::HashSet::new();
    let mut pts = Vec::new();
    while pts.len() < n {
        let mut c = [0i64; 3];
        for x in c.iter_mut().take(d) {
            // Skewed coordinates so fan-in varies.
            *x = (rng.gen::<f64>().powi(2) * domain as f64) as i64;
        }
        if seen.insert(c) {
            pts.push(Point { c, v: rng.gen_range(0..50) });
        }
    }
    Dataset { name: format!("syn{d}"), d, pts }
}

const ALL: Tails = Tails { group: true, semi: true, quantile: true };

fn check_layered(ds: &Dataset, cap: Option<usize>, branching: usize) {
    let mut sam: DryRunSam<Cell> = DryRunSam::new(AccessPolicy::MULTI_WRITE);
    let cfg = LayeredConfig { branching, tails: ALL, fanin_cap: cap };
    let mut l = Layered::build(ds, cfg, &mut sam).unwrap();
    let mut ctx = Ctx::new(branching).unwrap();
    for kind in [Kind::Group, Kind::Semi, Kind::Quantile] {
        let budget = l.analytic_budget(kind);
        for q in random_queries(ds, 150, 11) {
            let want = oracle(ds, &q);
            let r0 = reads(&sam);
            let got = l.query(&mut ctx, &mut sam, &q, kind).unwrap();
            let real = reads(&sam) - r0;
            assert!(real <= budget, "{kind:?}: {real} reads exceed budget {budget}");
            dummy_reads(&mut sam, budget - real).unwrap();
            assert_eq!(reads(&sam) - r0, budget, "padded trace length must be constant");
            match kind {
                Kind::Group => assert_eq!(got.g, want.group, "{q:?}"),
                Kind::Semi => assert_eq!(got.s, want.semi, "{q:?}"),
                Kind::Quantile => assert_eq!(got.median, want.median, "{q:?}"),
            }
        }
    }
}

#[test]
fn layered_1d_2d_3d_raw_and_capped() {
    for (d, n, dom) in [(1, 300, 2000), (2, 400, 60), (2, 300, 2000), (3, 300, 20), (3, 250, 400)] {
        let ds = synthetic(d, n, dom, d as u64 * 31 + n as u64);
        check_layered(&ds, None, 30);
        check_layered(&ds, None, 3);
        check_layered(&ds, Some(4), 30);
        check_layered(&ds, None, 4); // small fanout exercises deep alias trees
    }
}

#[test]
fn sparq_matches_oracle_within_schedule() {
    for (d, n, dom) in [(1, 300, 2000), (2, 400, 60), (3, 300, 20)] {
        let ds = synthetic(d, n, dom, 99 + d as u64);
        let mut sam: DryRunSam<Cell> = DryRunSam::new(AccessPolicy::MULTI_WRITE);
        let mut s = Sparq::build(&ds, 30, &mut sam).unwrap();
        let mut ctx = Ctx::new(30).unwrap();
        let budget = s.analytic_budget();
        for q in random_queries(&ds, 200, 5) {
            let want = oracle(&ds, &q);
            let r0 = reads(&sam);
            let got = s.query(&mut ctx, &mut sam, &q).unwrap();
            assert!(reads(&sam) - r0 <= budget);
            assert_eq!(got.g, want.group);
            assert_eq!(got.s, want.semi);
        }
        assert_eq!(s.level_overflows, 0, "schedule lemma violated");
    }
}
