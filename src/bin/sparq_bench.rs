//! Benchmark: build in a dry-run SAM (no crypto), optionally install the
//! snapshot into an encrypted Path OSAM+ tree, then run random range queries.
//!
//! Every query performs only its real reads and is then padded once, at the
//! end, to a public total (the whole trace length is what the server sees).
//!
//! ```text
//! sparq_bench --dataset datasets/cali-1024x1024.pts [--scheme layered|sparq|both]
//!     [--kinds group,semi,quantile] [--queries 100] [--seed 1] [--values hash16|dataset]
//!     [--crypto] [--block-size 128|256|512|1024] [--branching 30] [--fanin-cap K]
//!     [--pad analytic|none|N] [--stash-size 40] [--check] [--csv results.csv]
//! ```

use osam_plus::StashSize;
use sam_model::pointer::FixedSizeCodec;
use sam_model::pointer::{BalancedCell, BalancedCellValueCodec, ValueCodec};
use sam_model::{AccessPolicy, AccessStrategy, DryRunSam, PathOsamSam, SamError, SnapshotBlock};
use efficient_agro::ctx::{dummy_reads, reads, Cell, Ctx, Res, Sam};
use efficient_agro::data::{oracle, random_queries, Dataset, Query, Truth, ValueMode};
use efficient_agro::layered::{Kind, Layered, LayeredConfig, Tails};
use efficient_agro::obj::ObjCodec;
use efficient_agro::sparq::Sparq;
use std::io::Write;
use std::time::Instant;

#[derive(Clone, Debug)]
struct Config {
    dataset: String,
    schemes: Vec<String>,
    kinds: Vec<Kind>,
    queries: usize,
    seed: u64,
    values: ValueMode,
    crypto: bool,
    block_size: usize,
    branching: usize,
    fanin_cap: Option<usize>,
    pad: Pad,
    stash_size: StashSize,
    check: bool,
    csv: Option<String>,
    /// Report the largest encoded cells in a dry-run (copies the whole store).
    cell_sizes: bool,
}

#[derive(Clone, Copy, Debug)]
enum Pad {
    Analytic,
    None,
    Fixed(u64),
}

fn parse() -> Result<Config, String> {
    let mut c = Config {
        dataset: String::new(),
        schemes: vec!["layered".into(), "sparq".into()],
        kinds: vec![Kind::Group, Kind::Semi, Kind::Quantile],
        queries: 100,
        seed: 1,
        values: ValueMode::Hash16,
        crypto: false,
        block_size: 256,
        branching: 30,
        fanin_cap: None,
        pad: Pad::Analytic,
        stash_size: 40,
        check: false,
        cell_sizes: false,
        csv: None,
    };
    let args: Vec<String> = std::env::args().skip(1).collect();
    let mut i = 0;
    let val = |i: &mut usize| -> Result<String, String> {
        *i += 1;
        args.get(*i).cloned().ok_or_else(|| format!("missing value for {}", args[*i - 1]))
    };
    while i < args.len() {
        match args[i].as_str() {
            "--dataset" => c.dataset = val(&mut i)?,
            "--scheme" => {
                c.schemes = match val(&mut i)?.as_str() {
                    "both" => vec!["layered".into(), "sparq".into()],
                    s @ ("layered" | "sparq") => vec![s.into()],
                    s => return Err(format!("unknown scheme {s}")),
                }
            }
            "--kinds" => {
                c.kinds = val(&mut i)?.split(',').map(Kind::parse).collect::<Result<_, _>>()?;
            }
            "--queries" => c.queries = val(&mut i)?.parse().map_err(|_| "bad --queries")?,
            "--seed" => c.seed = val(&mut i)?.parse().map_err(|_| "bad --seed")?,
            "--values" => c.values = ValueMode::parse(&val(&mut i)?)?,
            "--crypto" => c.crypto = true,
            "--block-size" => c.block_size = val(&mut i)?.parse().map_err(|_| "bad --block-size")?,
            "--branching" => c.branching = val(&mut i)?.parse().map_err(|_| "bad --branching")?,
            "--fanin-cap" => c.fanin_cap = Some(val(&mut i)?.parse().map_err(|_| "bad --fanin-cap")?),
            "--pad" => {
                c.pad = match val(&mut i)?.as_str() {
                    "analytic" => Pad::Analytic,
                    "none" => Pad::None,
                    n => Pad::Fixed(n.parse().map_err(|_| "bad --pad")?),
                }
            }
            "--stash-size" => c.stash_size = val(&mut i)?.parse().map_err(|_| "bad --stash-size")?,
            "--check" => c.check = true,
            "--cell-sizes" => c.cell_sizes = true,
            "--csv" => c.csv = Some(val(&mut i)?),
            "-h" | "--help" => return Err("see the header of src/bin/sparq_bench.rs".into()),
            a => return Err(format!("unknown argument {a}")),
        }
        i += 1;
    }
    if c.dataset.is_empty() {
        return Err("--dataset is required".into());
    }
    Ok(c)
}

/// What one scheme needs from the benchmark loop.
trait Scheme {
    fn label(&self) -> String;
    /// Query kinds this scheme runs (one query may answer several).
    fn kinds(&self, requested: &[Kind]) -> Vec<Kind>;
    fn budget(&self, kind: Kind) -> u64;
    fn run<S: Sam>(&mut self, ctx: &mut Ctx, sam: &mut S, q: &Query, kind: Kind) -> Res<(Truth, u64)>;
    fn check(&self, kind: Kind, got: &Truth, want: &Truth) -> bool;
}

impl Scheme for Layered {
    fn label(&self) -> String {
        "layered".into()
    }
    fn kinds(&self, requested: &[Kind]) -> Vec<Kind> {
        requested.to_vec()
    }
    fn budget(&self, kind: Kind) -> u64 {
        self.analytic_budget(kind)
    }
    fn run<S: Sam>(&mut self, ctx: &mut Ctx, sam: &mut S, q: &Query, kind: Kind) -> Res<(Truth, u64)> {
        let o = self.query(ctx, sam, q, kind)?;
        Ok((
            Truth {
                group: o.g,
                semi: o.s,
                median: o.median,
            },
            o.rounds,
        ))
    }
    fn check(&self, kind: Kind, got: &Truth, want: &Truth) -> bool {
        match kind {
            Kind::Group => got.group == want.group,
            Kind::Semi => got.semi == want.semi,
            Kind::Quantile => got.median == want.median,
        }
    }
}

impl Scheme for Sparq {
    fn label(&self) -> String {
        "sparq".into()
    }
    fn kinds(&self, requested: &[Kind]) -> Vec<Kind> {
        // One traversal answers COUNT/SUM/SUMSQ and MIN/MAX together.
        if requested.iter().any(|k| *k != Kind::Quantile) {
            vec![Kind::Group]
        } else {
            vec![]
        }
    }
    fn budget(&self, _kind: Kind) -> u64 {
        self.analytic_budget()
    }
    fn run<S: Sam>(&mut self, ctx: &mut Ctx, sam: &mut S, q: &Query, _kind: Kind) -> Res<(Truth, u64)> {
        let o = self.query(ctx, sam, q)?;
        Ok((
            Truth {
                group: o.g,
                semi: o.s,
                median: None,
            },
            o.rounds,
        ))
    }
    fn check(&self, _kind: Kind, got: &Truth, want: &Truth) -> bool {
        got.group == want.group && got.semi == want.semi
    }
}

#[derive(Default)]
struct Series {
    real: Vec<u64>,
    /// Real reads + writes per query: the quantity OSAM+ leaks (a write
    /// reads a dummy path).
    ops: Vec<u64>,
    rounds: Vec<u64>,
    nanos: Vec<u128>,
    writes: u64,
    mismatches: usize,
    overflows: usize,
    padded_total: u64,
}

fn stats(v: &[u64]) -> (f64, f64, u64, u64) {
    let n = v.len().max(1) as f64;
    let mean = v.iter().sum::<u64>() as f64 / n;
    let var = v.iter().map(|&x| (x as f64 - mean).powi(2)).sum::<f64>() / n;
    (mean, var, *v.iter().min().unwrap_or(&0), *v.iter().max().unwrap_or(&0))
}

fn window_max(v: &[u64], k: usize) -> f64 {
    if v.len() < k {
        return f64::NAN;
    }
    v.windows(k).map(|w| w.iter().sum::<u64>()).max().unwrap() as f64 / k as f64
}

#[allow(clippy::too_many_arguments)]
fn run_queries<X: Scheme, S: Sam>(
    cfg: &Config,
    scheme: &mut X,
    sam: &mut S,
    queries: &[Query],
    truths: &[Option<Truth>],
    mode: &str,
    out: &mut dyn Write,
    csv: &mut Option<std::fs::File>,
    build_line: &str,
) -> Res<()> {
    let mut ctx = Ctx::new(cfg.branching)?;
    for kind in scheme.kinds(&cfg.kinds) {
        let budget = match cfg.pad {
            Pad::Analytic => Some(scheme.budget(kind)),
            Pad::None => None,
            Pad::Fixed(n) => Some(n),
        };
        let mut s = Series::default();
        sam.reset_stash_maximum();
        for (q, truth) in queries.iter().zip(truths) {
            let r0 = reads(sam);
            let w0 = sam.stats().operations.writes;
            let clock = Instant::now();
            let (got, rounds) = scheme.run(&mut ctx, sam, q, kind)?;
            let real = reads(sam) - r0;
            if let Some(t) = budget {
                if real > t {
                    s.overflows += 1;
                } else {
                    dummy_reads(sam, t - real)?;
                }
            }
            s.nanos.push(clock.elapsed().as_nanos());
            s.writes += sam.stats().operations.writes - w0;
            s.padded_total += reads(sam) - r0;
            s.real.push(real);
            s.ops.push(real + (sam.stats().operations.writes - w0));
            s.rounds.push(rounds);
            if let Some(t) = truth {
                if !scheme.check(kind, &got, t) {
                    s.mismatches += 1;
                    if s.mismatches <= 3 {
                        eprintln!("MISMATCH {} {:?}: got {:?} want {:?}", scheme.label(), kind, got, t);
                    }
                }
            }
        }
        let (mean, var, min, max) = stats(&s.real);
        let (rmean, _, _, rmax) = stats(&s.rounds);
        let (omean, ovar, omin, omax) = stats(&s.ops);
        let ms = s.nanos.iter().sum::<u128>() as f64 / 1e6 / s.real.len().max(1) as f64;
        let stash = sam.stats().stash.map(|x| x.maximum).unwrap_or(0);
        let label = if scheme.label() == "sparq" { "group+semi".to_string() } else { kind.label().to_string() };
        let line = format!(
            "result scheme={} kind={} mode={} queries={} budget={} real_mean={:.1} real_var={:.1} real_min={} real_max={} \
             win10_max={:.1} win100_max={:.1} padded_per_query={:.1} overflows={} rounds_mean={:.1} rounds_max={} \
             writes_per_query={:.1} ops_mean={:.1} ops_var={:.1} ops_min={} ops_max={} ops_win100_max={:.1} \
             ms_per_query={:.3} max_stash={} checked={} mismatches={}",
            scheme.label(),
            label,
            mode,
            s.real.len(),
            budget.map(|b| b.to_string()).unwrap_or("none".into()),
            mean,
            var,
            min,
            max,
            window_max(&s.real, 10),
            window_max(&s.real, 100),
            s.padded_total as f64 / s.real.len().max(1) as f64,
            s.overflows,
            rmean,
            rmax,
            s.writes as f64 / s.real.len().max(1) as f64,
            omean,
            ovar,
            omin,
            omax,
            window_max(&s.ops, 100),
            ms,
            stash,
            truths.iter().filter(|t| t.is_some()).count(),
            s.mismatches
        );
        writeln!(out, "{line}").ok();
        if let Some(f) = csv.as_mut() {
            let mut kv = std::collections::BTreeMap::new();
            for l in [build_line, line.as_str()] {
                for p in l.split_whitespace().skip(1) {
                    if let Some((k, v)) = p.split_once('=') {
                        let k = if l.starts_with("build") { format!("build_{k}") } else { k.to_string() };
                        kv.insert(k, v.to_string());
                    }
                }
            }
            let row: Vec<String> = CSV_KEYS.iter().map(|k| kv.get(*k).cloned().unwrap_or_default()).collect();
            writeln!(f, "{},{}", cfg_csv(cfg), row.join(",")).ok();
        }
    }
    Ok(())
}

const CSV_KEYS: &[&str] = &[
    "scheme", "kind", "mode", "queries", "budget", "real_mean", "real_var", "real_min", "real_max",
    "win10_max", "win100_max", "padded_per_query", "overflows", "rounds_mean", "rounds_max",
    "writes_per_query", "ops_mean", "ops_var", "ops_min", "ops_max", "ops_win100_max", "ms_per_query", "max_stash", "checked", "mismatches", "build_objects",
    "build_entries", "build_extra", "build_build_s", "build_max_fanin", "build_heights", "build_sigma",
];

fn cfg_csv(c: &Config) -> String {
    format!(
        "{},{},{},{},{},{},{}",
        c.dataset,
        c.values_label(),
        c.block_size,
        c.branching,
        c.fanin_cap.map(|k| k.to_string()).unwrap_or("raw".into()),
        c.seed,
        c.crypto
    )
}

impl Config {
    fn values_label(&self) -> &'static str {
        match self.values {
            ValueMode::Dataset => "dataset",
            ValueMode::Hash16 => "hash16",
        }
    }
}

fn install_and_run<X: Scheme, const B: usize>(
    cfg: &Config,
    scheme: &mut X,
    dry: DryRunSam<Cell>,
    queries: &[Query],
    truths: &[Option<Truth>],
    out: &mut dyn Write,
    csv: &mut Option<std::fs::File>,
    build_line: &str,
) -> Res<()> {
    let clock = Instant::now();
    let snapshot = dry.snapshot();
    drop(dry);
    cell_report(cfg, &scheme.label(), &snapshot.blocks, out)?;
    let allocated = snapshot.next_identifier.saturating_sub(1).max(2);
    let capacity = allocated
        .checked_next_power_of_two()
        .ok_or_else(|| SamError::Backend("capacity overflow".into()))?;
    let codec = FixedSizeCodec::new(BalancedCellValueCodec::new(ObjCodec::new(cfg.branching)));
    let mut enc = PathOsamSam::<Cell, _, B, 4, 1>::from_snapshot(
        snapshot,
        capacity,
        cfg.stash_size,
        true,
        AccessPolicy::MULTI_WRITE,
        AccessStrategy::MULTI_WRITE_RARY,
        true,
        cfg.seed,
        codec,
    )?
    .with_policy_checks(false);
    let install_s = clock.elapsed().as_secs_f64();
    let levels = capacity.ilog2(); // bucket levels = height + 1 = log2(capacity)
    writeln!(
        out,
        "install scheme={} capacity={} block_bytes={} bucket_levels={} path_bytes={} install_s={:.2} install_max_stash={}",
        scheme.label(),
        capacity,
        B,
        levels,
        levels as usize * 4 * B,
        install_s,
        enc.build_max_stash_occupancy()
    )
    .ok();
    run_queries(cfg, scheme, &mut enc, queries, truths, "crypto", out, csv, build_line)
}

/// Largest encoded cell of each kind (with the 4-byte block envelope). Fails
/// if any cell would not fit the configured block, so a dry-run already
/// proves that the encrypted install will work.
fn cell_report(cfg: &Config, label: &str, blocks: &[SnapshotBlock<Cell>], out: &mut dyn Write) -> Res<()> {
    let codec = BalancedCellValueCodec::new(ObjCodec::new(cfg.branching));
    let (mut root, mut node, mut kids) = (0usize, 0usize, 0usize);
    let mut buf = Vec::new();
    for block in blocks {
        buf.clear();
        codec.encode_value(&block.value, &mut buf)?;
        let n = buf.len() + 4;
        match &block.value {
            BalancedCell::Root { .. } => root = root.max(n),
            BalancedCell::Node { .. } => node = node.max(n),
            BalancedCell::Kids { .. } => kids = kids.max(n),
        }
    }
    writeln!(out, "cells scheme={label} block_bytes={} max_root={root} max_node={node} max_kids={kids}", cfg.block_size).ok();
    let worst = root.max(node).max(kids);
    if worst > cfg.block_size {
        return Err(SamError::Backend(format!(
            "a {label} cell needs {worst} bytes but blocks are {} bytes (lower --branching or raise --block-size)",
            cfg.block_size
        )));
    }
    Ok(())
}

fn dispatch<X: Scheme>(
    cfg: &Config,
    scheme: &mut X,
    dry: DryRunSam<Cell>,
    queries: &[Query],
    truths: &[Option<Truth>],
    out: &mut dyn Write,
    csv: &mut Option<std::fs::File>,
    build_line: &str,
) -> Res<()> {
    {
        // Server storage of the encrypted tree (Z = 4 blocks per bucket), as
        // the crypto install sizes it: capacity = next power of two >= the
        // addresses allocated during the build.
        let allocated = sam_model::SingleAccessMachine::stats(&dry).operations.allocations.max(2);
        let capacity = allocated.next_power_of_two();
        let levels = capacity.ilog2() as u64;
        writeln!(
            out,
            "storage scheme={} allocated={} capacity={} bucket_levels={} block_bytes={} server_bytes={} path_bytes={}",
            scheme.label(),
            allocated,
            capacity,
            levels,
            cfg.block_size,
            capacity * 4 * cfg.block_size as u64,
            levels * 4 * cfg.block_size as u64
        )
        .ok();
    }
    if cfg.cell_sizes && !cfg.crypto {
        cell_report(cfg, &scheme.label(), &dry.snapshot().blocks, out)?;
    }
    if !cfg.crypto {
        let mut dry = dry;
        return run_queries(cfg, scheme, &mut dry, queries, truths, "dry-run", out, csv, build_line);
    }
    match cfg.block_size {
        128 => install_and_run::<X, 128>(cfg, scheme, dry, queries, truths, out, csv, build_line),
        256 => install_and_run::<X, 256>(cfg, scheme, dry, queries, truths, out, csv, build_line),
        512 => install_and_run::<X, 512>(cfg, scheme, dry, queries, truths, out, csv, build_line),
        1024 => install_and_run::<X, 1024>(cfg, scheme, dry, queries, truths, out, csv, build_line),
        _ => Err(SamError::InvalidParameter("block size must be 128, 256, 512 or 1024")),
    }
}

fn main() {
    let cfg = match parse() {
        Ok(c) => c,
        Err(e) => {
            eprintln!("error: {e}");
            std::process::exit(2);
        }
    };
    if let Err(e) = real_main(&cfg) {
        eprintln!("error: {e}");
        std::process::exit(1);
    }
}

fn real_main(cfg: &Config) -> Result<(), Box<dyn std::error::Error>> {
    let ds = Dataset::load(&cfg.dataset, cfg.values)?;
    let stdout = std::io::stdout();
    let mut out = stdout.lock();
    let dims: Vec<String> = (0..ds.d).map(|i| ds.distinct(i).to_string()).collect();
    writeln!(
        out,
        "config dataset={} d={} N={} n={} values={} queries={} seed={} crypto={} block_size={} branching={} fanin_cap={} pad={:?}",
        ds.name,
        ds.d,
        ds.pts.len(),
        dims.join("x"),
        cfg.values_label(),
        cfg.queries,
        cfg.seed,
        cfg.crypto,
        cfg.block_size,
        cfg.branching,
        cfg.fanin_cap.map(|k| k.to_string()).unwrap_or("raw".into()),
        cfg.pad
    )?;
    let queries = random_queries(&ds, cfg.queries, cfg.seed);
    let truths: Vec<Option<Truth>> = queries
        .iter()
        .map(|q| cfg.check.then(|| oracle(&ds, q)))
        .collect();
    let mut csv = match &cfg.csv {
        None => None,
        Some(p) => {
            let new = !std::path::Path::new(p).exists();
            let mut f = std::fs::OpenOptions::new().create(true).append(true).open(p)?;
            if new {
                writeln!(f, "dataset,values,block_size,branching,fanin_cap,seed,crypto,{}", CSV_KEYS.join(","))?;
            }
            Some(f)
        }
    };
    for scheme in &cfg.schemes {
        let mut dry: DryRunSam<Cell> = DryRunSam::new(AccessPolicy::MULTI_WRITE);
        let clock = Instant::now();
        if scheme == "layered" {
            let tails = Tails {
                group: cfg.kinds.contains(&Kind::Group),
                semi: cfg.kinds.contains(&Kind::Semi),
                quantile: cfg.kinds.contains(&Kind::Quantile),
            };
            let lcfg = LayeredConfig {
                branching: cfg.branching,
                tails,
                fanin_cap: cfg.fanin_cap,
            };
            let mut l = Layered::build(&ds, lcfg, &mut dry)?;
            let st = l.stats.clone();
            let line = format!(
                "build scheme=layered objects={} entries={} extra=copies:{}/seg:{}/wav:{}/bst:{} build_s={:.2} \
                 max_fanin=casc:{}/dims:{}/wav:{} heights=skel:{}/bst:{}/seg:{}/wav:{} sigma={}",
                dry.live_blocks(),
                st.entries,
                st.copies,
                st.seg_nodes,
                st.wav_entries,
                st.bst_nodes,
                clock.elapsed().as_secs_f64(),
                st.max_fanin_casc,
                st.max_fanin_dim[..ds.d.saturating_sub(1)].iter().map(|h| h.to_string()).collect::<Vec<_>>().join("x"),
                st.max_fanin_wav,
                st.heights[..ds.d.saturating_sub(1)].iter().map(|h| h.to_string()).collect::<Vec<_>>().join("x"),
                st.bst_height,
                st.seg_height,
                st.wav_height,
                st.sigma
            );
            writeln!(out, "{line}")?;
            for k in &cfg.kinds {
                writeln!(out, "budget scheme=layered kind={} analytic_reads={}", k.label(), l.analytic_budget(*k))?;
            }
            dispatch(cfg, &mut l, dry, &queries, &truths, &mut out, &mut csv, &line)?;
        } else {
            let mut s = Sparq::build(&ds, cfg.branching, &mut dry)?;
            let line = format!(
                "build scheme=sparq objects={} entries={} extra=none build_s={:.2} max_fanin=1 heights={}",
                dry.live_blocks(),
                s.stats.nodes,
                clock.elapsed().as_secs_f64(),
                s.stats.heights[..ds.d].iter().map(|h| h.to_string()).collect::<Vec<_>>().join("x")
            );
            writeln!(out, "{line}")?;
            writeln!(out, "budget scheme=sparq kind=group+semi analytic_reads={}", s.analytic_budget())?;
            dispatch(cfg, &mut s, dry, &queries, &truths, &mut out, &mut csv, &line)?;
            writeln!(
                out,
                "sparq_schedule max_level_ratio={:.3} level_overflows={}",
                s.max_level_ratio, s.level_overflows
            )?;
        }
    }
    Ok(())
}
