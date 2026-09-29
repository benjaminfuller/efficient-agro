# efficient-agro

Oblivious multidimensional range aggregation on the Rust Path OSAM+ backend.
The crate depends on the local checkout `../rust_osam_plus` (sam-model
and osam_plus) through path dependencies, because the balanced r-ary pointer
lives there (`crates/sam-model/src/pointer/balanced.rs`). Two schemes, built on
the same backend so they can be compared directly:

* **SPARQ** (baseline): the paper's nested multidimensional segment tree,
  traversed breadth-first. Its schedule follows the schedule lemma: dimensions
  are processed in sequence, and level `l` of dimension `t` may touch at most
  `P_{t-1} * min(2^l, 4)` nodes, where `P_t = prod_{j<=t} (2 H_j - 1)`.
* **Layered** (cascaded rank tree): a skeleton over dimensions `0..d-2` whose
  nodes each own a *list* (the distinct last-dimension values of their
  points). Entries cascade to their predecessors in the children's lists.
  Three pluggable tails hang off the lists:
  * group (COUNT, SUM, SUMSQ, hence AVG and STD): prefix aggregates are inlined
    in the entries, and the canonical siblings' values sit in the parent entry,
    so a 2D query costs only the two boundary paths;
  * semigroup (MIN, MAX): a position segment tree walked bottom-up from the
    two boundary leaves;
  * quantile (MEDIAN): a wavelet tree over the value ranks, descended in
    lockstep across all canonical lists.

Every object sits behind a **balanced r-ary multi-write pointer**
(`sam_model::pointer::BalancedPointers`, described below), so a shared list entry with fan-in `f` costs
exactly `ceil(log_b f) + 1` reads per dereference, whatever happened before.
That makes every budget below a proven worst case with raw fan-in. The
`--fanin-cap K` sampling option is still there, but it is no longer needed.

## The balanced r-ary pointer (`sam_model::pointer::balanced`)

The aliases of a value are the leaves of the **complete b-ary tree fixed by
the alias count f**. All leaves are at depth `H = ceil(log_b f)`, and leaf `j`
follows the base-b digits of `j`. There are three cell types:

* `Root { value, count, height, down }`
* `Node { parent, group, index }`: one per alias and per internal node. It
  holds the parent's whole child group, as `(node, down)` entries.
* `Kids { group }`: the **downward pointer**, the child list of an internal
  node.

Each member's cell is determined by its parent and its group, so a
dereference can rewrite every off-path member without reading it.

* **Dereference:** read the H + 1 cells up the path, give the path fresh
  addresses, and rewrite each path node's group and child list in place.
  There is no restructuring, packing or splaying. It costs exactly H + 1
  reads (fewer only when the walk meets a root that is open in the client
  cache), and at most b + 1 writes per read.
* **Copy:** insert at position f, walking down through the child lists.
  **Delete:** move the last alias into the hole, then contract the root if it
  is left with one child. Both keep the tree complete and cost exactly 2H
  reads (padded), or 1 when H = 0.
* The backend implements sam-model's `SmartPointerBackend` and
  `CacheablePointerBackend`, so the library's `CachedPointers` client cache,
  the dry-run SAM and encrypted Path OSAM+ (`MULTI_WRITE_RARY`) all work
  unchanged. The cell codec (`BalancedCellValueCodec`) stores 4-byte
  addresses, so a node cell is 6 + 8b bytes. **The defaults are 256 B blocks
  and b = 30**, the largest fanout whose node cell (246 B) fits a 256 B block.
  (128 B blocks allow b = 14; 1 KB blocks allow b = 126.)
* Tests (`cargo test --lib balanced` in `crates/sam-model`): exact read counts for dereferences
  (b = 2..30, f = 1..952); 3,000 random copies, deletes, puts and gets per
  fanout, with exact costs and every alias checked; copies that grow the tree
  while its root is open in the client cache; and a random mix run through
  `CachedPointers`.

## Build and run

Rust >= 1.89 (the backend's MSRV). The first build fetches the backend from
GitHub.

```bash
cargo build --release        # needs ../rust_osam_plus next to this repo
python3 tools/convert_datasets.py datasets/original datasets   # once: pickles -> .pts (datasets/ is not tracked)
./target/release/sparq_bench --dataset datasets/cali-1024x1024.pts --check
./target/release/sparq_bench --dataset datasets/cali-1024x1024.pts --crypto --fanin-cap 30
./run_all.sh            # dry-run read counts: every dataset and tail, 1000 queries
./run_all.sh --crypto   # encrypted timing/stash: small datasets, group + SPARQ
                        # (+ semi on spitz and nh), 100 queries, 4 runs in parallel
                        # (overrides at the top of run_all.sh: QUERIES, KINDS,
                        # EXTRA_SEMI, LARGE, DATASETS, JOBS)
cargo test --release    # randomized 1D/2D/3D correctness + constant-trace tests
```

`datasets/original/` is a copy of `~/Research/agroramsse/datasets`. The
converter mirrors `agroramsse/utils.py`: duplicates are removed, and every
point gets value 1 unless the dataset stores values (cali). Queries follow
`agroramsse`'s `random_range`: `lo` is uniform over the data range of each
dimension, then `hi` is uniform in `[lo, max]`.

**Values.** `--values hash16` (the default) gives every point a
deterministic pseudo-random 16-bit value, so SUM, STD, MIN, MAX and MEDIAN are
non-trivial on every dataset. `--values dataset` uses the stored values, which
are mostly 1. Read counts don't depend on the values, except the median tail,
whose height is `ceil(log2 sigma)`.

## Crypto mode ("dry-run build, encrypted queries")

The structure is always built in `DryRunSam` (no cryptography). With
`--crypto`, the snapshot is installed into `PathOsamSam<.., B, Z=4, P=1>` with
encryption on, the r-ary access strategy (`read_multi_paths`, local writes)
and ordered eviction, exactly as `oblivious_graph_bench --crypto` does. The
queries then run encrypted. Block size `--block-size 128|256|512|1024`,
default 256. Objects use a compact encoding: zigzag varints for integers and
pointer identifiers, and a one-byte tag for the MIN/MAX identity. The largest
object cell is 96 B (a 3D list entry), SPARQ nodes are 50 B, and a pointer node
cell is 6 + 8b bytes (246 B at b = 30). Every build checks the largest encoded
cell of each kind against the block size (`cells ...` line), and fails if one
does not fit. This is automatic in `--crypto` runs; use `--cell-sizes` for a
dry-run. Each run also prints the encrypted tree's exact server storage
(`storage ... server_bytes=`).

## Obliviousness: padding the whole trace

With the r-ary access strategy (`MULTI_WRITE_RARY`), writes are local. They go
to the stash, and the P eviction paths that each read downloads drain them
(P = 1 here, much less than b). The server sees only reads, each of 1 + P
paths, plus public flushes during bulk builds. (Standard OSAM+, where every
write reads a dummy path, would leak reads + writes instead; the `ops_*`
columns report that total.) Each query performs its real reads and is then
padded **once, at the end**, with reads of fresh, never-written addresses, up
to a public total:

* `--pad analytic` (default): a worst-case read budget computed from the
  structure's public shape (tree heights, maximum fan-ins, sigma). It is
  printed as `budget ... analytic_reads=`. Any query that exceeds it is
  counted in `overflows=`.
* `--pad none` reports the real read counts only; `--pad N` uses a fixed total.

Each `result` line reports the real reads per query (`real_mean`,
`real_var`, `real_min`, `real_max`) and the worst average over windows of 10
and 100 consecutive queries (`win10_max`, `win100_max`). These feed the
trace-variance study. `rounds_*` is the longest chain of dependent reads, i.e.
the round trips if independent reads of a level were batched, excluding the
final padding batch. `ops_mean`, `ops_var`, `ops_max` and `ops_win100_max` give the
real reads + writes per query.

**Budgets are worst-case bounds.** `Layered::deref_bound(f)` is the balanced
pointer's exact cost, `ceil(log_b f) + 1`. The per-query budget charges every
cascade step of skeleton dimension j at the largest fan-in a dimension-j step
can reach (`max_fanin=casc:.../dims:...`). So no query can overflow, and
`overflows=` stays at 0. (The splay-based `RaryPointer` used before had no
such bound; see the theory note.)

## Results with the balanced pointer (b = 30, 256 B blocks, dry-run, 500 queries, raw fan-in, all answers checked)

Proven per-query read budget, with the largest real read count in
parentheses (`results/b30/`):

| dataset | SPARQ (group+semi) | layered group | layered semi | layered median | group vs SPARQ |
|---|---|---|---|---|---|
| amazon-books 1D | 51 | 30 (29) | 78 (77) | 82 (78) | 1.7x |
| spitz 2D | 780 | 144 (87) | 1056 (290) | 1232 (344) | 5.4x |
| cali 2D | 780 | 144 (90) | 1056 (427) | 2112 (448) | 5.4x |
| gowalla 2D-50K | 1036 | 172 (97) | 1460 (384) | 2812 (499) | 6.0x |
| gowalla 2D-100K | 1326 | 238 (116) | 1790 (535) | — | 5.6x |
| nh 3D | 3059 | 1122 (340) | 6834 (850) | 11538 (1250) | 2.7x |
| gowalla 3D | 4181 | 1338 (216) | 7890 (1222) | — | 3.1x |

There were no overflows and no wrong answers.

**Storage and bandwidth** (group tail only; both schemes at 256 B blocks). Storage
is the encrypted tree's `server_bytes`. Bandwidth per query uses the SPARQ
paper's accounting, padded reads x tree levels x block size. The "paper"
columns are the SPARQ paper's Path ORAM table (256 B blocks).

| dataset | storage: SPARQ paper / our SPARQ / layered | bandwidth per query: SPARQ paper / Demertzis (paper) / layered |
|---|---|---|
| Books 1D | 17 MB / 17 MB / 34 MB | 171 KB / 190 KB / 115 KB |
| Spitz 2D | 69 MB / 34 MB / 34 MB | 3.4 MB / 408 KB / 553 KB |
| cali 2D | 560 MB / 268 MB / 268 MB | 4.2 MB / 408 KB / 664 KB |
| gowalla 2D-50K | 1.1 GB / 537 MB / 537 MB | 5.7 MB / 1.1 MB / 837 KB |
| gowalla 2D-100K | 2.2 GB / 1.07 GB / 1.07 GB | 8.1 MB / 1.4 MB / 1.22 MB |
| nh 3D | 278 MB / 134 MB / 268 MB | 17 MB / 665 KB / 5.2 MB |
| gowalla 3D | 1.1 GB / 537 MB / 537 MB | 18 MB / 816 KB / 6.5 MB |

Encrypted Path OSAM+ (20 queries, 256 B blocks; `results/b30/*.crypto.log`):
all answers correct, no overflows. On amazon-books the maximum stash was
3,619 blocks (1,224 for SPARQ). On spitz, a group query took 0.73 s
against 2.65 s for SPARQ (stash 1,756 vs 1,249). The median tails of
gowalla 3D and gowalla 2D-100K need more than the VM's 6 GB. The b = 64 / 1 KB
runs are in `results/b64/`.

## Earlier results with the splay-based `RaryPointer` (dry-run, 500 queries, b = 30)

Padded reads per query, with the (mean / max) real reads in parentheses. Every
run checked every answer against a brute-force scan: 0 mismatches, 0 overflows.

| dataset | SPARQ (group+semi) | layered group, raw fan-in | layered group, cap 30 | layered semi, cap 30 | layered median, cap 30 |
|---|---|---|---|---|---|
| amazon-books 1D | 51 (43/51) | 30 | 30 (27/29) | 78 (68/77) | 82 (70/78) |
| spitz 2D | 780 (66/172) | 344 | 104 (50/86) | 976 (97/288) | 1152 (119/343) |
| cali 2D | 780 (115/269) | 424 | 104 (57/88) | 976 (173/428) | 5112 (190/449) |
| gowalla 2D-50K | 1036 (59/171) | 508 | 124 (40/94) | 1364 (82/372) | — |
| gowalla 2D-100K | 1326 (89/239) | 550 | 134 (52/111) | 1582 (119/532) | — |
| nh 3D | 3059 (98/335) | 3090 | 786 (102/339) | 6162 (161/850) | 10866 (218/1250) |
| gowalla 3D | 4181 (63/525) | 3538 | 898 (54/216) | 7058 (119/1220) | 55442 (137/1505) |

Stored objects (the SAM blocks, including r-ary alias nodes) with only the
group tail, compared with SPARQ's nodes: amazon 23.5K vs 15.7K, spitz 23.5K vs
20.3K, cali 174K vs 182K, gowalla 2D-50K 322K vs 371K, gowalla 2D-100K 837K vs
904K, nh 126K vs 127K, gowalla 3D 274K vs 410K.

Takeaways:

* The group tail pads to 7.5–10x fewer reads than SPARQ in 2D and 3.9–4.7x fewer
  in 3D, at equal or smaller storage. Raw fan-in loses most of that gain to the
  loose alias-depth bound; cap 30 costs at most 0.1% extra entries.
* Mean real reads are similar for both schemes. The difference is the
  worst case: SPARQ's padded schedule is 7–18x its mean in 2D and 31–66x in
  3D, while the group tail's is 2–3x in 2D and 8–17x in 3D.
* The semigroup and median tails pad to more reads than SPARQ's combined
  group+semi schedule, and the median tail's storage is large (the wavelet
  entries dominate: 14.9M for gowalla 3D). For MIN/MAX in more than one
  dimension, SPARQ is still the better choice.
* The analytic budgets are tight in 2D (the worst real group query uses 80–85%
  of it) but loose in 3D (24–43%). A finer count of the 3D sub-walks would
  lower them.

### Reads + writes (only relevant under standard OSAM+)

If the scheme ran on standard OSAM+, where every write costs a dummy path, the
leaked total would be R + W. The layered scheme's r-ary alias rebuilds write
2–3x as much as they read, while SPARQ writes once per read. With a fan-in cap
K <= b, (1 + (K+1)/2) x the read budget is then a provable R + W budget.
Measurements are in `results/ops/`: the `ops_*` fields, and the b/K sweep in
`sweep.sh`. On the r-ary backend used here, the read budgets above are the
relevant ones.

Encrypted runs were verified on amazon-books (20 queries, all correct). On
this 2-core VM, install costs about 2 ms per object and a query about
130 ms, and the stash holds about 1,500 blocks with P = 1. The 1M-point
datasets need the server.

## Layout

```
src/agg.rs        Group (COUNT/SUM/SUMSQ) and Semi (MIN/MAX) algebras
src/data.rs       .pts loader, value modes, query generator, brute-force oracle
src/obj.rs        SAM objects (entries, BST, segment, wavelet, SPARQ nodes) + block codec
src/ctx.rs        open/share/release handles over sam-model's CachedPointers; padding
src/layered.rs    cascaded rank tree: plain build, sampling, materialization, queries, budgets
src/sparq.rs      SPARQ baseline and its schedule
src/probe.rs      r-ary alias-depth probes (src/bin/probe.rs, probe2.rs)
src/bin/sparq_bench.rs   benchmark driver (key=value lines, --csv output)
tests/correctness.rs     randomized tests vs the oracle, constant padded trace
tools/convert_datasets.py   agroramsse datasets -> .pts
run_all.sh        sweep script
results/          logs and CSVs from the runs above
```

## Trace-variance study (theory/)

`theory/trace_variance.tex` (and the compiled `.pdf`) collects the results so far on how much r-ary alias
trees stretch an OSAM+ trace. `src/aliasmodel.rs` is a standalone model of the
`RaryPointer` rebuild rule. It matches the library read-for-read on 150,000
dereferences (`validate_model`) and has an optional randomized layout.

* Proved, and checked on 648,000 dereferences (`lemmas`): at most b items
  are ever pending in a rebuild, and each dereference deepens any alias by at
  most one level. Objects with f <= b always cost exactly 2 reads.
* Proved: against an oblivious adversary with randomized layout, the total
  trace exceeds its expectation by at most sqrt(s/2 * sum_o r_o^2) except with
  probability e^-s. The ranges r_o are computed from the schedule.
* Conditional (bounded single-step influence): an excess of O(c sqrt(sK)).
  Supported by `concentration`, where sd(T_k)/sqrt(k) stays in 0.6–1.3, and by
  `influence2`.
* A simple rank potential does not certify an amortized bound for the
  deterministic r-ary layout (`potential`).
* r-ary path halving (Section 7, `access_halving` in `src/aliasmodel.rs`): an r-ary
  splay rule with a proven deterministic bound, reads <= 2 log2 f + 1 + Phi - Phi'.
  It keeps every property of the current rule: at most b items pending, depth
  growth <= 1, and f <= b costs exactly 2 reads. Checked on 180K dereferences
  (`halving`). Against schedules and a greedy adaptive adversary
  (`halving_adv`), the amortized cost stays at about 2 log_b f + 1, which is
  stated as a conjecture. It beats the current rule on mean and worst-case
  reads, except on bursts.
* The model also counts writes, matched against the library (`validate_model`).
* Binary multi-write pointer (`potential_binary`): proved a deterministic
  access lemma, reads <= 2 log2 f + 1 + Phi - Phi'. Any k dereferences, even
  chosen adaptively, then cost at most k(2 log2 f + 1) + 2f. It holds on 126,000 library
  dereferences with at least 1 read of slack on every step.
