//! Cascaded rank tree ("layered range tree") with three pluggable tails.
//!
//! Skeleton (dimensions 0..d-2): nested balanced trees over the distinct
//! coordinate values. Every tree node owns one *list*: the distinct
//! last-dimension values of its points, with a -inf sentinel at rank 0. Each
//! list entry cascades to the entries of its children's lists that hold the
//! predecessor of its value (fractional cascading). Every list has exactly
//! one parent list, so the lists form a tree; a list entry is shared by every
//! parent entry whose predecessor it is (its fan-in), and those parents hold
//! r-ary multi-write aliases of it.
//!
//! Tails, attached to every last-dimension list and entered at the two
//! cascaded positions (e1 excluded, e2 included):
//! * group (COUNT/SUM/SUMSQ): prefix aggregates inlined in the entries, and
//!   in the parent entries for canonical siblings, so they cost no reads;
//! * semigroup (MIN/MAX): a position segment tree over the list, walked
//!   bottom-up from the two leaves (each internal node has fan-in <= 2);
//! * quantile (MEDIAN): a wavelet tree over the value ranks of the list's
//!   points, descended in lockstep across all canonical lists.
//!
//! Queries do only their real reads; the trace is padded once, at the end,
//! to a public total (see [`Layered::analytic_budget`]).

use crate::agg::{Group, Semi};
use crate::ctx::{Ctx, Res, Sam, H};
use crate::data::{Dataset, Query};
use crate::obj::{Bst, Casc, Entry, Kids, Obj, Ptr, Seg, Wav};
use sam_model::pointer::{BalancedPointer, BalancedPointers};
use sam_model::pointer::CachedPointer;
use sam_model::SamError;
use std::collections::HashMap;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Kind {
    Group,
    Semi,
    Quantile,
}

impl Kind {
    pub fn label(self) -> &'static str {
        match self {
            Kind::Group => "group",
            Kind::Semi => "semi",
            Kind::Quantile => "quantile",
        }
    }
    pub fn parse(s: &str) -> Result<Self, String> {
        match s {
            "group" => Ok(Kind::Group),
            "semi" => Ok(Kind::Semi),
            "quantile" | "median" => Ok(Kind::Quantile),
            _ => Err(format!("unknown query kind {s}")),
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct Tails {
    pub group: bool,
    pub semi: bool,
    pub quantile: bool,
}

#[derive(Clone, Debug)]
pub struct LayeredConfig {
    pub branching: usize,
    pub tails: Tails,
    /// Sampling cap on list-entry fan-in (None = raw fan-in).
    pub fanin_cap: Option<usize>,
}

#[derive(Clone, Copy, Debug)]
struct CascSpec {
    key: i64,
    kids: Option<(usize, usize)>,
}

#[derive(Clone, Debug, Default)]
struct PList {
    vals: Vec<i64>,
    own: Vec<Group>,
    semi: Vec<Semi>,
    /// Quantile tail: value ranks of the list's points in (last coord, value)
    /// order.
    seq: Vec<u32>,
    casc: [Option<CascSpec>; 2],
    parent: Option<usize>,
    fanin: Vec<u32>,
}

impl PList {
    fn kids(&self) -> Vec<usize> {
        let mut out = Vec::new();
        for c in self.casc.iter().flatten() {
            if let Some((l, r)) = c.kids {
                out.push(l);
                out.push(r);
            }
        }
        out
    }
    fn pred(&self, y: i64) -> usize {
        self.vals.partition_point(|&v| v <= y) - 1
    }
}

/// Build-time statistics and public schedule parameters.
#[derive(Clone, Debug, Default)]
pub struct LayeredStats {
    pub lists: usize,
    pub entries: usize,
    pub copies: usize,
    pub seg_nodes: usize,
    pub wav_entries: usize,
    pub bst_nodes: usize,
    /// Maximum tree height (edges) per skeleton dimension.
    pub heights: [usize; 2],
    pub bst_height: usize,
    /// Maximum log2 of a segment-tree leaf count.
    pub seg_height: usize,
    pub max_fanin_casc: u32,
    /// Largest alias count of an entry that a dimension-j cascade step can
    /// reach (the entries of lists that are dimension-j children).
    pub max_fanin_dim: [u32; 2],
    pub max_fanin_wav: u32,
    pub sigma: usize,
    pub wav_height: usize,
}

pub struct Layered {
    pub d: usize,
    pub cfg: LayeredConfig,
    pub bst_root: Ptr,
    pub alphabet: Vec<i64>,
    pub stats: LayeredStats,
}

// ---------------------------------------------------------------------------
// Plain (plaintext) construction.
// ---------------------------------------------------------------------------

struct Plain<'a> {
    ds: &'a Dataset,
    d: usize,
    lists: Vec<PList>,
    vrank: HashMap<i64, u32>,
    need_seq: bool,
    heights: [usize; 2],
}

impl Plain<'_> {
    fn make_list(&mut self, mut idx: Vec<u32>) -> usize {
        let z = self.d - 1;
        let pts = &self.ds.pts;
        idx.sort_unstable_by_key(|&i| (pts[i as usize].c[z], pts[i as usize].v));
        let mut l = PList {
            vals: vec![i64::MIN],
            own: vec![Group::ZERO],
            semi: vec![Semi::IDENTITY],
            ..Default::default()
        };
        let mut acc = Group::ZERO;
        for &i in &idx {
            let p = &pts[i as usize];
            acc = acc.add(&Group::of(p.v));
            if *l.vals.last().unwrap() == p.c[z] {
                *l.own.last_mut().unwrap() = acc;
                let s = l.semi.last_mut().unwrap();
                *s = s.merge(&Semi::of(p.v));
            } else {
                l.vals.push(p.c[z]);
                l.own.push(acc);
                l.semi.push(Semi::of(p.v));
            }
            if self.need_seq {
                l.seq.push(self.vrank[&p.v]);
            }
        }
        self.lists.push(l);
        self.lists.len() - 1
    }

    fn build_struct(&mut self, t: usize, mut idx: Vec<u32>) -> usize {
        if t == self.d - 1 {
            return self.make_list(idx);
        }
        let pts = &self.ds.pts;
        idx.sort_unstable_by_key(|&i| pts[i as usize].c[t]);
        let mut values = Vec::new();
        let mut starts = Vec::new();
        for (k, &i) in idx.iter().enumerate() {
            let v = pts[i as usize].c[t];
            if values.last() != Some(&v) {
                values.push(v);
                starts.push(k);
            }
        }
        starts.push(idx.len());
        self.build_node(t, &idx, &values, &starts, 0, values.len(), 0)
    }

    #[allow(clippy::too_many_arguments)]
    fn build_node(
        &mut self,
        t: usize,
        idx: &[u32],
        values: &[i64],
        starts: &[usize],
        lo: usize,
        hi: usize,
        depth: usize,
    ) -> usize {
        let r = self.build_struct(t + 1, idx[starts[lo]..starts[hi]].to_vec());
        self.heights[t] = self.heights[t].max(depth);
        if hi - lo == 1 {
            self.lists[r].casc[t] = Some(CascSpec {
                key: values[lo],
                kids: None,
            });
        } else {
            let mid = lo + (hi - lo) / 2;
            let l = self.build_node(t, idx, values, starts, lo, mid, depth + 1);
            let rr = self.build_node(t, idx, values, starts, mid, hi, depth + 1);
            self.lists[r].casc[t] = Some(CascSpec {
                key: values[mid],
                kids: Some((l, rr)),
            });
            self.lists[l].parent = Some(r);
            self.lists[rr].parent = Some(r);
        }
        r
    }

    /// Fractional-cascading sampling: insert copies of parent values so no
    /// run of more than `k - 1` parent entries maps to one child entry.
    fn augment(&mut self, root: usize, k: usize) -> usize {
        let mut copies = 0;
        let mut stack = vec![root];
        while let Some(p) = stack.pop() {
            for c in self.lists[p].kids() {
                let pv = self.lists[p].vals.clone();
                let child = &self.lists[c];
                let (mut vals, mut own, mut semi) = (Vec::new(), Vec::new(), Vec::new());
                let mut ci = 0;
                let mut run = 0;
                for &y in &pv {
                    if ci < child.vals.len() && child.vals[ci] == y {
                        vals.push(y);
                        own.push(child.own[ci]);
                        semi.push(child.semi[ci]);
                        ci += 1;
                        run = 0;
                    } else {
                        run += 1;
                        if run == k {
                            vals.push(y);
                            own.push(*own.last().unwrap());
                            semi.push(Semi::IDENTITY);
                            copies += 1;
                            run = 0;
                        }
                    }
                }
                assert_eq!(ci, child.vals.len(), "child values must be a subset of the parent's");
                let child = &mut self.lists[c];
                child.vals = vals;
                child.own = own;
                child.semi = semi;
                stack.push(c);
            }
        }
        copies
    }

    fn compute_fanin(&mut self, root: usize) {
        for id in 0..self.lists.len() {
            let n = self.lists[id].vals.len();
            self.lists[id].fanin = vec![0; n];
        }
        for p in 0..self.lists.len() {
            for c in self.lists[p].kids() {
                let pv = self.lists[p].vals.clone();
                let child = &mut self.lists[c];
                for &y in &pv {
                    let r = child.vals.partition_point(|&v| v <= y) - 1;
                    child.fanin[r] += 1;
                }
            }
        }
        for f in self.lists[root].fanin.iter_mut() {
            *f = 1; // one BST leaf per root-list entry
        }
    }
}

// ---------------------------------------------------------------------------
// Materialization into the (dry-run) SAM.
// ---------------------------------------------------------------------------

fn install<S: Sam>(sam: &mut S, value: Obj, count: u32, b: usize) -> Res<Vec<BalancedPointer>> {
    BalancedPointers::new(b)?.install(sam, value, count.max(1) as usize)
}

fn wrap(p: BalancedPointer) -> Ptr {
    CachedPointer::from_raw(p)
}

struct Mat<'a, S> {
    sam: &'a mut S,
    lists: &'a [PList],
    d: usize,
    b: usize,
    tails: Tails,
    hs: usize,
    pools: Vec<Option<Vec<Vec<BalancedPointer>>>>,
    stats: LayeredStats,
}

impl<S: Sam> Mat<'_, S> {
    fn materialize(&mut self, id: usize) -> Res<()> {
        // Post-order over the list tree (children first). Iterative to keep
        // the stack shallow on deep trees.
        let mut order = Vec::new();
        let mut stack = vec![(id, false)];
        while let Some((l, done)) = stack.pop() {
            if done {
                order.push(l);
            } else {
                stack.push((l, true));
                for c in self.lists[l].kids() {
                    stack.push((c, false));
                }
            }
        }
        for l in order {
            self.materialize_list(l)?;
        }
        Ok(())
    }

    fn materialize_list(&mut self, id: usize) -> Res<()> {
        let lists = self.lists;
        let list = &lists[id];
        let n = list.vals.len();
        self.stats.lists += 1;
        self.stats.entries += n;

        // Semigroup tail: position segment tree, created top-down so every
        // node can be given an alias of its parent.
        let mut seg_leaf: Vec<Option<BalancedPointer>> = vec![None; n];
        if self.tails.semi && n >= 2 {
            let m = n.next_power_of_two();
            let lg = m.trailing_zeros() as usize;
            self.stats.seg_height = self.stats.seg_height.max(lg);
            let mut agg = vec![Semi::IDENTITY; 2 * m];
            for i in 0..n {
                agg[m + i] = list.semi[i];
            }
            for h in (1..m).rev() {
                agg[h] = agg[2 * h].merge(&agg[2 * h + 1]);
            }
            // A heap node exists iff its leftmost leaf is a real position.
            let exists = |h: usize| {
                let depth = (usize::BITS - 1 - h.leading_zeros()) as usize;
                ((h << (lg - depth)) - m) < n
            };
            let mut slot: Vec<Option<BalancedPointer>> = vec![None; 2 * m];
            for h in 1..m {
                if !exists(h) {
                    continue;
                }
                let kids: Vec<usize> = [2 * h, 2 * h + 1].into_iter().filter(|&c| exists(c)).collect();
                let value = Obj::Seg(Box::new(Seg {
                    left: agg[2 * h],
                    right: agg[2 * h + 1],
                    parent: slot[h].take().map(wrap),
                    is_left: h % 2 == 0,
                }));
                let aliases = install(self.sam, value, kids.len() as u32, self.b)?;
                self.stats.seg_nodes += 1;
                for (c, a) in kids.into_iter().zip(aliases) {
                    slot[c] = Some(a);
                }
            }
            for (i, leaf) in seg_leaf.iter_mut().enumerate() {
                *leaf = slot[m + i].take();
            }
        }

        // Quantile tail: wavelet position entries reachable from the entries.
        let mut wav_root: HashMap<u32, Vec<BalancedPointer>> = HashMap::new();
        if self.tails.quantile {
            let mut fin: HashMap<u32, u32> = HashMap::new();
            for g in &list.own {
                *fin.entry(g.cnt as u32).or_default() += 1;
            }
            let mut positions: Vec<(u32, u32)> = fin.into_iter().collect();
            positions.sort_unstable();
            let seq = list.seq.clone();
            wav_root = self.build_wav(0, &seq, &positions)?;
        }

        // Entries.
        let mut pool: Vec<Vec<BalancedPointer>> = Vec::with_capacity(n);
        for k in 0..n {
            let y = list.vals[k];
            let mut casc: [Option<Casc>; 2] = [None, None];
            for (j, spec) in list.casc.iter().enumerate() {
                let Some(spec) = spec else { continue };
                let kids = match spec.kids {
                    None => None,
                    Some((lc, rc)) => {
                        let (ll, rl) = (&lists[lc], &lists[rc]);
                        let (rank_l, rank_r) = (ll.pred(y), rl.pred(y));
                        let pre = (self.tails.group && j + 2 == self.d)
                            .then(|| (ll.own[rank_l], rl.own[rank_r]));
                        let l = self.pools[lc].as_mut().unwrap()[rank_l]
                            .pop()
                            .ok_or(SamError::Backend("cascade alias pool exhausted".into()))?;
                        let r = self.pools[rc].as_mut().unwrap()[rank_r]
                            .pop()
                            .ok_or(SamError::Backend("cascade alias pool exhausted".into()))?;
                        Some(Kids {
                            l: wrap(l),
                            r: wrap(r),
                            rank_l: rank_l as u32,
                            rank_r: rank_r as u32,
                            pre,
                        })
                    }
                };
                casc[j] = Some(Casc { key: spec.key, kids });
            }
            let wav = if self.tails.quantile {
                let cnt = list.own[k].cnt as u32;
                Some(wrap(wav_root.get_mut(&cnt).and_then(|v| v.pop()).ok_or(
                    SamError::Backend("wavelet root alias pool exhausted".into()),
                )?))
            } else {
                None
            };
            let value = Obj::Entry(Box::new(Entry {
                rank: k as u32,
                own: list.own[k],
                semi: list.semi[k],
                casc,
                seg: seg_leaf[k].take().map(wrap),
                wav,
            }));
            pool.push(install(self.sam, value, list.fanin[k], self.b)?);
        }
        // Children's pools are fully consumed now.
        for c in list.kids() {
            if let Some(p) = self.pools[c].take() {
                debug_assert!(p.iter().all(|v| v.is_empty()));
            }
        }
        self.pools[id] = Some(pool);
        Ok(())
    }

    /// Builds the wavelet position entries of one node reachable from
    /// `positions` (position, fan-in); returns their aliases by position.
    fn build_wav(
        &mut self,
        level: usize,
        seq: &[u32],
        positions: &[(u32, u32)],
    ) -> Res<HashMap<u32, Vec<BalancedPointer>>> {
        let shift = self.hs - 1 - level;
        let mut zeros = Vec::with_capacity(seq.len() + 1);
        zeros.push(0u32);
        let (mut left, mut right) = (Vec::new(), Vec::new());
        for &s in seq {
            let bit = (s >> shift) & 1;
            zeros.push(zeros.last().unwrap() + (1 - bit));
            if bit == 0 {
                left.push(s);
            } else {
                right.push(s);
            }
        }
        let last = level + 1 == self.hs;
        let mut lpool = HashMap::new();
        let mut rpool = HashMap::new();
        if !last {
            let mut lf: HashMap<u32, u32> = HashMap::new();
            let mut rf: HashMap<u32, u32> = HashMap::new();
            for &(i, _) in positions {
                let z = zeros[i as usize];
                *lf.entry(z).or_default() += 1;
                *rf.entry(i - z).or_default() += 1;
            }
            if !left.is_empty() {
                let mut lp: Vec<_> = lf.into_iter().collect();
                lp.sort_unstable();
                lpool = self.build_wav(level + 1, &left, &lp)?;
            }
            if !right.is_empty() {
                let mut rp: Vec<_> = rf.into_iter().collect();
                rp.sort_unstable();
                rpool = self.build_wav(level + 1, &right, &rp)?;
            }
        }
        let mut out = HashMap::with_capacity(positions.len());
        for &(i, fin) in positions {
            let z = zeros[i as usize];
            let o = i - z;
            let l = lpool.get_mut(&z).and_then(|v: &mut Vec<BalancedPointer>| v.pop()).map(wrap);
            let r = rpool.get_mut(&o).and_then(|v: &mut Vec<BalancedPointer>| v.pop()).map(wrap);
            self.stats.max_fanin_wav = self.stats.max_fanin_wav.max(fin);
            let value = Obj::Wav(Box::new(Wav { zeros: z, ones: o, l, r }));
            out.insert(i, install(self.sam, value, fin, self.b)?);
            self.stats.wav_entries += 1;
        }
        Ok(out)
    }

    fn build_bst(&mut self, vals: &[i64], pool: &mut [Vec<BalancedPointer>], lo: usize, hi: usize, depth: usize) -> Res<BalancedPointer> {
        self.stats.bst_nodes += 1;
        self.stats.bst_height = self.stats.bst_height.max(depth);
        let value = if hi - lo == 1 {
            Obj::Bst(Box::new(Bst {
                key: vals[lo],
                l: None,
                r: None,
                leaf_rank: lo as u32,
                entry: Some(wrap(pool[lo].pop().ok_or(SamError::Backend("root pool".into()))?)),
            }))
        } else {
            let mid = lo + (hi - lo) / 2;
            let l = self.build_bst(vals, pool, lo, mid, depth + 1)?;
            let r = self.build_bst(vals, pool, mid, hi, depth + 1)?;
            Obj::Bst(Box::new(Bst {
                key: vals[mid],
                l: Some(wrap(l)),
                r: Some(wrap(r)),
                leaf_rank: 0,
                entry: None,
            }))
        };
        Ok(install(self.sam, value, 1, self.b)?.remove(0))
    }
}

fn ceil_log2(x: usize) -> usize {
    if x <= 1 {
        0
    } else {
        (usize::BITS - (x - 1).leading_zeros()) as usize
    }
}

impl Layered {
    /// Builds the structure into `sam` (normally a `DryRunSam`).
    pub fn build<S: Sam>(ds: &Dataset, cfg: LayeredConfig, sam: &mut S) -> Res<Self> {
        let d = ds.d;
        let mut alphabet: Vec<i64> = ds.pts.iter().map(|p| p.v).collect();
        alphabet.sort_unstable();
        alphabet.dedup();
        let vrank: HashMap<i64, u32> = alphabet.iter().enumerate().map(|(i, &v)| (v, i as u32)).collect();
        let mut plain = Plain {
            ds,
            d,
            lists: Vec::new(),
            vrank,
            need_seq: cfg.tails.quantile,
            heights: [0; 2],
        };
        let root = plain.build_struct(0, (0..ds.pts.len() as u32).collect());
        let copies = match cfg.fanin_cap {
            Some(k) if d >= 2 => plain.augment(root, k),
            _ => 0,
        };
        plain.compute_fanin(root);
        let max_fanin_casc = plain
            .lists
            .iter()
            .enumerate()
            .filter(|(i, _)| *i != root)
            .flat_map(|(_, l)| l.fanin.iter().copied())
            .max()
            .unwrap_or(1);
        let mut max_fanin_dim = [1u32; 2];
        for l in &plain.lists {
            for (j, c) in l.casc.iter().enumerate() {
                if let Some((a, b)) = c.as_ref().and_then(|c| c.kids) {
                    for k in [a, b] {
                        let m = plain.lists[k].fanin.iter().copied().max().unwrap_or(1);
                        max_fanin_dim[j] = max_fanin_dim[j].max(m);
                    }
                }
            }
        }
        let hs = ceil_log2(alphabet.len()).max(1);
        let lists = std::mem::take(&mut plain.lists);
        let mut mat = Mat {
            sam,
            lists: &lists,
            d,
            b: cfg.branching,
            tails: cfg.tails,
            hs,
            pools: vec![None; lists.len()],
            stats: LayeredStats::default(),
        };
        mat.materialize(root)?;
        let mut pool = mat.pools[root].take().unwrap();
        let root_vals = lists[root].vals.clone();
        let bst = mat.build_bst(&root_vals, &mut pool, 0, root_vals.len(), 0)?;
        let mut stats = mat.stats;
        stats.copies = copies;
        stats.heights = plain.heights;
        stats.max_fanin_casc = max_fanin_casc;
        stats.max_fanin_dim = max_fanin_dim;
        stats.sigma = alphabet.len();
        stats.wav_height = hs;
        Ok(Self {
            d,
            cfg,
            bst_root: wrap(bst),
            alphabet,
            stats,
        })
    }

    /// Worst-case reads for one dereference of an object with `f` aliases:
    /// exactly `ceil(log_b f) + 1` for the balanced pointer (`sam_model::pointer::BalancedPointers`),
    /// whatever the history, and fewer only on a cache hit.
    pub fn deref_bound(&self, f: u32) -> u64 {
        sam_model::pointer::balanced_deref_reads(f.max(1), self.cfg.branching)
    }

    /// Public per-query read budget for `kind`, from the structure's shape.
    pub fn analytic_budget(&self, kind: Kind) -> u64 {
        let st = &self.stats;
        let rcs = [self.deref_bound(st.max_fanin_dim[0]), self.deref_bound(st.max_fanin_dim[1])];
        let rw = self.deref_bound(st.max_fanin_wav);
        let rs = self.deref_bound(2);
        let bst = 2 * (st.bst_height as u64 + 1) + 2;
        let tail = match kind {
            Kind::Group => 0,
            Kind::Semi => 2 * (st.seg_height.saturating_sub(1) as u64) * rs,
            Kind::Quantile => 2 * rw,
        };
        // (reads, canonical tails) of a walk of skeleton dimension j, whose
        // cascade steps each cost at most rcs[j].
        fn walk(l: &Layered, kind: Kind, j: usize, rcs: [u64; 2], tail: u64) -> (u64, u64) {
            let h = l.stats.heights[j] as u64;
            let rc = rcs[j];
            let (sib, sib_c, leaf, leaf_c) = if j + 2 == l.d {
                let s = if kind == Kind::Group { 0 } else { 2 * rc + tail };
                (s, 1, tail, 1)
            } else {
                let (w, c) = walk(l, kind, j + 1, rcs, tail);
                (2 * rc + w, c, w, c)
            };
            (4 * h * rc + 2 * h * sib + 2 * leaf, 2 * h * sib_c + 2 * leaf_c)
        }
        let (skel, canon) = if self.d == 1 { (tail, 1) } else { walk(self, kind, 0, rcs, tail) };
        let descent = if kind == Kind::Quantile {
            (st.wav_height as u64 - 1) * canon * 2 * rw
        } else {
            0
        };
        bst + skel + descent
    }

    // -----------------------------------------------------------------------
    // Queries.
    // -----------------------------------------------------------------------

    pub fn query<S: Sam>(&mut self, ctx: &mut Ctx, sam: &mut S, q: &Query, kind: Kind) -> Res<QOut> {
        let d = self.d;
        let (e1, e2) = self.bst_search(ctx, sam, q.lo[d - 1].saturating_sub(1), q.hi[d - 1])?;
        let mut st = QState {
            q,
            d,
            kind,
            g: Group::ZERO,
            s: Semi::IDENTITY,
            canon: Vec::new(),
            tmax: e1.t.max(e2.t),
        };
        if d == 1 {
            st.tail(ctx, sam, e1, e2)?;
        } else {
            st.walk(ctx, sam, 0, e1, e2)?;
        }
        let median = if kind == Kind::Quantile {
            st.descend(ctx, sam, self.stats.wav_height, &self.alphabet)?
        } else {
            None
        };
        if ctx.open_values() != 0 {
            return Err(SamError::Backend(format!("{} objects left open", ctx.open_values())));
        }
        Ok(QOut {
            g: st.g,
            s: st.s,
            median,
            rounds: st.tmax,
        })
    }

    /// Two lockstep predecessor searches: largest value <= q1 and <= q2.
    fn bst_search<S: Sam>(&mut self, ctx: &mut Ctx, sam: &mut S, q1: i64, q2: i64) -> Res<(H, H)> {
        let a = ctx.open(sam, &mut self.bst_root, 0)?;
        let b = ctx.share(&a)?;
        let (mut a, mut b, mut together) = (a, b, true);
        loop {
            let (ka, la) = {
                let n = ctx.val(&a).bst();
                (n.key, n.l.is_none())
            };
            let lb = ctx.val(&b).bst().l.is_none();
            if la && lb {
                break;
            }
            if together {
                let (ra, rb) = (q1 >= ka, q2 >= ka);
                if ra == rb {
                    let n = ctx.open_field(sam, &a, bst_sel(ra))?;
                    ctx.close(sam, a)?;
                    ctx.close(sam, b)?;
                    b = ctx.share(&n)?;
                    a = n;
                } else {
                    let na = ctx.open_field(sam, &a, bst_sel(ra))?;
                    let nb = ctx.open_field(sam, &b, bst_sel(rb))?;
                    ctx.close(sam, a)?;
                    ctx.close(sam, b)?;
                    a = na;
                    b = nb;
                    together = false;
                }
            } else {
                if !la {
                    let na = ctx.open_field(sam, &a, bst_sel(q1 >= ka))?;
                    ctx.close(sam, a)?;
                    a = na;
                }
                if !lb {
                    let kb = ctx.val(&b).bst().key;
                    let nb = ctx.open_field(sam, &b, bst_sel(q2 >= kb))?;
                    ctx.close(sam, b)?;
                    b = nb;
                }
            }
        }
        let e1 = ctx.open_field(sam, &a, |o| o.bst_mut().entry.as_mut())?;
        let e2 = if together {
            ctx.share(&e1)?
        } else {
            ctx.open_field(sam, &b, |o| o.bst_mut().entry.as_mut())?
        };
        ctx.close(sam, a)?;
        ctx.close(sam, b)?;
        Ok((e1, e2))
    }
}

fn bst_sel(right: bool) -> impl for<'a> FnOnce(&'a mut Obj) -> Option<&'a mut Ptr> {
    move |o: &mut Obj| {
        let n = o.bst_mut();
        if right {
            n.r.as_mut()
        } else {
            n.l.as_mut()
        }
    }
}

fn casc_sel(j: usize, right: bool) -> impl for<'a> FnOnce(&'a mut Obj) -> Option<&'a mut Ptr> {
    move |o: &mut Obj| {
        o.entry_mut().casc[j]
            .as_mut()
            .and_then(|c| c.kids.as_mut())
            .map(|k| if right { &mut k.r } else { &mut k.l })
    }
}

fn seg_sel() -> impl for<'a> FnOnce(&'a mut Obj) -> Option<&'a mut Ptr> {
    |o: &mut Obj| match o {
        Obj::Entry(e) => e.seg.as_mut(),
        Obj::Seg(s) => s.parent.as_mut(),
        _ => None,
    }
}

fn wav_sel(right: bool) -> impl for<'a> FnOnce(&'a mut Obj) -> Option<&'a mut Ptr> {
    move |o: &mut Obj| match o {
        Obj::Entry(e) => e.wav.as_mut(),
        Obj::Wav(w) => {
            if right {
                w.r.as_mut()
            } else {
                w.l.as_mut()
            }
        }
        _ => None,
    }
}

#[derive(Clone, Debug)]
pub struct QOut {
    pub g: Group,
    pub s: Semi,
    pub median: Option<i64>,
    /// Longest chain of dependent reads (round trips with same-level
    /// batching), excluding the final padding batch.
    pub rounds: u64,
}

struct QState<'a> {
    q: &'a Query,
    d: usize,
    kind: Kind,
    g: Group,
    s: Semi,
    canon: Vec<(H, H)>,
    tmax: u64,
}

impl QState<'_> {
    fn close<S: Sam>(&mut self, ctx: &mut Ctx, sam: &mut S, h: H) -> Res<()> {
        self.tmax = self.tmax.max(h.t);
        ctx.close(sam, h)
    }

    fn close2<S: Sam>(&mut self, ctx: &mut Ctx, sam: &mut S, p: (H, H)) -> Res<()> {
        self.close(ctx, sam, p.0)?;
        self.close(ctx, sam, p.1)
    }

    fn entry<'c>(ctx: &'c Ctx, h: &H) -> &'c Entry {
        ctx.val(h).entry()
    }

    /// Cascades both entries of a pair into child `right`; None if the two
    /// positions coincide there (the range is empty below).
    fn step<S: Sam>(&mut self, ctx: &mut Ctx, sam: &mut S, j: usize, p: &(H, H), right: bool) -> Res<Option<(H, H)>> {
        let rank = |h: &H| {
            let k = Self::entry(ctx, h).casc[j].as_ref().unwrap().kids.as_ref().unwrap();
            if right {
                k.rank_r
            } else {
                k.rank_l
            }
        };
        if rank(&p.0) == rank(&p.1) {
            return Ok(None);
        }
        let a = ctx.open_field(sam, &p.0, casc_sel(j, right))?;
        let b = ctx.open_field(sam, &p.1, casc_sel(j, right))?;
        Ok(Some((a, b)))
    }

    fn casc_key(ctx: &Ctx, h: &H, j: usize) -> (i64, bool) {
        let c = Self::entry(ctx, h).casc[j].as_ref().expect("skeleton entry has cascade data");
        (c.key, c.kids.is_none())
    }

    /// Aggregates the points of a dimension-j structure inside the query;
    /// (e1, e2) are entries of its root list. Consumes the handles.
    fn walk<S: Sam>(&mut self, ctx: &mut Ctx, sam: &mut S, j: usize, e1: H, e2: H) -> Res<()> {
        if Self::entry(ctx, &e1).rank == Self::entry(ctx, &e2).rank {
            return self.close2(ctx, sam, (e1, e2));
        }
        let (x1, x2) = (self.q.lo[j], self.q.hi[j]);
        let mut cur = (e1, e2);
        loop {
            let (key, leaf) = Self::casc_key(ctx, &cur.0, j);
            if leaf {
                if x1 <= key && key <= x2 {
                    return self.include(ctx, sam, j, cur);
                }
                return self.close2(ctx, sam, cur);
            }
            let (ra, rb) = (x1 >= key, x2 >= key);
            if ra == rb {
                let next = self.step(ctx, sam, j, &cur, ra)?;
                self.close2(ctx, sam, cur)?;
                match next {
                    Some(n) => cur = n,
                    None => return Ok(()),
                }
                continue;
            }
            if ra && !rb {
                return self.close2(ctx, sam, cur); // empty range in dimension j
            }
            let a = self.step(ctx, sam, j, &cur, false)?;
            let b = self.step(ctx, sam, j, &cur, true)?;
            self.close2(ctx, sam, cur)?;
            self.side(ctx, sam, j, a, true)?;
            return self.side(ctx, sam, j, b, false);
        }
    }

    /// One boundary path below the split: `a_side` follows x1 (siblings on
    /// its right are canonical), otherwise x2 (siblings on its left).
    fn side<S: Sam>(&mut self, ctx: &mut Ctx, sam: &mut S, j: usize, mut cur: Option<(H, H)>, a_side: bool) -> Res<()> {
        let (x1, x2) = (self.q.lo[j], self.q.hi[j]);
        while let Some(p) = cur {
            let (key, leaf) = Self::casc_key(ctx, &p.0, j);
            if leaf {
                let inc = if a_side { key >= x1 } else { key <= x2 };
                if inc {
                    return self.include(ctx, sam, j, p);
                }
                return self.close2(ctx, sam, p);
            }
            let right = if a_side { x1 >= key } else { x2 >= key };
            if a_side && !right {
                self.sibling(ctx, sam, j, &p, true)?;
            }
            if !a_side && right {
                self.sibling(ctx, sam, j, &p, false)?;
            }
            let next = self.step(ctx, sam, j, &p, right)?;
            self.close2(ctx, sam, p)?;
            cur = next;
        }
        Ok(())
    }

    /// The child `right` of the pair's node is fully inside dimension j.
    fn sibling<S: Sam>(&mut self, ctx: &mut Ctx, sam: &mut S, j: usize, p: &(H, H), right: bool) -> Res<()> {
        if j + 2 == self.d && self.kind == Kind::Group {
            let pre = |h: &H| {
                let (l, r) = Self::entry(ctx, h).casc[j].as_ref().unwrap().kids.as_ref().unwrap().pre.unwrap();
                if right {
                    r
                } else {
                    l
                }
            };
            self.g = self.g.add(&pre(&p.1).sub(&pre(&p.0)));
            return Ok(());
        }
        if let Some(s) = self.step(ctx, sam, j, p, right)? {
            if j + 2 == self.d {
                self.tail(ctx, sam, s.0, s.1)?;
            } else {
                self.walk(ctx, sam, j + 1, s.0, s.1)?;
            }
        }
        Ok(())
    }

    /// A leaf of the dimension-j tree whose value is inside the query.
    fn include<S: Sam>(&mut self, ctx: &mut Ctx, sam: &mut S, j: usize, p: (H, H)) -> Res<()> {
        if j + 2 == self.d {
            self.tail(ctx, sam, p.0, p.1)
        } else {
            self.walk(ctx, sam, j + 1, p.0, p.1)
        }
    }

    /// Last-dimension list entries e1 (excluded) and e2 (included).
    fn tail<S: Sam>(&mut self, ctx: &mut Ctx, sam: &mut S, e1: H, e2: H) -> Res<()> {
        let (r1, r2) = (Self::entry(ctx, &e1).rank, Self::entry(ctx, &e2).rank);
        if r1 >= r2 {
            return self.close2(ctx, sam, (e1, e2));
        }
        match self.kind {
            Kind::Group => {
                let g = Self::entry(ctx, &e2).own.sub(&Self::entry(ctx, &e1).own);
                self.g = self.g.add(&g);
                self.close2(ctx, sam, (e1, e2))
            }
            Kind::Semi => {
                self.s = self.s.merge(&Self::entry(ctx, &e2).semi);
                let (lr, rr) = (r1 as u64, r2 as u64);
                let (mut l, mut r) = (e1, e2);
                let mut k = 0;
                while (lr >> (k + 1)) != (rr >> (k + 1)) {
                    let pl = ctx.open_field(sam, &l, seg_sel())?;
                    let pr = ctx.open_field(sam, &r, seg_sel())?;
                    if (lr >> k) & 1 == 0 {
                        self.s = self.s.merge(&ctx.val(&pl).seg().right);
                    }
                    if (rr >> k) & 1 == 1 {
                        self.s = self.s.merge(&ctx.val(&pr).seg().left);
                    }
                    self.close2(ctx, sam, (l, r))?;
                    l = pl;
                    r = pr;
                    k += 1;
                }
                self.close2(ctx, sam, (l, r))
            }
            Kind::Quantile => {
                if Self::entry(ctx, &e1).own.cnt == Self::entry(ctx, &e2).own.cnt {
                    return self.close2(ctx, sam, (e1, e2)); // sampled copies only
                }
                let w1 = ctx.open_field(sam, &e1, wav_sel(false))?;
                let w2 = ctx.open_field(sam, &e2, wav_sel(false))?;
                self.close2(ctx, sam, (e1, e2))?;
                self.canon.push((w1, w2));
                Ok(())
            }
        }
    }

    /// Lockstep wavelet descent for the lower median of all canonical ranges.
    fn descend<S: Sam>(&mut self, ctx: &mut Ctx, sam: &mut S, hs: usize, alphabet: &[i64]) -> Res<Option<i64>> {
        let pos = |ctx: &Ctx, h: &H| {
            let w = ctx.val(h).wav();
            (w.zeros, w.ones)
        };
        let total: u64 = self
            .canon
            .iter()
            .map(|(a, b)| {
                let (z1, o1) = pos(ctx, a);
                let (z2, o2) = pos(ctx, b);
                ((z2 + o2) - (z1 + o1)) as u64
            })
            .sum();
        if total == 0 {
            for p in std::mem::take(&mut self.canon) {
                self.close2(ctx, sam, p)?;
            }
            return Ok(None);
        }
        let mut k = (total - 1) / 2;
        let mut sym = 0usize;
        for level in 0..hs {
            let z: u64 = self
                .canon
                .iter()
                .map(|(a, b)| (pos(ctx, b).0 - pos(ctx, a).0) as u64)
                .sum();
            let bit = if k < z {
                0
            } else {
                k -= z;
                1
            };
            sym = (sym << 1) | bit;
            if level + 1 == hs {
                break;
            }
            let t = self.canon.iter().map(|(a, b)| a.t.max(b.t)).max().unwrap_or(0);
            let mut next = Vec::new();
            for (a, b) in std::mem::take(&mut self.canon) {
                let (pa, pb) = (pos(ctx, &a), pos(ctx, &b));
                let (ia, ib) = if bit == 0 { (pa.0, pb.0) } else { (pa.1, pb.1) };
                if ia == ib {
                    self.close2(ctx, sam, (a, b))?;
                    continue;
                }
                let na = ctx.open_field_at(sam, &a, t, wav_sel(bit == 1))?;
                let nb = ctx.open_field_at(sam, &b, t, wav_sel(bit == 1))?;
                self.close2(ctx, sam, (a, b))?;
                next.push((na, nb));
            }
            self.canon = next;
        }
        for p in std::mem::take(&mut self.canon) {
            self.close2(ctx, sam, p)?;
        }
        Ok(alphabet.get(sym).copied())
    }
}
