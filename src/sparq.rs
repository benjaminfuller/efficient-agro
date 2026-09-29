//! SPARQ baseline: nested multidimensional segment tree over the distinct
//! coordinate values, traversed breadth-first on OSAM+.
//!
//! Every node has exactly one parent (fan-in 1). Following the paper's
//! schedule lemma, a query processes the dimensions in sequence: fully covered
//! nodes of dimension i are held open until the dimension-i trees are done,
//! then their next-dimension roots start together. The public budget is the
//! lemma's schedule, sum_i P_{i-1} * sum_l min(2^l, 4) with P_i the product of
//! the per-dimension bounds on returned nodes, 2 H_j - 1.

use crate::agg::{Group, Semi};
use crate::ctx::{Ctx, Res, Sam, H};
use crate::data::{Dataset, Query};
use crate::obj::{Obj, Ptr, StNode};
use sam_model::pointer::{BalancedPointer, BalancedPointers};
use sam_model::pointer::CachedPointer;
use sam_model::SamError;

#[derive(Clone, Debug, Default)]
pub struct SparqStats {
    pub nodes: usize,
    pub heights: [usize; 3],
}

pub struct Sparq {
    pub d: usize,
    pub root: Ptr,
    pub stats: SparqStats,
    /// Largest (real accesses / schedule cap) seen at any level.
    pub max_level_ratio: f64,
    pub level_overflows: u64,
}

pub struct SparqOut {
    pub g: Group,
    pub s: Semi,
    pub rounds: u64,
}

struct Builder<'a, S> {
    ds: &'a Dataset,
    sam: &'a mut S,
    b: usize,
    stats: SparqStats,
}

impl<S: Sam> Builder<'_, S> {
    fn tree(&mut self, t: usize, mut idx: Vec<u32>) -> Res<BalancedPointer> {
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
        self.node(t, &idx, &values, &starts, 0, values.len(), 0)
    }

    #[allow(clippy::too_many_arguments)]
    fn node(&mut self, t: usize, idx: &[u32], values: &[i64], starts: &[usize], lo: usize, hi: usize, depth: usize) -> Res<BalancedPointer> {
        self.stats.heights[t] = self.stats.heights[t].max(depth);
        let slice = &idx[starts[lo]..starts[hi]];
        let (mut g, mut s) = (Group::ZERO, Semi::IDENTITY);
        for &i in slice {
            let v = self.ds.pts[i as usize].v;
            g = g.add(&Group::of(v));
            s = s.merge(&Semi::of(v));
        }
        let next = if t + 1 < self.ds.d {
            Some(CachedPointer::from_raw(self.tree(t + 1, slice.to_vec())?))
        } else {
            None
        };
        let (l, r) = if hi - lo > 1 {
            let mid = lo + (hi - lo) / 2;
            let l = self.node(t, idx, values, starts, lo, mid, depth + 1)?;
            let r = self.node(t, idx, values, starts, mid, hi, depth + 1)?;
            (Some(CachedPointer::from_raw(l)), Some(CachedPointer::from_raw(r)))
        } else {
            (None, None)
        };
        self.stats.nodes += 1;
        let value = Obj::St(Box::new(StNode {
            lo: values[lo],
            hi: values[hi - 1],
            g,
            s,
            l,
            r,
            next,
        }));
        Ok(BalancedPointers::new(self.b)?.install(self.sam, value, 1)?.remove(0))
    }
}

fn returned_bound(h: usize) -> u64 {
    if h == 0 {
        1
    } else {
        2 * h as u64 - 1
    }
}

impl Sparq {
    pub fn build<S: Sam>(ds: &Dataset, branching: usize, sam: &mut S) -> Res<Self> {
        let mut b = Builder {
            ds,
            sam,
            b: branching,
            stats: SparqStats::default(),
        };
        let root = b.tree(0, (0..ds.pts.len() as u32).collect())?;
        Ok(Self {
            d: ds.d,
            root: CachedPointer::from_raw(root),
            stats: b.stats,
            max_level_ratio: 0.0,
            level_overflows: 0,
        })
    }

    /// Accesses allowed at level `l` of dimension `t` by the schedule.
    pub fn cap(&self, t: usize, l: usize) -> u64 {
        let p: u64 = (0..t).map(|j| returned_bound(self.stats.heights[j])).product();
        p * (1u64 << l.min(2))
    }

    pub fn analytic_budget(&self) -> u64 {
        (0..self.d)
            .map(|t| (0..=self.stats.heights[t]).map(|l| self.cap(t, l)).sum::<u64>())
            .sum()
    }

    pub fn query<S: Sam>(&mut self, ctx: &mut Ctx, sam: &mut S, q: &Query) -> Res<SparqOut> {
        let (mut g, mut s) = (Group::ZERO, Semi::IDENTITY);
        let mut tmax = 0;
        let mut frontier: Vec<H> = vec![ctx.open(sam, &mut self.root, 0)?];
        for t in 0..self.d {
            let mut deferred: Vec<H> = Vec::new();
            let mut level = 0;
            while !frontier.is_empty() {
                let cap = self.cap(t, level);
                let used = frontier.len() as u64;
                self.max_level_ratio = self.max_level_ratio.max(used as f64 / cap as f64);
                if used > cap {
                    self.level_overflows += 1;
                }
                let mut next = Vec::new();
                for h in frontier.drain(..) {
                    let (lo, hi, has_kids, ng, ns) = {
                        let n = ctx.val(&h).st();
                        (n.lo, n.hi, n.l.is_some(), n.g, n.s)
                    };
                    tmax = tmax.max(h.t);
                    if hi < q.lo[t] || lo > q.hi[t] {
                        ctx.close(sam, h)?;
                    } else if q.lo[t] <= lo && hi <= q.hi[t] {
                        if t + 1 == self.d {
                            g = g.add(&ng);
                            s = s.merge(&ns);
                            ctx.close(sam, h)?;
                        } else {
                            deferred.push(h);
                        }
                    } else {
                        if !has_kids {
                            return Err(SamError::Backend("partial leaf".into()));
                        }
                        next.push(ctx.open_field(sam, &h, |o| o.st_mut().l.as_mut())?);
                        next.push(ctx.open_field(sam, &h, |o| o.st_mut().r.as_mut())?);
                        ctx.close(sam, h)?;
                    }
                }
                frontier = next;
                level += 1;
            }
            for h in deferred {
                let n = ctx.open_field(sam, &h, |o| o.st_mut().next.as_mut())?;
                ctx.close(sam, h)?;
                frontier.push(n);
            }
        }
        if ctx.open_values() != 0 {
            return Err(SamError::Backend("SPARQ left objects open".into()));
        }
        Ok(SparqOut { g, s, rounds: tmax })
    }
}
