//! Standalone model of the r-ary multi-write alias tree of sam-model
//! (`RaryPointer`), faithful to its rebuild rule, with an optional
//! randomized layout.
//!
//! Tree: leaves are aliases; internal nodes have at most `b` ordered
//! children; the root holds the value. Dereferencing the alias at depth D
//! (edges to the root) reads the D + 1 cells on its path. The path is then
//! dissolved and rebuilt bottom-up: the fresh leaf plus the sibling blocks of
//! every path cell are accumulated level by level; whenever more than `b`
//! items are pending they are split into `ceil(n/b)` balanced chunks, each
//! under a new node; at the root, packing repeats until at most `b` remain.
//! Sibling blocks come from a virtual balanced binary split of the `b` slots
//! (`RaryPointer::sibling_blocks`), smallest block first.
//!
//! `randomized = true` shuffles the pending members before every chunking and
//! before they are placed under the new root, and the initial alias-to-leaf
//! assignment; everything else is unchanged.

use rand::{rngs::StdRng, seq::SliceRandom, SeedableRng};

const NONE: usize = usize::MAX;

#[derive(Clone, Debug)]
struct N {
    parent: usize,
    children: Vec<usize>,
}

#[derive(Clone)]
pub struct AliasTree {
    b: usize,
    nodes: Vec<N>,
    free: Vec<usize>,
    root: usize,
    leaf_of: Vec<usize>,
    randomized: bool,
    rng: StdRng,
    /// Writes of the most recent `access`: one per member of every rebuilt
    /// sibling group (each member's cell is rewritten with its new parent),
    /// plus the write-back of the value to the new root.
    last_writes: u64,
}

/// Mirrors `RaryPointer::sibling_blocks(index, b, b)`: slot blocks from the
/// leaf upward.
fn sibling_blocks(index: usize, b: usize) -> Vec<(usize, usize)> {
    let mut top_down = Vec::new();
    let (mut lo, mut hi) = (0, b);
    while hi - lo > 1 {
        let mid = (lo + hi).div_ceil(2);
        if index < mid {
            top_down.push((mid, hi));
            hi = mid;
        } else {
            top_down.push((lo, mid));
            lo = mid;
        }
    }
    top_down.reverse();
    top_down
}

impl AliasTree {
    /// `RaryPointer::install(value, f, b)`: a balanced tree, aliases numbered
    /// in depth-first order.
    pub fn install(f: usize, b: usize, randomized: bool, seed: u64) -> Self {
        let mut t = Self {
            b,
            nodes: Vec::new(),
            free: Vec::new(),
            root: 0,
            leaf_of: Vec::new(),
            randomized,
            rng: StdRng::seed_from_u64(seed),
            last_writes: 0,
        };
        t.root = t.alloc(NONE);
        let mut leaves = Vec::new();
        if f == 1 {
            // A single pointer is the root itself (0 edges).
            leaves.push(t.root);
        } else {
            t.build(t.root, f, &mut leaves);
        }
        if randomized {
            leaves.shuffle(&mut t.rng);
        }
        t.leaf_of = leaves;
        t
    }

    fn alloc(&mut self, parent: usize) -> usize {
        let n = N { parent, children: Vec::new() };
        if let Some(i) = self.free.pop() {
            self.nodes[i] = n;
            i
        } else {
            self.nodes.push(n);
            self.nodes.len() - 1
        }
    }

    fn release(&mut self, i: usize) {
        self.nodes[i].children.clear();
        self.nodes[i].parent = NONE;
        self.free.push(i);
    }

    fn set_children(&mut self, parent: usize, kids: Vec<usize>) {
        for &k in &kids {
            self.nodes[k].parent = parent;
        }
        self.nodes[parent].children = kids;
    }

    fn build(&mut self, parent: usize, count: usize, leaves: &mut Vec<usize>) {
        let b = self.b;
        if count <= b {
            let kids: Vec<usize> = (0..count).map(|_| self.alloc(parent)).collect();
            leaves.extend(&kids);
            self.set_children(parent, kids);
            return;
        }
        let base = count / b;
        let rem = count % b;
        let kids: Vec<usize> = (0..b).map(|_| self.alloc(parent)).collect();
        self.set_children(parent, kids.clone());
        for (i, k) in kids.into_iter().enumerate() {
            let quota = base + usize::from(i < rem);
            if quota == 1 {
                leaves.push(k);
            } else {
                self.build(k, quota, leaves);
            }
        }
    }

    /// Replaces the random source (used for common-random-number coupling:
    /// one independent stream per step).
    pub fn reseed(&mut self, seed: u64) {
        self.rng = StdRng::seed_from_u64(seed);
    }

    /// Depth of every alias (a coarse fingerprint of the layout).
    pub fn depths(&self) -> Vec<usize> {
        (0..self.aliases()).map(|a| self.depth(a)).collect()
    }

    pub fn aliases(&self) -> usize {
        self.leaf_of.len()
    }

    pub fn depth(&self, alias: usize) -> usize {
        let mut d = 0;
        let mut c = self.leaf_of[alias];
        while c != self.root {
            c = self.nodes[c].parent;
            d += 1;
        }
        d
    }

    /// Height of the tree (maximum alias depth).
    pub fn height(&self) -> usize {
        (0..self.aliases()).map(|a| self.depth(a)).max().unwrap_or(0)
    }

    fn sibling_sets(&self, node: usize) -> Vec<Vec<usize>> {
        let p = self.nodes[node].parent;
        let group = &self.nodes[p].children;
        let index = group.iter().position(|&c| c == node).unwrap();
        sibling_blocks(index, self.b)
            .into_iter()
            .filter_map(|(lo, hi)| {
                let m: Vec<usize> = (lo..hi).filter_map(|s| group.get(s).copied()).collect();
                (!m.is_empty()).then_some(m)
            })
            .collect()
    }

    fn pack(&mut self, mut members: Vec<usize>) -> Vec<Vec<usize>> {
        let b = self.b;
        if members.len() <= b {
            return vec![members];
        }
        if self.randomized {
            members.shuffle(&mut self.rng);
        }
        let chunks = members.len().div_ceil(b);
        let base = members.len() / chunks;
        let rem = members.len() % chunks;
        let mut out = Vec::new();
        let mut it = members.into_iter();
        for i in 0..chunks {
            let size = base + usize::from(i < rem);
            let chunk: Vec<usize> = it.by_ref().take(size).collect();
            if chunk.len() == 1 {
                out.push(chunk);
            } else {
                let p = self.alloc(NONE);
                self.last_writes += chunk.len() as u64;
                self.set_children(p, chunk);
                out.push(vec![p]);
            }
        }
        out
    }

    /// Writes performed by the most recent [`Self::access`].
    pub fn last_writes(&self) -> u64 {
        self.last_writes
    }

    /// Dereferences `alias`; returns the number of reads.
    pub fn access(&mut self, alias: usize) -> u64 {
        let leaf = self.leaf_of[alias];
        // The value's write-back to the (new) root.
        self.last_writes = 1;
        if leaf == self.root {
            // Sole pointer: one read of the root, rewritten in place.
            return 1;
        }
        let d = self.depth(alias);
        let fresh = self.alloc(NONE);
        self.leaf_of[alias] = fresh;
        let mut pending: Vec<Vec<usize>> = vec![vec![fresh]];
        pending.extend(self.sibling_sets(leaf));
        let mut current = self.nodes[leaf].parent;
        self.release(leaf);
        loop {
            if current == self.root {
                let mut members: Vec<usize> = pending.concat();
                while members.len() > self.b {
                    members = self.pack(members).concat();
                }
                if self.randomized {
                    members.shuffle(&mut self.rng);
                }
                let root = self.alloc(NONE);
                self.last_writes += members.len() as u64;
                self.set_children(root, members);
                self.release(current);
                self.root = root;
                break;
            }
            pending.extend(self.sibling_sets(current));
            let members: Vec<usize> = pending.concat();
            pending = self.pack(members);
            let up = self.nodes[current].parent;
            self.release(current);
            current = up;
        }
        d as u64 + 1
    }

    // ---------------------------------------------------------------
    // r-ary path halving with collapse (candidate rule with a provable
    // amortized bound). The carried set K is the new child set of the
    // current path node p. Two path nodes are read per step.
    // ---------------------------------------------------------------

    /// Number of aliases below `v`.
    pub fn size_of(&self, v: usize) -> usize {
        let ch = &self.nodes[v].children;
        if ch.is_empty() { 1 } else { ch.iter().map(|&c| self.size_of(c)).sum() }
    }

    /// Potential: sum over internal nodes of log2(size).
    pub fn potential(&self) -> f64 {
        fn rec(t: &AliasTree, v: usize) -> (usize, f64) {
            let ch = &t.nodes[v].children;
            if ch.is_empty() { return (1, 0.0); }
            let (mut s, mut phi) = (0, 0.0);
            for &c in ch { let (a, b) = rec(t, c); s += a; phi += b; }
            (s, phi + (s as f64).log2())
        }
        if self.nodes[self.root].children.is_empty() { 0.0 } else { rec(self, self.root).1 }
    }

    fn new_node(&mut self, kids: Vec<usize>) -> usize {
        let v = self.alloc(NONE);
        self.last_writes += kids.len() as u64;
        self.set_children(v, kids);
        v
    }

    fn siblings(&self, v: usize) -> Vec<usize> {
        let p = self.nodes[v].parent;
        self.nodes[p].children.iter().copied().filter(|&c| c != v).collect()
    }

    /// Dereferences `alias` with r-ary path halving; returns the reads.
    pub fn access_halving(&mut self, alias: usize) -> u64 {
        let b = self.b;
        self.last_writes = 1;
        let leaf = self.leaf_of[alias];
        if leaf == self.root { return 1; }
        let fresh = self.alloc(NONE);
        self.leaf_of[alias] = fresh;
        let mut p = self.nodes[leaf].parent;
        let mut k: Vec<usize> = self.nodes[p].children.iter().map(|&c| if c == leaf { fresh } else { c }).collect();
        self.release(leaf);
        let mut reads = 1u64;
        loop {
            reads += 1; // read p
            if p == self.root {
                let r = self.new_node(k);
                self.release(p);
                self.root = r;
                return reads;
            }
            let pp = self.nodes[p].parent;
            let s1 = self.siblings(p);
            reads += 1; // read pp
            if pp == self.root {
                let kids = if k.len() + s1.len() <= b { let mut v = k; v.extend(s1); v } else {
                    let l = self.new_node(k); let mut v = vec![l]; v.extend(s1); v };
                let r = self.new_node(kids);
                self.release(p); self.release(pp);
                self.root = r;
                return reads;
            }
            let ppp = self.nodes[pp].parent;
            let mut s2 = self.siblings(pp);
            let next: Vec<usize> = if k.len() + s1.len() + s2.len() <= b {
                let mut v = k; v.extend(s1); v.extend(s2); v            // collapse both levels
            } else if 1 + s1.len() + s2.len() <= b {
                let l = self.new_node(k); let mut v = vec![l]; v.extend(s1); v.extend(s2); v
            } else {
                let l = self.new_node(k);
                let nd = (s2.len() + 2).saturating_sub(b);         // 0 or 1
                s2.sort_by_key(|&c| self.size_of(c));               // push down the smallest
                let d: Vec<usize> = s2.drain(..nd).collect();
                let mut rset = s1; rset.extend(d);
                let mut v = vec![l];
                if rset.len() == 1 { v.push(rset[0]); } else { let r = self.new_node(rset); v.push(r); }
                v.extend(s2); v
            };
            assert!(next.len() <= b, "carried set exceeds b");
            self.release(p); self.release(pp);
            // ppp keeps its address until it is read in the next step; its
            // child list is replaced by `next`.
            self.set_children(ppp, next);
            k = self.nodes[ppp].children.clone();
            p = ppp;
        }
    }
}
