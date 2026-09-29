//! Query-time access layer over sam-model's client-side cache.
//!
//! [`CachedPointers`] gives explicit handles: dereferencing a pointer field of
//! an open object keeps the owner open (its pointer field is refreshed in
//! place), and releasing a handle writes the object back. This lets a query
//! keep many objects open across lockstep levels while obeying SAM+ rules.
//!
//! Every handle carries a logical time: the number of sequential reads on
//! the longest dependency chain that produced it. With all independent reads
//! of a level batched into one request, a query's round trips equal the
//! maximum handle time (plus one batch for the final padding reads).

use crate::obj::{Obj, Ptr};
use sam_model::pointer::{BalancedCell, BalancedPointers};
use sam_model::pointer::{CacheObject, CachedPointers};
use sam_model::{MemoryClass, SamError, SingleAccessMachine};

pub type Cell = BalancedCell<Obj>;
pub type Res<T> = Result<T, SamError>;

pub trait Sam: SingleAccessMachine<Cell> {}
impl<T: SingleAccessMachine<Cell>> Sam for T {}

#[derive(Clone, Copy, Debug)]
pub struct H {
    pub o: CacheObject,
    pub t: u64,
}

pub struct Ctx {
    pub cache: CachedPointers<Obj, BalancedPointers>,
}

pub fn reads<S: Sam>(sam: &S) -> u64 {
    sam.stats().operations.reads
}

impl Ctx {
    pub fn new(branching_factor: usize) -> Res<Self> {
        Ok(Self {
            cache: CachedPointers::new(BalancedPointers::new(branching_factor)?),
        })
    }

    /// Dereferences a client-held pointer.
    pub fn open<S: Sam>(&mut self, sam: &mut S, p: &mut Ptr, t: u64) -> Res<H> {
        let r0 = reads(sam);
        let o = self
            .cache
            .get(sam, p)?
            .ok_or(SamError::InvalidPointerCell("dereferenced an empty pointer"))?;
        Ok(H {
            o,
            t: t + (reads(sam) - r0),
        })
    }

    /// Dereferences a pointer field of an open object, starting at time `t`.
    pub fn open_field_at<S, F>(&mut self, sam: &mut S, owner: &H, t: u64, select: F) -> Res<H>
    where
        S: Sam,
        F: for<'a> FnOnce(&'a mut Obj) -> Option<&'a mut Ptr>,
    {
        let r0 = reads(sam);
        let o = self
            .cache
            .get_pointer_field(sam, &owner.o, select)?
            .ok_or(SamError::InvalidPointerCell("dereferenced an absent field"))?;
        Ok(H {
            o,
            t: t.max(owner.t) + (reads(sam) - r0),
        })
    }

    pub fn open_field<S, F>(&mut self, sam: &mut S, owner: &H, select: F) -> Res<H>
    where
        S: Sam,
        F: for<'a> FnOnce(&'a mut Obj) -> Option<&'a mut Ptr>,
    {
        self.open_field_at(sam, owner, owner.t, select)
    }

    pub fn val(&self, h: &H) -> &Obj {
        self.cache.value(&h.o).expect("handle is live")
    }

    /// A second reference to an open object (no reads).
    pub fn share(&mut self, h: &H) -> Res<H> {
        Ok(H {
            o: self.cache.clone_object(h.o)?,
            t: h.t,
        })
    }

    pub fn close<S: Sam>(&mut self, sam: &mut S, h: H) -> Res<()> {
        let mut o = h.o;
        self.cache.release(sam, &mut o)
    }

    pub fn open_values(&self) -> usize {
        self.cache.cached_values()
    }
}

/// Pads the trace with reads of fresh, never-written addresses (valid in
/// SAM+; indistinguishable from real reads in Path OSAM+).
pub fn dummy_reads<S: Sam>(sam: &mut S, n: u64) -> Res<()> {
    for _ in 0..n {
        let a = sam.alloc(MemoryClass::Oblivious, "Padding");
        let _ = sam.read(a, "Padding")?;
    }
    Ok(())
}
