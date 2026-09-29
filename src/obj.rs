//! Objects stored in the SAM and their fixed-size block codec.
//!
//! Every object lives behind a balanced r-ary multi-write smart pointer
//! (`BalancedCell<Obj>` cells, see `sam_model::pointer::BalancedPointers`). Pointer fields persist as stable logical
//! identifiers (0 = absent), exactly as the sam-model graph codecs do.

use crate::agg::{Group, Semi};
use sam_model::pointer::BalancedPointer;
use sam_model::pointer::{CachedPointer, ValueCodec};
use sam_model::{Address, SamError};

pub type Ptr = CachedPointer<BalancedPointer>;

/// Cascade data for one skeleton dimension.
#[derive(Clone, Debug)]
pub struct Casc {
    /// Internal tree node: split key (go right iff q >= key).
    /// Leaf: the leaf's coordinate value.
    pub key: i64,
    pub kids: Option<Kids>,
}

#[derive(Clone, Debug)]
pub struct Kids {
    /// Aliases of the entries this entry cascades to in each child list.
    pub l: Ptr,
    pub r: Ptr,
    /// Ranks of those entries in their (augmented) lists.
    pub rank_l: u32,
    pub rank_r: u32,
    /// Children's prefix aggregates at those ranks (group tail, last skeleton
    /// dimension only), so canonical siblings cost no reads.
    pub pre: Option<(Group, Group)>,
}

/// One entry of a last-dimension list.
#[derive(Clone, Debug)]
pub struct Entry {
    /// Index in its (augmented) list; 0 is the -inf sentinel.
    pub rank: u32,
    /// Prefix aggregate of the list's points with last coordinate <= this
    /// entry's value. `own.cnt` is also this entry's wavelet-root position.
    pub own: Group,
    /// MIN/MAX of the points whose last coordinate equals this entry's value
    /// (identity for the sentinel and for sampled copies).
    pub semi: Semi,
    /// Cascade data per skeleton dimension (index = dimension).
    pub casc: [Option<Casc>; 2],
    /// Semigroup tail: alias of this leaf's parent in the position segment tree.
    pub seg: Option<Ptr>,
    /// Quantile tail: alias of the wavelet-root position entry at `own.cnt`.
    pub wav: Option<Ptr>,
}

/// Node of the predecessor-search tree over the global root list.
#[derive(Clone, Debug)]
pub struct Bst {
    pub key: i64,
    pub l: Option<Ptr>,
    pub r: Option<Ptr>,
    /// Leaf only: rank of the root-list entry and an alias to it.
    pub leaf_rank: u32,
    pub entry: Option<Ptr>,
}

/// Internal node of a semigroup-tail position segment tree.
#[derive(Clone, Debug)]
pub struct Seg {
    pub left: Semi,
    pub right: Semi,
    pub parent: Option<Ptr>,
    /// True if this node is its parent's left child.
    pub is_left: bool,
}

/// Wavelet-tree position entry: rank counts at one position of a node.
#[derive(Clone, Debug)]
pub struct Wav {
    pub zeros: u32,
    pub ones: u32,
    pub l: Option<Ptr>,
    pub r: Option<Ptr>,
}

/// SPARQ baseline: node of a nested multidimensional segment tree.
#[derive(Clone, Debug)]
pub struct StNode {
    pub lo: i64,
    pub hi: i64,
    pub g: Group,
    pub s: Semi,
    pub l: Option<Ptr>,
    pub r: Option<Ptr>,
    pub next: Option<Ptr>,
}

#[derive(Clone, Debug)]
pub enum Obj {
    Entry(Box<Entry>),
    Bst(Box<Bst>),
    Seg(Box<Seg>),
    Wav(Box<Wav>),
    St(Box<StNode>),
}

impl Obj {
    pub fn entry(&self) -> &Entry {
        match self {
            Obj::Entry(e) => e,
            _ => panic!("expected Entry, found {:?}", self.tag()),
        }
    }
    pub fn entry_mut(&mut self) -> &mut Entry {
        match self {
            Obj::Entry(e) => e,
            _ => panic!("expected Entry"),
        }
    }
    pub fn bst(&self) -> &Bst {
        match self {
            Obj::Bst(e) => e,
            _ => panic!("expected Bst"),
        }
    }
    pub fn bst_mut(&mut self) -> &mut Bst {
        match self {
            Obj::Bst(e) => e,
            _ => panic!("expected Bst"),
        }
    }
    pub fn seg(&self) -> &Seg {
        match self {
            Obj::Seg(e) => e,
            _ => panic!("expected Seg"),
        }
    }
    pub fn seg_mut(&mut self) -> &mut Seg {
        match self {
            Obj::Seg(e) => e,
            _ => panic!("expected Seg"),
        }
    }
    pub fn wav(&self) -> &Wav {
        match self {
            Obj::Wav(e) => e,
            _ => panic!("expected Wav"),
        }
    }
    pub fn wav_mut(&mut self) -> &mut Wav {
        match self {
            Obj::Wav(e) => e,
            _ => panic!("expected Wav"),
        }
    }
    pub fn st(&self) -> &StNode {
        match self {
            Obj::St(e) => e,
            _ => panic!("expected St"),
        }
    }
    pub fn st_mut(&mut self) -> &mut StNode {
        match self {
            Obj::St(e) => e,
            _ => panic!("expected St"),
        }
    }
    fn tag(&self) -> &'static str {
        match self {
            Obj::Entry(_) => "Entry",
            Obj::Bst(_) => "Bst",
            Obj::Seg(_) => "Seg",
            Obj::Wav(_) => "Wav",
            Obj::St(_) => "St",
        }
    }
}

/// Byte codec for [`Obj`] (wrapped by `BalancedCellValueCodec` and `FixedSizeCodec`).
#[derive(Clone, Debug)]
pub struct ObjCodec {
    _branching: usize,
}

impl ObjCodec {
    pub fn new(branching_factor: usize) -> Self {
        Self {
            _branching: branching_factor,
        }
    }

    fn put_ptr(&self, p: &Option<Ptr>, out: &mut Vec<u8>) -> Result<(), SamError> {
        let id = match p {
            None => 0,
            Some(p) => {
                let raw = p
                    .raw()
                    .ok_or(SamError::InvalidPointerCell("deleted pointer field"))?;
                match raw.head() {
                    Some(Address::Oblivious(id)) if id != 0 => id,
                    _ => return Err(SamError::InvalidPointerCell("pointer without an oblivious head")),
                }
            }
        };
        put_var(out, id);
        Ok(())
    }

    fn get_ptr(&self, input: &mut &[u8]) -> Result<Option<Ptr>, SamError> {
        let id = get_var(input)?;
        if id == 0 {
            return Ok(None);
        }
        Ok(Some(CachedPointer::from_raw(BalancedPointer::from_head(Address::Oblivious(id)))))
    }

    fn put_req(&self, p: &Ptr, out: &mut Vec<u8>) -> Result<(), SamError> {
        self.put_ptr(&Some(p.clone()), out)
    }

    fn get_req(&self, input: &mut &[u8]) -> Result<Ptr, SamError> {
        self.get_ptr(input)?
            .ok_or(SamError::InvalidPointerCell("missing required pointer"))
    }
}

fn get_bytes<'a>(input: &mut &'a [u8], n: usize) -> Result<&'a [u8], SamError> {
    if input.len() < n {
        return Err(SamError::InvalidPointerCell("truncated object"));
    }
    let (a, b) = input.split_at(n);
    *input = b;
    Ok(a)
}
/// LEB128 varint of an unsigned value.
fn put_var(out: &mut Vec<u8>, mut x: u64) {
    loop {
        let byte = (x & 0x7f) as u8;
        x >>= 7;
        if x == 0 {
            out.push(byte);
            return;
        }
        out.push(byte | 0x80);
    }
}
fn get_var(input: &mut &[u8]) -> Result<u64, SamError> {
    let mut x = 0u64;
    for shift in (0..64).step_by(7) {
        let byte = get_u8(input)?;
        x |= ((byte & 0x7f) as u64) << shift;
        if byte & 0x80 == 0 {
            return Ok(x);
        }
    }
    Err(SamError::InvalidPointerCell("varint too long"))
}
/// Zigzag varint of a signed value (small magnitudes take few bytes).
fn put_int(out: &mut Vec<u8>, x: i64) {
    put_var(out, ((x << 1) ^ (x >> 63)) as u64);
}
fn get_i64(input: &mut &[u8]) -> Result<i64, SamError> {
    let z = get_var(input)?;
    Ok(((z >> 1) as i64) ^ -((z & 1) as i64))
}
fn get_u32(input: &mut &[u8]) -> Result<u32, SamError> {
    u32::try_from(get_i64(input)?).map_err(|_| SamError::InvalidPointerCell("u32 field out of range"))
}
fn get_u8(input: &mut &[u8]) -> Result<u8, SamError> {
    Ok(get_bytes(input, 1)?[0])
}
fn put_group(g: &Group, out: &mut Vec<u8>) {
    put_int(out, g.cnt);
    put_int(out, g.sum);
    put_int(out, g.sq);
}
fn get_group(input: &mut &[u8]) -> Result<Group, SamError> {
    Ok(Group {
        cnt: get_i64(input)?,
        sum: get_i64(input)?,
        sq: get_i64(input)?,
    })
}
/// The MIN/MAX identity (i64::MAX, i64::MIN) is common (sentinels, sampled
/// copies) and would take 20 bytes as varints, so it is a one-byte tag.
fn put_semi(s: &Semi, out: &mut Vec<u8>) {
    if *s == Semi::IDENTITY {
        out.push(0);
    } else {
        out.push(1);
        put_int(out, s.min);
        put_int(out, s.max);
    }
}
fn get_semi(input: &mut &[u8]) -> Result<Semi, SamError> {
    Ok(match get_u8(input)? {
        0 => Semi::IDENTITY,
        _ => Semi {
            min: get_i64(input)?,
            max: get_i64(input)?,
        },
    })
}

impl ValueCodec<Obj> for ObjCodec {
    fn encode_value(&self, value: &Obj, out: &mut Vec<u8>) -> Result<(), SamError> {
        match value {
            Obj::Entry(e) => {
                out.push(1);
                put_int(out, e.rank as i64);
                put_group(&e.own, out);
                put_semi(&e.semi, out);
                let mut flags = 0u8;
                for (j, c) in e.casc.iter().enumerate() {
                    if let Some(c) = c {
                        flags |= 1 << (2 * j);
                        if let Some(k) = &c.kids {
                            flags |= 2 << (2 * j);
                            if k.pre.is_some() {
                                flags |= 16 << j;
                            }
                        }
                    }
                }
                out.push(flags);
                for c in e.casc.iter().flatten() {
                    put_int(out, c.key as i64);
                    if let Some(k) = &c.kids {
                        self.put_req(&k.l, out)?;
                        self.put_req(&k.r, out)?;
                        put_int(out, k.rank_l as i64);
                        put_int(out, k.rank_r as i64);
                        if let Some((a, b)) = &k.pre {
                            put_group(a, out);
                            put_group(b, out);
                        }
                    }
                }
                self.put_ptr(&e.seg, out)?;
                self.put_ptr(&e.wav, out)?;
            }
            Obj::Bst(b) => {
                out.push(2);
                put_int(out, b.key as i64);
                self.put_ptr(&b.l, out)?;
                self.put_ptr(&b.r, out)?;
                put_int(out, b.leaf_rank as i64);
                self.put_ptr(&b.entry, out)?;
            }
            Obj::Seg(s) => {
                out.push(3);
                put_semi(&s.left, out);
                put_semi(&s.right, out);
                self.put_ptr(&s.parent, out)?;
                out.push(u8::from(s.is_left));
            }
            Obj::Wav(w) => {
                out.push(4);
                put_int(out, w.zeros as i64);
                put_int(out, w.ones as i64);
                self.put_ptr(&w.l, out)?;
                self.put_ptr(&w.r, out)?;
            }
            Obj::St(s) => {
                out.push(5);
                put_int(out, s.lo as i64);
                put_int(out, s.hi as i64);
                put_group(&s.g, out);
                put_semi(&s.s, out);
                self.put_ptr(&s.l, out)?;
                self.put_ptr(&s.r, out)?;
                self.put_ptr(&s.next, out)?;
            }
        }
        Ok(())
    }

    fn decode_value(&self, input: &mut &[u8]) -> Result<Obj, SamError> {
        Ok(match get_u8(input)? {
            1 => {
                let rank = get_u32(input)?;
                let own = get_group(input)?;
                let semi = get_semi(input)?;
                let flags = get_u8(input)?;
                let mut casc: [Option<Casc>; 2] = [None, None];
                for (j, slot) in casc.iter_mut().enumerate() {
                    if flags & (1 << (2 * j)) == 0 {
                        continue;
                    }
                    let key = get_i64(input)?;
                    let kids = if flags & (2 << (2 * j)) != 0 {
                        let l = self.get_req(input)?;
                        let r = self.get_req(input)?;
                        let rank_l = get_u32(input)?;
                        let rank_r = get_u32(input)?;
                        let pre = if flags & (16 << j) != 0 {
                            Some((get_group(input)?, get_group(input)?))
                        } else {
                            None
                        };
                        Some(Kids {
                            l,
                            r,
                            rank_l,
                            rank_r,
                            pre,
                        })
                    } else {
                        None
                    };
                    *slot = Some(Casc { key, kids });
                }
                let seg = self.get_ptr(input)?;
                let wav = self.get_ptr(input)?;
                Obj::Entry(Box::new(Entry {
                    rank,
                    own,
                    semi,
                    casc,
                    seg,
                    wav,
                }))
            }
            2 => Obj::Bst(Box::new(Bst {
                key: get_i64(input)?,
                l: self.get_ptr(input)?,
                r: self.get_ptr(input)?,
                leaf_rank: get_u32(input)?,
                entry: self.get_ptr(input)?,
            })),
            3 => Obj::Seg(Box::new(Seg {
                left: get_semi(input)?,
                right: get_semi(input)?,
                parent: self.get_ptr(input)?,
                is_left: get_u8(input)? != 0,
            })),
            4 => Obj::Wav(Box::new(Wav {
                zeros: get_u32(input)?,
                ones: get_u32(input)?,
                l: self.get_ptr(input)?,
                r: self.get_ptr(input)?,
            })),
            5 => Obj::St(Box::new(StNode {
                lo: get_i64(input)?,
                hi: get_i64(input)?,
                g: get_group(input)?,
                s: get_semi(input)?,
                l: self.get_ptr(input)?,
                r: self.get_ptr(input)?,
                next: self.get_ptr(input)?,
            })),
            t => return Err(SamError::InvalidPointerCell(if t == 0 { "zero tag" } else { "unknown object tag" })),
        })
    }
}
