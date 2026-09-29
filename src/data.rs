//! Datasets, queries, and brute-force oracles.

use crate::agg::{Group, Semi};
use rand::{rngs::StdRng, Rng, SeedableRng};
use std::io::{BufRead, BufReader};

/// One record: up to three coordinates and an integer value.
#[derive(Clone, Copy, Debug)]
pub struct Point {
    pub c: [i64; 3],
    pub v: i64,
}

#[derive(Clone, Debug)]
pub struct Dataset {
    pub name: String,
    pub d: usize,
    pub pts: Vec<Point>,
}

/// How record values are chosen.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ValueMode {
    /// Values stored in the dataset (1 for every point except cali).
    Dataset,
    /// Deterministic pseudo-random 16-bit values derived from the coordinates,
    /// so SUM/STD/MIN/MAX/MEDIAN are non-trivial on every dataset.
    Hash16,
}

impl ValueMode {
    pub fn parse(s: &str) -> Result<Self, String> {
        match s {
            "dataset" => Ok(Self::Dataset),
            "hash16" => Ok(Self::Hash16),
            _ => Err(format!("unknown value mode {s}")),
        }
    }
}

fn splitmix(mut x: u64) -> u64 {
    x = x.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let mut z = x;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

impl Dataset {
    /// Loads the `.pts` format written by `tools/convert_datasets.py`.
    pub fn load(path: &str, values: ValueMode) -> Result<Self, String> {
        let file = std::fs::File::open(path).map_err(|e| format!("{path}: {e}"))?;
        let mut lines = BufReader::new(file).lines();
        let header = lines.next().ok_or("empty file")?.map_err(|e| e.to_string())?;
        let mut it = header.split_whitespace();
        let d: usize = it.next().ok_or("bad header")?.parse().map_err(|_| "bad d")?;
        let n: usize = it.next().ok_or("bad header")?.parse().map_err(|_| "bad N")?;
        if !(1..=3).contains(&d) {
            return Err(format!("unsupported dimension {d}"));
        }
        let mut pts = Vec::with_capacity(n);
        for line in lines {
            let line = line.map_err(|e| e.to_string())?;
            if line.trim().is_empty() {
                continue;
            }
            let nums: Vec<i64> = line
                .split_whitespace()
                .map(|x| x.parse::<i64>().map_err(|_| format!("bad number {x}")))
                .collect::<Result<_, _>>()?;
            if nums.len() != d + 1 {
                return Err(format!("expected {} columns, got {}", d + 1, nums.len()));
            }
            let mut c = [0i64; 3];
            c[..d].copy_from_slice(&nums[..d]);
            let v = match values {
                ValueMode::Dataset => nums[d],
                ValueMode::Hash16 => {
                    let mut h = 0u64;
                    for x in &c[..d] {
                        h = splitmix(h ^ (*x as u64));
                    }
                    (h & 0xFFFF) as i64
                }
            };
            pts.push(Point { c, v });
        }
        if pts.len() != n {
            return Err(format!("header says {n} points, file has {}", pts.len()));
        }
        let name = std::path::Path::new(path)
            .file_stem()
            .map(|s| s.to_string_lossy().into_owned())
            .unwrap_or_default();
        Ok(Self { name, d, pts })
    }

    pub fn bounds(&self) -> Vec<(i64, i64)> {
        (0..self.d)
            .map(|i| {
                let lo = self.pts.iter().map(|p| p.c[i]).min().unwrap();
                let hi = self.pts.iter().map(|p| p.c[i]).max().unwrap();
                (lo, hi)
            })
            .collect()
    }

    pub fn distinct(&self, dim: usize) -> usize {
        let mut v: Vec<i64> = self.pts.iter().map(|p| p.c[dim]).collect();
        v.sort_unstable();
        v.dedup();
        v.len()
    }
}

/// A closed hyper-rectangle `[lo_i, hi_i]` per dimension.
#[derive(Clone, Debug)]
pub struct Query {
    pub lo: [i64; 3],
    pub hi: [i64; 3],
}

impl Query {
    pub fn contains(&self, p: &Point, d: usize) -> bool {
        (0..d).all(|i| self.lo[i] <= p.c[i] && p.c[i] <= self.hi[i])
    }
}

/// agroramsse's `random_range`: lo uniform in the dimension's data range,
/// hi uniform in [lo, max].
pub fn random_queries(ds: &Dataset, count: usize, seed: u64) -> Vec<Query> {
    let mut rng = StdRng::seed_from_u64(seed);
    let bounds = ds.bounds();
    (0..count)
        .map(|_| {
            let mut q = Query { lo: [0; 3], hi: [0; 3] };
            for (i, (a, b)) in bounds.iter().enumerate() {
                let lo = rng.gen_range(*a..=*b);
                let hi = rng.gen_range(lo..=*b);
                q.lo[i] = lo;
                q.hi[i] = hi;
            }
            q
        })
        .collect()
}

/// Exact answers computed by scanning the dataset.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Truth {
    pub group: Group,
    pub semi: Semi,
    /// Lower median (element of rank floor((cnt-1)/2)), None if empty.
    pub median: Option<i64>,
}

pub fn oracle(ds: &Dataset, q: &Query) -> Truth {
    let mut g = Group::ZERO;
    let mut s = Semi::IDENTITY;
    let mut vals = Vec::new();
    for p in &ds.pts {
        if q.contains(p, ds.d) {
            g = g.add(&Group::of(p.v));
            s = s.merge(&Semi::of(p.v));
            vals.push(p.v);
        }
    }
    let median = if vals.is_empty() {
        None
    } else {
        vals.sort_unstable();
        Some(vals[(vals.len() - 1) / 2])
    };
    Truth { group: g, semi: s, median }
}
