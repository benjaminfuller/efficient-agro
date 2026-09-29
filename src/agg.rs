//! Aggregate algebras.
//!
//! * [`Group`]: COUNT, SUM and sum of squares, an abelian group, so prefix
//!   differences are valid (AVG and STD are derived client-side).
//! * [`Semi`]: MIN and MAX, a semigroup with identity (no inverses).

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Group {
    pub cnt: i64,
    pub sum: i64,
    pub sq: i64,
}

impl Group {
    pub const ZERO: Self = Self { cnt: 0, sum: 0, sq: 0 };

    pub fn of(v: i64) -> Self {
        Self { cnt: 1, sum: v, sq: v * v }
    }

    pub fn add(&self, o: &Self) -> Self {
        Self {
            cnt: self.cnt + o.cnt,
            sum: self.sum + o.sum,
            sq: self.sq + o.sq,
        }
    }

    pub fn sub(&self, o: &Self) -> Self {
        Self {
            cnt: self.cnt - o.cnt,
            sum: self.sum - o.sum,
            sq: self.sq - o.sq,
        }
    }

    pub fn avg(&self) -> Option<f64> {
        (self.cnt > 0).then(|| self.sum as f64 / self.cnt as f64)
    }

    /// Population standard deviation.
    pub fn std(&self) -> Option<f64> {
        (self.cnt > 0).then(|| {
            let n = self.cnt as f64;
            let m = self.sum as f64 / n;
            (self.sq as f64 / n - m * m).max(0.0).sqrt()
        })
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Semi {
    pub min: i64,
    pub max: i64,
}

impl Default for Semi {
    fn default() -> Self {
        Self::IDENTITY
    }
}

impl Semi {
    pub const IDENTITY: Self = Self {
        min: i64::MAX,
        max: i64::MIN,
    };

    pub fn of(v: i64) -> Self {
        Self { min: v, max: v }
    }

    pub fn merge(&self, o: &Self) -> Self {
        Self {
            min: self.min.min(o.min),
            max: self.max.max(o.max),
        }
    }

    pub fn is_empty(&self) -> bool {
        self.min > self.max
    }
}
