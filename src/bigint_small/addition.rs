use super::CheckedUnsignedAbs::{Negative, Positive};
use super::Sign::{Minus, NoSign, Plus};
use super::{BigIntSmall, UnsignedAbs};

use crate::{IsizePromotion, UsizePromotion};

use core::cmp::Ordering::{Equal, Greater, Less};
use core::iter::Sum;
use core::mem;
use core::ops::{Add, AddAssign};
use num_traits::{CheckedAdd, Zero};

// We want to forward to BigUint::add, but it's not clear how that will go until
// we compare both sign and magnitude.  So we duplicate this body for every
// val/ref combination, deferring that decision to BigUint's own forwarding.
macro_rules! bigint_add {
    ($a:expr, $a_owned:expr, $a_data:expr, $b:expr, $b_owned:expr, $b_data:expr) => {
        match ($a.sign, $b.sign) {
            (_, NoSign) => $a_owned,
            (NoSign, _) => $b_owned,
            // same sign => keep the sign with the sum of magnitudes
            (Plus, Plus) | (Minus, Minus) => BigIntSmall::from_biguint($a.sign, $a_data + $b_data),
            // opposite signs => keep the sign of the larger with the difference of magnitudes
            (Plus, Minus) | (Minus, Plus) => match $a.data.cmp(&$b.data) {
                Less => BigIntSmall::from_biguint($b.sign, $b_data - $a_data),
                Greater => BigIntSmall::from_biguint($a.sign, $a_data - $b_data),
                Equal => Zero::zero(),
            },
        }
    };
}

impl<'a, 'b> Add<&'b BigIntSmall> for &'a BigIntSmall {
    type Output = BigIntSmall;

    #[inline]
    fn add(self, other: &BigIntSmall) -> BigIntSmall {
        bigint_add!(
            self,
            self.clone(),
            &self.data,
            other,
            other.clone(),
            &other.data
        )
    }
}

impl<'a> Add<BigIntSmall> for &'a BigIntSmall {
    type Output = BigIntSmall;

    #[inline]
    fn add(self, other: BigIntSmall) -> BigIntSmall {
        bigint_add!(self, self.clone(), &self.data, other, other, other.data)
    }
}

impl<'a> Add<&'a BigIntSmall> for BigIntSmall {
    type Output = BigIntSmall;

    #[inline]
    fn add(self, other: &BigIntSmall) -> BigIntSmall {
        bigint_add!(self, self, self.data, other, other.clone(), &other.data)
    }
}

impl Add<BigIntSmall> for BigIntSmall {
    type Output = BigIntSmall;

    #[inline]
    fn add(self, other: BigIntSmall) -> BigIntSmall {
        bigint_add!(self, self, self.data, other, other, other.data)
    }
}

impl<'a> AddAssign<&'a BigIntSmall> for BigIntSmall {
    #[inline]
    fn add_assign(&mut self, other: &BigIntSmall) {
        let n = mem::replace(self, BigIntSmall::zero());
        *self = n + other;
    }
}
forward_val_assign!(impl AddAssign for BigIntSmall, add_assign);

promote_all_scalars!(impl Add for BigIntSmall, add);
promote_all_scalars_assign!(impl AddAssign for BigIntSmall, add_assign);
forward_all_scalar_binop_to_val_val_commutative!(impl Add<u32> for BigIntSmall, add);
forward_all_scalar_binop_to_val_val_commutative!(impl Add<u64> for BigIntSmall, add);
forward_all_scalar_binop_to_val_val_commutative!(impl Add<u128> for BigIntSmall, add);

impl Add<u32> for BigIntSmall {
    type Output = BigIntSmall;

    #[inline]
    fn add(self, other: u32) -> BigIntSmall {
        match self.sign {
            NoSign => From::from(other),
            Plus => BigIntSmall::from(self.data + other),
            Minus => match self.data.cmp(&From::from(other)) {
                Equal => Zero::zero(),
                Less => BigIntSmall::from(other - self.data),
                Greater => -BigIntSmall::from(self.data - other),
            },
        }
    }
}

impl AddAssign<u32> for BigIntSmall {
    #[inline]
    fn add_assign(&mut self, other: u32) {
        let n = mem::replace(self, BigIntSmall::zero());
        *self = n + other;
    }
}

impl Add<u64> for BigIntSmall {
    type Output = BigIntSmall;

    #[inline]
    fn add(self, other: u64) -> BigIntSmall {
        match self.sign {
            NoSign => From::from(other),
            Plus => BigIntSmall::from(self.data + other),
            Minus => match self.data.cmp(&From::from(other)) {
                Equal => Zero::zero(),
                Less => BigIntSmall::from(other - self.data),
                Greater => -BigIntSmall::from(self.data - other),
            },
        }
    }
}

impl AddAssign<u64> for BigIntSmall {
    #[inline]
    fn add_assign(&mut self, other: u64) {
        let n = mem::replace(self, BigIntSmall::zero());
        *self = n + other;
    }
}

impl Add<u128> for BigIntSmall {
    type Output = BigIntSmall;

    #[inline]
    fn add(self, other: u128) -> BigIntSmall {
        match self.sign {
            NoSign => BigIntSmall::from(other),
            Plus => BigIntSmall::from(self.data + other),
            Minus => match self.data.cmp(&From::from(other)) {
                Equal => BigIntSmall::zero(),
                Less => BigIntSmall::from(other - self.data),
                Greater => -BigIntSmall::from(self.data - other),
            },
        }
    }
}
impl AddAssign<u128> for BigIntSmall {
    #[inline]
    fn add_assign(&mut self, other: u128) {
        let n = mem::replace(self, BigIntSmall::zero());
        *self = n + other;
    }
}

forward_all_scalar_binop_to_val_val_commutative!(impl Add<i32> for BigIntSmall, add);
forward_all_scalar_binop_to_val_val_commutative!(impl Add<i64> for BigIntSmall, add);
forward_all_scalar_binop_to_val_val_commutative!(impl Add<i128> for BigIntSmall, add);

impl Add<i32> for BigIntSmall {
    type Output = BigIntSmall;

    #[inline]
    fn add(self, other: i32) -> BigIntSmall {
        match other.checked_uabs() {
            Positive(u) => self + u,
            Negative(u) => self - u,
        }
    }
}
impl AddAssign<i32> for BigIntSmall {
    #[inline]
    fn add_assign(&mut self, other: i32) {
        match other.checked_uabs() {
            Positive(u) => *self += u,
            Negative(u) => *self -= u,
        }
    }
}

impl Add<i64> for BigIntSmall {
    type Output = BigIntSmall;

    #[inline]
    fn add(self, other: i64) -> BigIntSmall {
        match other.checked_uabs() {
            Positive(u) => self + u,
            Negative(u) => self - u,
        }
    }
}
impl AddAssign<i64> for BigIntSmall {
    #[inline]
    fn add_assign(&mut self, other: i64) {
        match other.checked_uabs() {
            Positive(u) => *self += u,
            Negative(u) => *self -= u,
        }
    }
}

impl Add<i128> for BigIntSmall {
    type Output = BigIntSmall;

    #[inline]
    fn add(self, other: i128) -> BigIntSmall {
        match other.checked_uabs() {
            Positive(u) => self + u,
            Negative(u) => self - u,
        }
    }
}
impl AddAssign<i128> for BigIntSmall {
    #[inline]
    fn add_assign(&mut self, other: i128) {
        match other.checked_uabs() {
            Positive(u) => *self += u,
            Negative(u) => *self -= u,
        }
    }
}

impl CheckedAdd for BigIntSmall {
    #[inline]
    fn checked_add(&self, v: &BigIntSmall) -> Option<BigIntSmall> {
        Some(self.add(v))
    }
}

impl_sum_iter_type!(BigIntSmall);
