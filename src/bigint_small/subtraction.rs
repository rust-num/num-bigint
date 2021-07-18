use super::CheckedUnsignedAbs::{Negative, Positive};
use super::Sign::{Minus, NoSign, Plus};
use super::{BigIntSmall, BigUint, UnsignedAbs};

use crate::{IsizePromotion, UsizePromotion};

use core::cmp::Ordering::{Equal, Greater, Less};
use core::mem;
use core::ops::{Sub, SubAssign};
use num_traits::{CheckedSub, Zero};

// We want to forward to BigUint::sub, but it's not clear how that will go until
// we compare both sign and magnitude.  So we duplicate this body for every
// val/ref combination, deferring that decision to BigUint's own forwarding.
macro_rules! bigint_sub {
    ($a:expr, $a_owned:expr, $a_data:expr, $b:expr, $b_owned:expr, $b_data:expr) => {
        match ($a.sign(), $b.sign()) {
            (_, NoSign) => $a_owned,
            (NoSign, _) => -$b_owned,
            // opposite signs => keep the sign of the left with the sum of magnitudes
            (Plus, Minus) | (Minus, Plus) => {
                BigIntSmall::from_biguint($a.sign(), $a_data + $b_data)
            }
            // same sign => keep or toggle the sign of the left with the difference of magnitudes
            (Plus, Plus) | (Minus, Minus) => match $a.data().cmp(&$b.data()) {
                Less => BigIntSmall::from_biguint(-$a.sign(), $b_data - $a_data),
                Greater => BigIntSmall::from_biguint($a.sign(), $a_data - $b_data),
                Equal => Zero::zero(),
            },
        }
    };
}

impl<'a, 'b> Sub<&'b BigIntSmall> for &'a BigIntSmall {
    type Output = BigIntSmall;

    #[inline]
    fn sub(self, other: &BigIntSmall) -> BigIntSmall {
        bigint_sub!(
            self,
            self.clone(),
            &self.data() as &BigUint,
            other,
            other.clone(),
            &other.data() as &BigUint
        )
    }
}

impl<'a> Sub<BigIntSmall> for &'a BigIntSmall {
    type Output = BigIntSmall;

    #[inline]
    fn sub(self, other: BigIntSmall) -> BigIntSmall {
        bigint_sub!(
            self,
            self.clone(),
            &self.data() as &BigUint,
            other,
            other,
            other.to_biguint_unchecked()
        )
    }
}

impl<'a> Sub<&'a BigIntSmall> for BigIntSmall {
    type Output = BigIntSmall;

    #[inline]
    fn sub(self, other: &BigIntSmall) -> BigIntSmall {
        bigint_sub!(
            self,
            self,
            self.to_biguint_unchecked(),
            other,
            other.clone(),
            &other.data() as &BigUint
        )
    }
}

impl Sub<BigIntSmall> for BigIntSmall {
    type Output = BigIntSmall;

    #[inline]
    fn sub(self, other: BigIntSmall) -> BigIntSmall {
        bigint_sub!(
            self,
            self,
            self.to_biguint_unchecked(),
            other,
            other,
            other.to_biguint_unchecked()
        )
    }
}

impl<'a> SubAssign<&'a BigIntSmall> for BigIntSmall {
    #[inline]
    fn sub_assign(&mut self, other: &BigIntSmall) {
        let n = mem::replace(self, BigIntSmall::zero());
        *self = n - other;
    }
}
forward_val_assign!(impl SubAssign for BigIntSmall, sub_assign);

promote_all_scalars!(impl Sub for BigIntSmall, sub);
promote_all_scalars_assign!(impl SubAssign for BigIntSmall, sub_assign);
forward_all_scalar_binop_to_val_val!(impl Sub<u32> for BigIntSmall, sub);
forward_all_scalar_binop_to_val_val!(impl Sub<u64> for BigIntSmall, sub);
forward_all_scalar_binop_to_val_val!(impl Sub<u128> for BigIntSmall, sub);

impl Sub<u32> for BigIntSmall {
    type Output = BigIntSmall;

    #[inline]
    fn sub(self, other: u32) -> BigIntSmall {
        match self.sign() {
            NoSign => -BigIntSmall::from(other),
            Minus => -BigIntSmall::from(self.to_biguint_unchecked() + other),
            Plus => match self.data().cmp(&BigIntSmall::from(other).data()) {
                Equal => Zero::zero(),
                Greater => BigIntSmall::from(self.to_biguint_unchecked() - other),
                Less => -BigIntSmall::from(other - self.to_biguint_unchecked()),
            },
        }
    }
}
impl SubAssign<u32> for BigIntSmall {
    #[inline]
    fn sub_assign(&mut self, other: u32) {
        let n = mem::replace(self, BigIntSmall::zero());
        *self = n - other;
    }
}

impl Sub<BigIntSmall> for u32 {
    type Output = BigIntSmall;

    #[inline]
    fn sub(self, other: BigIntSmall) -> BigIntSmall {
        -(other - self)
    }
}

impl Sub<BigIntSmall> for u64 {
    type Output = BigIntSmall;

    #[inline]
    fn sub(self, other: BigIntSmall) -> BigIntSmall {
        -(other - self)
    }
}

impl Sub<BigIntSmall> for u128 {
    type Output = BigIntSmall;

    #[inline]
    fn sub(self, other: BigIntSmall) -> BigIntSmall {
        -(other - self)
    }
}

impl Sub<u64> for BigIntSmall {
    type Output = BigIntSmall;

    #[inline]
    fn sub(self, other: u64) -> BigIntSmall {
        match self.sign() {
            NoSign => -BigIntSmall::from(other),
            Minus => -BigIntSmall::from(self.to_biguint_unchecked() + other),
            Plus => match self.data().cmp(&BigIntSmall::from(other).data()) {
                Equal => Zero::zero(),
                Greater => BigIntSmall::from(self.to_biguint_unchecked() - other),
                Less => -BigIntSmall::from(other - self.to_biguint_unchecked()),
            },
        }
    }
}

impl SubAssign<u64> for BigIntSmall {
    #[inline]
    fn sub_assign(&mut self, other: u64) {
        let n = mem::replace(self, BigIntSmall::zero());
        *self = n - other;
    }
}

impl Sub<u128> for BigIntSmall {
    type Output = BigIntSmall;

    #[inline]
    fn sub(self, other: u128) -> BigIntSmall {
        match self.sign() {
            NoSign => -BigIntSmall::from(other),
            Minus => -BigIntSmall::from(self.to_biguint_unchecked() + other),
            Plus => match self.data().cmp(&BigIntSmall::from(other).data()) {
                Equal => Zero::zero(),
                Greater => BigIntSmall::from(self.to_biguint_unchecked() - other),
                Less => -BigIntSmall::from(other - self.to_biguint_unchecked()),
            },
        }
    }
}

impl SubAssign<u128> for BigIntSmall {
    #[inline]
    fn sub_assign(&mut self, other: u128) {
        let n = mem::replace(self, BigIntSmall::zero());
        *self = n - other;
    }
}

forward_all_scalar_binop_to_val_val!(impl Sub<i32> for BigIntSmall, sub);
forward_all_scalar_binop_to_val_val!(impl Sub<i64> for BigIntSmall, sub);
forward_all_scalar_binop_to_val_val!(impl Sub<i128> for BigIntSmall, sub);

impl Sub<i32> for BigIntSmall {
    type Output = BigIntSmall;

    #[inline]
    fn sub(self, other: i32) -> BigIntSmall {
        match other.checked_uabs() {
            Positive(u) => self - u,
            Negative(u) => self + u,
        }
    }
}
impl SubAssign<i32> for BigIntSmall {
    #[inline]
    fn sub_assign(&mut self, other: i32) {
        match other.checked_uabs() {
            Positive(u) => *self -= u,
            Negative(u) => *self += u,
        }
    }
}

impl Sub<BigIntSmall> for i32 {
    type Output = BigIntSmall;

    #[inline]
    fn sub(self, other: BigIntSmall) -> BigIntSmall {
        match self.checked_uabs() {
            Positive(u) => u - other,
            Negative(u) => -other - u,
        }
    }
}

impl Sub<i64> for BigIntSmall {
    type Output = BigIntSmall;

    #[inline]
    fn sub(self, other: i64) -> BigIntSmall {
        match other.checked_uabs() {
            Positive(u) => self - u,
            Negative(u) => self + u,
        }
    }
}
impl SubAssign<i64> for BigIntSmall {
    #[inline]
    fn sub_assign(&mut self, other: i64) {
        match other.checked_uabs() {
            Positive(u) => *self -= u,
            Negative(u) => *self += u,
        }
    }
}

impl Sub<BigIntSmall> for i64 {
    type Output = BigIntSmall;

    #[inline]
    fn sub(self, other: BigIntSmall) -> BigIntSmall {
        match self.checked_uabs() {
            Positive(u) => u - other,
            Negative(u) => -other - u,
        }
    }
}

impl Sub<i128> for BigIntSmall {
    type Output = BigIntSmall;

    #[inline]
    fn sub(self, other: i128) -> BigIntSmall {
        match other.checked_uabs() {
            Positive(u) => self - u,
            Negative(u) => self + u,
        }
    }
}

impl SubAssign<i128> for BigIntSmall {
    #[inline]
    fn sub_assign(&mut self, other: i128) {
        match other.checked_uabs() {
            Positive(u) => *self -= u,
            Negative(u) => *self += u,
        }
    }
}

impl Sub<BigIntSmall> for i128 {
    type Output = BigIntSmall;

    #[inline]
    fn sub(self, other: BigIntSmall) -> BigIntSmall {
        match self.checked_uabs() {
            Positive(u) => u - other,
            Negative(u) => -other - u,
        }
    }
}

impl CheckedSub for BigIntSmall {
    #[inline]
    fn checked_sub(&self, v: &BigIntSmall) -> Option<BigIntSmall> {
        Some(self.sub(v))
    }
}
