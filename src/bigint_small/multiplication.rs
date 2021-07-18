use super::CheckedUnsignedAbs::{Negative, Positive};
use super::Sign::{self, Minus, NoSign, Plus};
use super::{BigIntSmall, UnsignedAbs};

use crate::{IsizePromotion, UsizePromotion};

use core::iter::Product;
use core::ops::{Mul, MulAssign};
use num_traits::{CheckedMul, One, Zero};

impl Mul<Sign> for Sign {
    type Output = Sign;

    #[inline]
    fn mul(self, other: Sign) -> Sign {
        match (self, other) {
            (NoSign, _) | (_, NoSign) => NoSign,
            (Plus, Plus) | (Minus, Minus) => Plus,
            (Plus, Minus) | (Minus, Plus) => Minus,
        }
    }
}

macro_rules! impl_mul {
    ($(impl<$($a:lifetime),*> Mul<$Other:ty> for $Self:ty;)*) => {$(
        impl<$($a),*> Mul<$Other> for $Self {
            type Output = BigIntSmall;

            #[inline]
            fn mul(self, other: $Other) -> BigIntSmall {
                // automatically match value/ref
                let BigIntSmall { data: x, .. } = self;
                let BigIntSmall { data: y, .. } = other;
                BigIntSmall::from_biguint(self.sign * other.sign, x * y)
            }
        }
    )*}
}
impl_mul! {
    impl<> Mul<BigIntSmall> for BigIntSmall;
    impl<'b> Mul<&'b BigIntSmall> for BigIntSmall;
    impl<'a> Mul<BigIntSmall> for &'a BigIntSmall;
    impl<'a, 'b> Mul<&'b BigIntSmall> for &'a BigIntSmall;
}

macro_rules! impl_mul_assign {
    ($(impl<$($a:lifetime),*> MulAssign<$Other:ty> for BigInt;)*) => {$(
        impl<$($a),*> MulAssign<$Other> for BigIntSmall {
            #[inline]
            fn mul_assign(&mut self, other: $Other) {
                // automatically match value/ref
                let BigIntSmall { data: y, .. } = other;
                self.data *= y;
                if self.data.is_zero() {
                    self.sign = NoSign;
                } else {
                    self.sign = self.sign * other.sign;
                }
            }
        }
    )*}
}
impl_mul_assign! {
    impl<> MulAssign<BigIntSmall> for BigInt;
    impl<'a> MulAssign<&'a BigIntSmall> for BigInt;
}

promote_all_scalars!(impl Mul for BigIntSmall, mul);
promote_all_scalars_assign!(impl MulAssign for BigIntSmall, mul_assign);
forward_all_scalar_binop_to_val_val_commutative!(impl Mul<u32> for BigIntSmall, mul);
forward_all_scalar_binop_to_val_val_commutative!(impl Mul<u64> for BigIntSmall, mul);
forward_all_scalar_binop_to_val_val_commutative!(impl Mul<u128> for BigIntSmall, mul);

impl Mul<u32> for BigIntSmall {
    type Output = BigIntSmall;

    #[inline]
    fn mul(self, other: u32) -> BigIntSmall {
        BigIntSmall::from_biguint(self.sign, self.data * other)
    }
}

impl MulAssign<u32> for BigIntSmall {
    #[inline]
    fn mul_assign(&mut self, other: u32) {
        self.data *= other;
        if self.data.is_zero() {
            self.sign = NoSign;
        }
    }
}

impl Mul<u64> for BigIntSmall {
    type Output = BigIntSmall;

    #[inline]
    fn mul(self, other: u64) -> BigIntSmall {
        BigIntSmall::from_biguint(self.sign, self.data * other)
    }
}

impl MulAssign<u64> for BigIntSmall {
    #[inline]
    fn mul_assign(&mut self, other: u64) {
        self.data *= other;
        if self.data.is_zero() {
            self.sign = NoSign;
        }
    }
}

impl Mul<u128> for BigIntSmall {
    type Output = BigIntSmall;

    #[inline]
    fn mul(self, other: u128) -> BigIntSmall {
        BigIntSmall::from_biguint(self.sign, self.data * other)
    }
}

impl MulAssign<u128> for BigIntSmall {
    #[inline]
    fn mul_assign(&mut self, other: u128) {
        self.data *= other;
        if self.data.is_zero() {
            self.sign = NoSign;
        }
    }
}

forward_all_scalar_binop_to_val_val_commutative!(impl Mul<i32> for BigIntSmall, mul);
forward_all_scalar_binop_to_val_val_commutative!(impl Mul<i64> for BigIntSmall, mul);
forward_all_scalar_binop_to_val_val_commutative!(impl Mul<i128> for BigIntSmall, mul);

impl Mul<i32> for BigIntSmall {
    type Output = BigIntSmall;

    #[inline]
    fn mul(self, other: i32) -> BigIntSmall {
        match other.checked_uabs() {
            Positive(u) => self * u,
            Negative(u) => -self * u,
        }
    }
}

impl MulAssign<i32> for BigIntSmall {
    #[inline]
    fn mul_assign(&mut self, other: i32) {
        match other.checked_uabs() {
            Positive(u) => *self *= u,
            Negative(u) => {
                self.sign = -self.sign;
                self.data *= u;
            }
        }
    }
}

impl Mul<i64> for BigIntSmall {
    type Output = BigIntSmall;

    #[inline]
    fn mul(self, other: i64) -> BigIntSmall {
        match other.checked_uabs() {
            Positive(u) => self * u,
            Negative(u) => -self * u,
        }
    }
}

impl MulAssign<i64> for BigIntSmall {
    #[inline]
    fn mul_assign(&mut self, other: i64) {
        match other.checked_uabs() {
            Positive(u) => *self *= u,
            Negative(u) => {
                self.sign = -self.sign;
                self.data *= u;
            }
        }
    }
}

impl Mul<i128> for BigIntSmall {
    type Output = BigIntSmall;

    #[inline]
    fn mul(self, other: i128) -> BigIntSmall {
        match other.checked_uabs() {
            Positive(u) => self * u,
            Negative(u) => -self * u,
        }
    }
}

impl MulAssign<i128> for BigIntSmall {
    #[inline]
    fn mul_assign(&mut self, other: i128) {
        match other.checked_uabs() {
            Positive(u) => *self *= u,
            Negative(u) => {
                self.sign = -self.sign;
                self.data *= u;
            }
        }
    }
}

impl CheckedMul for BigIntSmall {
    #[inline]
    fn checked_mul(&self, v: &BigIntSmall) -> Option<BigIntSmall> {
        Some(self.mul(v))
    }
}

impl_product_iter_type!(BigIntSmall);
