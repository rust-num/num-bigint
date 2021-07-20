use super::CheckedUnsignedAbs::{Negative, Positive};
use super::{BigIntSmall, BigUint, UnsignedAbs};

use crate::{IsizePromotion, UsizePromotion};

use core::iter::Product;
use core::ops::{Mul, MulAssign};
use num_traits::{CheckedMul, One};

macro_rules! impl_mul {
    ($(impl<$($a:lifetime),*> Mul<$Other:ty> for $Self:ty;)*) => {$(
        impl<$($a),*> Mul<$Other> for $Self {
            type Output = BigIntSmall;

            #[inline]
            fn mul(self, other: $Other) -> BigIntSmall {
                use BigIntSmall::*;
                match (&self, &other) {
                    (PlusSmall(a), PlusSmall(b)) => {
                        let c = *a as u128 * *b as u128;
                        return BigIntSmall::from(c);
                    }
                    (MinusSmall(a), MinusSmall(b)) => {
                        let c = *a as u128 * *b as u128;
                        return -BigIntSmall::from(c);
                    }
                    _ => {
                        // dbg!(&self, &other);
                    }
                }
                // automatically match value/ref
                let self_data = self.data();
                let other_data = other.data();
                // BigUint multiplication doesn't reuse buffers so it doesn't
                // hurt performance to always use references.
                BigIntSmall::from_biguint(self.sign() * other.sign(), &self_data as &BigUint * &other_data as &BigUint)
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
    ($(impl<$($a:lifetime),*> MulAssign<$Other:ty> for BigIntSmall;)*) => {$(
        impl<$($a),*> MulAssign<$Other> for BigIntSmall {
            #[inline]
            fn mul_assign(&mut self, other: $Other) {
                // automatically match value/ref
                let sign = self.sign() * other.sign();
                let uint = &self.data() as &BigUint * &other.data() as &BigUint;
                *self = BigIntSmall::from_biguint( sign, uint)

                // let BigIntSmall { data: y, .. } = other;
                // self.data *= y;
                // if self.data.is_zero() {
                //     self.sign = NoSign;
                // } else {
                //     self.sign = self.sign * other.sign;
                // }
            }
        }
    )*}
}
impl_mul_assign! {
    impl<> MulAssign<BigIntSmall> for BigIntSmall;
    impl<'a> MulAssign<&'a BigIntSmall> for BigIntSmall;
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
        BigIntSmall::from_biguint(self.sign(), &self.data() as &BigUint * other)
    }
}

impl MulAssign<u32> for BigIntSmall {
    #[inline]
    fn mul_assign(&mut self, other: u32) {
        let owned = self.take();
        let (sign, mut uint) = owned.into_parts();
        uint *= other;
        *self = BigIntSmall::from_biguint(sign, uint)
        // self.data *= other;
        // if self.data.is_zero() {
        //     self.sign = NoSign;
        // }
    }
}

impl Mul<u64> for BigIntSmall {
    type Output = BigIntSmall;

    #[inline]
    fn mul(self, other: u64) -> BigIntSmall {
        BigIntSmall::from_biguint(self.sign(), &self.data() as &BigUint * other)
    }
}

impl MulAssign<u64> for BigIntSmall {
    #[inline]
    fn mul_assign(&mut self, other: u64) {
        let owned = self.take();
        let (sign, mut uint) = owned.into_parts();
        uint *= other;
        *self = BigIntSmall::from_biguint(sign, uint)
        // self.data *= other;
        // if self.data.is_zero() {
        //     self.sign = NoSign;
        // }
    }
}

impl Mul<u128> for BigIntSmall {
    type Output = BigIntSmall;

    #[inline]
    fn mul(self, other: u128) -> BigIntSmall {
        BigIntSmall::from_biguint(self.sign(), &self.data() as &BigUint * other)
    }
}

impl MulAssign<u128> for BigIntSmall {
    #[inline]
    fn mul_assign(&mut self, other: u128) {
        let owned = self.take();
        let (sign, mut uint) = owned.into_parts();
        uint *= other;
        *self = BigIntSmall::from_biguint(sign, uint)
        // self.data *= other;
        // if self.data.is_zero() {
        //     self.sign = NoSign;
        // }
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
                // *self.mut_sign() = -self.sign();
                // self.data *= u;
                let owned = self.take();
                let (sign, mut uint) = owned.into_parts();
                uint *= u;
                *self = BigIntSmall::from_biguint(-sign, uint)
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
                // *self.mut_sign() = -self.sign();
                // self.data *= u;
                let owned = self.take();
                let (sign, mut uint) = owned.into_parts();
                uint *= u;
                *self = BigIntSmall::from_biguint(-sign, uint)
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
                // *self.mut_sign() = -self.sign();
                // self.data *= u;
                let owned = self.take();
                let (sign, mut uint) = owned.into_parts();
                uint *= u;
                *self = BigIntSmall::from_biguint(-sign, uint)
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
