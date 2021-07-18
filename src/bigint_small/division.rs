use super::CheckedUnsignedAbs::{Negative, Positive};
use super::{BigIntSmall, BigUint, UnsignedAbs};

use crate::{IsizePromotion, UsizePromotion};

use core::ops::{Div, DivAssign, Rem, RemAssign};
use num_integer::Integer;
use num_traits::{CheckedDiv, ToPrimitive, Zero};

forward_all_binop_to_ref_ref!(impl Div for BigIntSmall, div);

impl<'a, 'b> Div<&'b BigIntSmall> for &'a BigIntSmall {
    type Output = BigIntSmall;

    #[inline]
    fn div(self, other: &BigIntSmall) -> BigIntSmall {
        let (q, _) = self.div_rem(other);
        q
    }
}

impl<'a> DivAssign<&'a BigIntSmall> for BigIntSmall {
    #[inline]
    fn div_assign(&mut self, other: &BigIntSmall) {
        *self = &*self / other;
    }
}
forward_val_assign!(impl DivAssign for BigIntSmall, div_assign);

promote_all_scalars!(impl Div for BigIntSmall, div);
promote_all_scalars_assign!(impl DivAssign for BigIntSmall, div_assign);
forward_all_scalar_binop_to_val_val!(impl Div<u32> for BigIntSmall, div);
forward_all_scalar_binop_to_val_val!(impl Div<u64> for BigIntSmall, div);
forward_all_scalar_binop_to_val_val!(impl Div<u128> for BigIntSmall, div);

impl Div<u32> for BigIntSmall {
    type Output = BigIntSmall;

    #[inline]
    fn div(self, other: u32) -> BigIntSmall {
        BigIntSmall::from_biguint(self.sign(), &self.data() as &BigUint / other)
    }
}

impl DivAssign<u32> for BigIntSmall {
    #[inline]
    fn div_assign(&mut self, other: u32) {
        let owned = std::mem::replace(self, BigIntSmall::zero());
        let (sign, mut uint) = owned.into_parts();
        uint /= other;
        *self = BigIntSmall::from_biguint(sign, uint)
        // self.data /= other;
        // if self.data.is_zero() {
        //     *self.mut_sign() = NoSign;
        // }
    }
}

impl Div<BigIntSmall> for u32 {
    type Output = BigIntSmall;

    #[inline]
    fn div(self, other: BigIntSmall) -> BigIntSmall {
        BigIntSmall::from_biguint(other.sign, self / other.data)
    }
}

impl Div<u64> for BigIntSmall {
    type Output = BigIntSmall;

    #[inline]
    fn div(self, other: u64) -> BigIntSmall {
        BigIntSmall::from_biguint(self.sign(), &self.data() as &BigUint / other)
    }
}

impl DivAssign<u64> for BigIntSmall {
    #[inline]
    fn div_assign(&mut self, other: u64) {
        let owned = std::mem::replace(self, BigIntSmall::zero());
        let (sign, mut uint) = owned.into_parts();
        uint /= other;
        *self = BigIntSmall::from_biguint(sign, uint)
        // self.data /= other;
        // if self.data.is_zero() {
        //     self.sign = NoSign;
        // }
    }
}

impl Div<BigIntSmall> for u64 {
    type Output = BigIntSmall;

    #[inline]
    fn div(self, other: BigIntSmall) -> BigIntSmall {
        BigIntSmall::from_biguint(other.sign, self / other.data)
    }
}

impl Div<u128> for BigIntSmall {
    type Output = BigIntSmall;

    #[inline]
    fn div(self, other: u128) -> BigIntSmall {
        BigIntSmall::from_biguint(self.sign(), &self.data() as &BigUint / other)
    }
}

impl DivAssign<u128> for BigIntSmall {
    #[inline]
    fn div_assign(&mut self, other: u128) {
        let owned = std::mem::replace(self, BigIntSmall::zero());
        let (sign, mut uint) = owned.into_parts();
        uint /= other;
        *self = BigIntSmall::from_biguint(sign, uint)
        // self.data /= other;
        // if self.data.is_zero() {
        //     self.sign = NoSign;
        // }
    }
}

impl Div<BigIntSmall> for u128 {
    type Output = BigIntSmall;

    #[inline]
    fn div(self, other: BigIntSmall) -> BigIntSmall {
        BigIntSmall::from_biguint(other.sign(), self / &other.data() as &BigUint)
    }
}

forward_all_scalar_binop_to_val_val!(impl Div<i32> for BigIntSmall, div);
forward_all_scalar_binop_to_val_val!(impl Div<i64> for BigIntSmall, div);
forward_all_scalar_binop_to_val_val!(impl Div<i128> for BigIntSmall, div);

impl Div<i32> for BigIntSmall {
    type Output = BigIntSmall;

    #[inline]
    fn div(self, other: i32) -> BigIntSmall {
        match other.checked_uabs() {
            Positive(u) => self / u,
            Negative(u) => -self / u,
        }
    }
}

impl DivAssign<i32> for BigIntSmall {
    #[inline]
    fn div_assign(&mut self, other: i32) {
        match other.checked_uabs() {
            Positive(u) => *self /= u,
            Negative(u) => {
                *self.mut_sign() = -self.sign();
                *self /= u;
            }
        }
    }
}

impl Div<BigIntSmall> for i32 {
    type Output = BigIntSmall;

    #[inline]
    fn div(self, other: BigIntSmall) -> BigIntSmall {
        match self.checked_uabs() {
            Positive(u) => u / other,
            Negative(u) => u / -other,
        }
    }
}

impl Div<i64> for BigIntSmall {
    type Output = BigIntSmall;

    #[inline]
    fn div(self, other: i64) -> BigIntSmall {
        match other.checked_uabs() {
            Positive(u) => self / u,
            Negative(u) => -self / u,
        }
    }
}

impl DivAssign<i64> for BigIntSmall {
    #[inline]
    fn div_assign(&mut self, other: i64) {
        match other.checked_uabs() {
            Positive(u) => *self /= u,
            Negative(u) => {
                *self.mut_sign() = -self.sign();
                *self /= u;
            }
        }
    }
}

impl Div<BigIntSmall> for i64 {
    type Output = BigIntSmall;

    #[inline]
    fn div(self, other: BigIntSmall) -> BigIntSmall {
        match self.checked_uabs() {
            Positive(u) => u / other,
            Negative(u) => u / -other,
        }
    }
}

impl Div<i128> for BigIntSmall {
    type Output = BigIntSmall;

    #[inline]
    fn div(self, other: i128) -> BigIntSmall {
        match other.checked_uabs() {
            Positive(u) => self / u,
            Negative(u) => -self / u,
        }
    }
}

impl DivAssign<i128> for BigIntSmall {
    #[inline]
    fn div_assign(&mut self, other: i128) {
        match other.checked_uabs() {
            Positive(u) => *self /= u,
            Negative(u) => {
                *self.mut_sign() = -self.sign();
                *self /= u;
            }
        }
    }
}

impl Div<BigIntSmall> for i128 {
    type Output = BigIntSmall;

    #[inline]
    fn div(self, other: BigIntSmall) -> BigIntSmall {
        match self.checked_uabs() {
            Positive(u) => u / other,
            Negative(u) => u / -other,
        }
    }
}

forward_all_binop_to_ref_ref!(impl Rem for BigIntSmall, rem);

impl<'a, 'b> Rem<&'b BigIntSmall> for &'a BigIntSmall {
    type Output = BigIntSmall;

    #[inline]
    fn rem(self, other: &BigIntSmall) -> BigIntSmall {
        if let Some(other) = other.to_u32() {
            self % other
        } else if let Some(other) = other.to_i32() {
            self % other
        } else {
            let (_, r) = self.div_rem(other);
            r
        }
    }
}

impl<'a> RemAssign<&'a BigIntSmall> for BigIntSmall {
    #[inline]
    fn rem_assign(&mut self, other: &BigIntSmall) {
        *self = &*self % other;
    }
}
forward_val_assign!(impl RemAssign for BigIntSmall, rem_assign);

promote_all_scalars!(impl Rem for BigIntSmall, rem);
promote_all_scalars_assign!(impl RemAssign for BigIntSmall, rem_assign);
forward_all_scalar_binop_to_val_val!(impl Rem<u32> for BigIntSmall, rem);
forward_all_scalar_binop_to_val_val!(impl Rem<u64> for BigIntSmall, rem);
forward_all_scalar_binop_to_val_val!(impl Rem<u128> for BigIntSmall, rem);

impl Rem<u32> for BigIntSmall {
    type Output = BigIntSmall;

    #[inline]
    fn rem(self, other: u32) -> BigIntSmall {
        BigIntSmall::from_biguint(self.sign(), &self.data() as &BigUint % other)
    }
}

impl RemAssign<u32> for BigIntSmall {
    #[inline]
    fn rem_assign(&mut self, other: u32) {
        let owned = std::mem::replace(self, BigIntSmall::zero());
        let (sign, mut uint) = owned.into_parts();
        uint %= other;
        *self = BigIntSmall::from_biguint(sign, uint)
        // self.data %= other;
        // if self.data.is_zero() {
        //     self.sign = NoSign;
        // }
    }
}

impl Rem<BigIntSmall> for u32 {
    type Output = BigIntSmall;

    #[inline]
    fn rem(self, other: BigIntSmall) -> BigIntSmall {
        BigIntSmall::from(self % &other.data() as &BigUint)
    }
}

impl Rem<u64> for BigIntSmall {
    type Output = BigIntSmall;

    #[inline]
    fn rem(self, other: u64) -> BigIntSmall {
        BigIntSmall::from_biguint(self.sign(), &self.data() as &BigUint % other)
    }
}

impl RemAssign<u64> for BigIntSmall {
    #[inline]
    fn rem_assign(&mut self, other: u64) {
        let owned = std::mem::replace(self, BigIntSmall::zero());
        let (sign, mut uint) = owned.into_parts();
        uint %= other;
        *self = BigIntSmall::from_biguint(sign, uint)
        // self.data %= other;
        // if self.data.is_zero() {
        //     self.sign = NoSign;
        // }
    }
}

impl Rem<BigIntSmall> for u64 {
    type Output = BigIntSmall;

    #[inline]
    fn rem(self, other: BigIntSmall) -> BigIntSmall {
        BigIntSmall::from(self % &other.data() as &BigUint)
    }
}

impl Rem<u128> for BigIntSmall {
    type Output = BigIntSmall;

    #[inline]
    fn rem(self, other: u128) -> BigIntSmall {
        BigIntSmall::from_biguint(self.sign(), &self.data() as &BigUint % other)
    }
}

impl RemAssign<u128> for BigIntSmall {
    #[inline]
    fn rem_assign(&mut self, other: u128) {
        let owned = std::mem::replace(self, BigIntSmall::zero());
        let (sign, mut uint) = owned.into_parts();
        uint %= other;
        *self = BigIntSmall::from_biguint(sign, uint)
        // self.data %= other;
        // if self.data.is_zero() {
        //     self.sign = NoSign;
        // }
    }
}

impl Rem<BigIntSmall> for u128 {
    type Output = BigIntSmall;

    #[inline]
    fn rem(self, other: BigIntSmall) -> BigIntSmall {
        BigIntSmall::from(self % other.data)
    }
}

forward_all_scalar_binop_to_val_val!(impl Rem<i32> for BigIntSmall, rem);
forward_all_scalar_binop_to_val_val!(impl Rem<i64> for BigIntSmall, rem);
forward_all_scalar_binop_to_val_val!(impl Rem<i128> for BigIntSmall, rem);

impl Rem<i32> for BigIntSmall {
    type Output = BigIntSmall;

    #[inline]
    fn rem(self, other: i32) -> BigIntSmall {
        self % other.uabs()
    }
}

impl RemAssign<i32> for BigIntSmall {
    #[inline]
    fn rem_assign(&mut self, other: i32) {
        *self %= other.uabs();
    }
}

impl Rem<BigIntSmall> for i32 {
    type Output = BigIntSmall;

    #[inline]
    fn rem(self, other: BigIntSmall) -> BigIntSmall {
        match self.checked_uabs() {
            Positive(u) => u % other,
            Negative(u) => -(u % other),
        }
    }
}

impl Rem<i64> for BigIntSmall {
    type Output = BigIntSmall;

    #[inline]
    fn rem(self, other: i64) -> BigIntSmall {
        self % other.uabs()
    }
}

impl RemAssign<i64> for BigIntSmall {
    #[inline]
    fn rem_assign(&mut self, other: i64) {
        *self %= other.uabs();
    }
}

impl Rem<BigIntSmall> for i64 {
    type Output = BigIntSmall;

    #[inline]
    fn rem(self, other: BigIntSmall) -> BigIntSmall {
        match self.checked_uabs() {
            Positive(u) => u % other,
            Negative(u) => -(u % other),
        }
    }
}

impl Rem<i128> for BigIntSmall {
    type Output = BigIntSmall;

    #[inline]
    fn rem(self, other: i128) -> BigIntSmall {
        self % other.uabs()
    }
}

impl RemAssign<i128> for BigIntSmall {
    #[inline]
    fn rem_assign(&mut self, other: i128) {
        *self %= other.uabs();
    }
}

impl Rem<BigIntSmall> for i128 {
    type Output = BigIntSmall;

    #[inline]
    fn rem(self, other: BigIntSmall) -> BigIntSmall {
        match self.checked_uabs() {
            Positive(u) => u % other,
            Negative(u) => -(u % other),
        }
    }
}

impl CheckedDiv for BigIntSmall {
    #[inline]
    fn checked_div(&self, v: &BigIntSmall) -> Option<BigIntSmall> {
        if v.is_zero() {
            return None;
        }
        Some(self.div(v))
    }
}
