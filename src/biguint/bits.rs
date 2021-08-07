use super::{BigUint, IntDigits};

use crate::big_digit::{self, BigDigit};
use crate::{IsizePromotion, UsizePromotion};

use core::ops::{BitAnd, BitAndAssign, BitOr, BitOrAssign, BitXor, BitXorAssign};
use num_traits::Zero;

forward_val_val_binop!(impl BitAnd for BigUint, bitand);
forward_ref_val_binop!(impl BitAnd for BigUint, bitand);

// do not use forward_ref_ref_binop_commutative! for bitand so that we can
// clone the smaller value rather than the larger, avoiding over-allocation
impl<'a, 'b> BitAnd<&'b BigUint> for &'a BigUint {
    type Output = BigUint;

    #[inline]
    fn bitand(self, other: &BigUint) -> BigUint {
        // forward to val-ref, choosing the smaller to clone
        if self.data.len() <= other.data.len() {
            self.clone() & other
        } else {
            other.clone() & self
        }
    }
}

forward_val_assign!(impl BitAndAssign for BigUint, bitand_assign);

impl<'a> BitAnd<&'a BigUint> for BigUint {
    type Output = BigUint;

    #[inline]
    fn bitand(mut self, other: &BigUint) -> BigUint {
        self &= other;
        self
    }
}
impl<'a> BitAndAssign<&'a BigUint> for BigUint {
    #[inline]
    fn bitand_assign(&mut self, other: &BigUint) {
        for (ai, &bi) in self.data.iter_mut().zip(other.data.iter()) {
            *ai &= bi;
        }
        self.data.truncate(other.data.len());
        self.normalize();
    }
}

forward_all_binop_to_val_ref_commutative!(impl BitOr for BigUint, bitor);
forward_val_assign!(impl BitOrAssign for BigUint, bitor_assign);

impl<'a> BitOr<&'a BigUint> for BigUint {
    type Output = BigUint;

    fn bitor(mut self, other: &BigUint) -> BigUint {
        self |= other;
        self
    }
}
impl<'a> BitOrAssign<&'a BigUint> for BigUint {
    #[inline]
    fn bitor_assign(&mut self, other: &BigUint) {
        for (ai, &bi) in self.data.iter_mut().zip(other.data.iter()) {
            *ai |= bi;
        }
        if other.data.len() > self.data.len() {
            let extra = &other.data[self.data.len()..];
            self.data.extend(extra.iter().cloned());
        }
    }
}

forward_all_binop_to_val_ref_commutative!(impl BitXor for BigUint, bitxor);
forward_val_assign!(impl BitXorAssign for BigUint, bitxor_assign);

impl<'a> BitXor<&'a BigUint> for BigUint {
    type Output = BigUint;

    fn bitxor(mut self, other: &BigUint) -> BigUint {
        self ^= other;
        self
    }
}
impl<'a> BitXorAssign<&'a BigUint> for BigUint {
    #[inline]
    fn bitxor_assign(&mut self, other: &BigUint) {
        for (ai, &bi) in self.data.iter_mut().zip(other.data.iter()) {
            *ai ^= bi;
        }
        if other.data.len() > self.data.len() {
            let extra = &other.data[self.data.len()..];
            self.data.extend(extra.iter().cloned());
        }
        self.normalize();
    }
}

promote_all_scalars!(impl BitAnd for BigUint, bitand);
promote_all_scalars_assign!(impl BitAndAssign for BigUint, bitand_assign);
forward_all_scalar_binop_to_val_val_commutative!(impl BitAnd<u32> for BigUint, bitand);
forward_all_scalar_binop_to_val_val_commutative!(impl BitAnd<u64> for BigUint, bitand);
forward_all_scalar_binop_to_val_val_commutative!(impl BitAnd<u128> for BigUint, bitand);

impl BitAnd<u32> for BigUint {
    type Output = BigUint;

    fn bitand(mut self, rhs: u32) -> Self::Output {
        self &= rhs;
        self
    }
}

impl BitAndAssign<u32> for BigUint {
    fn bitand_assign(&mut self, rhs: u32) {
        if !self.is_zero() {
            self.data[0] &= rhs as BigDigit;
            self.data.drain(1..);
        }
    }
}

impl BitAnd<u64> for BigUint {
    type Output = BigUint;

    fn bitand(mut self, rhs: u64) -> Self::Output {
        self &= rhs;
        self
    }
}

#[cfg(u64_digit)]
impl BitAndAssign<u64> for BigUint {
    fn bitand_assign(&mut self, rhs: u64) {
        if !self.is_zero() {
            self.data[0] &= rhs as BigDigit;
            self.data.drain(1..);
        }
    }
}

#[cfg(not(u64_digit))]
impl BitAndAssign<u64> for BigUint {
    fn bitand_assign(&mut self, rhs: u64) {
        if !self.is_zero() {
            self.data[0] &= rhs as BigDigit;
            if self.data.len() > 1 {
                self.data[1] &= (rhs >> big_digit::BITS) as BigDigit;
                self.data.drain(2..);
            }
        }
    }
}

impl BitAnd<u128> for BigUint {
    type Output = BigUint;

    fn bitand(mut self, rhs: u128) -> Self::Output {
        self &= rhs;
        self
    }
}

#[cfg(u64_digit)]
impl BitAndAssign<u128> for BigUint {
    fn bitand_assign(&mut self, rhs: u128) {
        if !self.is_zero() {
            self.data[0] &= rhs as BigDigit;
            if self.data.len() > 1 {
                self.data[1] &= (rhs >> big_digit::BITS) as BigDigit;
                self.data.drain(2..);
            }
        }
    }
}

#[cfg(not(u64_digit))]
impl BitAndAssign<u128> for BigUint {
    fn bitand_assign(&mut self, rhs: u128) {
        match self.data.len() {
            0 => {}
            1 => self.data[0] &= rhs as BigDigit,
            2 => {
                self.data[0] &= rhs as BigDigit;
                self.data[1] &= (rhs >> big_digit::BITS) as BigDigit;
            }
            3 => {
                self.data[0] &= rhs as BigDigit;
                self.data[1] &= (rhs >> big_digit::BITS) as BigDigit;
                self.data[2] &= (rhs >> (big_digit::BITS * 2)) as BigDigit;
            }
            _ => {
                self.data[0] &= rhs as BigDigit;
                self.data[1] &= (rhs >> big_digit::BITS) as BigDigit;
                self.data[2] &= (rhs >> (big_digit::BITS * 2)) as BigDigit;
                self.data[3] &= (rhs >> (big_digit::BITS * 3)) as BigDigit;
                self.data.drain(4..);
            }
        }
    }
}

// Implementation note: The signed bitwise variants work because `i* as u*`
// produces numbers with identical bit patterns, thus bitwise operations will
// return identical results. The only semantic difference lies in the leading
// zeroes, which (conceptually) are all 0 if positive, and all one when
// negative. As such, this will only truncate the leading bits when positive,
// since x & 0 == 0.

forward_all_scalar_binop_to_val_val_commutative!(impl BitAnd<i32> for BigUint, bitand);
forward_all_scalar_binop_to_val_val_commutative!(impl BitAnd<i64> for BigUint, bitand);
forward_all_scalar_binop_to_val_val_commutative!(impl BitAnd<i128> for BigUint, bitand);

impl BitAnd<i32> for BigUint {
    type Output = BigUint;

    fn bitand(mut self, rhs: i32) -> Self::Output {
        self &= rhs;
        self
    }
}

impl BitAndAssign<i32> for BigUint {
    fn bitand_assign(&mut self, rhs: i32) {
        if !self.is_zero() {
            self.data[0] &= rhs as BigDigit;
            if rhs >= 0 {
                self.data.drain(1..);
            }
        }
    }
}

impl BitAnd<i64> for BigUint {
    type Output = BigUint;

    fn bitand(mut self, rhs: i64) -> Self::Output {
        self &= rhs;
        self
    }
}

#[cfg(u64_digit)]
impl BitAndAssign<i64> for BigUint {
    fn bitand_assign(&mut self, rhs: i64) {
        if !self.is_zero() {
            self.data[0] &= rhs as BigDigit;
            if rhs >= 0 {
                self.data.drain(1..);
            }
        }
    }
}

#[cfg(not(u64_digit))]
impl BitAndAssign<i64> for BigUint {
    fn bitand_assign(&mut self, rhs: i64) {
        if !self.is_zero() {
            self.data[0] &= rhs as BigDigit;
            if self.data.len() > 1 {
                self.data[1] &= (rhs >> big_digit::BITS) as BigDigit;
                if rhs >= 0 {
                    self.data.drain(2..);
                }
            }
        }
    }
}

impl BitAnd<i128> for BigUint {
    type Output = BigUint;

    fn bitand(mut self, rhs: i128) -> Self::Output {
        self &= rhs;
        self
    }
}

#[cfg(u64_digit)]
impl BitAndAssign<i128> for BigUint {
    fn bitand_assign(&mut self, rhs: i128) {
        if !self.is_zero() {
            self.data[0] &= rhs as BigDigit;
            if self.data.len() > 1 {
                self.data[1] &= (rhs >> big_digit::BITS) as BigDigit;
                if rhs >= 0 {
                    self.data.drain(2..);
                }
            }
        }
    }
}

#[cfg(not(u64_digit))]
impl BitAndAssign<i128> for BigUint {
    fn bitand_assign(&mut self, rhs: i128) {
        match self.data.len() {
            0 => {}
            1 => self.data[0] &= rhs as BigDigit,
            2 => {
                self.data[0] &= rhs as BigDigit;
                self.data[1] &= (rhs >> big_digit::BITS) as BigDigit;
            }
            3 => {
                self.data[0] &= rhs as BigDigit;
                self.data[1] &= (rhs >> big_digit::BITS) as BigDigit;
                self.data[2] &= (rhs >> (big_digit::BITS * 2)) as BigDigit;
            }
            _ => {
                self.data[0] &= rhs as BigDigit;
                self.data[1] &= (rhs >> big_digit::BITS) as BigDigit;
                self.data[2] &= (rhs >> (big_digit::BITS * 2)) as BigDigit;
                self.data[3] &= (rhs >> (big_digit::BITS * 3)) as BigDigit;
                if rhs >= 0 {
                    self.data.drain(4..);
                }
            }
        }
    }
}
