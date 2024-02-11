use super::{BigUint, IntDigits};

use crate::big_digit::{self, BigDigit};
use crate::UsizePromotion;

#[cfg(not(u64_digit))]
use crate::std_alloc::Vec;

use core::ops::{BitAnd, BitAndAssign, BitOr, BitOrAssign, BitXor, BitXorAssign};
use num_traits::Zero;

forward_val_val_binop!(impl BitAnd for BigUint, bitand);
forward_ref_val_binop!(impl BitAnd for BigUint, bitand);

// do not use forward_ref_ref_binop_commutative! for bitand so that we can
// clone the smaller value rather than the larger, avoiding over-allocation
impl BitAnd<&BigUint> for &BigUint {
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

impl BitAnd<&BigUint> for BigUint {
    type Output = BigUint;

    #[inline]
    fn bitand(mut self, other: &BigUint) -> BigUint {
        self &= other;
        self
    }
}
impl BitAndAssign<&BigUint> for BigUint {
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

impl BitOr<&BigUint> for BigUint {
    type Output = BigUint;

    fn bitor(mut self, other: &BigUint) -> BigUint {
        self |= other;
        self
    }
}
impl BitOrAssign<&BigUint> for BigUint {
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

impl BitXor<&BigUint> for BigUint {
    type Output = BigUint;

    fn bitxor(mut self, other: &BigUint) -> BigUint {
        self ^= other;
        self
    }
}
impl BitXorAssign<&BigUint> for BigUint {
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

promote_unsigned_scalars!(impl BitAnd for BigUint, bitand);
promote_unsigned_scalars_assign!(impl BitAndAssign for BigUint, bitand_assign);
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

// Implementation note: Bitwise or (and xor) are not implemented for signed
// types because there is no reasonable value for the result to be if rhs is
// negative.

promote_unsigned_scalars!(impl BitOr for BigUint, bitor);
promote_unsigned_scalars_assign!(impl BitOrAssign for BigUint, bitor_assign);
forward_all_scalar_binop_to_val_val_commutative!(impl BitOr<u32> for BigUint, bitor);
forward_all_scalar_binop_to_val_val_commutative!(impl BitOr<u64> for BigUint, bitor);
forward_all_scalar_binop_to_val_val_commutative!(impl BitOr<u128> for BigUint, bitor);

impl BitOr<u32> for BigUint {
    type Output = BigUint;

    fn bitor(mut self, rhs: u32) -> Self::Output {
        self |= rhs;
        self
    }
}

impl BitOrAssign<u32> for BigUint {
    fn bitor_assign(&mut self, rhs: u32) {
        if !self.is_zero() {
            self.data[0] |= rhs as BigDigit;
        } else {
            *self = rhs.into();
        }
    }
}

impl BitOr<u64> for BigUint {
    type Output = BigUint;

    fn bitor(mut self, rhs: u64) -> Self::Output {
        self |= rhs;
        self
    }
}

#[cfg(u64_digit)]
impl BitOrAssign<u64> for BigUint {
    fn bitor_assign(&mut self, rhs: u64) {
        if !self.is_zero() {
            self.data[0] |= rhs;
        } else {
            self.data.push(rhs);
        }
    }
}

#[cfg(not(u64_digit))]
impl BitOrAssign<u64> for BigUint {
    fn bitor_assign(&mut self, rhs: u64) {
        match self.data.len() {
            0 => *self = rhs.into(),
            1 => {
                self.data[0] |= rhs as BigDigit;
                if rhs > big_digit::MAX {
                    self.data.push((rhs >> big_digit::BITS) as BigDigit);
                }
            }
            _ => {
                self.data[0] |= rhs as BigDigit;
                self.data[1] |= (rhs >> big_digit::BITS) as u32;
            }
        }
    }
}

impl BitOr<u128> for BigUint {
    type Output = BigUint;

    fn bitor(mut self, rhs: u128) -> Self::Output {
        self |= rhs;
        self
    }
}

#[cfg(u64_digit)]
impl BitOrAssign<u128> for BigUint {
    fn bitor_assign(&mut self, rhs: u128) {
        if !self.is_zero() {
            self.data[0] |= rhs as BigDigit;
            if self.data.len() > 1 {
                self.data[1] |= (rhs >> big_digit::BITS) as BigDigit;
            } else if rhs > big_digit::MAX as u128 {
                self.data.push((rhs >> big_digit::BITS) as BigDigit);
            }
        } else {
            *self = rhs.into();
        }
    }
}

#[inline]
#[cfg(not(u64_digit))]
fn push_nonzero<T: Zero + Copy>(data: &mut Vec<T>, to_add: &[T]) {
    for i in to_add {
        if i.is_zero() {
            return;
        } else {
            data.push(*i);
        }
    }
}

#[cfg(not(u64_digit))]
impl BitOrAssign<u128> for BigUint {
    fn bitor_assign(&mut self, rhs: u128) {
        let a = rhs as BigDigit;
        let b = (rhs >> big_digit::BITS) as BigDigit;
        let c = (rhs >> (big_digit::BITS * 2)) as BigDigit;
        let d = (rhs >> (big_digit::BITS * 2)) as BigDigit;
        match self.data.len() {
            0 => *self = rhs.into(),
            1 => {
                self.data[0] &= a;
                push_nonzero(&mut self.data, &[b, c, d]);
            }
            2 => {
                self.data[0] &= a;
                self.data[1] &= b;
                push_nonzero(&mut self.data, &[c, d]);
            }
            3 => {
                self.data[0] &= a;
                self.data[1] &= b;
                self.data[2] &= c;
                push_nonzero(&mut self.data, &[d]);
            }
            _ => {
                self.data[0] &= a;
                self.data[1] &= b;
                self.data[2] &= c;
                self.data[3] &= d;
            }
        }
    }
}

promote_unsigned_scalars!(impl BitXor for BigUint, bitxor);
promote_unsigned_scalars_assign!(impl BitXorAssign for BigUint, bitxor_assign);
forward_all_scalar_binop_to_val_val_commutative!(impl BitXor<u32> for BigUint, bitxor);
forward_all_scalar_binop_to_val_val_commutative!(impl BitXor<u64> for BigUint, bitxor);
forward_all_scalar_binop_to_val_val_commutative!(impl BitXor<u128> for BigUint, bitxor);

impl BitXor<u32> for BigUint {
    type Output = BigUint;

    fn bitxor(mut self, rhs: u32) -> Self::Output {
        self ^= rhs;
        self
    }
}

impl BitXorAssign<u32> for BigUint {
    fn bitxor_assign(&mut self, rhs: u32) {
        if !self.is_zero() {
            self.data[0] ^= rhs as BigDigit;
        } else {
            *self = rhs.into();
        }
    }
}

impl BitXor<u64> for BigUint {
    type Output = BigUint;

    fn bitxor(mut self, rhs: u64) -> Self::Output {
        self ^= rhs;
        self
    }
}

#[cfg(u64_digit)]
impl BitXorAssign<u64> for BigUint {
    fn bitxor_assign(&mut self, rhs: u64) {
        if !self.is_zero() {
            self.data[0] ^= rhs;
        } else {
            self.data.push(rhs);
        }
    }
}

#[cfg(not(u64_digit))]
impl BitXorAssign<u64> for BigUint {
    fn bitxor_assign(&mut self, rhs: u64) {
        match self.data.len() {
            0 => *self = rhs.into(),
            1 => {
                self.data[0] ^= rhs as BigDigit;
                if rhs > big_digit::MAX {
                    self.data.push((rhs >> big_digit::BITS) as BigDigit);
                }
            }
            _ => {
                self.data[0] ^= rhs as BigDigit;
                self.data[1] ^= (rhs >> big_digit::BITS) as u32;
            }
        }
    }
}

impl BitXor<u128> for BigUint {
    type Output = BigUint;

    fn bitxor(mut self, rhs: u128) -> Self::Output {
        self ^= rhs;
        self
    }
}

#[cfg(u64_digit)]
impl BitXorAssign<u128> for BigUint {
    fn bitxor_assign(&mut self, rhs: u128) {
        if !self.is_zero() {
            self.data[0] ^= rhs as BigDigit;
            if self.data.len() > 1 {
                self.data[1] ^= (rhs >> big_digit::BITS) as BigDigit;
            } else if rhs > big_digit::MAX as u128 {
                self.data.push((rhs >> big_digit::BITS) as BigDigit);
            }
        } else {
            *self = rhs.into();
        }
    }
}

#[cfg(not(u64_digit))]
impl BitXorAssign<u128> for BigUint {
    fn bitxor_assign(&mut self, rhs: u128) {
        let a = rhs as BigDigit;
        let b = (rhs >> big_digit::BITS) as BigDigit;
        let c = (rhs >> (big_digit::BITS * 2)) as BigDigit;
        let d = (rhs >> (big_digit::BITS * 2)) as BigDigit;
        match self.data.len() {
            0 => *self = rhs.into(),
            1 => {
                self.data[0] ^= a;
                push_nonzero(&mut self.data, &[b, c, d]);
            }
            2 => {
                self.data[0] ^= a;
                self.data[1] ^= b;
                push_nonzero(&mut self.data, &[c, d]);
            }
            3 => {
                self.data[0] ^= a;
                self.data[1] ^= b;
                self.data[2] ^= c;
                push_nonzero(&mut self.data, &[d]);
            }
            _ => {
                self.data[0] ^= a;
                self.data[1] ^= b;
                self.data[2] ^= c;
                self.data[3] ^= d;
            }
        }
    }
}
