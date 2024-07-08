use super::{BigUint, IntDigits};
#[cfg(target_arch = "x86_64")]
use std::arch::asm;

use crate::big_digit::{self, BigDigit};
use crate::UsizePromotion;

use core::iter::Sum;
use core::ops::{Add, AddAssign};
use num_traits::CheckedAdd;

#[cfg(target_arch = "x86_64")]
use core::arch::x86_64 as arch;

#[cfg(target_arch = "x86")]
use core::arch::x86 as arch;

// Add with carry:
#[cfg(target_arch = "x86_64")]
cfg_64!(
    #[inline]
    fn adc(carry: u8, a: u64, b: u64, out: &mut u64) -> u8 {
        // Safety: There are absolutely no safety concerns with calling `_addcarry_u64`.
        // It's just unsafe for API consistency with other intrinsics.
        unsafe { arch::_addcarry_u64(carry, a, b, out) }
    }
);

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
cfg_32!(
    #[inline]
    fn adc(carry: u8, a: u32, b: u32, out: &mut u32) -> u8 {
        // Safety: There are absolutely no safety concerns with calling `_addcarry_u32`.
        // It's just unsafe for API consistency with other intrinsics.
        unsafe { arch::_addcarry_u32(carry, a, b, out) }
    }
);

// fallback for environments where we don't have an addcarry intrinsic
// (copied from the standard library's `carrying_add`)
#[cfg(not(any(target_arch = "x86", target_arch = "x86_64")))]
#[inline]
fn adc(carry: u8, lhs: BigDigit, rhs: BigDigit, out: &mut BigDigit) -> u8 {
    let (a, b) = lhs.overflowing_add(rhs);
    let (c, d) = a.overflowing_add(carry as BigDigit);
    *out = c;
    u8::from(b || d)
}

/// Performs a part of the addition. Returns a tuple containing the carry state
/// and the number of integers that were added
///
/// By using as many registers as possible, we treat digits 5 by 5
#[cfg(target_arch = "x86_64")]
unsafe fn schoolbook_add_assign_x86_64(
    lhs: *mut u64,
    rhs: *const u64,
    mut size: usize,
) -> (bool, usize) {
    size /= 5;
    if size == 0 {
        return (false, 0);
    }

    let mut c: u8;
    let mut idx = 0;

    asm!(
        // Clear the carry flag
        "clc",

        "3:",

        // Copy a in registers
        "mov {a_tmp1}, qword ptr [{a} + 8*{idx}]",
        "mov {a_tmp2}, qword ptr [{a} + 8*{idx} + 8]",
        "mov {a_tmp3}, qword ptr [{a} + 8*{idx} + 16]",
        "mov {a_tmp4}, qword ptr [{a} + 8*{idx} + 24]",
        "mov {a_tmp5}, qword ptr [{a} + 8*{idx} + 32]",

        // Copy b in registers
        "mov {b_tmp1}, qword ptr [{b} + 8*{idx}]",
        "mov {b_tmp2}, qword ptr [{b} + 8*{idx} + 8]",
        "mov {b_tmp3}, qword ptr [{b} + 8*{idx} + 16]",
        "mov {b_tmp4}, qword ptr [{b} + 8*{idx} + 24]",
        "mov {b_tmp5}, qword ptr [{b} + 8*{idx} + 32]",

        // Perform the addition
        "adc {a_tmp1}, {b_tmp1}",
        "adc {a_tmp2}, {b_tmp2}",
        "adc {a_tmp3}, {b_tmp3}",
        "adc {a_tmp4}, {b_tmp4}",
        "adc {a_tmp5}, {b_tmp5}",

        // Copy the return values
        "mov qword ptr [{a} + 8*{idx}], {a_tmp1}",
        "mov qword ptr [{a} + 8*{idx} + 8], {a_tmp2}",
        "mov qword ptr [{a} + 8*{idx} + 16], {a_tmp3}",
        "mov qword ptr [{a} + 8*{idx} + 24], {a_tmp4}",
        "mov qword ptr [{a} + 8*{idx} + 32], {a_tmp5}",

        // Increment loop counter
        // `inc` and `dec` aren't modifying carry flag
        "inc {idx}",
        "inc {idx}",
        "inc {idx}",
        "inc {idx}",
        "inc {idx}",
        "dec {size}",
        "jnz 3b",

        // Output carry flag and clear
        "setc {c}",
        "clc",

        size = in(reg) size,
        a = in(reg) lhs,
        b = in(reg) rhs,
        c = lateout(reg_byte) c,
        idx = inout(reg) idx,

        a_tmp1 = out(reg) _,
        a_tmp2 = out(reg) _,
        a_tmp3 = out(reg) _,
        a_tmp4 = out(reg) _,
        a_tmp5 = out(reg) _,

        b_tmp1 = out(reg) _,
        b_tmp2 = out(reg) _,
        b_tmp3 = out(reg) _,
        b_tmp4 = out(reg) _,
        b_tmp5 = out(reg) _,

        options(nostack),
    );

    (c > 0, idx)
}

/// Two argument addition of raw slices, `a += b`, returning the carry.
///
/// This is used when the data `Vec` might need to resize to push a non-zero carry, so we perform
/// the addition first hoping that it will fit.
///
/// The caller _must_ ensure that `a` is at least as long as `b`.
#[inline]
pub(super) fn __add2(a: &mut [BigDigit], b: &[BigDigit]) -> BigDigit {
    debug_assert!(a.len() >= b.len());

    let (a_lo, a_hi) = a.split_at_mut(b.len());

    // On x86_64 machine, perform most of the addition via inline assembly
    #[cfg(target_arch = "x86_64")]
    let (c, done) = unsafe { schoolbook_add_assign_x86_64(a_lo.as_mut_ptr(), b.as_ptr(), b.len()) };
    #[cfg(not(target_arch = "x86_64"))]
    let (c, done) = (false, 0);

    let mut carry = c as u8;

    for (a, b) in a_lo[done..].iter_mut().zip(b[done..].iter()) {
        carry = adc(carry, *a, *b, a);
    }

    if carry != 0 {
        for a in a_hi {
            carry = adc(carry, *a, 0, a);
            if carry == 0 {
                break;
            }
        }
    }

    carry as BigDigit
}

/// Two argument addition of raw slices:
/// a += b
///
/// The caller _must_ ensure that a is big enough to store the result - typically this means
/// resizing a to max(a.len(), b.len()) + 1, to fit a possible carry.
pub(super) fn add2(a: &mut [BigDigit], b: &[BigDigit]) {
    let carry = __add2(a, b);

    debug_assert!(carry == 0);
}

forward_all_binop_to_val_ref_commutative!(impl Add for BigUint, add);
forward_val_assign!(impl AddAssign for BigUint, add_assign);

impl Add<&BigUint> for BigUint {
    type Output = BigUint;

    fn add(mut self, other: &BigUint) -> BigUint {
        self += other;
        self
    }
}
impl AddAssign<&BigUint> for BigUint {
    #[inline]
    fn add_assign(&mut self, other: &BigUint) {
        let self_len = self.data.len();
        let carry = if self_len < other.data.len() {
            let lo_carry = __add2(&mut self.data[..], &other.data[..self_len]);
            self.data.extend_from_slice(&other.data[self_len..]);
            __add2(&mut self.data[self_len..], &[lo_carry])
        } else {
            __add2(&mut self.data[..], &other.data[..])
        };
        if carry != 0 {
            self.data.push(carry);
        }
    }
}

promote_unsigned_scalars!(impl Add for BigUint, add);
promote_unsigned_scalars_assign!(impl AddAssign for BigUint, add_assign);
forward_all_scalar_binop_to_val_val_commutative!(impl Add<u32> for BigUint, add);
forward_all_scalar_binop_to_val_val_commutative!(impl Add<u64> for BigUint, add);
forward_all_scalar_binop_to_val_val_commutative!(impl Add<u128> for BigUint, add);

impl Add<u32> for BigUint {
    type Output = BigUint;

    #[inline]
    fn add(mut self, other: u32) -> BigUint {
        self += other;
        self
    }
}

impl AddAssign<u32> for BigUint {
    #[inline]
    fn add_assign(&mut self, other: u32) {
        if other != 0 {
            if self.data.is_empty() {
                self.data.push(0);
            }

            let carry = __add2(&mut self.data, &[other as BigDigit]);
            if carry != 0 {
                self.data.push(carry);
            }
        }
    }
}

impl Add<u64> for BigUint {
    type Output = BigUint;

    #[inline]
    fn add(mut self, other: u64) -> BigUint {
        self += other;
        self
    }
}

impl AddAssign<u64> for BigUint {
    cfg_digit!(
        #[inline]
        fn add_assign(&mut self, other: u64) {
            let (hi, lo) = big_digit::from_doublebigdigit(other);
            if hi == 0 {
                *self += lo;
            } else {
                while self.data.len() < 2 {
                    self.data.push(0);
                }

                let carry = __add2(&mut self.data, &[lo, hi]);
                if carry != 0 {
                    self.data.push(carry);
                }
            }
        }

        #[inline]
        fn add_assign(&mut self, other: u64) {
            if other != 0 {
                if self.data.is_empty() {
                    self.data.push(0);
                }

                let carry = __add2(&mut self.data, &[other as BigDigit]);
                if carry != 0 {
                    self.data.push(carry);
                }
            }
        }
    );
}

impl Add<u128> for BigUint {
    type Output = BigUint;

    #[inline]
    fn add(mut self, other: u128) -> BigUint {
        self += other;
        self
    }
}

impl AddAssign<u128> for BigUint {
    cfg_digit!(
        #[inline]
        fn add_assign(&mut self, other: u128) {
            if other <= u128::from(u64::MAX) {
                *self += other as u64
            } else {
                let (a, b, c, d) = super::u32_from_u128(other);
                let carry = if a > 0 {
                    while self.data.len() < 4 {
                        self.data.push(0);
                    }
                    __add2(&mut self.data, &[d, c, b, a])
                } else {
                    debug_assert!(b > 0);
                    while self.data.len() < 3 {
                        self.data.push(0);
                    }
                    __add2(&mut self.data, &[d, c, b])
                };

                if carry != 0 {
                    self.data.push(carry);
                }
            }
        }

        #[inline]
        fn add_assign(&mut self, other: u128) {
            let (hi, lo) = big_digit::from_doublebigdigit(other);
            if hi == 0 {
                *self += lo;
            } else {
                while self.data.len() < 2 {
                    self.data.push(0);
                }

                let carry = __add2(&mut self.data, &[lo, hi]);
                if carry != 0 {
                    self.data.push(carry);
                }
            }
        }
    );
}

impl CheckedAdd for BigUint {
    #[inline]
    fn checked_add(&self, v: &BigUint) -> Option<BigUint> {
        Some(self.add(v))
    }
}

impl_sum_iter_type!(BigUint);
