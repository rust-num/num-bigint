use crate::std_alloc::{Cow, Vec};
use core::cmp::Ordering::{self, Equal, Greater, Less};
use core::iter::repeat;
use core::mem;
use num_traits::{One, PrimInt, Zero};

use crate::biguint::biguint_from_vec;
use crate::biguint::BigUint;

use crate::big_digit::{self, BigDigit, DoubleBigDigit};

use super::addition::__add2;

/// Divide a two digit numerator by a one digit divisor, returns quotient and remainder:
///
/// Note: the caller must ensure that both the quotient and remainder will fit into a single digit.
/// This is _not_ true for an arbitrary numerator/denominator.
///
/// (This function also matches what the x86 divide instruction does).
#[inline]
fn div_wide(hi: BigDigit, lo: BigDigit, divisor: BigDigit) -> (BigDigit, BigDigit) {
    debug_assert!(hi < divisor);

    let lhs = big_digit::to_doublebigdigit(hi, lo);
    let rhs = DoubleBigDigit::from(divisor);
    ((lhs / rhs) as BigDigit, (lhs % rhs) as BigDigit)
}

/// For small divisors, we can divide without promoting to `DoubleBigDigit` by
/// using half-size pieces of digit, like long-division.
#[inline]
fn div_half(rem: BigDigit, digit: BigDigit, divisor: BigDigit) -> (BigDigit, BigDigit) {
    use crate::big_digit::{HALF, HALF_BITS};
    use num_integer::Integer;

    debug_assert!(rem < divisor && divisor <= HALF);
    let (hi, rem) = ((rem << HALF_BITS) | (digit >> HALF_BITS)).div_rem(&divisor);
    let (lo, rem) = ((rem << HALF_BITS) | (digit & HALF)).div_rem(&divisor);
    ((hi << HALF_BITS) | lo, rem)
}

#[inline]
pub(crate) fn div_rem_digit(mut a: BigUint, b: BigDigit) -> (BigUint, BigDigit) {
    let mut rem = 0;

    if b <= big_digit::HALF {
        for d in a.data.iter_mut().rev() {
            let (q, r) = div_half(rem, *d, b);
            *d = q;
            rem = r;
        }
    } else {
        for d in a.data.iter_mut().rev() {
            let (q, r) = div_wide(rem, *d, b);
            *d = q;
            rem = r;
        }
    }

    (a.normalized(), rem)
}

#[inline]
pub(crate) fn rem_digit(a: &BigUint, b: BigDigit) -> BigDigit {
    let mut rem = 0;

    if b <= big_digit::HALF {
        for &digit in a.data.iter().rev() {
            let (_, r) = div_half(rem, digit, b);
            rem = r;
        }
    } else {
        for &digit in a.data.iter().rev() {
            let (_, r) = div_wide(rem, digit, b);
            rem = r;
        }
    }

    rem
}

/// Subtract a multiple.
/// a -= b * c
/// Returns a borrow (if a < b then borrow > 0).
fn sub_mul_digit_same_len(a: &mut [BigDigit], b: &[BigDigit], c: BigDigit) -> BigDigit {
    debug_assert!(a.len() == b.len());

    // carry is between -big_digit::MAX and 0, so to avoid overflow we store
    // offset_carry = carry + big_digit::MAX
    let mut offset_carry = big_digit::MAX;

    for (x, y) in a.iter_mut().zip(b) {
        // We want to calculate sum = x - y * c + carry.
        // sum >= -(big_digit::MAX * big_digit::MAX) - big_digit::MAX
        // sum <= big_digit::MAX
        // Offsetting sum by (big_digit::MAX << big_digit::BITS) puts it in DoubleBigDigit range.
        let offset_sum = big_digit::to_doublebigdigit(big_digit::MAX, *x)
            - big_digit::MAX as DoubleBigDigit
            + offset_carry as DoubleBigDigit
            - *y as DoubleBigDigit * c as DoubleBigDigit;

        let (new_offset_carry, new_x) = big_digit::from_doublebigdigit(offset_sum);
        offset_carry = new_offset_carry;
        *x = new_x;
    }

    // Return the borrow.
    big_digit::MAX - offset_carry
}

pub(crate) fn div_rem(mut u: BigUint, mut d: BigUint) -> (BigUint, BigUint) {
    if d.is_zero() {
        panic!("attempt to divide by zero")
    }
    if u.is_zero() {
        return (Zero::zero(), Zero::zero());
    }

    if d.data.len() == 1 {
        if d.data == [1] {
            return (u, Zero::zero());
        }
        let (div, rem) = div_rem_digit(u, d.data[0]);
        // reuse d
        d.data.clear();
        d += rem;
        return (div, d);
    }

    // Required or the q_len calculation below can underflow:
    match u.cmp(&d) {
        Less => return (Zero::zero(), u),
        Equal => {
            u.set_one();
            return (u, Zero::zero());
        }
        Greater => {} // Do nothing
    }

    // This algorithm is from Knuth, TAOCP vol 2 section 4.3, algorithm D:
    //
    // First, normalize the arguments so the highest bit in the highest digit of the divisor is
    // set: the main loop uses the highest digit of the divisor for generating guesses, so we
    // want it to be the largest number we can efficiently divide by.
    //
    let shift = d.data.last().unwrap().leading_zeros() as usize;

    let (q, r) = if shift == 0 {
        // no need to clone d
        div_rem_core(u, &d)
    } else {
        div_rem_core(u << shift, &(d << shift))
    };
    // renormalize the remainder
    (q, r >> shift)
}

pub(crate) fn div_rem_ref(u: &BigUint, d: &BigUint) -> (BigUint, BigUint) {
    if d.is_zero() {
        panic!("attempt to divide by zero")
    }
    if u.is_zero() {
        return (Zero::zero(), Zero::zero());
    }

    if d.data.len() == 1 {
        if d.data == [1] {
            return (u.clone(), Zero::zero());
        }

        let (div, rem) = div_rem_digit(u.clone(), d.data[0]);
        return (div, rem.into());
    }

    // Required or the q_len calculation below can underflow:
    match u.cmp(d) {
        Less => return (Zero::zero(), u.clone()),
        Equal => return (One::one(), Zero::zero()),
        Greater => {} // Do nothing
    }

    // This algorithm is from Knuth, TAOCP vol 2 section 4.3, algorithm D:
    //
    // First, normalize the arguments so the highest bit in the highest digit of the divisor is
    // set: the main loop uses the highest digit of the divisor for generating guesses, so we
    // want it to be the largest number we can efficiently divide by.
    //
    let shift = d.data.last().unwrap().leading_zeros() as usize;

    let (q, r) = if shift == 0 {
        // no need to clone d
        div_rem_core(u.clone(), d)
    } else {
        div_rem_core(u << shift, &(d << shift))
    };
    // renormalize the remainder
    (q, r >> shift)
}

/// An implementation of the base division algorithm.
/// Knuth, TAOCP vol 2 section 4.3.1, algorithm D, with an improvement from exercises 19-21.
fn div_rem_core(mut a: BigUint, b: &BigUint) -> (BigUint, BigUint) {
    debug_assert!(
        a.data.len() >= b.data.len()
            && b.data.len() > 1
            && b.data.last().unwrap().leading_zeros() == 0
    );

    // The algorithm works by incrementally calculating "guesses", q0, for the next digit of the
    // quotient. Once we have any number q0 such that (q0 << j) * b <= a, we can set
    //
    //     q += q0 << j
    //     a -= (q0 << j) * b
    //
    // and then iterate until a < b. Then, (q, a) will be our desired quotient and remainder.
    //
    // q0, our guess, is calculated by dividing the last three digits of a by the last two digits of
    // b - this will give us a guess that is close to the actual quotient, but is possibly greater.
    // It can only be greater by 1 and only in rare cases, with probability at most
    // 2^-(big_digit::BITS-1) for random a, see TAOCP 4.3.1 exercise 21.
    //
    // If the quotient turns out to be too large, we adjust it by 1:
    // q -= 1 << j
    // a += b << j

    // a0 stores an additional extra most significant digit of the dividend, not stored in a.
    let mut a0 = 0;

    // [b1, b0] are the two most significant digits of the divisor. They never change.
    let b0 = *b.data.last().unwrap();
    let b1 = b.data[b.data.len() - 2];

    let q_len = a.data.len() - b.data.len() + 1;
    let mut q = BigUint {
        data: vec![0; q_len],
    };

    for j in (0..q_len).rev() {
        debug_assert!(a.data.len() == b.data.len() + j);

        let a1 = *a.data.last().unwrap();
        let a2 = a.data[a.data.len() - 2];

        // The first q0 estimate is [a1,a0] / b0. It will never be too small, it may be too large
        // by at most 2.
        let (mut q0, mut r) = if a0 < b0 {
            let (q0, r) = div_wide(a0, a1, b0);
            (q0, r as DoubleBigDigit)
        } else {
            debug_assert!(a0 == b0);
            // Avoid overflowing q0, we know the quotient fits in BigDigit.
            // [a1,a0] = b0 * (1<<BITS - 1) + (a0 + a1)
            (big_digit::MAX, a0 as DoubleBigDigit + a1 as DoubleBigDigit)
        };

        // r = [a1,a0] - q0 * b0
        //
        // Now we want to compute a more precise estimate [a2,a1,a0] / [b1,b0] which can only be
        // less or equal to the current q0.
        //
        // q0 is too large if:
        // [a2,a1,a0] < q0 * [b1,b0]
        // (r << BITS) + a2 < q0 * b1
        while r <= big_digit::MAX as DoubleBigDigit
            && big_digit::to_doublebigdigit(r as BigDigit, a2)
                < q0 as DoubleBigDigit * b1 as DoubleBigDigit
        {
            q0 -= 1;
            r += b0 as DoubleBigDigit;
        }

        // q0 is now either the correct quotient digit, or in rare cases 1 too large.
        // Subtract (q0 << j) from a. This may overflow, in which case we will have to correct.

        let mut borrow = sub_mul_digit_same_len(&mut a.data[j..], &b.data, q0);
        if borrow > a0 {
            // q0 is too large. We need to add back one multiple of b.
            q0 -= 1;
            borrow -= __add2(&mut a.data[j..], &b.data);
        }
        // The top digit of a, stored in a0, has now been zeroed.
        debug_assert!(borrow == a0);

        q.data[j] = q0;

        // Pop off the next top digit of a.
        a0 = a.data.pop().unwrap();
    }

    a.data.push(a0);
    a.normalize();

    debug_assert!(a < *b);

    (q.normalized(), a)
}

/// Find last set bit
/// fls(0) == 0, fls(u32::MAX) == 32
pub(crate) fn fls<T: PrimInt>(v: T) -> u8 {
    mem::size_of::<T>() as u8 * 8 - v.leading_zeros() as u8
}

pub(crate) fn ilog2<T: PrimInt>(v: T) -> u8 {
    fls(v) - 1
}

#[inline]
pub(crate) fn biguint_shl<T: PrimInt>(n: Cow<'_, BigUint>, shift: T) -> BigUint {
    if shift < T::zero() {
        panic!("attempt to shift left with negative");
    }
    if n.is_zero() {
        return n.into_owned();
    }
    let bits = T::from(big_digit::BITS).unwrap();
    let digits = (shift / bits).to_usize().expect("capacity overflow");
    let shift = (shift % bits).to_u8().unwrap();
    biguint_shl2(n, digits, shift)
}

fn biguint_shl2(n: Cow<'_, BigUint>, digits: usize, shift: u8) -> BigUint {
    let mut data = match digits {
        0 => n.into_owned().data,
        _ => {
            let len = digits.saturating_add(n.data.len() + 1);
            let mut data = Vec::with_capacity(len);
            data.extend(repeat(0).take(digits));
            data.extend(n.data.iter());
            data
        }
    };

    if shift > 0 {
        let mut carry = 0;
        let carry_shift = big_digit::BITS as u8 - shift;
        for elem in data[digits..].iter_mut() {
            let new_carry = *elem >> carry_shift;
            *elem = (*elem << shift) | carry;
            carry = new_carry;
        }
        if carry != 0 {
            data.push(carry);
        }
    }

    biguint_from_vec(data)
}

#[inline]
pub(crate) fn biguint_shr<T: PrimInt>(n: Cow<'_, BigUint>, shift: T) -> BigUint {
    if shift < T::zero() {
        panic!("attempt to shift right with negative");
    }
    if n.is_zero() {
        return n.into_owned();
    }
    let bits = T::from(big_digit::BITS).unwrap();
    let digits = (shift / bits).to_usize().unwrap_or(core::usize::MAX);
    let shift = (shift % bits).to_u8().unwrap();
    biguint_shr2(n, digits, shift)
}

fn biguint_shr2(n: Cow<'_, BigUint>, digits: usize, shift: u8) -> BigUint {
    if digits >= n.data.len() {
        let mut n = n.into_owned();
        n.set_zero();
        return n;
    }
    let mut data = match n {
        Cow::Borrowed(n) => n.data[digits..].to_vec(),
        Cow::Owned(mut n) => {
            n.data.drain(..digits);
            n.data
        }
    };

    if shift > 0 {
        let mut borrow = 0;
        let borrow_shift = big_digit::BITS as u8 - shift;
        for elem in data.iter_mut().rev() {
            let new_borrow = *elem << borrow_shift;
            *elem = (*elem >> shift) | borrow;
            borrow = new_borrow;
        }
    }

    biguint_from_vec(data)
}

pub(crate) fn cmp_slice(a: &[BigDigit], b: &[BigDigit]) -> Ordering {
    debug_assert!(a.last() != Some(&0));
    debug_assert!(b.last() != Some(&0));

    match Ord::cmp(&a.len(), &b.len()) {
        Equal => Iterator::cmp(a.iter().rev(), b.iter().rev()),
        other => other,
    }
}
