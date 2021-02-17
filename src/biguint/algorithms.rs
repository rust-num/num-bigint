use core::cmp::Ordering::{self, Equal};
use core::mem;
use num_traits::PrimInt;

use crate::big_digit::BigDigit;

/// Find last set bit
/// fls(0) == 0, fls(u32::MAX) == 32
pub(crate) fn fls<T: PrimInt>(v: T) -> u8 {
    mem::size_of::<T>() as u8 * 8 - v.leading_zeros() as u8
}

pub(crate) fn ilog2<T: PrimInt>(v: T) -> u8 {
    fls(v) - 1
}

pub(crate) fn cmp_slice(a: &[BigDigit], b: &[BigDigit]) -> Ordering {
    debug_assert!(a.last() != Some(&0));
    debug_assert!(b.last() != Some(&0));

    match Ord::cmp(&a.len(), &b.len()) {
        Equal => Iterator::cmp(a.iter().rev(), b.iter().rev()),
        other => other,
    }
}
