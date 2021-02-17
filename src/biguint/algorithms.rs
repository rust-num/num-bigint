use core::cmp::Ordering::{self, Equal};

use crate::big_digit::BigDigit;

pub(crate) fn cmp_slice(a: &[BigDigit], b: &[BigDigit]) -> Ordering {
    debug_assert!(a.last() != Some(&0));
    debug_assert!(b.last() != Some(&0));

    match Ord::cmp(&a.len(), &b.len()) {
        Equal => Iterator::cmp(a.iter().rev(), b.iter().rev()),
        other => other,
    }
}
