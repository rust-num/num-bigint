use core::iter::FusedIterator;

#[cfg(not(u64_digit))]
use super::u32_chunk_to_u64;

/// An iterator of `u32` digits representation of a `BigUint` or `BigInt`,
/// ordered most significant digit first.
pub struct U32DigitsBe<'a> {
    #[cfg(u64_digit)]
    data: &'a [u64],
    #[cfg(u64_digit)]
    next_is_hi: bool,

    #[cfg(not(u64_digit))]
    it: core::slice::Iter<'a, u32>,
}

#[cfg(u64_digit)]
impl<'a> U32DigitsBe<'a> {
    #[inline]
    pub(super) fn new(data: &'a [u64]) -> Self {
        let last_hi_is_zero = data
            .last()
            .map(|&last| {
                let last_hi = (last >> 32) as u32;
                last_hi == 0
            })
            .unwrap_or(false);
        U32DigitsBe {
            data,
            next_is_hi: !last_hi_is_zero,
        }
    }
}

#[cfg(u64_digit)]
impl Iterator for U32DigitsBe<'_> {
    type Item = u32;
    #[inline]
    fn next(&mut self) -> Option<u32> {
        match self.data.split_last() {
            Some((&first, data)) => {
                let next_is_hi = self.next_is_hi;
                self.next_is_hi = !next_is_hi;
                if next_is_hi {
                    Some((first >> 32) as u32)
                } else {
                    self.data = data;
                    Some(first as u32)
                }
            }
            None => None,
        }
    }

    #[inline]
    fn size_hint(&self) -> (usize, Option<usize>) {
        let len = self.len();
        (len, Some(len))
    }

    #[inline]
    fn last(self) -> Option<u32> {
        self.data.first().map(|&last| last as u32)
    }

    #[inline]
    fn count(self) -> usize {
        self.len()
    }
}

#[cfg(u64_digit)]
impl ExactSizeIterator for U32DigitsBe<'_> {
    #[inline]
    fn len(&self) -> usize {
        self.data.len() * 2 - usize::from(!self.next_is_hi)
    }
}

#[cfg(not(u64_digit))]
impl<'a> U32DigitsBe<'a> {
    #[inline]
    pub(super) fn new(data: &'a [u32]) -> Self {
        Self { it: data.iter() }
    }
}

#[cfg(not(u64_digit))]
impl Iterator for U32DigitsBe<'_> {
    type Item = u32;
    #[inline]
    fn next(&mut self) -> Option<u32> {
        self.it.next_back().cloned()
    }

    #[inline]
    fn size_hint(&self) -> (usize, Option<usize>) {
        self.it.size_hint()
    }

    #[inline]
    fn nth(&mut self, n: usize) -> Option<u32> {
        self.it.nth_back(n).cloned()
    }

    #[inline]
    fn last(mut self) -> Option<u32> {
        self.it.next().cloned()
    }

    #[inline]
    fn count(self) -> usize {
        self.it.count()
    }
}

#[cfg(not(u64_digit))]
impl ExactSizeIterator for U32DigitsBe<'_> {
    #[inline]
    fn len(&self) -> usize {
        self.it.len()
    }
}

impl FusedIterator for U32DigitsBe<'_> {}

/// An iterator of `u64` digits representation of a `BigUint` or `BigInt`,
/// ordered most significant digit first.
pub struct U64DigitsBe<'a> {
    #[cfg(not(u64_digit))]
    it: core::slice::Chunks<'a, u32>,

    #[cfg(u64_digit)]
    it: core::slice::Iter<'a, u64>,
}

#[cfg(not(u64_digit))]
impl<'a> U64DigitsBe<'a> {
    #[inline]
    pub(super) fn new(data: &'a [u32]) -> Self {
        U64DigitsBe { it: data.chunks(2) }
    }
}

#[cfg(not(u64_digit))]
impl Iterator for U64DigitsBe<'_> {
    type Item = u64;
    #[inline]
    fn next(&mut self) -> Option<u64> {
        self.it.next_back().map(u32_chunk_to_u64)
    }

    #[inline]
    fn size_hint(&self) -> (usize, Option<usize>) {
        let len = self.len();
        (len, Some(len))
    }

    #[inline]
    fn last(mut self) -> Option<u64> {
        self.it.next().map(u32_chunk_to_u64)
    }

    #[inline]
    fn count(self) -> usize {
        self.len()
    }
}

#[cfg(u64_digit)]
impl<'a> U64DigitsBe<'a> {
    #[inline]
    pub(super) fn new(data: &'a [u64]) -> Self {
        Self { it: data.iter() }
    }
}

#[cfg(u64_digit)]
impl Iterator for U64DigitsBe<'_> {
    type Item = u64;
    #[inline]
    fn next(&mut self) -> Option<u64> {
        self.it.next_back().cloned()
    }

    #[inline]
    fn size_hint(&self) -> (usize, Option<usize>) {
        self.it.size_hint()
    }

    #[inline]
    fn nth(&mut self, n: usize) -> Option<u64> {
        self.it.nth_back(n).cloned()
    }

    #[inline]
    fn last(mut self) -> Option<u64> {
        self.it.next().cloned()
    }

    #[inline]
    fn count(self) -> usize {
        self.it.count()
    }
}

impl ExactSizeIterator for U64DigitsBe<'_> {
    #[inline]
    fn len(&self) -> usize {
        self.it.len()
    }
}

impl FusedIterator for U64DigitsBe<'_> {}

#[test]
fn test_iter_u32_digits() {
    let n = super::BigUint::from(5u8);
    let mut it = n.iter_u32_digits_be();
    assert_eq!(it.len(), 1);
    assert_eq!(it.next(), Some(5));
    assert_eq!(it.len(), 0);
    assert_eq!(it.next(), None);
    assert_eq!(it.len(), 0);
    assert_eq!(it.next(), None);

    let n = super::BigUint::from(112500000000u64);
    let mut it = n.iter_u32_digits_be();
    assert_eq!(it.len(), 2);
    assert_eq!(it.next(), Some(26));
    assert_eq!(it.len(), 1);
    assert_eq!(it.next(), Some(830850304));
    assert_eq!(it.len(), 0);
    assert_eq!(it.next(), None);
}

#[test]
fn test_iter_u64_digits() {
    let n = super::BigUint::from(5u8);
    let mut it = n.iter_u64_digits_be();
    assert_eq!(it.len(), 1);
    assert_eq!(it.next(), Some(5));
    assert_eq!(it.len(), 0);
    assert_eq!(it.next(), None);
    assert_eq!(it.len(), 0);
    assert_eq!(it.next(), None);

    let n = super::BigUint::from(18_446_744_073_709_551_616u128);
    let mut it = n.iter_u64_digits_be();
    assert_eq!(it.len(), 2);
    assert_eq!(it.next(), Some(1));
    assert_eq!(it.len(), 1);
    assert_eq!(it.next(), Some(0));
    assert_eq!(it.len(), 0);
    assert_eq!(it.next(), None);
}
