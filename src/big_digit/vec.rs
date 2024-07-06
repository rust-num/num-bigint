#![cfg(not(feature = "inline"))]

use super::BigDigit;
use alloc::vec::Vec;

// #[derive(Debug)]
#[repr(transparent)]
pub(crate) struct BigDigits {
    vec: Vec<BigDigit>,
}

impl BigDigits {
    pub(crate) const ZERO: Self = Self { vec: Vec::new() };

    #[inline]
    pub(crate) fn from_digit(x: BigDigit) -> Self {
        if x == 0 {
            Self::ZERO
        } else {
            Self { vec: vec![x] }
        }
    }

    #[inline]
    pub(crate) fn from_slice(slice: &[BigDigit]) -> Self {
        Self {
            vec: slice.to_vec(),
        }
    }

    #[inline]
    pub(crate) fn from_vec(vec: Vec<BigDigit>) -> Self {
        Self { vec }
    }

    #[inline]
    pub(crate) fn clear(&mut self) {
        self.vec.clear();
    }

    #[inline]
    pub(crate) fn push(&mut self, y: BigDigit) {
        self.vec.push(y);
    }

    #[inline]
    pub(crate) fn pop(&mut self) -> Option<BigDigit> {
        self.vec.pop()
    }

    #[inline]
    pub(crate) fn last(&self) -> Option<&BigDigit> {
        self.vec.last()
    }

    #[inline]
    pub(crate) fn len(&self) -> usize {
        self.vec.len()
    }

    #[inline]
    pub(crate) fn is_empty(&self) -> bool {
        self.vec.is_empty()
    }

    #[inline]
    pub(crate) fn capacity(&self) -> usize {
        self.vec.capacity()
    }

    #[inline]
    pub(crate) fn shrink(&mut self) {
        let xs = &mut self.vec;
        if xs.len() < xs.capacity() / 2 {
            xs.shrink_to(xs.len() + 1);
        }
    }

    /// Returns `true` if the most-significant digit (if any) is nonzero.
    #[inline]
    pub(crate) fn is_normal(&self) -> bool {
        !matches!(*self.vec, [.., 0])
    }

    /// Strips off trailing zero bigdigits - most algorithms require
    /// the most significant digit in the number to be nonzero.
    #[inline]
    pub(crate) fn normalize(&mut self) {
        let xs = &mut self.vec;
        if let [.., 0] = **xs {
            let len = xs.iter().rposition(|&d| d != 0).map_or(0, |i| i + 1);
            xs.truncate(len);
        }
        self.shrink()
    }

    #[inline]
    pub(crate) fn truncate(&mut self, len: usize) {
        self.vec.truncate(len);
    }

    #[inline]
    pub(crate) fn drain_front(&mut self, len: usize) {
        self.vec.drain(..len);
    }

    #[inline]
    pub(crate) fn resize(&mut self, len: usize, value: BigDigit) {
        self.vec.resize(len, value);
    }

    #[inline]
    pub(crate) fn extend_from_slice(&mut self, ys: &[BigDigit]) {
        self.vec.extend_from_slice(ys);
    }

    #[inline]
    pub(crate) fn extend<I>(&mut self, iter: I)
    where
        I: ExactSizeIterator<Item = BigDigit>,
    {
        self.vec.extend(iter);
    }
}

impl Clone for BigDigits {
    #[inline]
    fn clone(&self) -> Self {
        Self {
            vec: self.vec.clone(),
        }
    }

    #[inline]
    fn clone_from(&mut self, source: &Self) {
        self.vec.clone_from(&source.vec);
    }
}

impl core::ops::Deref for BigDigits {
    type Target = [BigDigit];

    #[inline]
    fn deref(&self) -> &Self::Target {
        &self.vec
    }
}

impl core::ops::DerefMut for BigDigits {
    #[inline]
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.vec
    }
}
