use core::iter::Step;

use num_traits::CheckedSub;

use crate::{BigInt, BigUint};

impl Step for BigUint {
    fn steps_between(start: &Self, end: &Self) -> Option<usize> {
        end.checked_sub(start)?.try_into().ok()
    }

    fn forward_checked(start: Self, count: usize) -> Option<Self> {
        Some(Self::forward(start, count))
    }

    fn backward_checked(start: Self, count: usize) -> Option<Self> {
        start.checked_sub(&count.into())
    }

    fn forward(start: Self, count: usize) -> Self {
        start + count
    }
}

impl Step for BigInt {
    fn steps_between(start: &Self, end: &Self) -> Option<usize> {
        (end - start).try_into().ok()
    }

    fn forward_checked(start: Self, count: usize) -> Option<Self> {
        Some(Self::forward(start, count))
    }

    fn backward_checked(start: Self, count: usize) -> Option<Self> {
        Some(Self::backward(start, count))
    }

    fn forward(start: Self, count: usize) -> Self {
        start + count
    }

    fn backward(start: Self, count: usize) -> Self {
        start - count
    }
}

#[cfg(test)]
mod tests {
    use itertools::assert_equal;
    use num_traits::Zero;

    use crate::{BigInt, BigUint};

    #[test]
    fn ranges() {
        assert_equal(
            BigUint::zero()..BigUint::from(6u8),
            (0u8..6).map(BigUint::from),
        );
        assert_equal(BigInt::from(-6)..BigInt::from(6), (-6..6).map(BigInt::from));
    }
}
