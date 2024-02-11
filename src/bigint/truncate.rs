use num_traits::{WrappingNeg, Zero};

use super::{BigUint, BigInt, Sign};
use crate::biguint::TruncateFrom;

impl BigInt {
    /// Returns the input reduced modulo 2 to the power of the size of the target type.
    /// This is analogous to the behavior of `as` conversion on integral types.
    #[inline]
    pub fn truncate<T: TruncateFrom<BigUint> + Zero + WrappingNeg>(&self) -> T {
        match self.sign {
            Sign::Minus => self.data.truncate::<T>().wrapping_neg(),
            Sign::NoSign => T::zero(),
            Sign::Plus => self.data.truncate(),
        }
    }
}
