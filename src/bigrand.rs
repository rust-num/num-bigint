//! Randomization of big integers
#![cfg(feature = "rand_core")]

use crate::BigInt;
use crate::BigUint;
use crate::Sign::*;

use crate::big_digit::BigDigit;
use crate::biguint::biguint_from_vec;

use num_traits::Zero;

#[cfg(feature = "rand")]
pub use self::distr::{RandomBits, UniformBigInt, UniformBigUint};

/// A trait for sampling random big integers.
///
/// The `rand` feature must be enabled to use this. See crate-level documentation for details.
pub trait RandBigInt {
    /// Generate a random [`BigUint`] of the given bit size.
    fn random_biguint(&mut self, bit_size: u64) -> BigUint;

    /// Generate a random [ BigInt`] of the given bit size.
    fn random_bigint(&mut self, bit_size: u64) -> BigInt;

    /// Generate a random [`BigUint`] less than the given bound. Fails
    /// when the bound is zero.
    fn random_biguint_below(&mut self, bound: &BigUint) -> BigUint;

    /// Generate a random [`BigUint`] within the given range. The lower
    /// bound is inclusive; the upper bound is exclusive. Fails when
    /// the upper bound is not greater than the lower bound.
    fn random_biguint_range(&mut self, lbound: &BigUint, ubound: &BigUint) -> BigUint;

    /// Generate a random [`BigInt`] within the given range. The lower
    /// bound is inclusive; the upper bound is exclusive. Fails when
    /// the upper bound is not greater than the lower bound.
    fn random_bigint_range(&mut self, lbound: &BigInt, ubound: &BigInt) -> BigInt;

    #[deprecated(since = "0.5.0", note = "Renamed to `random_biguint`")]
    fn gen_biguint(&mut self, bit_size: u64) -> BigUint {
        self.random_biguint(bit_size)
    }

    #[deprecated(since = "0.5.0", note = "Renamed to `random_bigint`")]
    fn gen_bigint(&mut self, bit_size: u64) -> BigInt {
        self.random_bigint(bit_size)
    }

    #[deprecated(since = "0.5.0", note = "Renamed to `random_biguint_below`")]
    fn gen_biguint_below(&mut self, bound: &BigUint) -> BigUint {
        self.random_biguint_below(bound)
    }

    #[deprecated(since = "0.5.0", note = "Renamed to `random_biguint_range`")]
    fn gen_biguint_range(&mut self, lbound: &BigUint, ubound: &BigUint) -> BigUint {
        self.random_biguint_range(lbound, ubound)
    }

    #[deprecated(since = "0.5.0", note = "Renamed to `random_bigint_range`")]
    fn gen_bigint_range(&mut self, lbound: &BigInt, ubound: &BigInt) -> BigInt {
        self.random_bigint_range(lbound, ubound)
    }
}

impl<R: rand_core::Rng + ?Sized> RandBigInt for R {
    fn random_biguint(&mut self, bit_size: u64) -> BigUint {
        let bytes = usize::try_from(bit_size.div_ceil(8)).expect("capacity overflow");
        let mask = (1 << (bit_size % 8)) - 1;
        let digits = usize::try_from(bit_size.div_ceil(BigDigit::BITS.into())).unwrap();
        debug_assert_eq!(digits, bytes.div_ceil(size_of::<BigDigit>()));
        let mut data = vec![0 as BigDigit; digits];
        {
            let ptr = data.as_mut_ptr().cast::<u8>();
            let slice = unsafe { core::slice::from_raw_parts_mut(ptr, bytes) };
            self.fill_bytes(slice);
            if mask != 0 {
                slice[bytes - 1] &= mask;
            }
        }
        #[cfg(not(target_endian = "little"))]
        for digit in &mut data {
            *digit = BigDigit::from_le(*digit);
        }
        biguint_from_vec(data)
    }

    fn random_bigint(&mut self, bit_size: u64) -> BigInt {
        loop {
            // Generate a random BigUint...
            let biguint = self.random_biguint(bit_size);
            // ...and then randomly assign it a Sign...
            let negative = (self.next_u32() as i32) < 0;
            let sign = if biguint.is_zero() {
                // ...except that if the BigUint is zero, we need to try
                // again with probability 0.5. This is because otherwise,
                // the probability of generating a zero BigInt would be
                // double that of any other number.
                if negative {
                    continue;
                }
                NoSign
            } else if negative {
                Minus
            } else {
                Plus
            };
            return BigInt::from_biguint(sign, biguint);
        }
    }

    fn random_biguint_below(&mut self, bound: &BigUint) -> BigUint {
        assert!(!bound.is_zero());
        let bits = bound.bits();
        if bound.trailing_zeros() == Some(bits - 1) {
            // The exclusive bound is trivially met in this case
            return self.random_biguint(bits - 1);
        }
        loop {
            let n = self.random_biguint(bits);
            if n < *bound {
                return n;
            }
        }
    }

    fn random_biguint_range(&mut self, lbound: &BigUint, ubound: &BigUint) -> BigUint {
        assert!(*lbound < *ubound);
        if lbound.is_zero() {
            self.random_biguint_below(ubound)
        } else {
            lbound + self.random_biguint_below(&(ubound - lbound))
        }
    }

    fn random_bigint_range(&mut self, lbound: &BigInt, ubound: &BigInt) -> BigInt {
        assert!(*lbound < *ubound);
        if lbound.is_zero() {
            BigInt::from(self.random_biguint_below(ubound.magnitude()))
        } else if ubound.is_zero() {
            lbound + BigInt::from(self.random_biguint_below(lbound.magnitude()))
        } else {
            let delta = ubound - lbound;
            lbound + BigInt::from(self.random_biguint_below(delta.magnitude()))
        }
    }
}

#[cfg(feature = "rand")]
mod distr {
    use super::{BigInt, BigUint, RandBigInt};

    use rand::distr::uniform::{Error, SampleBorrow, SampleUniform, UniformSampler};
    use rand::distr::Distribution;
    use rand_core::Rng;

    /// The back-end implementing rand's [`UniformSampler`] for [`BigUint`].
    #[derive(Clone, Debug)]
    pub struct UniformBigUint {
        base: BigUint,
        len: BigUint,
    }

    impl UniformSampler for UniformBigUint {
        type X = BigUint;

        #[inline]
        fn new<B1, B2>(low_b: B1, high_b: B2) -> Result<Self, Error>
        where
            B1: SampleBorrow<Self::X> + Sized,
            B2: SampleBorrow<Self::X> + Sized,
        {
            let low = low_b.borrow();
            let high = high_b.borrow();
            if low < high {
                Ok(UniformBigUint {
                    len: high - low,
                    base: low.clone(),
                })
            } else {
                Err(Error::EmptyRange)
            }
        }

        #[inline]
        fn new_inclusive<B1, B2>(low_b: B1, high_b: B2) -> Result<Self, Error>
        where
            B1: SampleBorrow<Self::X> + Sized,
            B2: SampleBorrow<Self::X> + Sized,
        {
            let low = low_b.borrow();
            let high = high_b.borrow();
            if low <= high {
                Ok(UniformBigUint {
                    len: high - low + 1u32,
                    base: low.clone(),
                })
            } else {
                Err(Error::EmptyRange)
            }
        }

        #[inline]
        fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> Self::X {
            &self.base + rng.random_biguint_below(&self.len)
        }

        #[inline]
        fn sample_single<R: Rng + ?Sized, B1, B2>(
            low: B1,
            high: B2,
            rng: &mut R,
        ) -> Result<Self::X, Error>
        where
            B1: SampleBorrow<Self::X> + Sized,
            B2: SampleBorrow<Self::X> + Sized,
        {
            let low = low.borrow();
            let high = high.borrow();
            if low < high {
                Ok(rng.random_biguint_range(low, high))
            } else {
                Err(Error::EmptyRange)
            }
        }

        #[inline]
        fn sample_single_inclusive<R: Rng + ?Sized, B1, B2>(
            low: B1,
            high: B2,
            rng: &mut R,
        ) -> Result<Self::X, Error>
        where
            B1: SampleBorrow<Self::X> + Sized,
            B2: SampleBorrow<Self::X> + Sized,
        {
            let low = low.borrow();
            let high = high.borrow();
            if low <= high {
                Ok(rng.random_biguint_range(low, &(high + 1u32)))
            } else {
                Err(Error::EmptyRange)
            }
        }
    }

    impl SampleUniform for BigUint {
        type Sampler = UniformBigUint;
    }

    /// The back-end implementing rand's [`UniformSampler`] for [`BigInt`].
    #[derive(Clone, Debug)]
    pub struct UniformBigInt {
        base: BigInt,
        len: BigUint,
    }

    impl UniformSampler for UniformBigInt {
        type X = BigInt;

        #[inline]
        fn new<B1, B2>(low_b: B1, high_b: B2) -> Result<Self, Error>
        where
            B1: SampleBorrow<Self::X> + Sized,
            B2: SampleBorrow<Self::X> + Sized,
        {
            let low = low_b.borrow();
            let high = high_b.borrow();
            if low < high {
                Ok(UniformBigInt {
                    len: (high - low).into_parts().1,
                    base: low.clone(),
                })
            } else {
                Err(Error::EmptyRange)
            }
        }

        #[inline]
        fn new_inclusive<B1, B2>(low_b: B1, high_b: B2) -> Result<Self, Error>
        where
            B1: SampleBorrow<Self::X> + Sized,
            B2: SampleBorrow<Self::X> + Sized,
        {
            let low = low_b.borrow();
            let high = high_b.borrow();
            if low <= high {
                Ok(UniformBigInt {
                    len: (high - low).into_parts().1 + 1u32,
                    base: low.clone(),
                })
            } else {
                Err(Error::EmptyRange)
            }
        }

        #[inline]
        fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> Self::X {
            &self.base + BigInt::from(rng.random_biguint_below(&self.len))
        }

        #[inline]
        fn sample_single<R: Rng + ?Sized, B1, B2>(
            low: B1,
            high: B2,
            rng: &mut R,
        ) -> Result<Self::X, Error>
        where
            B1: SampleBorrow<Self::X> + Sized,
            B2: SampleBorrow<Self::X> + Sized,
        {
            let low = low.borrow();
            let high = high.borrow();
            if low < high {
                Ok(rng.random_bigint_range(low, high))
            } else {
                Err(Error::EmptyRange)
            }
        }

        #[inline]
        fn sample_single_inclusive<R: Rng + ?Sized, B1, B2>(
            low: B1,
            high: B2,
            rng: &mut R,
        ) -> Result<Self::X, Error>
        where
            B1: SampleBorrow<Self::X> + Sized,
            B2: SampleBorrow<Self::X> + Sized,
        {
            let low = low.borrow();
            let high = high.borrow();
            if low <= high {
                Ok(rng.random_bigint_range(low, &(high + 1u32)))
            } else {
                Err(Error::EmptyRange)
            }
        }
    }

    impl SampleUniform for BigInt {
        type Sampler = UniformBigInt;
    }

    /// A random distribution for [`BigUint`] and [`BigInt`] values of a particular bit size.
    ///
    /// The `rand` feature must be enabled to use this. See crate-level documentation for details.
    #[derive(Clone, Copy, Debug)]
    pub struct RandomBits {
        bits: u64,
    }

    impl RandomBits {
        #[inline]
        pub fn new(bits: u64) -> RandomBits {
            RandomBits { bits }
        }
    }

    impl Distribution<BigUint> for RandomBits {
        #[inline]
        fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> BigUint {
            rng.random_biguint(self.bits)
        }
    }

    impl Distribution<BigInt> for RandomBits {
        #[inline]
        fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> BigInt {
            rng.random_bigint(self.bits)
        }
    }
}
