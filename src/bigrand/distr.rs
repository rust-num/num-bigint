//! Random distributions of big integers

use super::{rand, BigRng};
use crate::{BigInt, BigUint, RandomBits, UniformBigInt, UniformBigUint};

use rand::distr::uniform::{Error, SampleBorrow, SampleUniform, UniformSampler};
use rand::distr::Distribution;
use rand::Rng;

/// Registers [`UniformBigUint`] as the back-end enabling
/// [`Uniform`][rand::distr::Uniform] for [`BigUint`].
impl SampleUniform for BigUint {
    type Sampler = UniformBigUint;
}

/// Enables [`Uniform`][rand::distr::Uniform] for [`BigUint`].
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

/// Registers [`UniformBigInt`] as the back-end enabling
/// [`Uniform`][rand::distr::Uniform] for [`BigInt`].
impl SampleUniform for BigInt {
    type Sampler = UniformBigInt;
}

/// Enables [`Uniform`][rand::distr::Uniform] for [`BigInt`].
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
