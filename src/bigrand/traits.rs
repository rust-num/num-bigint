//! Randomization of big integers

use super::Rng;
use crate::Sign::*;
use crate::{BigInt, BigUint};

use crate::big_digit::BigDigit;
use crate::biguint::biguint_from_vec;

use core::mem::size_of;
use num_integer::Integer;
use num_traits::Zero;

/// A trait extending [`Rng`] for sampling random big integers.
pub trait BigRng: Rng {
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

impl<R: Rng + ?Sized> BigRng for R {
    fn random_biguint(&mut self, bit_size: u64) -> BigUint {
        let mask = (1 << (bit_size % 8)) - 1;
        let bytes = Integer::div_ceil(&bit_size, &8);
        let bytes = usize::try_from(bytes).expect("capacity overflow");

        let digits = Integer::div_ceil(&bit_size, &BigDigit::BITS.into());
        let digits = usize::try_from(digits).unwrap();
        debug_assert_eq!(digits, Integer::div_ceil(&bytes, &size_of::<BigDigit>()));

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
