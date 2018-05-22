//! Randomization of big integers

// Some of the tests of non-RNG-based functionality are randomized using the
// RNG-based functionality, so the RNG-based functionality needs to be enabled
// for tests.
use rand::Rng;

use BigInt;
use BigUint;
use Sign::*;
use big_digit::BigDigit;

use traits::Zero;
use integer::Integer;

pub trait RandBigInt {
    /// Generate a random `BigUint` of the given bit size.
    fn gen_biguint(&mut self, bit_size: usize) -> BigUint;

    /// Generate a random BigInt of the given bit size.
    fn gen_bigint(&mut self, bit_size: usize) -> BigInt;

    /// Generate a random `BigUint` less than the given bound. Fails
    /// when the bound is zero.
    fn gen_biguint_below(&mut self, bound: &BigUint) -> BigUint;

    /// Generate a random `BigUint` within the given range. The lower
    /// bound is inclusive; the upper bound is exclusive. Fails when
    /// the upper bound is not greater than the lower bound.
    fn gen_biguint_range(&mut self, lbound: &BigUint, ubound: &BigUint) -> BigUint;

    /// Generate a random `BigInt` within the given range. The lower
    /// bound is inclusive; the upper bound is exclusive. Fails when
    /// the upper bound is not greater than the lower bound.
    fn gen_bigint_range(&mut self, lbound: &BigInt, ubound: &BigInt) -> BigInt;
}

#[cfg(any(feature = "rand", test))]
impl<R: Rng> RandBigInt for R {
    fn gen_biguint(&mut self, bit_size: usize) -> BigUint {
        use super::big_digit::BITS;
        let (digits, rem) = bit_size.div_rem(&BITS);
        let mut data = Vec::with_capacity(digits + 1);
        for _ in 0..digits {
            data.push(self.gen());
        }
        if rem > 0 {
            let final_digit: BigDigit = self.gen();
            data.push(final_digit >> (BITS - rem));
        }
        BigUint::new(data)
    }

    fn gen_bigint(&mut self, bit_size: usize) -> BigInt {
        loop {
            // Generate a random BigUint...
            let biguint = self.gen_biguint(bit_size);
            // ...and then randomly assign it a Sign...
            let sign = if biguint.is_zero() {
                // ...except that if the BigUint is zero, we need to try
                // again with probability 0.5. This is because otherwise,
                // the probability of generating a zero BigInt would be
                // double that of any other number.
                if self.gen() {
                    continue;
                } else {
                    NoSign
                }
            } else if self.gen() {
                Plus
            } else {
                Minus
            };
            return BigInt::from_biguint(sign, biguint);
        }
    }

    fn gen_biguint_below(&mut self, bound: &BigUint) -> BigUint {
        assert!(!bound.is_zero());
        let bits = bound.bits();
        loop {
            let n = self.gen_biguint(bits);
            if n < *bound {
                return n;
            }
        }
    }

    fn gen_biguint_range(&mut self, lbound: &BigUint, ubound: &BigUint) -> BigUint {
        assert!(*lbound < *ubound);
        return lbound + self.gen_biguint_below(&(ubound - lbound));
    }

    fn gen_bigint_range(&mut self, lbound: &BigInt, ubound: &BigInt) -> BigInt {
        assert!(*lbound < *ubound);
        let delta = (ubound - lbound).to_biguint().unwrap();
        return lbound + BigInt::from(self.gen_biguint_below(&delta))
    }
}
