//! Example generating prime numbers in a few large bit sizes.

#![cfg_attr(not(feature = "rand"), allow(unused))]

extern crate num_bigint;
extern crate num_integer;
extern crate num_traits;
extern crate rand;

use num_bigint::{BigUint, RandBigInt};
use num_integer::Integer;
use num_traits::{One, ToPrimitive, Zero};
use rand::{Rng, SeedableRng, XorShiftRng};
use std::env;
use std::u64;

#[cfg(not(feature = "rand"))]
fn main() {
    println!("This example requires `--features rand`");
}

#[cfg(feature = "rand")]
fn main() {
    // Use a seeded Rng for benchmarking purposes.
    let seeded = env::var_os("SEEDED_PRIMES").is_some();
    let (mut thread_rng, mut seeded_rng);
    let mut rng: &mut Rng = if seeded {
        seeded_rng = XorShiftRng::from_seed([4, 3, 2, 1]);
        &mut seeded_rng
    } else {
        thread_rng = rand::thread_rng();
        &mut thread_rng
    };

    for &bits in &[256, 512, 1024, 2048] {
        let min = BigUint::one() << (bits - 1);
        for i in 1usize.. {
            let mut n = &min + rng.gen_biguint_below(&min);
            if n.is_even() {
                n += 1u32;
            }
            if is_prime(&n) {
                println!("Found a {}-bit prime in {} attempts:", bits, i);
                println!("{}", n);
                println!("");
                break;
            }
        }
    }
}

/// Attempt to determine whether n is prime.  This algorithm is exact for
/// numbers less than 2^64, and highly probabilistic for other numbers.
#[cfg(feature = "rand")]
fn is_prime(n: &BigUint) -> bool {
    // guarantees of primality given by Sloane's [A014233](https://oeis.org/A014233)
    const A014233: [(u8, u64); 12] = [
        (2, 2_047),
        (3, 1_373_653),
        (5, 25_326_001),
        (7, 3_215_031_751),
        (11, 2_152_302_898_747),
        (13, 3_474_749_660_383),
        (17, 341_550_071_728_321),
        (19, 341_550_071_728_321),
        (23, 3_825_123_056_546_413_051),
        (29, 3_825_123_056_546_413_051),
        (31, 3_825_123_056_546_413_051),
        (37, u64::MAX),
    ];

    let n64 = n.to_u64();
    if n.is_even() {
        return n64 == Some(2);
    } else if n64 == Some(1) {
        return false;
    }

    // try some simple divisors
    for &(p, _) in &A014233[1..] {
        if n64 == Some(u64::from(p)) {
            return true;
        }
        if (n % p).is_zero() {
            return false;
        }
    }

    // try series of primes as witnesses
    for &(p, u) in &A014233 {
        if witness(&BigUint::from(p), n) {
            return false;
        }
        if let Some(n) = n64 {
            if n < u || u == u64::MAX {
                return true;
            }
        }
    }

    // Now we're into the "probably prime" arena...
    // Generate a few witnesses in the range `[2, n - 2]`
    let mut rng = rand::thread_rng();
    let max = n - 3u32;
    for _ in 0..10 {
        let w = rng.gen_biguint_below(&max) + 2u32;
        if witness(&w, n) {
            return false;
        }
    }
    true
}

/// Test a possible witness to the compositeness of n.
///
/// Computes a**(n-1) (mod n), and checks if the result gives evidence that
/// n is composite.  Returning false means that n is likely prime, but not
/// necessarily.
#[cfg(feature = "rand")]
fn witness(a: &BigUint, n: &BigUint) -> bool {
    // let n-1 = (2**t)*u, where t ≥ 1 and u is odd
    // TODO `trailing_zeros` would make this trivial
    let n_m1 = n - 1u32;
    let (t, u) = if n_m1.is_even() {
        let mut t = 1usize;
        let mut u = n_m1.clone() >> 1;
        while u.is_even() {
            t += 1;
            u >>= 1;
        }
        (t, u)
    } else {
        (0, n_m1.clone())
    };

    let mut x = a.modpow(&u, n);
    if x.is_one() || x == n_m1 {
        return false;
    }

    for _ in 1..t {
        x = &x * &x % n;
        if x.is_one() {
            return true;
        }
        if x == n_m1 {
            return false;
        }
    }

    true
}
