//! Test randomization of `BigUint` and `BigInt`
//!
//! This test is in a completely separate crate so `rand::rng()`
//! can be available without "infecting" the rest of the build with
//! `rand`'s default features, especially not `rand/std`.

#![cfg(test)]

use num_bigint::BigRng09 as BigRng;

#[path = "../../big_rand_0_10/src/bigint.rs"]
mod bigint;
#[path = "../../big_rand_0_10/src/biguint.rs"]
mod biguint;
#[path = "../../big_rand_0_10/src/torture.rs"]
mod torture;
