//! Test randomization of `BigUint` and `BigInt`
//!
//! This test is in a completely separate crate so `rand::rng()`
//! can be available without "infecting" the rest of the build with
//! `rand`'s default features, especially not `rand/std`.

#![cfg(test)]

use num_bigint::BigRng010 as BigRng;

mod bigint;
mod biguint;
mod torture;
