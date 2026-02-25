// Copyright 2013-2014 The Rust Project Developers. See the COPYRIGHT
// file at the top-level directory of this distribution and at
// http://rust-lang.org/COPYRIGHT.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

//! Big Integer Types for Rust
//!
//! * A [`BigUint`] is unsigned and represented as a vector of digits.
//! * A [`BigInt`] is signed and is a combination of [`BigUint`] and [`Sign`].
//!
//! Common numerical operations are overloaded, so we can treat them
//! the same way we treat other numbers.
//!
//! ## Example
//!
//! ```rust
//! # fn main() {
//! use num_bigint::BigUint;
//! use num_traits::One;
//!
//! // Calculate large fibonacci numbers.
//! fn fib(n: usize) -> BigUint {
//!     let mut f0 = BigUint::ZERO;
//!     let mut f1 = BigUint::one();
//!     for _ in 0..n {
//!         let f2 = f0 + &f1;
//!         f0 = f1;
//!         f1 = f2;
//!     }
//!     f0
//! }
//!
//! // This is a very large number.
//! println!("fib(1000) = {}", fib(1000));
//! # }
//! ```
//!
//! It's easy to generate large random numbers:
//!
#![cfg_attr(feature = "rand_0_10", doc = "```")]
#![cfg_attr(not(feature = "rand_0_10"), doc = "```ignore")]
//! # use ::rand_0_10 as rand;
//! use rand::{SeedableRng, rngs::SmallRng};
//! use num_bigint::{BigInt, BigRng010};
//!
//! // This seed is just for demonstration, but in most cases
//! // you'll probably want a non-deterministic `rng`.
//! let mut rng = SmallRng::seed_from_u64(42);
//! let a = rng.random_bigint(1000);
//!
//! let low = BigInt::from(-10000);
//! let high = BigInt::from(10000);
//! let b = rng.random_bigint_range(&low, &high);
//!
//! // Probably an even larger number.
//! println!("{}", a * b);
#![doc = "```"]
//!
//! See the "Features" section for instructions for enabling random number generation.
//!
//! ## Features
//!
//! The `std` crate feature is enabled by default, which enables [`std::error::Error`]
//! implementations and some internal use of floating point approximations. This can be disabled by
//! depending on `num-bigint` with `default-features = false`. Either way, the `alloc` crate is
//! always required for heap allocation of the `BigInt`/`BigUint` digits.
//!
//! ### Random Generation
//!
//! `num-bigint` supports the generation of random big integers when either of the `rand_0_9` or
//! `rand_0_10` features are enabled. The [`BigRng09`] and [`BigRng010`] traits provide extension
//! methods for any `rand_core` RNG of their respective version, while the structs [`RandomBits`],
//! [`UniformBigInt`], and [`UniformBigUint`] fulfill further functionality for random
//! distributions in `rand::distr`.
//!
//! For example, using `rand v0.10` in your `Cargo.toml` may look like this:
//!
//! ```toml
//! rand = "0.10"
//! num-bigint = { version = "0.5", features = ["rand_0_10"] }
//! ```
//!
//! Note that you must use the same version of `rand` as the feature you enable in `num-bigint`.
//! It's also fine for multiple versions to be enabled at once -- the random-distribution structs
//! will be shared with trait implementations for each `rand` feature that is enabled, while the
//! `BigRng` traits are distinct.
//!
//! You can instead use `rand_core_0_9` or `rand_core_0_10` for a more restricted subset, with
//! *only* the `BigRng` traits.
//!
//! ### Arbitrary Big Integers
//!
//! `num-bigint` supports `arbitrary` and `quickcheck` features to implement
//! [`arbitrary::Arbitrary`] and [`quickcheck::Arbitrary`], respectively, for both `BigInt` and
//! `BigUint`. These are useful for fuzzing and other forms of randomized testing.
//!
//! ### Serialization
//!
//! The `serde` feature adds implementations of [`Serialize`][serde::Serialize] and
//! [`Deserialize`][serde::Deserialize] for both `BigInt` and `BigUint`. Their serialized data is
//! generated portably, regardless of platform differences like the internal digit size.
//!
//!
//! ## Compatibility
//!
//! The `num-bigint` crate is tested for rustc 1.60 and greater.

#![cfg_attr(docsrs, feature(doc_cfg))]
#![warn(rust_2018_idioms)]
#![no_std]

#[macro_use]
extern crate alloc;

#[cfg(feature = "std")]
extern crate std;

use core::fmt;

#[macro_use]
mod macros;

mod bigint;
mod bigrand;
mod biguint;

#[cfg(target_pointer_width = "32")]
type UsizePromotion = u32;
#[cfg(target_pointer_width = "64")]
type UsizePromotion = u64;

#[cfg(target_pointer_width = "32")]
type IsizePromotion = i32;
#[cfg(target_pointer_width = "64")]
type IsizePromotion = i64;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ParseBigIntError {
    kind: BigIntErrorKind,
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum BigIntErrorKind {
    Empty,
    InvalidDigit,
}

impl ParseBigIntError {
    fn __description(&self) -> &str {
        use crate::BigIntErrorKind::*;
        match self.kind {
            Empty => "cannot parse integer from empty string",
            InvalidDigit => "invalid digit found in string",
        }
    }

    fn empty() -> Self {
        ParseBigIntError {
            kind: BigIntErrorKind::Empty,
        }
    }

    fn invalid() -> Self {
        ParseBigIntError {
            kind: BigIntErrorKind::InvalidDigit,
        }
    }
}

impl fmt::Display for ParseBigIntError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.__description().fmt(f)
    }
}

#[cfg(feature = "std")]
#[cfg_attr(docsrs, doc(cfg(feature = "std")))]
impl std::error::Error for ParseBigIntError {
    fn description(&self) -> &str {
        self.__description()
    }
}

/// The error type returned when a checked conversion regarding big integer fails.
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub struct TryFromBigIntError<T> {
    original: T,
}

impl<T> TryFromBigIntError<T> {
    fn new(original: T) -> Self {
        TryFromBigIntError { original }
    }

    fn __description(&self) -> &str {
        "out of range conversion regarding big integer attempted"
    }

    /// Extract the original value, if available. The value will be available
    /// if the type before conversion was either [`BigInt`] or [`BigUint`].
    pub fn into_original(self) -> T {
        self.original
    }
}

#[cfg(feature = "std")]
#[cfg_attr(docsrs, doc(cfg(feature = "std")))]
impl<T> std::error::Error for TryFromBigIntError<T>
where
    T: fmt::Debug,
{
    fn description(&self) -> &str {
        self.__description()
    }
}

impl<T> fmt::Display for TryFromBigIntError<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.__description().fmt(f)
    }
}

pub use crate::biguint::BigUint;
pub use crate::biguint::ToBigUint;
pub use crate::biguint::U32Digits;
pub use crate::biguint::U64Digits;

pub use crate::bigint::BigInt;
pub use crate::bigint::Sign;
pub use crate::bigint::ToBigInt;

#[cfg(feature = "rand_core_0_10")]
pub use crate::bigrand::BigRng010;

#[cfg(feature = "rand_core_0_9")]
pub use crate::bigrand::BigRng09;

#[cfg(any(feature = "rand_0_10", feature = "rand_0_9"))]
pub use crate::bigrand::{RandomBits, UniformBigInt, UniformBigUint};

mod big_digit {
    // A [`BigDigit`] is a [`BigUint`]'s composing element.
    cfg_digit!(
        pub(crate) type BigDigit = u32;
        pub(crate) type BigDigit = u64;
    );

    // A [`DoubleBigDigit`] is the internal type used to do the computations.  Its
    // size is the double of the size of [`BigDigit`].
    cfg_digit!(
        pub(crate) type DoubleBigDigit = u64;
        pub(crate) type DoubleBigDigit = u128;
    );

    pub(crate) const BITS: u8 = BigDigit::BITS as u8;
    pub(crate) const HALF_BITS: u8 = BITS / 2;
    pub(crate) const HALF: BigDigit = (1 << HALF_BITS) - 1;

    pub(crate) const MAX: BigDigit = BigDigit::MAX;
    const LO_MASK: DoubleBigDigit = MAX as DoubleBigDigit;

    #[inline]
    fn get_hi(n: DoubleBigDigit) -> BigDigit {
        (n >> BITS) as BigDigit
    }
    #[inline]
    fn get_lo(n: DoubleBigDigit) -> BigDigit {
        (n & LO_MASK) as BigDigit
    }

    /// Split one [`DoubleBigDigit`] into two [`BigDigit`]s.
    #[inline]
    pub(crate) fn from_doublebigdigit(n: DoubleBigDigit) -> (BigDigit, BigDigit) {
        (get_hi(n), get_lo(n))
    }

    /// Join two [`BigDigit`]s into one [`DoubleBigDigit`].
    #[inline]
    pub(crate) fn to_doublebigdigit(hi: BigDigit, lo: BigDigit) -> DoubleBigDigit {
        DoubleBigDigit::from(lo) | (DoubleBigDigit::from(hi) << BITS)
    }
}
