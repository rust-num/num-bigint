//! Randomization of big integers

#[cfg(any(feature = "rand_0_10", feature = "rand_0_9"))]
mod structs {
    use crate::{BigInt, BigUint};

    /// The back-end implementing rand's `UniformSampler` for [`BigUint`].
    #[derive(Clone, Debug)]
    pub struct UniformBigUint {
        pub(super) base: BigUint,
        pub(super) len: BigUint,
    }

    /// The back-end implementing rand's `UniformSampler` for [`BigInt`].
    #[derive(Clone, Debug)]
    pub struct UniformBigInt {
        pub(super) base: BigInt,
        pub(super) len: BigUint,
    }

    /// A random distribution for [`BigUint`] and [`BigInt`] values of a particular bit size.
    #[derive(Clone, Copy, Debug)]
    pub struct RandomBits {
        pub(super) bits: u64,
    }

    impl RandomBits {
        #[inline]
        pub fn new(bits: u64) -> RandomBits {
            RandomBits { bits }
        }
    }
}

#[cfg(any(feature = "rand_0_10", feature = "rand_0_9"))]
pub use structs::{RandomBits, UniformBigInt, UniformBigUint};

#[cfg(feature = "rand_core_0_10")]
#[path = "bigrand"]
mod impl_0_10 {
    use rand_core_0_10::Rng;

    mod traits;

    pub use traits::BigRng;

    #[cfg(feature = "rand_0_10")]
    use rand_0_10 as rand;

    #[cfg(feature = "rand_0_10")]
    #[cfg_attr(docsrs, doc(cfg(feature = "rand_0_10")))]
    mod distr;
}

#[cfg(feature = "rand_core_0_10")]
pub use impl_0_10::BigRng as BigRng010;

#[cfg(feature = "rand_core_0_9")]
#[path = "bigrand"]
mod impl_0_9 {
    use rand_core_0_9::RngCore as Rng;

    mod traits;

    pub use traits::BigRng;

    #[cfg(feature = "rand_0_9")]
    use rand_0_9 as rand;

    #[cfg(feature = "rand_0_9")]
    #[cfg_attr(docsrs, doc(cfg(feature = "rand_0_9")))]
    mod distr;
}

#[cfg(feature = "rand_core_0_9")]
pub use impl_0_9::BigRng as BigRng09;
