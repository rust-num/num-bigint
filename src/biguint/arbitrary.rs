#![cfg(any(feature = "quickcheck", feature = "arbitrary", feature = "proptest"))]

use super::{biguint_from_vec, BigUint};

use crate::big_digit::BigDigit;
#[cfg(feature = "quickcheck")]
use alloc::boxed::Box;
use alloc::vec::Vec;

#[cfg(feature = "quickcheck")]
#[cfg_attr(docsrs, doc(cfg(feature = "quickcheck")))]
impl quickcheck::Arbitrary for BigUint {
    fn arbitrary(g: &mut quickcheck::Gen) -> Self {
        // Use arbitrary from Vec
        biguint_from_vec(Vec::<BigDigit>::arbitrary(g))
    }

    fn shrink(&self) -> Box<dyn Iterator<Item = Self>> {
        // Use shrinker from Vec
        Box::new(self.data.to_vec().shrink().map(biguint_from_vec))
    }
}

#[cfg(feature = "arbitrary")]
#[cfg_attr(docsrs, doc(cfg(feature = "arbitrary")))]
impl arbitrary::Arbitrary<'_> for BigUint {
    fn arbitrary(u: &mut arbitrary::Unstructured<'_>) -> arbitrary::Result<Self> {
        Ok(biguint_from_vec(Vec::<BigDigit>::arbitrary(u)?))
    }

    fn arbitrary_take_rest(u: arbitrary::Unstructured<'_>) -> arbitrary::Result<Self> {
        Ok(biguint_from_vec(Vec::<BigDigit>::arbitrary_take_rest(u)?))
    }

    fn size_hint(depth: usize) -> (usize, Option<usize>) {
        Vec::<BigDigit>::size_hint(depth)
    }
}

#[cfg(feature = "proptest")]
#[cfg_attr(docsrs, doc(cfg(feature = "proptest")))]
impl proptest::arbitrary::Arbitrary for BigUint {
    type Parameters = ();
    type Strategy = proptest::strategy::BoxedStrategy<Self>;

    fn arbitrary_with(_args: Self::Parameters) -> Self::Strategy {
        use proptest::prelude::*;

        let bigint_strategy = any::<Vec<BigDigit>>().prop_map(biguint_from_vec);
        bigint_strategy.boxed()
    }
}

// These are just smoke tests for the arbitrary implementations to ensure they can be called.
#[cfg(test)]
mod tests {
    use super::*;

    #[cfg(feature = "quickcheck")]
    #[test]
    fn test_quickcheck_arbitrary() {
        let mut gen = quickcheck::Gen::new(10);
        let bigint = <BigUint as quickcheck::Arbitrary>::arbitrary(&mut gen);
        assert_eq!(bigint, bigint);
    }

    #[cfg(feature = "arbitrary")]
    #[test]
    fn test_arbitrary_arbitrary() {
        let mut u = arbitrary::Unstructured::new(&[1, 2, 3]);
        let bigint = <BigUint as arbitrary::Arbitrary>::arbitrary(&mut u);
        assert_eq!(bigint, bigint);
    }

    #[cfg(feature = "proptest")]
    proptest::proptest! {
        #[test]
        fn test_proptest_arbitrary(bigint: BigUint) {
            proptest::prop_assert_eq!(&bigint, &bigint);
        }
    }
}
