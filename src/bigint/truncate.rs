use super::{BigInt, Sign};

macro_rules! impl_truncate_bigint {
    ($T:ty, $truncate_N:ident) => {
        /// Returns the input reduced modulo 2 to the power of the size of the target type.
        /// This is analogous to the behavior of `as` conversion on integral types.
        #[inline]
        pub fn $truncate_N(&self) -> $T {
            match self.sign {
                Sign::Minus => self.data.$truncate_N().wrapping_neg(),
                Sign::NoSign => 0,
                Sign::Plus => self.data.$truncate_N(),
            }
        }
    };
}

impl BigInt {
    impl_truncate_bigint!(u128, truncate_u128);
    impl_truncate_bigint!(usize, truncate_usize);
    impl_truncate_bigint!(u64, truncate_u64);
    impl_truncate_bigint!(u32, truncate_u32);
    impl_truncate_bigint!(u16, truncate_u16);
    impl_truncate_bigint!(u8, truncate_u8);
    impl_truncate_bigint!(i128, truncate_i128);
    impl_truncate_bigint!(isize, truncate_isize);
    impl_truncate_bigint!(i64, truncate_i64);
    impl_truncate_bigint!(i32, truncate_i32);
    impl_truncate_bigint!(i16, truncate_i16);
    impl_truncate_bigint!(i8, truncate_i8);
}
