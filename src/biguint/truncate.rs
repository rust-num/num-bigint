use super::{big_digit, BigUint};

macro_rules! impl_truncate_large {
    ($T:ty, $truncate_N:ident, $size:expr) => {
        /// Returns the input reduced modulo 2 to the power of the size of the target type.
        /// This is analogous to the behavior of `as` conversion on integral types.
        #[inline]
        pub fn $truncate_N(&self) -> $T {
            let mut ret = 0;
            let mut bits = 0;

            for i in self.data.iter() {
                if bits >= $size {
                    break;
                }

                ret += <$T>::from(*i) << bits;
                bits += big_digit::BITS;
            }
            ret
        }
    };
}

macro_rules! impl_truncate_from {
    ($truncate_from:ident: $($T:ty, $truncate_to:ident);* $(;)*) => {
        $(
            #[inline]
            pub fn $truncate_to(&self) -> $T {
                self.$truncate_from() as $T
            }
        )*
    };
}

impl BigUint {
    #[cfg(not(u64_digit))]
    impl_truncate_large!(u32, truncate_u32, 32);
    impl_truncate_large!(u64, truncate_u64, 64);
    impl_truncate_large!(u128, truncate_u128, 128);

    #[cfg(not(u64_digit))]
    impl_truncate_from! {
        truncate_u32:
        u16, truncate_u16;
        u8, truncate_u8;
    }

    #[cfg(u64_digit)]
    impl_truncate_from! {
        truncate_u64:
        u32, truncate_u32;
        u16, truncate_u16;
        u8, truncate_u8;
    }

    #[cfg(target_pointer_width = "64")]
    impl_truncate_from!(truncate_u64: usize, truncate_usize);
    #[cfg(target_pointer_width = "32")]
    impl_truncate_from!(truncate_u32: usize, truncate_usize);
    #[cfg(target_pointer_width = "16")]
    impl_truncate_from!(truncate_u16: usize, truncate_usize);

    impl_truncate_from!(truncate_u128: i128, truncate_i128);
    impl_truncate_from!(truncate_usize: isize, truncate_isize);
    impl_truncate_from!(truncate_u64: i64, truncate_i64);
    impl_truncate_from!(truncate_u32: i32, truncate_i32);
    impl_truncate_from!(truncate_u16: i16, truncate_i16);
    impl_truncate_from!(truncate_u8: i8, truncate_i8);
}
