use super::{big_digit, BigUint};

pub trait TruncateFrom<S> {
    /// Returns the input reduced modulo 2 to the power of the size of the target type.
    /// This is analogous to the behavior of `as` conversion on integral types.
    fn truncate_from(_: &S) -> Self;
}

impl BigUint {
    /// Returns the input reduced modulo 2 to the power of the size of the target type.
    /// This is analogous to the behavior of `as` conversion on integral types.
    #[inline]
    pub fn truncate<T: TruncateFrom<Self>>(&self) -> T {
        TruncateFrom::truncate_from(self)
    }
}

macro_rules! impl_truncate_large {
    ($T:ty, $size:expr) => {
        impl TruncateFrom<BigUint> for $T {
            #[inline]
            fn truncate_from(n: &BigUint) -> Self {
                let mut ret = 0;
                let mut bits = 0;

                for i in n.data.iter() {
                    if bits >= $size {
                        break;
                    }

                    ret += Self::from(*i) << bits;
                    bits += big_digit::BITS;
                }
                ret
            }
        }
    };
}

macro_rules! impl_truncate_from {
    ($from:ty: $($T:ty),* $(,)*) => {
        $(
            impl TruncateFrom<BigUint> for $T {
                #[inline]
                fn truncate_from(n: &BigUint) -> Self {
                    n.truncate::<$from>() as Self
                }
            }
        )*
    };
}

#[cfg(not(u64_digit))]
impl_truncate_large!(u32, 32);
impl_truncate_large!(u64, 64);
impl_truncate_large!(u128, 128);

#[cfg(not(u64_digit))]
impl_truncate_from! { u32: u16, u8 }

#[cfg(u64_digit)]
impl_truncate_from! { u64: u32, u16, u8 }

#[cfg(target_pointer_width = "64")]
impl_truncate_from!(u64: usize);
#[cfg(target_pointer_width = "32")]
impl_truncate_from!(u32: usize);
#[cfg(target_pointer_width = "16")]
impl_truncate_from!(u16: usize);

impl_truncate_from!(u128: i128);
impl_truncate_from!(usize: isize);
impl_truncate_from!(u64: i64);
impl_truncate_from!(u32: i32);
impl_truncate_from!(u16: i16);
impl_truncate_from!(u8: i8);
