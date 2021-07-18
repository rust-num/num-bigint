use super::BigIntSmall;
use super::Sign::NoSign;

use core::ops::{Shl, ShlAssign, Shr, ShrAssign};
use num_traits::{PrimInt, Signed, Zero};

macro_rules! impl_shift {
    (@ref $Shx:ident :: $shx:ident, $ShxAssign:ident :: $shx_assign:ident, $rhs:ty) => {
        impl<'b> $Shx<&'b $rhs> for BigIntSmall {
            type Output = BigIntSmall;

            #[inline]
            fn $shx(self, rhs: &'b $rhs) -> BigIntSmall {
                $Shx::$shx(self, *rhs)
            }
        }
        impl<'a, 'b> $Shx<&'b $rhs> for &'a BigIntSmall {
            type Output = BigIntSmall;

            #[inline]
            fn $shx(self, rhs: &'b $rhs) -> BigIntSmall {
                $Shx::$shx(self, *rhs)
            }
        }
        impl<'b> $ShxAssign<&'b $rhs> for BigIntSmall {
            #[inline]
            fn $shx_assign(&mut self, rhs: &'b $rhs) {
                $ShxAssign::$shx_assign(self, *rhs);
            }
        }
    };
    ($($rhs:ty),+) => {$(
        impl Shl<$rhs> for BigIntSmall {
            type Output = BigIntSmall;

            #[inline]
            fn shl(self, rhs: $rhs) -> BigIntSmall {
                BigIntSmall::from_biguint(self.sign, self.data << rhs)
            }
        }
        impl<'a> Shl<$rhs> for &'a BigIntSmall {
            type Output = BigIntSmall;

            #[inline]
            fn shl(self, rhs: $rhs) -> BigIntSmall {
                BigIntSmall::from_biguint(self.sign, &self.data << rhs)
            }
        }
        impl ShlAssign<$rhs> for BigIntSmall {
            #[inline]
            fn shl_assign(&mut self, rhs: $rhs) {
                self.data <<= rhs
            }
        }
        impl_shift! { @ref Shl::shl, ShlAssign::shl_assign, $rhs }

        impl Shr<$rhs> for BigIntSmall {
            type Output = BigIntSmall;

            #[inline]
            fn shr(self, rhs: $rhs) -> BigIntSmall {
                let round_down = shr_round_down(&self, rhs);
                let data = self.data >> rhs;
                let data = if round_down { data + 1u8 } else { data };
                BigIntSmall::from_biguint(self.sign, data)
            }
        }
        impl<'a> Shr<$rhs> for &'a BigIntSmall {
            type Output = BigIntSmall;

            #[inline]
            fn shr(self, rhs: $rhs) -> BigIntSmall {
                let round_down = shr_round_down(self, rhs);
                let data = &self.data >> rhs;
                let data = if round_down { data + 1u8 } else { data };
                BigIntSmall::from_biguint(self.sign, data)
            }
        }
        impl ShrAssign<$rhs> for BigIntSmall {
            #[inline]
            fn shr_assign(&mut self, rhs: $rhs) {
                let round_down = shr_round_down(self, rhs);
                self.data >>= rhs;
                if round_down {
                    self.data += 1u8;
                } else if self.data.is_zero() {
                    self.sign = NoSign;
                }
            }
        }
        impl_shift! { @ref Shr::shr, ShrAssign::shr_assign, $rhs }
    )*};
}

impl_shift! { u8, u16, u32, u64, u128, usize }
impl_shift! { i8, i16, i32, i64, i128, isize }

// Negative values need a rounding adjustment if there are any ones in the
// bits that are getting shifted out.
fn shr_round_down<T: PrimInt>(i: &BigIntSmall, shift: T) -> bool {
    if i.is_negative() {
        let zeros = i.trailing_zeros().expect("negative values are non-zero");
        shift > T::zero() && shift.to_u64().map(|shift| zeros < shift).unwrap_or(true)
    } else {
        false
    }
}
