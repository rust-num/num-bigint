// `Add`/`Sub` ops may flip from `BigInt` to its `BigUint` magnitude
#![allow(clippy::suspicious_arithmetic_impl)]

use crate::std_alloc::{String, Vec};
// use crate::BigInt;
use crate::Sign;
use core::cmp::Ordering::{self, Equal};
use core::default::Default;
use core::fmt;
use core::hash;
use core::ops::{Neg, Not};
use core::str;
use core::{i128, u128};
use core::{i64, u64};
use std::borrow::Cow;

use num_integer::{Integer, Roots};
use num_traits::{Num, One, Pow, Signed, Zero};

use smallvec::SmallVec;

use self::Sign::{Minus, NoSign, Plus};

use crate::big_digit::BigDigit;
use crate::biguint::to_str_radix_reversed;
use crate::biguint::{BigUint, IntDigits, U32Digits, U64Digits};

mod addition;
mod division;
mod multiplication;
mod subtraction;

mod bits;
mod convert;
mod power;
mod shift;

/// A big signed integer type.
// pub struct BigIntSmall {
//     sign_i: Sign,
//     data_i: BigUint,
// }
pub enum BigIntSmall {
    MinusBig(BigUint),
    MinusMedium([BigDigit; 3]),
    MinusSmall(BigDigit),
    PlusSmall(BigDigit),
    PlusMedium([BigDigit; 3]),
    PlusBig(BigUint),
}
use BigIntSmall::*;

// Note: derived `Clone` doesn't specialize `clone_from`,
// but we want to keep the allocation in `data`.
impl Clone for BigIntSmall {
    #[inline]
    fn clone(&self) -> Self {
        match self {
            MinusBig(big) => MinusBig(big.clone()),
            &MinusMedium(medium) => MinusMedium(medium),
            &MinusSmall(small) => MinusSmall(small),
            &PlusSmall(small) => PlusSmall(small),
            &PlusMedium(medium) => PlusMedium(medium),
            PlusBig(big) => PlusBig(big.clone()),
        }
        // BigIntSmall {
        //     sign_i: self.sign(),
        //     data_i: self.data_i.clone(),
        // }
    }

    // #[inline]
    // fn clone_from(&mut self, other: &Self) {
    //     *self.mut_sign() = other.sign();
    //     self.data_i.clone_from(&other.data_i);
    // }
}

impl hash::Hash for BigIntSmall {
    #[inline]
    fn hash<H: hash::Hasher>(&self, state: &mut H) {
        debug_assert!((self.sign() != NoSign) ^ self.data().is_zero());
        self.sign().hash(state);
        if self.sign() != NoSign {
            self.data().hash(state);
        }
    }
}

impl PartialEq for BigIntSmall {
    #[inline]
    fn eq(&self, other: &BigIntSmall) -> bool {
        debug_assert!((self.sign() != NoSign) ^ self.data().is_zero());
        debug_assert!((other.sign() != NoSign) ^ other.data().is_zero());
        self.sign() == other.sign() && (self.sign() == NoSign || self.data() == other.data())
    }
}

impl Eq for BigIntSmall {}

impl PartialOrd for BigIntSmall {
    #[inline]
    fn partial_cmp(&self, other: &BigIntSmall) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for BigIntSmall {
    #[inline]
    fn cmp(&self, other: &BigIntSmall) -> Ordering {
        debug_assert!((self.sign() != NoSign) ^ self.data().is_zero());
        debug_assert!((other.sign() != NoSign) ^ other.data().is_zero());
        let scmp = self.sign().cmp(&other.sign());
        if scmp != Equal {
            return scmp;
        }

        match self.sign() {
            NoSign => Equal,
            Plus => self.data().cmp(&other.data()),
            Minus => other.data().cmp(&self.data()),
        }
    }
}

impl Default for BigIntSmall {
    #[inline]
    fn default() -> BigIntSmall {
        Zero::zero()
    }
}

impl fmt::Debug for BigIntSmall {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt::Display::fmt(self, f)
    }
}

impl fmt::Display for BigIntSmall {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.pad_integral(!self.is_negative(), "", &self.data().to_str_radix(10))
    }
}

impl fmt::Binary for BigIntSmall {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.pad_integral(!self.is_negative(), "0b", &self.data().to_str_radix(2))
    }
}

impl fmt::Octal for BigIntSmall {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.pad_integral(!self.is_negative(), "0o", &self.data().to_str_radix(8))
    }
}

impl fmt::LowerHex for BigIntSmall {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.pad_integral(!self.is_negative(), "0x", &self.data().to_str_radix(16))
    }
}

impl fmt::UpperHex for BigIntSmall {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut s = self.data().to_str_radix(16);
        s.make_ascii_uppercase();
        f.pad_integral(!self.is_negative(), "0x", &s)
    }
}

// !-2 = !...f fe = ...0 01 = +1
// !-1 = !...f ff = ...0 00 =  0
// ! 0 = !...0 00 = ...f ff = -1
// !+1 = !...0 01 = ...f fe = -2
impl Not for BigIntSmall {
    type Output = BigIntSmall;

    fn not(self) -> BigIntSmall {
        let (sign, mut uint) = self.into_parts();
        match sign {
            NoSign | Plus => {
                uint += 1u32;
                BigIntSmall::from_biguint(Minus, uint)
            }
            Minus => {
                uint -= 1u32;
                BigIntSmall::from_biguint(Plus, uint)
            }
        }
        // self
    }
}

impl<'a> Not for &'a BigIntSmall {
    type Output = BigIntSmall;

    fn not(self) -> BigIntSmall {
        match self.sign() {
            NoSign => -BigIntSmall::one(),
            Plus => -BigIntSmall::from(&self.data() as &BigUint + 1u32),
            Minus => BigIntSmall::from(&self.data() as &BigUint - 1u32),
        }
    }
}

impl Zero for BigIntSmall {
    #[inline]
    fn zero() -> BigIntSmall {
        BigIntSmall::PlusBig(BigUint::zero())
        // BigIntSmall::from(BigUint::zero())
    }

    #[inline]
    fn set_zero(&mut self) {
        *self = Self::zero()
        // self.data.set_zero();
        // *self.mut_sign() = NoSign;
    }

    #[inline]
    fn is_zero(&self) -> bool {
        self.sign() == NoSign
    }
}

impl One for BigIntSmall {
    #[inline]
    fn one() -> BigIntSmall {
        BigIntSmall::PlusBig(BigUint::one())
        // BigIntSmall::from(BigUint::one())
    }

    #[inline]
    fn set_one(&mut self) {
        *self = Self::one();
        // self.data.set_one();
        // *self.mut_sign() = Plus;
    }

    #[inline]
    fn is_one(&self) -> bool {
        self.sign() == Plus && self.data().is_one()
    }
}

impl Signed for BigIntSmall {
    #[inline]
    fn abs(&self) -> BigIntSmall {
        match self.sign() {
            Plus | NoSign => self.clone(),
            Minus => BigIntSmall::from(self.data().into_owned()),
        }
    }

    #[inline]
    fn abs_sub(&self, other: &BigIntSmall) -> BigIntSmall {
        if *self <= *other {
            Zero::zero()
        } else {
            self - other
        }
    }

    #[inline]
    fn signum(&self) -> BigIntSmall {
        match self.sign() {
            Plus => BigIntSmall::one(),
            Minus => -BigIntSmall::one(),
            NoSign => BigIntSmall::zero(),
        }
    }

    #[inline]
    fn is_positive(&self) -> bool {
        self.sign() == Plus
    }

    #[inline]
    fn is_negative(&self) -> bool {
        self.sign() == Minus
    }
}

trait UnsignedAbs {
    type Unsigned;

    /// A convenience method for getting the absolute value of a signed primitive as unsigned
    /// See also `unsigned_abs`: https://github.com/rust-lang/rust/issues/74913
    fn uabs(self) -> Self::Unsigned;

    fn checked_uabs(self) -> CheckedUnsignedAbs<Self::Unsigned>;
}

enum CheckedUnsignedAbs<T> {
    Positive(T),
    Negative(T),
}
use self::CheckedUnsignedAbs::{Negative, Positive};

macro_rules! impl_unsigned_abs {
    ($Signed:ty, $Unsigned:ty) => {
        impl UnsignedAbs for $Signed {
            type Unsigned = $Unsigned;

            #[inline]
            fn uabs(self) -> $Unsigned {
                self.wrapping_abs() as $Unsigned
            }

            #[inline]
            fn checked_uabs(self) -> CheckedUnsignedAbs<Self::Unsigned> {
                if self >= 0 {
                    Positive(self as $Unsigned)
                } else {
                    Negative(self.wrapping_neg() as $Unsigned)
                }
            }
        }
    };
}
impl_unsigned_abs!(i8, u8);
impl_unsigned_abs!(i16, u16);
impl_unsigned_abs!(i32, u32);
impl_unsigned_abs!(i64, u64);
impl_unsigned_abs!(i128, u128);
impl_unsigned_abs!(isize, usize);

impl Neg for BigIntSmall {
    type Output = BigIntSmall;

    #[inline]
    fn neg(mut self) -> BigIntSmall {
        self.toggle_sign();
        self
    }
}

impl<'a> Neg for &'a BigIntSmall {
    type Output = BigIntSmall;

    #[inline]
    fn neg(self) -> BigIntSmall {
        self.clone().neg()
    }
}

impl Integer for BigIntSmall {
    #[inline]
    fn div_rem(&self, other: &BigIntSmall) -> (BigIntSmall, BigIntSmall) {
        // r.sign == self.sign
        let (d_ui, r_ui) = self.data().div_rem(&other.data());
        let d = BigIntSmall::from_biguint(self.sign(), d_ui);
        let r = BigIntSmall::from_biguint(self.sign(), r_ui);
        if other.is_negative() {
            (-d, r)
        } else {
            (d, r)
        }
    }

    #[inline]
    fn div_floor(&self, other: &BigIntSmall) -> BigIntSmall {
        let (d_ui, m) = self.data().div_mod_floor(&other.data());
        let d = BigIntSmall::from(d_ui);
        match (self.sign(), other.sign()) {
            (Plus, Plus) | (NoSign, Plus) | (Minus, Minus) => d,
            (Plus, Minus) | (NoSign, Minus) | (Minus, Plus) => {
                if m.is_zero() {
                    -d
                } else {
                    -d - 1u32
                }
            }
            (_, NoSign) => unreachable!(),
        }
    }

    #[inline]
    fn mod_floor(&self, other: &BigIntSmall) -> BigIntSmall {
        // m.sign == other.sign
        let m_ui = self.data().mod_floor(&other.data());
        let m = BigIntSmall::from_biguint(other.sign(), m_ui);
        match (self.sign(), other.sign()) {
            (Plus, Plus) | (NoSign, Plus) | (Minus, Minus) => m,
            (Plus, Minus) | (NoSign, Minus) | (Minus, Plus) => {
                if m.is_zero() {
                    m
                } else {
                    other - m
                }
            }
            (_, NoSign) => unreachable!(),
        }
    }

    fn div_mod_floor(&self, other: &BigIntSmall) -> (BigIntSmall, BigIntSmall) {
        // m.sign == other.sign
        let (d_ui, m_ui) = self.data().div_mod_floor(&other.data());
        let d = BigIntSmall::from(d_ui);
        let m = BigIntSmall::from_biguint(other.sign(), m_ui);
        match (self.sign(), other.sign()) {
            (Plus, Plus) | (NoSign, Plus) | (Minus, Minus) => (d, m),
            (Plus, Minus) | (NoSign, Minus) | (Minus, Plus) => {
                if m.is_zero() {
                    (-d, m)
                } else {
                    (-d - 1u32, other - m)
                }
            }
            (_, NoSign) => unreachable!(),
        }
    }

    #[inline]
    fn div_ceil(&self, other: &Self) -> Self {
        let (d_ui, m) = self.data().div_mod_floor(&other.data());
        let d = BigIntSmall::from(d_ui);
        match (self.sign(), other.sign()) {
            (Plus, Minus) | (NoSign, Minus) | (Minus, Plus) => -d,
            (Plus, Plus) | (NoSign, Plus) | (Minus, Minus) => {
                if m.is_zero() {
                    d
                } else {
                    d + 1u32
                }
            }
            (_, NoSign) => unreachable!(),
        }
    }

    /// Calculates the Greatest Common Divisor (GCD) of the number and `other`.
    ///
    /// The result is always positive.
    #[inline]
    fn gcd(&self, other: &BigIntSmall) -> BigIntSmall {
        match (self, other) {
            (PlusSmall(a), PlusSmall(b)) => PlusSmall(a.gcd(b)),
            (PlusSmall(a), MinusSmall(b)) => PlusSmall(a.gcd(b)),
            (MinusSmall(a), PlusSmall(b)) => PlusSmall(a.gcd(b)),
            (MinusSmall(a), MinusSmall(b)) => PlusSmall(a.gcd(b)),
            _ => BigIntSmall::from(self.data().gcd(&other.data())),
        }
    }

    /// Calculates the Lowest Common Multiple (LCM) of the number and `other`.
    #[inline]
    fn lcm(&self, other: &BigIntSmall) -> BigIntSmall {
        BigIntSmall::from(self.data().lcm(&other.data()))
    }

    /// Calculates the Greatest Common Divisor (GCD) and
    /// Lowest Common Multiple (LCM) together.
    #[inline]
    fn gcd_lcm(&self, other: &BigIntSmall) -> (BigIntSmall, BigIntSmall) {
        let (gcd, lcm) = self.data().gcd_lcm(&other.data());
        (BigIntSmall::from(gcd), BigIntSmall::from(lcm))
    }

    /// Greatest common divisor, least common multiple, and Bézout coefficients.
    #[inline]
    fn extended_gcd_lcm(
        &self,
        other: &BigIntSmall,
    ) -> (num_integer::ExtendedGcd<BigIntSmall>, BigIntSmall) {
        let egcd = self.extended_gcd(other);
        let lcm = if egcd.gcd.is_zero() {
            BigIntSmall::zero()
        } else {
            BigIntSmall::from(
                &self.data() as &BigUint / &egcd.gcd.data() as &BigUint * &other.data() as &BigUint,
            )
        };
        (egcd, lcm)
    }

    /// Deprecated, use `is_multiple_of` instead.
    #[inline]
    fn divides(&self, other: &BigIntSmall) -> bool {
        self.is_multiple_of(other)
    }

    /// Returns `true` if the number is a multiple of `other`.
    #[inline]
    fn is_multiple_of(&self, other: &BigIntSmall) -> bool {
        self.data().is_multiple_of(&other.data())
    }

    /// Returns `true` if the number is divisible by `2`.
    #[inline]
    fn is_even(&self) -> bool {
        self.data().is_even()
    }

    /// Returns `true` if the number is not divisible by `2`.
    #[inline]
    fn is_odd(&self) -> bool {
        self.data().is_odd()
    }

    /// Rounds up to nearest multiple of argument.
    #[inline]
    fn next_multiple_of(&self, other: &Self) -> Self {
        let m = self.mod_floor(other);
        if m.is_zero() {
            self.clone()
        } else {
            self + (other - m)
        }
    }
    /// Rounds down to nearest multiple of argument.
    #[inline]
    fn prev_multiple_of(&self, other: &Self) -> Self {
        self - self.mod_floor(other)
    }
}

impl Roots for BigIntSmall {
    fn nth_root(&self, n: u32) -> Self {
        assert!(
            !(self.is_negative() && n.is_even()),
            "root of degree {} is imaginary",
            n
        );

        BigIntSmall::from_biguint(self.sign(), self.data().nth_root(n))
    }

    fn sqrt(&self) -> Self {
        assert!(!self.is_negative(), "square root is imaginary");

        BigIntSmall::from_biguint(self.sign(), self.data().sqrt())
    }

    fn cbrt(&self) -> Self {
        BigIntSmall::from_biguint(self.sign(), self.data().cbrt())
    }
}

impl IntDigits for BigIntSmall {
    #[inline]
    fn digits(&self) -> &[BigDigit] {
        // Can easily be written for compact formats
        match self {
            MinusBig(big) => big.digits(),
            MinusMedium(digits) => digits,
            MinusSmall(digit) => std::slice::from_ref(digit),
            PlusSmall(digit) => std::slice::from_ref(digit),
            PlusMedium(digits) => digits,
            PlusBig(big) => big.digits(),
        }
        // self.data_i.digits()
    }

    #[inline]
    fn digits_mut(&mut self) -> &mut SmallVec<[BigDigit; BigUint::INLINED]> {
        todo!()
        // // Must convert from compact format first.
        // match self {
        //     MinusBig(big) => big.digits_mut(),
        //     MinusMedium(digits) => {
        //         let uint = BigUint::from_digits(digits);
        //         *self = MinusBig(uint);
        //         self.digits_mut()
        //     }
        //     MinusSmall(digit) => {
        //         let uint = BigUint::from(*digit);
        //         *self = MinusBig(uint);
        //         self.digits_mut()
        //     }
        //     PlusSmall(digit) => {
        //         let uint = BigUint::from(*digit);
        //         *self = MinusBig(uint);
        //         self.digits_mut()
        //     }
        //     PlusMedium(digits) => {
        //         let uint = BigUint::from_digits(digits);
        //         *self = PlusBig(uint);
        //         self.digits_mut()
        //     }
        //     PlusBig(big) => big.digits_mut(),
        // }
        // self.data_i.digits_mut()
    }

    #[inline]
    fn normalize(&mut self) {
        match self {
            MinusBig(big) => {
                big.normalize();
                if big.is_zero() {
                    self.set_zero();
                }
            }
            MinusMedium([0, 0, 0]) => {
                self.set_zero();
            }
            MinusMedium(_) => {}
            MinusSmall(0) => {
                self.set_zero();
            }
            MinusSmall(_) => {}
            PlusSmall(_) => {}
            PlusMedium(_) => {}
            PlusBig(big) => {
                big.normalize();
                if big.is_zero() {
                    self.set_zero();
                }
            }
        }
        // self.data_i.normalize();
        // if self.data_i.is_zero() {
        //     *self.mut_sign() = NoSign;
        // }
    }
    #[inline]
    fn capacity(&self) -> usize {
        self.data().capacity()
    }
    #[inline]
    fn len(&self) -> usize {
        self.data().len()
    }
}

/// A generic trait for converting a value to a `BigInt`. This may return
/// `None` when converting from `f32` or `f64`, and will always succeed
/// when converting from any integer or unsigned primitive, or `BigUint`.
pub trait ToBigInt {
    /// Converts the value of `self` to a `BigInt`.
    fn to_bigint(&self) -> Option<BigIntSmall>;
}

impl BigIntSmall {
    /// Creates and initializes a BigInt.
    ///
    /// The base 2<sup>32</sup> digits are ordered least significant digit first.
    #[inline]
    pub fn new(sign: Sign, digits: Vec<u32>) -> BigIntSmall {
        BigIntSmall::from_biguint(sign, BigUint::new(digits))
    }

    /// Creates and initializes a `BigInt`.
    ///
    /// The base 2<sup>32</sup> digits are ordered least significant digit first.
    #[inline]
    pub fn from_biguint(mut sign: Sign, mut data: BigUint) -> BigIntSmall {
        if sign == NoSign {
            data.assign_from_slice(&[]);
        } else if data.is_zero() {
            sign = NoSign;
        }

        match sign {
            Minus => MinusBig(data),
            NoSign => PlusBig(data),
            Plus => PlusBig(data),
        }
        // BigIntSmall {
        //     sign_i: sign,
        //     data_i: data,
        // }
    }

    /// Creates and initializes a `BigInt`.
    ///
    /// The base 2<sup>32</sup> digits are ordered least significant digit first.
    #[inline]
    pub fn from_slice(sign: Sign, slice: &[u32]) -> BigIntSmall {
        BigIntSmall::from_biguint(sign, BigUint::from_slice(slice))
    }

    /// Reinitializes a `BigInt`.
    ///
    /// The base 2<sup>32</sup> digits are ordered least significant digit first.
    #[inline]
    pub fn assign_from_slice(&mut self, sign: Sign, slice: &[u32]) {
        *self = BigIntSmall::from_slice(sign, slice)
        // if sign == NoSign {
        //     self.set_zero();
        // } else {
        //     self.data_i.assign_from_slice(slice);
        //     *self.mut_sign() = if self.data_i.is_zero() { NoSign } else { sign };
        // }
    }

    /// Creates and initializes a `BigInt`.
    ///
    /// The bytes are in big-endian byte order.
    ///
    /// # Examples
    ///
    /// ```
    /// use num_bigint::{BigInt, Sign};
    ///
    /// assert_eq!(BigInt::from_bytes_be(Sign::Plus, b"A"),
    ///            BigInt::parse_bytes(b"65", 10).unwrap());
    /// assert_eq!(BigInt::from_bytes_be(Sign::Plus, b"AA"),
    ///            BigInt::parse_bytes(b"16705", 10).unwrap());
    /// assert_eq!(BigInt::from_bytes_be(Sign::Plus, b"AB"),
    ///            BigInt::parse_bytes(b"16706", 10).unwrap());
    /// assert_eq!(BigInt::from_bytes_be(Sign::Plus, b"Hello world!"),
    ///            BigInt::parse_bytes(b"22405534230753963835153736737", 10).unwrap());
    /// ```
    #[inline]
    pub fn from_bytes_be(sign: Sign, bytes: &[u8]) -> BigIntSmall {
        BigIntSmall::from_biguint(sign, BigUint::from_bytes_be(bytes))
    }

    /// Creates and initializes a `BigInt`.
    ///
    /// The bytes are in little-endian byte order.
    #[inline]
    pub fn from_bytes_le(sign: Sign, bytes: &[u8]) -> BigIntSmall {
        BigIntSmall::from_biguint(sign, BigUint::from_bytes_le(bytes))
    }

    /// Creates and initializes a `BigInt` from an array of bytes in
    /// two's complement binary representation.
    ///
    /// The digits are in big-endian base 2<sup>8</sup>.
    #[inline]
    pub fn from_signed_bytes_be(digits: &[u8]) -> BigIntSmall {
        convert::from_signed_bytes_be(digits)
    }

    /// Creates and initializes a `BigInt` from an array of bytes in two's complement.
    ///
    /// The digits are in little-endian base 2<sup>8</sup>.
    #[inline]
    pub fn from_signed_bytes_le(digits: &[u8]) -> BigIntSmall {
        convert::from_signed_bytes_le(digits)
    }

    /// Creates and initializes a `BigInt`.
    ///
    /// # Examples
    ///
    /// ```
    /// use num_bigint::{BigInt, ToBigInt};
    ///
    /// assert_eq!(BigInt::parse_bytes(b"1234", 10), ToBigInt::to_bigint(&1234));
    /// assert_eq!(BigInt::parse_bytes(b"ABCD", 16), ToBigInt::to_bigint(&0xABCD));
    /// assert_eq!(BigInt::parse_bytes(b"G", 16), None);
    /// ```
    #[inline]
    pub fn parse_bytes(buf: &[u8], radix: u32) -> Option<BigIntSmall> {
        let s = str::from_utf8(buf).ok()?;
        BigIntSmall::from_str_radix(s, radix).ok()
    }

    /// Creates and initializes a `BigInt`. Each u8 of the input slice is
    /// interpreted as one digit of the number
    /// and must therefore be less than `radix`.
    ///
    /// The bytes are in big-endian byte order.
    /// `radix` must be in the range `2...256`.
    ///
    /// # Examples
    ///
    /// ```
    /// use num_bigint::{BigInt, Sign};
    ///
    /// let inbase190 = vec![15, 33, 125, 12, 14];
    /// let a = BigInt::from_radix_be(Sign::Minus, &inbase190, 190).unwrap();
    /// assert_eq!(a.to_radix_be(190), (Sign:: Minus, inbase190));
    /// ```
    pub fn from_radix_be(sign: Sign, buf: &[u8], radix: u32) -> Option<BigIntSmall> {
        let u = BigUint::from_radix_be(buf, radix)?;
        Some(BigIntSmall::from_biguint(sign, u))
    }

    /// Creates and initializes a `BigInt`. Each u8 of the input slice is
    /// interpreted as one digit of the number
    /// and must therefore be less than `radix`.
    ///
    /// The bytes are in little-endian byte order.
    /// `radix` must be in the range `2...256`.
    ///
    /// # Examples
    ///
    /// ```
    /// use num_bigint::{BigInt, Sign};
    ///
    /// let inbase190 = vec![14, 12, 125, 33, 15];
    /// let a = BigInt::from_radix_be(Sign::Minus, &inbase190, 190).unwrap();
    /// assert_eq!(a.to_radix_be(190), (Sign::Minus, inbase190));
    /// ```
    pub fn from_radix_le(sign: Sign, buf: &[u8], radix: u32) -> Option<BigIntSmall> {
        let u = BigUint::from_radix_le(buf, radix)?;
        Some(BigIntSmall::from_biguint(sign, u))
    }

    /// Returns the sign and the byte representation of the `BigInt` in big-endian byte order.
    ///
    /// # Examples
    ///
    /// ```
    /// use num_bigint::{ToBigInt, Sign};
    ///
    /// let i = -1125.to_bigint().unwrap();
    /// assert_eq!(i.to_bytes_be(), (Sign::Minus, vec![4, 101]));
    /// ```
    #[inline]
    pub fn to_bytes_be(&self) -> (Sign, Vec<u8>) {
        (self.sign(), self.data().to_bytes_be())
    }

    /// Returns the sign and the byte representation of the `BigInt` in little-endian byte order.
    ///
    /// # Examples
    ///
    /// ```
    /// use num_bigint::{ToBigInt, Sign};
    ///
    /// let i = -1125.to_bigint().unwrap();
    /// assert_eq!(i.to_bytes_le(), (Sign::Minus, vec![101, 4]));
    /// ```
    #[inline]
    pub fn to_bytes_le(&self) -> (Sign, Vec<u8>) {
        (self.sign(), self.data().to_bytes_le())
    }

    /// Returns the sign and the `u32` digits representation of the `BigInt` ordered least
    /// significant digit first.
    ///
    /// # Examples
    ///
    /// ```
    /// use num_bigint::{BigInt, Sign};
    ///
    /// assert_eq!(BigInt::from(-1125).to_u32_digits(), (Sign::Minus, vec![1125]));
    /// assert_eq!(BigInt::from(4294967295u32).to_u32_digits(), (Sign::Plus, vec![4294967295]));
    /// assert_eq!(BigInt::from(4294967296u64).to_u32_digits(), (Sign::Plus, vec![0, 1]));
    /// assert_eq!(BigInt::from(-112500000000i64).to_u32_digits(), (Sign::Minus, vec![830850304, 26]));
    /// assert_eq!(BigInt::from(112500000000i64).to_u32_digits(), (Sign::Plus, vec![830850304, 26]));
    /// ```
    #[inline]
    pub fn to_u32_digits(&self) -> (Sign, Vec<u32>) {
        (self.sign(), self.data().to_u32_digits())
    }

    /// Returns the sign and the `u64` digits representation of the `BigInt` ordered least
    /// significant digit first.
    ///
    /// # Examples
    ///
    /// ```
    /// use num_bigint::{BigInt, Sign};
    ///
    /// assert_eq!(BigInt::from(-1125).to_u64_digits(), (Sign::Minus, vec![1125]));
    /// assert_eq!(BigInt::from(4294967295u32).to_u64_digits(), (Sign::Plus, vec![4294967295]));
    /// assert_eq!(BigInt::from(4294967296u64).to_u64_digits(), (Sign::Plus, vec![4294967296]));
    /// assert_eq!(BigInt::from(-112500000000i64).to_u64_digits(), (Sign::Minus, vec![112500000000]));
    /// assert_eq!(BigInt::from(112500000000i64).to_u64_digits(), (Sign::Plus, vec![112500000000]));
    /// assert_eq!(BigInt::from(1u128 << 64).to_u64_digits(), (Sign::Plus, vec![0, 1]));
    /// ```
    #[inline]
    pub fn to_u64_digits(&self) -> (Sign, Vec<u64>) {
        (self.sign(), self.data().to_u64_digits())
    }

    /// Returns an iterator of `u32` digits representation of the `BigInt` ordered least
    /// significant digit first.
    ///
    /// # Examples
    ///
    /// ```
    /// use num_bigint::BigInt;
    ///
    /// assert_eq!(BigInt::from(-1125).iter_u32_digits().collect::<Vec<u32>>(), vec![1125]);
    /// assert_eq!(BigInt::from(4294967295u32).iter_u32_digits().collect::<Vec<u32>>(), vec![4294967295]);
    /// assert_eq!(BigInt::from(4294967296u64).iter_u32_digits().collect::<Vec<u32>>(), vec![0, 1]);
    /// assert_eq!(BigInt::from(-112500000000i64).iter_u32_digits().collect::<Vec<u32>>(), vec![830850304, 26]);
    /// assert_eq!(BigInt::from(112500000000i64).iter_u32_digits().collect::<Vec<u32>>(), vec![830850304, 26]);
    /// ```
    #[inline]
    pub fn iter_u32_digits(&self) -> U32Digits<'_> {
        match self {
            MinusBig(big) => big.iter_u32_digits(),
            MinusMedium(digits) => U32Digits::new(digits),
            MinusSmall(small) => U32Digits::new(std::slice::from_ref(small)),
            PlusSmall(small) => U32Digits::new(std::slice::from_ref(small)),
            PlusMedium(digits) => U32Digits::new(digits),
            PlusBig(big) => big.iter_u32_digits(),
        }
        // self.data_i.iter_u32_digits()
    }

    /// Returns an iterator of `u64` digits representation of the `BigInt` ordered least
    /// significant digit first.
    ///
    /// # Examples
    ///
    /// ```
    /// use num_bigint::BigInt;
    ///
    /// assert_eq!(BigInt::from(-1125).iter_u64_digits().collect::<Vec<u64>>(), vec![1125u64]);
    /// assert_eq!(BigInt::from(4294967295u32).iter_u64_digits().collect::<Vec<u64>>(), vec![4294967295u64]);
    /// assert_eq!(BigInt::from(4294967296u64).iter_u64_digits().collect::<Vec<u64>>(), vec![4294967296u64]);
    /// assert_eq!(BigInt::from(-112500000000i64).iter_u64_digits().collect::<Vec<u64>>(), vec![112500000000u64]);
    /// assert_eq!(BigInt::from(112500000000i64).iter_u64_digits().collect::<Vec<u64>>(), vec![112500000000u64]);
    /// assert_eq!(BigInt::from(1u128 << 64).iter_u64_digits().collect::<Vec<u64>>(), vec![0, 1]);
    /// ```
    #[inline]
    pub fn iter_u64_digits(&self) -> U64Digits<'_> {
        match self {
            MinusBig(big) => big.iter_u64_digits(),
            MinusMedium(digits) => U64Digits::new(digits),
            MinusSmall(small) => U64Digits::new(std::slice::from_ref(small)),
            PlusSmall(small) => U64Digits::new(std::slice::from_ref(small)),
            PlusMedium(digits) => U64Digits::new(digits),
            PlusBig(big) => big.iter_u64_digits(),
        }
        // self.data_i.iter_u64_digits()
    }

    /// Returns the two's-complement byte representation of the `BigInt` in big-endian byte order.
    ///
    /// # Examples
    ///
    /// ```
    /// use num_bigint::ToBigInt;
    ///
    /// let i = -1125.to_bigint().unwrap();
    /// assert_eq!(i.to_signed_bytes_be(), vec![251, 155]);
    /// ```
    #[inline]
    pub fn to_signed_bytes_be(&self) -> Vec<u8> {
        convert::to_signed_bytes_be(self)
    }

    /// Returns the two's-complement byte representation of the `BigInt` in little-endian byte order.
    ///
    /// # Examples
    ///
    /// ```
    /// use num_bigint::ToBigInt;
    ///
    /// let i = -1125.to_bigint().unwrap();
    /// assert_eq!(i.to_signed_bytes_le(), vec![155, 251]);
    /// ```
    #[inline]
    pub fn to_signed_bytes_le(&self) -> Vec<u8> {
        convert::to_signed_bytes_le(self)
    }

    /// Returns the integer formatted as a string in the given radix.
    /// `radix` must be in the range `2...36`.
    ///
    /// # Examples
    ///
    /// ```
    /// use num_bigint::BigInt;
    ///
    /// let i = BigInt::parse_bytes(b"ff", 16).unwrap();
    /// assert_eq!(i.to_str_radix(16), "ff");
    /// ```
    #[inline]
    pub fn to_str_radix(&self, radix: u32) -> String {
        let mut v = to_str_radix_reversed(&self.data() as &BigUint, radix);

        if self.is_negative() {
            v.push(b'-');
        }

        v.reverse();
        unsafe { String::from_utf8_unchecked(v) }
    }

    /// Returns the integer in the requested base in big-endian digit order.
    /// The output is not given in a human readable alphabet but as a zero
    /// based u8 number.
    /// `radix` must be in the range `2...256`.
    ///
    /// # Examples
    ///
    /// ```
    /// use num_bigint::{BigInt, Sign};
    ///
    /// assert_eq!(BigInt::from(-0xFFFFi64).to_radix_be(159),
    ///            (Sign::Minus, vec![2, 94, 27]));
    /// // 0xFFFF = 65535 = 2*(159^2) + 94*159 + 27
    /// ```
    #[inline]
    pub fn to_radix_be(&self, radix: u32) -> (Sign, Vec<u8>) {
        (self.sign(), self.data().to_radix_be(radix))
    }

    /// Returns the integer in the requested base in little-endian digit order.
    /// The output is not given in a human readable alphabet but as a zero
    /// based u8 number.
    /// `radix` must be in the range `2...256`.
    ///
    /// # Examples
    ///
    /// ```
    /// use num_bigint::{BigInt, Sign};
    ///
    /// assert_eq!(BigInt::from(-0xFFFFi64).to_radix_le(159),
    ///            (Sign::Minus, vec![27, 94, 2]));
    /// // 0xFFFF = 65535 = 27 + 94*159 + 2*(159^2)
    /// ```
    #[inline]
    pub fn to_radix_le(&self, radix: u32) -> (Sign, Vec<u8>) {
        (self.sign(), self.data().to_radix_le(radix))
    }

    /// Returns the sign of the `BigInt` as a `Sign`.
    ///
    /// # Examples
    ///
    /// ```
    /// use num_bigint::{BigInt, Sign};
    /// use num_traits::Zero;
    ///
    /// assert_eq!(BigInt::from(1234).sign(), Sign::Plus);
    /// assert_eq!(BigInt::from(-4321).sign(), Sign::Minus);
    /// assert_eq!(BigInt::zero().sign(), Sign::NoSign);
    /// ```
    #[inline]
    pub fn sign(&self) -> Sign {
        match self {
            MinusBig(_) => Minus,
            MinusMedium(_) => Minus,
            MinusSmall(_) => Minus,
            PlusSmall(_) => Plus,
            PlusMedium(_) => Plus,
            PlusBig(big) => {
                if big.is_zero() {
                    NoSign
                } else {
                    Plus
                }
            }
        }
        // self.sign_i
    }

    // sign = -sign;
    fn toggle_sign(&mut self) {
        if self.sign() == NoSign {
            return;
        }
        match self.take() {
            MinusBig(big) => *self = PlusBig(big),
            MinusMedium(digits) => *self = PlusMedium(digits),
            MinusSmall(small) => *self = PlusSmall(small),
            PlusSmall(small) => *self = MinusSmall(small),
            PlusMedium(digits) => *self = MinusMedium(digits),
            PlusBig(big) => *self = MinusBig(big),
        }
    }

    fn set_sign(&mut self, sign: Sign) {
        if self.sign() == sign {
            return;
        }
        let uint = self.take().into_biguint();
        match sign {
            NoSign => self.set_zero(),
            Minus => *self = MinusBig(uint),
            Plus => *self = PlusBig(uint),
        }
    }

    // fn mut_sign(&mut self) -> &mut Sign {
    //     &mut self.sign_i
    // }

    fn data(&self) -> Cow<'_, BigUint> {
        match self {
            MinusBig(big) => Cow::Borrowed(big),
            MinusMedium(digits) => Cow::Owned(BigUint::from_digits(digits)),
            &MinusSmall(small) => Cow::Owned(BigUint::from(small)),
            &PlusSmall(small) => Cow::Owned(BigUint::from(small)),
            PlusMedium(digits) => Cow::Owned(BigUint::from_digits(digits)),
            PlusBig(big) => Cow::Borrowed(big),
        }
        // Cow::Borrowed(&self.data_i)
    }

    /// Returns the magnitude of the `BigInt` as a `BigUint`.
    ///
    /// # Examples
    ///
    /// ```
    /// use num_bigint::{BigInt, BigUint};
    /// use num_traits::Zero;
    ///
    /// assert_eq!(BigInt::from(1234).magnitude(), &BigUint::from(1234u32));
    /// assert_eq!(BigInt::from(-4321).magnitude(), &BigUint::from(4321u32));
    /// assert!(BigInt::zero().magnitude().is_zero());
    /// ```
    #[inline]
    pub fn magnitude(&self) -> &BigUint {
        panic!("There isn't always a BigUint inside a BigIntSmall")
        // &self.data_i
    }

    /// Convert this `BigInt` into its `Sign` and `BigUint` magnitude,
    /// the reverse of `BigInt::from_biguint`.
    ///
    /// # Examples
    ///
    /// ```
    /// use num_bigint::{BigInt, BigUint, Sign};
    /// use num_traits::Zero;
    ///
    /// assert_eq!(BigInt::from(1234).into_parts(), (Sign::Plus, BigUint::from(1234u32)));
    /// assert_eq!(BigInt::from(-4321).into_parts(), (Sign::Minus, BigUint::from(4321u32)));
    /// assert_eq!(BigInt::zero().into_parts(), (Sign::NoSign, BigUint::zero()));
    /// ```
    #[inline]
    pub fn into_parts(self) -> (Sign, BigUint) {
        (self.sign(), self.into_biguint())
    }

    pub fn into_biguint(self) -> BigUint {
        match self {
            MinusBig(big) => big,
            MinusMedium(digits) => BigUint::from_digits(&digits),
            MinusSmall(small) => BigUint::from(small),
            PlusSmall(small) => BigUint::from(small),
            PlusMedium(digits) => BigUint::from_digits(&digits),
            PlusBig(big) => big,
        }
        // self.data_i
    }

    pub fn biguint_mut(&mut self) -> &mut BigUint {
        match self {
            MinusBig(big) => big,
            MinusMedium(digits) => {
                let uint = BigUint::from_digits(digits);
                *self = MinusBig(uint);
                self.biguint_mut()
            }
            MinusSmall(small) => {
                let uint = BigUint::from(*small);
                *self = MinusBig(uint);
                self.biguint_mut()
            }
            PlusSmall(small) => {
                let uint = BigUint::from(*small);
                *self = PlusBig(uint);
                self.biguint_mut()
            }
            PlusMedium(digits) => {
                let uint = BigUint::from_digits(digits);
                *self = PlusBig(uint);
                self.biguint_mut()
            }
            PlusBig(big) => big,
        }
        // &mut self.data_i
    }

    fn take(&mut self) -> BigIntSmall {
        std::mem::replace(self, PlusBig(BigUint::zero()))
    }

    /// Determines the fewest bits necessary to express the `BigInt`,
    /// not including the sign.
    #[inline]
    pub fn bits(&self) -> u64 {
        self.data().bits()
    }

    /// Converts this `BigInt` into a `BigUint`, if it's not negative.
    #[inline]
    pub fn to_biguint(&self) -> Option<BigUint> {
        match self.sign() {
            Plus => Some(self.data().into_owned()),
            NoSign => Some(Zero::zero()),
            Minus => None,
        }
    }

    #[inline]
    pub fn checked_add(&self, v: &BigIntSmall) -> Option<BigIntSmall> {
        Some(self + v)
    }

    #[inline]
    pub fn checked_sub(&self, v: &BigIntSmall) -> Option<BigIntSmall> {
        Some(self - v)
    }

    #[inline]
    pub fn checked_mul(&self, v: &BigIntSmall) -> Option<BigIntSmall> {
        Some(self * v)
    }

    #[inline]
    pub fn checked_div(&self, v: &BigIntSmall) -> Option<BigIntSmall> {
        if v.is_zero() {
            return None;
        }
        Some(self / v)
    }

    /// Returns `self ^ exponent`.
    pub fn pow(&self, exponent: u32) -> Self {
        Pow::pow(self, exponent)
    }

    /// Returns `(self ^ exponent) mod modulus`
    ///
    /// Note that this rounds like `mod_floor`, not like the `%` operator,
    /// which makes a difference when given a negative `self` or `modulus`.
    /// The result will be in the interval `[0, modulus)` for `modulus > 0`,
    /// or in the interval `(modulus, 0]` for `modulus < 0`
    ///
    /// Panics if the exponent is negative or the modulus is zero.
    pub fn modpow(&self, exponent: &Self, modulus: &Self) -> Self {
        power::modpow(self, exponent, modulus)
    }

    /// Returns the truncated principal square root of `self` --
    /// see [Roots::sqrt](https://docs.rs/num-integer/0.1/num_integer/trait.Roots.html#method.sqrt).
    pub fn sqrt(&self) -> Self {
        Roots::sqrt(self)
    }

    /// Returns the truncated principal cube root of `self` --
    /// see [Roots::cbrt](https://docs.rs/num-integer/0.1/num_integer/trait.Roots.html#method.cbrt).
    pub fn cbrt(&self) -> Self {
        Roots::cbrt(self)
    }

    /// Returns the truncated principal `n`th root of `self` --
    /// See [Roots::nth_root](https://docs.rs/num-integer/0.1/num_integer/trait.Roots.html#tymethod.nth_root).
    pub fn nth_root(&self, n: u32) -> Self {
        Roots::nth_root(self, n)
    }

    /// Returns the number of least-significant bits that are zero,
    /// or `None` if the entire number is zero.
    pub fn trailing_zeros(&self) -> Option<u64> {
        self.data().trailing_zeros()
    }

    /// Returns whether the bit in position `bit` is set,
    /// using the two's complement for negative numbers
    pub fn bit(&self, bit: u64) -> bool {
        if self.is_negative() {
            // Let the binary representation of a number be
            //   ... 0  x 1 0 ... 0
            // Then the two's complement is
            //   ... 1 !x 1 0 ... 0
            // where !x is obtained from x by flipping each bit
            if bit >= u64::from(crate::big_digit::BITS) * self.len() as u64 {
                true
            } else {
                let trailing_zeros = self.data().trailing_zeros().unwrap();
                match Ord::cmp(&bit, &trailing_zeros) {
                    Ordering::Less => false,
                    Ordering::Equal => true,
                    Ordering::Greater => !self.data().bit(bit),
                }
            }
        } else {
            self.data().bit(bit)
        }
    }

    /// Sets or clears the bit in the given position,
    /// using the two's complement for negative numbers
    ///
    /// Note that setting/clearing a bit (for positive/negative numbers,
    /// respectively) greater than the current bit length, a reallocation
    /// may be needed to store the new digits
    pub fn set_bit(&mut self, bit: u64, value: bool) {
        match self.sign() {
            Sign::Plus => {
                let uint = self.biguint_mut();
                uint.set_bit(bit, value)
            }
            Sign::Minus => bits::set_negative_bit(self, bit, value),
            Sign::NoSign => {
                if value {
                    let uint = self.biguint_mut();
                    uint.set_bit(bit, true);
                    self.set_sign(Sign::Plus);
                } else {
                    // Clearing a bit for zero is a no-op
                }
            }
        }
        // The top bit may have been cleared, so normalize
        self.normalize();
    }
}
