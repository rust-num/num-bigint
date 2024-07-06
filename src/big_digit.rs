mod inline;
mod vec;

#[cfg(feature = "inline")]
pub(crate) use inline::BigDigits;

#[cfg(not(feature = "inline"))]
pub(crate) use vec::BigDigits;

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
