use super::BigUint;
use num_integer::Roots;
use num_traits::{One, ToPrimitive, Zero};

// sqrt, cbrt, and nth_root use Newton's method to compute
// principal root of a given degree for a given integer.
//
// Reference:
// Brent & Zimmermann, Modern Computer Arithmetic, v0.5.9, Algorithms 1.13 and 1.14

impl BigUint {
    /// Returns the truncated principal square root of `self` --
    /// see [`Roots::sqrt()`].
    pub fn sqrt(&self) -> Self {
        if self.is_zero() || self.is_one() {
            return self.clone();
        }

        // If we fit in `u64`, compute the root that way.
        if let Some(x) = self.to_u64() {
            return x.sqrt().into();
        }

        let bits = self.bits();
        let max_bits = bits / 2 + 1;

        #[cfg(feature = "std")]
        let guess = match self.to_f64() {
            Some(f) if f.is_finite() => {
                use num_traits::FromPrimitive;

                // We fit in `f64` (lossy), so get a better initial guess from that.
                BigUint::from_f64(f.sqrt()).unwrap()
            }
            _ => {
                // Try to guess by scaling down such that it does fit in `f64`.
                // With some (x * 2²ᵏ), its sqrt ≈ (√x * 2ᵏ)
                let extra_bits = bits - (f64::MAX_EXP as u64 - 1);
                let root_scale = (extra_bits + 1) / 2;
                let scale = root_scale * 2;
                (self >> scale).sqrt() << root_scale
            }
        };

        #[cfg(not(feature = "std"))]
        let guess = BigUint::ONE << max_bits;

        fixpoint(guess, max_bits, move |s| {
            let q = self / s;
            let t = s + q;
            t >> 1
        })
    }

    /// Returns the truncated principal cube root of `self` --
    /// see [`Roots::cbrt()`].
    pub fn cbrt(&self) -> Self {
        if self.is_zero() || self.is_one() {
            return self.clone();
        }

        // If we fit in `u64`, compute the root that way.
        if let Some(x) = self.to_u64() {
            return x.cbrt().into();
        }

        let bits = self.bits();
        let max_bits = bits / 3 + 1;

        #[cfg(feature = "std")]
        let guess = match self.to_f64() {
            Some(f) if f.is_finite() => {
                use num_traits::FromPrimitive;

                // We fit in `f64` (lossy), so get a better initial guess from that.
                BigUint::from_f64(f.cbrt()).unwrap()
            }
            _ => {
                // Try to guess by scaling down such that it does fit in `f64`.
                // With some (x * 2³ᵏ), its cbrt ≈ (∛x * 2ᵏ)
                let extra_bits = bits - (f64::MAX_EXP as u64 - 1);
                let root_scale = (extra_bits + 2) / 3;
                let scale = root_scale * 3;
                (self >> scale).cbrt() << root_scale
            }
        };

        #[cfg(not(feature = "std"))]
        let guess = BigUint::ONE << max_bits;

        fixpoint(guess, max_bits, move |s| {
            let q = self / (s * s);
            let t = (s << 1) + q;
            t / 3u32
        })
    }

    /// Returns the truncated principal `n`th root of `self` --
    /// see [`Roots::nth_root()`].
    pub fn nth_root(&self, n: u32) -> Self {
        assert!(n > 0, "root degree n must be at least 1");

        if self.is_zero() || self.is_one() {
            return self.clone();
        }

        match n {
            // Optimize for small n
            1 => return self.clone(),
            2 => return self.sqrt(),
            3 => return self.cbrt(),
            _ => (),
        }

        // The root of non-zero values less than 2ⁿ can only be 1.
        let bits = self.bits();
        let n64 = u64::from(n);
        if bits <= n64 {
            return BigUint::ONE;
        }

        // If we fit in `u64`, compute the root that way.
        if let Some(x) = self.to_u64() {
            return x.nth_root(n).into();
        }

        let max_bits = bits / n64 + 1;

        #[cfg(feature = "std")]
        let guess = match self.to_f64() {
            Some(f) if f.is_finite() => {
                use num_traits::FromPrimitive;

                // We fit in `f64` (lossy), so get a better initial guess from that.
                BigUint::from_f64((f.ln() / f64::from(n)).exp()).unwrap()
            }
            _ => {
                use num_integer::Integer;

                // Try to guess by scaling down such that it does fit in `f64`.
                // With some (x * 2ⁿᵏ), its nth root ≈ (ⁿ√x * 2ᵏ)

                let extra_bits = bits - (f64::MAX_EXP as u64 - 1);
                let root_scale = Integer::div_ceil(&extra_bits, &n64);
                let scale = root_scale * n64;
                if scale < bits && bits - scale > n64 {
                    (self >> scale).nth_root(n) << root_scale
                } else {
                    BigUint::ONE << max_bits
                }
            }
        };

        #[cfg(not(feature = "std"))]
        let guess = BigUint::ONE << max_bits;

        let n_min_1 = n - 1;
        fixpoint(guess, max_bits, move |s| {
            let q = self / s.pow(n_min_1);
            let t = n_min_1 * s + q;
            t / n
        })
    }
}

#[inline]
fn fixpoint<F>(mut x: BigUint, max_bits: u64, f: F) -> BigUint
where
    F: Fn(&BigUint) -> BigUint,
{
    let mut xn = f(&x);

    // If the value increased, then the initial guess must have been low.
    // Repeat until we reverse course.
    while x < xn {
        // Sometimes an increase will go way too far, especially with large
        // powers, and then take a long time to walk back.  We know an upper
        // bound based on bit size, so saturate on that.
        x = if xn.bits() > max_bits {
            BigUint::ONE << max_bits
        } else {
            xn
        };
        xn = f(&x);
    }

    // Now keep repeating while the estimate is decreasing.
    while x > xn {
        x = xn;
        xn = f(&x);
    }
    x
}

#[deny(unconditional_recursion)]
impl Roots for BigUint {
    #[inline]
    fn sqrt(&self) -> Self {
        self.sqrt()
    }

    #[inline]
    fn cbrt(&self) -> Self {
        self.cbrt()
    }

    #[inline]
    fn nth_root(&self, n: u32) -> Self {
        self.nth_root(n)
    }
}
