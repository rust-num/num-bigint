use num_integer::Integer;
use num_traits::ToPrimitive;

use crate::BigUint;

impl BigUint {
    pub fn true_div(&self, other: &BigUint) -> f64 {
        // Compute exponent, biased by 52 (or 53)
        let mut exponent =
            (self.bits() as i32) - (other.bits() as i32) - (f64::MANTISSA_DIGITS as i32);

        // Actual division, where the rhs is shifted, so that
        // the quotient must have 52 or 53 bits in length
        let rhs = match exponent {
            i if i > 0 => self >> i,
            i if i < 0 => self << -i,
            _ => self.clone(),
        };
        let (mut q, r) = rhs.div_rem(other);

        // rounding
        if r > (other >> 1) {
            q += 1u32;
        }

        // Get the exact 53 bits of mantissa and actual exponent
        let extra_bits = q.bits() as u32 - f64::MANTISSA_DIGITS;
        q >>= extra_bits;
        exponent += (f64::MANTISSA_DIGITS + extra_bits - 1) as i32;

        // Handle overflow or underflow
        if exponent > core::f64::MAX_EXP {
            return f64::INFINITY;
        } else if exponent < core::f64::MIN_EXP {
            return 0f64;
        }

        // Get actual mantissa as float and apply exponent
        let mantissa_as_float = q.to_f64().unwrap() * 2.0f64.powi(1 - f64::MANTISSA_DIGITS as i32);
        mantissa_as_float * 2.0f64.powi(exponent)
    }
}
