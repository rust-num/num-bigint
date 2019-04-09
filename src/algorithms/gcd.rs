use std::borrow::Cow;

use integer::Integer;
use num_traits::Zero;

use crate::big_digit::{BigDigit, DoubleBigDigit, BITS};
use crate::bigint::Sign::*;
use crate::bigint::{BigInt, ToBigInt};
use crate::biguint::{BigUint, IntDigits};

/// Uses the lehemer algorithm.
/// Based on https://github.com/golang/go/blob/master/src/math/big/int.go#L612
/// If `extended` is set, the Bezout coefficients are calculated, otherwise they are `None`.
pub fn extended_gcd(
    a_in: Cow<BigUint>,
    b_in: Cow<BigUint>,
    extended: bool,
) -> (BigInt, Option<BigInt>, Option<BigInt>) {
    if a_in.is_zero() && b_in.is_zero() {
        if extended {
            return (b_in.to_bigint().unwrap(), Some(0.into()), Some(0.into()));
        } else {
            return (b_in.to_bigint().unwrap(), None, None);
        }
    }

    if a_in.is_zero() {
        if extended {
            return (b_in.to_bigint().unwrap(), Some(0.into()), Some(1.into()));
        } else {
            return (b_in.to_bigint().unwrap(), None, None);
        }
    }

    if b_in.is_zero() {
        if extended {
            return (a_in.to_bigint().unwrap(), Some(1.into()), Some(0.into()));
        } else {
            return (a_in.to_bigint().unwrap(), None, None);
        }
    }

    let a_in = a_in.to_bigint().unwrap();
    let b_in = b_in.to_bigint().unwrap();

    let mut a = a_in.clone();
    let mut b = b_in.clone();

    // `ua` (`ub`) tracks how many times input `a_in` has beeen accumulated into `a` (`b`).
    let mut ua = if extended { Some(1.into()) } else { None };
    let mut ub = if extended { Some(0.into()) } else { None };

    // Ensure that a >= b
    if a < b {
        std::mem::swap(&mut a, &mut b);
        std::mem::swap(&mut ua, &mut ub);
    }

    let mut q: BigInt = 0.into();
    let mut r: BigInt = 0.into();
    let mut s: BigInt = 0.into();
    let mut t: BigInt = 0.into();

    while b.len() > 1 {
        // Attempt to calculate in single-precision using leading words of a and b.
        let (u0, u1, v0, v1, even) = lehmer_simulate(&a, &b);
        // multiprecision step
        if v0 != 0 {
            // Simulate the effect of the single-precision steps using cosequences.
            // a = u0 * a + v0 * b
            // b = u1 * a + v1 * b
            lehmer_update(
                &mut a, &mut b, &mut q, &mut r, &mut s, &mut t, u0, u1, v0, v1, even,
            );

            if extended {
                // ua = u0 * ua + v0 * ub
                // ub = u1 * ua + v1 * ub
                lehmer_update(
                    ua.as_mut().unwrap(),
                    ub.as_mut().unwrap(),
                    &mut q,
                    &mut r,
                    &mut s,
                    &mut t,
                    u0,
                    u1,
                    v0,
                    v1,
                    even,
                );
            }
        } else {
            // Single-digit calculations failed to simulate any quotients.
            euclid_udpate(
                &mut a, &mut b, &mut ua, &mut ub, &mut q, &mut r, &mut s, &mut t, extended,
            );
        }
    }

    if b.len() > 0 {
        // base case if b is a single digit
        if a.len() > 1 {
            // a is longer than a single word, so one update is needed
            euclid_udpate(
                &mut a, &mut b, &mut ua, &mut ub, &mut q, &mut r, &mut s, &mut t, extended,
            );
        }

        if b.len() > 0 {
            // a and b are both single word
            let mut a_word = a.digits()[0];
            let mut b_word = b.digits()[0];

            if extended {
                let mut ua_word: BigDigit = 1;
                let mut ub_word: BigDigit = 0;
                let mut va: BigDigit = 0;
                let mut vb: BigDigit = 1;
                let mut even = true;

                while b_word != 0 {
                    let q = a_word / b_word;
                    let r = a_word % b_word;
                    a_word = b_word;
                    b_word = r;

                    let k = ua_word.wrapping_add(q.wrapping_mul(ub_word));
                    ua_word = ub_word;
                    ub_word = k;

                    let k = va.wrapping_add(q.wrapping_mul(vb));
                    va = vb;
                    vb = k;
                    even = !even;
                }

                t.data.set_digit(ua_word);
                s.data.set_digit(va);
                
                println!("extended gcd: t pre = {:?}", t);
                println!("extended gcd: s pre = {:?}", s);

                t.sign = if even { Plus } else { Minus };
                s.sign = if even { Minus } else { Plus };

                if let Some(ua) = ua.as_mut() {
                    println!("extended gcd: enter ua = {:?}", ua);
                    println!("extended gcd: enter ua = t pre = {:?}", t);
                    println!("extended gcd: enter ua = s pre = {:?}", s);
                    t *= &*ua;
                    println!("extended gcd: enter ua = t = {:?}", t);
                    s *= ub.unwrap();
                    println!("extended gcd: enter ua = s = {:?}", s);

                    *ua = &t + &s;
                    println!("extended gcd: ua = {:?}", ua);
                }
            } else {
                while b_word != 0 {
                    let quotient = a_word % b_word;
                    a_word = b_word;
                    b_word = quotient;
                }
            }
            a.digits_mut()[0] = a_word;
        }
    }

    a.normalize();

    let y = if let Some(ref ua) = ua {
        // y = (z - a * x) / b
        Some((&a - (&a_in * ua)) / &b_in)
    } else {
        None
    };
    //println!("exteneded gcd: {:?}", (a, ua, y));
    (a, ua, y)
}

/// Attempts to simulate several Euclidean update steps using leading digits of `a` and `b`.
/// It returns `u0`, `u1`, `v0`, `v1` such that `a` and `b` can be updated as:
///     a = u0 * a + v0 * b
///     b = u1 * a + v1 * b
///
/// Requirements: `a >= b` and `b.len() > 1`.
/// Since we are calculating with full words to avoid overflow, `even` (the returned bool)
/// is used to track the sign of cosequences.
/// For even iterations: `u0, v1 >= 0 && u1, v0 <= 0`
/// For odd iterations: `u0, v1 <= && u1, v0 >= 0`
#[inline]
fn lehmer_simulate(a: &BigInt, b: &BigInt) -> (BigDigit, BigDigit, BigDigit, BigDigit, bool) {
    // m >= 2
    let m = b.len();
    // n >= m >= 2
    let n = a.len();

    debug_assert!(m >= 2);
    debug_assert!(n >= m);

    // extract the top word of bits from a and b
    let h = a.digits()[n - 1].leading_zeros();

    let mut a1: BigDigit = a.digits()[n - 1] << h
        | ((a.digits()[n - 2] as DoubleBigDigit) >> (BITS as u32 - h)) as BigDigit;

    // b may have implicit zero words in the high bits if the lengths differ
    let mut a2: BigDigit = if n == m {
        b.digits()[n - 1] << h
            | ((b.digits()[n - 2] as DoubleBigDigit) >> (BITS as u32 - h)) as BigDigit
    } else if n == m + 1 {
        ((b.digits()[n - 2] as DoubleBigDigit) >> (BITS as u32 - h)) as BigDigit
    } else {
        0
    };

    // odd, even tracking
    let mut even = false;

    let mut u0 = 0;
    let mut u1 = 1;
    let mut u2 = 0;

    let mut v0 = 0;
    let mut v1 = 0;
    let mut v2 = 1;

    // Calculate the quotient and cosequences using Collins' stoppting condition.
    while a2 >= v2 && a1.wrapping_sub(a2) >= v1 + v2 {
        let q = a1 / a2;
        let r = a1 % a2;

        a1 = a2;
        a2 = r;

        let k = u1 + q * u2;
        u0 = u1;
        u1 = u2;
        u2 = k;

        let k = v1 + q * v2;
        v0 = v1;
        v1 = v2;
        v2 = k;

        even = !even;
    }

    (u0, u1, v0, v1, even)
}

fn lehmer_update(
    a: &mut BigInt,
    b: &mut BigInt,
    q: &mut BigInt,
    r: &mut BigInt,
    s: &mut BigInt,
    t: &mut BigInt,
    u0: BigDigit,
    u1: BigDigit,
    v0: BigDigit,
    v1: BigDigit,
    even: bool,
) {
    t.data.set_digit(u0);
    s.data.set_digit(v0);
    if even {
        t.sign = Plus;
        s.sign = Minus
    } else {
        t.sign = Minus;
        s.sign = Plus;
    }

    *t *= &*a;
    *s *= &*b;

    r.data.set_digit(u1);
    q.data.set_digit(v1);
    if even {
        q.sign = Plus;
        r.sign = Minus
    } else {
        q.sign = Minus;
        r.sign = Plus;
    }

    *r *= &*a;
    *q *= &*b;

    *a = t + s;
    *b = r + q;
}

fn euclid_udpate(
    a: &mut BigInt,
    b: &mut BigInt,
    ua: &mut Option<BigInt>,
    ub: &mut Option<BigInt>,
    q: &mut BigInt,
    r: &mut BigInt,
    s: &mut BigInt,
    t: &mut BigInt,
    extended: bool,
) {
    let (q_new, r_new) = a.div_rem(b);
    *q = q_new;
    *r = r_new;

    std::mem::swap(a, b);
    std::mem::swap(b, r);

    if extended {
        // ua, ub = ub, ua - q * ub
        if let Some(ub) = ub.as_mut() {
            if let Some(ua) = ua.as_mut() {
                *t = ub.clone();
                *s = &*ub * &*q;
                *ub = &*ua - &*s;
                *ua = t.clone();
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use num_traits::FromPrimitive;

    #[cfg(feature = "rand")]
    use crate::bigrand::RandBigInt;
    #[cfg(feature = "rand")]
    use num_traits::{One, Zero};
    #[cfg(feature = "rand")]
    use rand::{SeedableRng, XorShiftRng};

    #[cfg(feature = "rand")]
    fn extended_gcd_euclid(a: Cow<BigUint>, b: Cow<BigUint>) -> (BigInt, BigInt, BigInt) {
        if a.is_zero() && b.is_zero() {
            return (0.into(), 0.into(), 0.into());
        }

        let (mut s, mut old_s) = (BigInt::zero(), BigInt::one());
        let (mut t, mut old_t) = (BigInt::one(), BigInt::zero());
        let (mut r, mut old_r) = (b.to_bigint().unwrap(), a.to_bigint().unwrap());

        while !r.is_zero() {
            let quotient = &old_r / &r;
            old_r = old_r - &quotient * &r;
            std::mem::swap(&mut old_r, &mut r);
            old_s = old_s - &quotient * &s;
            std::mem::swap(&mut old_s, &mut s);
            old_t = old_t - quotient * &t;
            std::mem::swap(&mut old_t, &mut t);
        }

        (old_r, old_s, old_t)
    }

    #[test]
    #[cfg(feature = "rand")]
    fn test_extended_gcd_assumptions() {
        let mut rng = XorShiftRng::from_seed([1u8; 16]);

        for i in 1usize..100 {
            for j in &[1usize, 64, 128] {
                println!("round {} - {}", i, j);
                let a = rng.gen_biguint(i * j);
                let b = rng.gen_biguint(i * j);
                let (q, s_k, t_k) = extended_gcd(Cow::Borrowed(&a), Cow::Borrowed(&b), true);

                let lhs = BigInt::from_biguint(Plus, a) * &s_k.unwrap();
                let rhs = BigInt::from_biguint(Plus, b) * &t_k.unwrap();

                assert_eq!(q.clone(), &lhs + &rhs, "{} = {} + {}", q, lhs, rhs);
            }
        }
    }

    #[test]
    fn test_extended_gcd_example() {
        // simple example for wikipedia
        let a = BigUint::from_u32(240).unwrap();
        let b = BigUint::from_u32(46).unwrap();
        let (q, s_k, t_k) = extended_gcd(Cow::Owned(a), Cow::Owned(b), true);

        assert_eq!(q, BigInt::from_i32(2).unwrap());
        assert_eq!(s_k.unwrap(), BigInt::from_i32(-9).unwrap());
        assert_eq!(t_k.unwrap(), BigInt::from_i32(47).unwrap());
    }

    #[test]
    fn test_extended_gcd_example_not_extended() {
        // simple example for wikipedia
        let a = BigUint::from_u32(240).unwrap();
        let b = BigUint::from_u32(46).unwrap();
        let (q, s_k, t_k) = extended_gcd(Cow::Owned(a), Cow::Owned(b), false);

        assert_eq!(q, BigInt::from_i32(2).unwrap());
        assert_eq!(s_k, None);
        assert_eq!(t_k, None);
    }

    #[test]
    #[cfg(feature = "rand")]
    fn test_gcd_lehmer_euclid_extended() {
        let mut rng = XorShiftRng::from_seed([1u8; 16]);

        for i in 1usize..80 {
            for j in &[1usize, 16, 24, 64, 128] {
                println!("round {} - {}", i, j);
                let a = rng.gen_biguint(i * j);
                let b = rng.gen_biguint(i * j);
                let (q, s_k, t_k) = extended_gcd(Cow::Borrowed(&a), Cow::Borrowed(&b), true);

                let expected = extended_gcd_euclid(Cow::Borrowed(&a), Cow::Borrowed(&b));
                assert_eq!(q, expected.0);
                assert_eq!(s_k.unwrap(), expected.1);
                assert_eq!(t_k.unwrap(), expected.2);
            }
        }
    }

    #[test]
    #[cfg(feature = "rand")]
    fn test_gcd_lehmer_euclid_not_extended() {
        let mut rng = XorShiftRng::from_seed([1u8; 16]);

        for i in 1usize..80 {
            for j in &[1usize, 16, 24, 64, 128] {
                println!("round {} - {}", i, j);
                let a = rng.gen_biguint(i * j);
                let b = rng.gen_biguint(i * j);
                let (q, s_k, t_k) = extended_gcd(Cow::Borrowed(&a), Cow::Borrowed(&b), false);

                let expected = extended_gcd_euclid(Cow::Borrowed(&a), Cow::Borrowed(&b));
                assert_eq!(
                    q, expected.0,
                    "gcd({}, {}) = {} != {}",
                    &a, &b, &q, expected.0
                );
                assert_eq!(s_k, None);
                assert_eq!(t_k, None);
            }
        }
    }
}
