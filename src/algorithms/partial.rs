// Partial Euclidean algorithm.
// Lehmer's version for computing GCD
// Ported from Flint-2.4.1, fmpz_xgcd_partial

use crate::big_digit::BITS as LIMB_BITS;
use crate::bigint::{BigInt, Sign};
use bigint::ToBigInt;
use integer::Integer;
use num_traits::{One, Pow, ToPrimitive, Zero};
use std::borrow::Cow;
use std::ops::Neg;

#[inline]
fn signed_shift(op: u64, shift: isize) -> u64 {
    if shift > 0 {
        return op << shift;
    }

    if shift <= -64 {
        return 0;
    }

    return op >> (-shift);
}

//Port to rust from c++
//https://github.com/Chia-Network/vdftrack1results/blob/b171a420ecfffb2d956979da7b1b183c4c88d2a0/akashnil/entry/vdf4.cpp#L82
// Return an approximation x of the large mpz_t op by an int64_t and the exponent e adjustment.
// We must have (x * 2^e) / op = constant approximately.
#[inline]
pub fn partial_bigint(op: &BigInt) -> (i64, isize) {
    //uint64_t size = mpz_size(op);
    let size = op.bits();
    //uint64_t last = mpz_getlimbn(op, size - 1);
    let last: f64 = (size % LIMB_BITS) as f64;
    //uint64_t ret;
    let mut ret: i64;

    let mut lg2: f64 = last.log2() + 1f64;

    let mut exp: isize = lg2.clone() as isize;

    ret = signed_shift(last as u64, 63isize - exp) as i64;

    if size > 1 {
        exp = exp + (size as isize - 1) * 64;
        // uint64_t prev = mpz_getlimbn(op, size - 2);
        let mut prev: u64 = (size % LIMB_BITS) as u64;
        ret += signed_shift(prev, (-1f64 - lg2) as isize) as i64;
    }

    //if (mpz_sgn(op) < 0) return - ((int64_t)ret);
    if op.sign() == Sign::Minus {
        return (-ret, exp);
    }

    return (ret, exp);
}

/// This function is an implementation of Lehmer extended GCD with early termination.
/// It terminates early when remainders fall below the specified bound.
/// The initial values r1 and r2 are treated as successive remainders in the Euclidean algorithm
/// and are replaced with the last two remainders computed. The values co1 and co2 are the last two
/// cofactors and satisfy the identity co2*r1 - co1*r2 == +/- r2_orig upon termination, where
/// r2_orig is the starting value of r2 supplied, and r1 and r2 are the final values.
pub fn partial_extended_gcd(
    r2_in: Cow<BigInt>,
    r1_in: Cow<BigInt>,
    bound: Cow<BigInt>,
) -> (BigInt, BigInt, BigInt, BigInt) {
    //returb (R2, R1, C1, C2)

    //Temp Computing helper variables
    let mut a2: isize = 0;
    let mut a1: isize = 0;
    let mut b2: isize = 0;
    let mut b1: isize = 0;
    let mut t: isize = 0;
    let mut rr2: isize = 0;
    let mut rr1: isize = 0;
    let mut qq: isize = 0;
    let mut bb: isize = 0;

    //
    let mut q: BigInt = Zero::zero();
    let mut r: BigInt = Zero::zero();

    let mut co1 = BigInt::one();
    co1.sign = co1.sign.neg();
    let mut co2 = BigInt::zero();

    let r1_in = r1_in.to_bigint().unwrap();
    let r2_in = r2_in.to_bigint().unwrap();

    let mut r1 = r1_in.clone();
    let mut r2 = r2_in.clone();

    //loop index
    let mut index = 0;

    while !r1.is_zero() && r1 > *bound {
        //get bits length
        let r2_bits = r2.bits();
        let r1_bits = r1.bits();
        let one_under_limb = LIMB_BITS - 1;

        let mut T = (r2_bits as isize) - (one_under_limb as isize);
        let mut T1 = (r1_bits as isize) - (one_under_limb as isize);

        //Bits
        if T < T1 {
            T = T1
        }
        if T < 0 {
            T = 0
        }

        t = T;

        //truncate for a positive number is same as floor
        //truncate for a negative number is same as ceil

        //r = R2 / (2 ^ T); truncate r
        let (d, m) = r2.div_mod_floor(&BigInt::from(2).pow(t as usize));

        //positive sign or no sign
        if d.sign() == Sign::Minus {
            //negative numbers we all to ceil

            if m.is_zero() {
                r = d;
            } else {
                r = d + BigInt::one();
            }
        } else {
            r = d;
        }

        rr2 = r.to_isize().unwrap();

        //r = R1 / (2 ^ T); truncate r
        let (d, m) = r1.div_mod_floor(&BigInt::from(2).pow(t as usize));

        //positive sign or no sign
        if d.sign() == Sign::Minus {
            //negative numbers we all to ceil
            println!("suckerssss 2");
            if m.is_zero() {
                r = d;
            } else {
                r = d + BigInt::one();
            }
        } else {
            r = d;
        }

        //positive sign or no sign

        rr1 = r.to_isize().unwrap();

        //r = R1 / (2 ^ T); truncate r
        r = bound.div_floor(&BigInt::from(2).pow(t as usize));

        //r = bound / (2 ^ T); truncate r
        let (d, m) = bound.div_mod_floor(&BigInt::from(2).pow(t as usize));

        //positive sign or no sign
        if d.sign() == Sign::Minus {
            //negative numbers we all to ceil

            if m.is_zero() {
                r = d;
            } else {
                r = d + BigInt::one();
            }
        } else {
            r = d;
        }

        bb = r.to_isize().unwrap(); //might need tobe isize

        //reset values
        a1 = 1;
        a2 = 0;
        b1 = 0;
        b2 = 1;

        //reset inner loop index
        index = 0;

        // Euclidean Step
        while rr1 != 0 && rr1 > bb {
            qq = rr2 / rr1;

            //t1
            t = rr2 - qq * rr1;
            rr2 = rr1;
            rr1 = t;

            //t2
            t = a2 - qq * a1;
            a2 = a1;
            a1 = t;

            //t3
            t = b2 - qq * b1;
            b2 = b1;
            b1 = t;

            //check if it is even or odd
            if index % 2 != 0 {
                //index & 1
                //its odd
                if rr1 < -b1 || rr2 - rr1 < a1 - a2 {
                    break;
                }
            } else {
                //its even
                if rr1 < -a1 || rr2 - rr1 < b1 - b2 {
                    break;
                }
            }

            //increment counter
            index += 1;
        }

        if index == 0 {
            // multiprecision step
            q = r2.div_floor(&r1);
            r2 = &r2 % &r1;
            std::mem::swap(&mut r2, &mut r1);
            co2 = &co2 - (&q * &co1);
            std::mem::swap(&mut co2, &mut co1);
        } else {
            // recombination
            r = &r2 * &b2;

            if a2 >= 0 {
                r = &r + &r1 * &a2;
            } else {
                r = &r - (&r1 * &-a2);
            }

            r1 = &r1 * &a1;
            if b1 >= 0 {
                r1 = &r1 + &r2 * &b1;
            } else {
                r1 = &r1 - (&r2 * &-b1);
            }

            r2 = r.clone();

            r = &co2 * &b2;

            if a2 >= 0 {
                r = &r + &co1 * &a2;
            } else {
                r = &r - (&co1 * &-a2);
            }

            co1 = &co1 * &a1;
            if b1 >= 0 {
                co1 = &co1 + &co2 * &b1;
            } else {
                co1 = &co1 - (&co2 * &-b1);
            }

            // C2 = r;
            co2 = r.clone();

            // make sure R1 is positive
            if r1.sign() == Sign::Minus {
                co1.sign = co1.sign.neg();
                r1.sign = r1.sign.neg();
            }
            // make sure R2 is positive
            if r2.sign() == Sign::Minus {
                co2.sign = co2.sign.neg();
                r2.sign = r2.sign.neg();
            }
        }
    }

    // make sure R2 is positive
    if r2.sign() == Sign::Minus {
        co1.negate_sign();
        co2.negate_sign();
        r2.negate_sign();
    }

    //return back
    (co2, co1, r2, r1)
}

#[cfg(test)]
mod test {

    #[cfg(feature = "rand")]
    use crate::bigrand::RandBigInt;
    #[cfg(feature = "rand")]
    use num_traits::Signed;
    #[cfg(feature = "rand")]
    use rand::{SeedableRng, XorShiftRng};

    #[test]
    #[cfg(feature = "rand")]
    fn test_partial_extended_gcd() {
        use super::*;

        let mut rng = XorShiftRng::from_seed([1u8; 16]);

        /* Test co2*r1 - co1*r2 = r2_orig */
        let mut co1 = BigInt::zero();
        let mut co2 = BigInt::zero();
        let mut f = BigInt::zero();
        let mut g = rng.gen_bigint(2000);
        //let mut g = BigInt::zero();
        let mut t1 = BigInt::zero();
        let mut t2 = BigInt::zero();
        let mut L = BigInt::zero();

        g += BigInt::one();
        f = BigInt::from_biguint(Sign::Plus, rng.gen_biguint_below(&g.to_biguint().unwrap()));
        L = rng.gen_bigint(1000);
        //println!("L: {:?}", L);

        t2 = g.clone();
        t2 = t2.abs();

        let (co2, co1, r2, r1) =
            partial_extended_gcd(Cow::Borrowed(&g), Cow::Borrowed(&f), Cow::Borrowed(&L));

        t1 = &co2 * &r1;
        t1 -= &co1 * &r2;
        t1 = t1.abs();
        // println!("------------------test_partial_extended_gcd");
        // println!("t1: {:?}", t1);
        // println!("t2: {:?}", t2);

        assert_eq!(&t1, &t2);
    }

}
