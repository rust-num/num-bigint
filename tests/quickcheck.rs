extern crate num_bigint;
extern crate num_integer;
extern crate num_traits;

#[cfg(test)]
#[cfg(feature = "quickcheck")]
#[macro_use]
pub extern crate quickcheck;

#[cfg(test)]
#[cfg(feature = "quickcheck")]
mod quickchecks {
    use num_bigint::{BigInt, BigUint};
    use num_traits::{Num, One, Pow, ToPrimitive};
    use quickcheck::TestResult;

    quickcheck! {
        fn quickcheck_add_primitive(a: u128, b: u128) -> TestResult {
            let a = a;
            let b = b;
            let a_big = BigUint::from(a);
            let b_big = BigUint::from(b);
            //maximum value of u64 means no overflow
            match a.checked_add(b) {
                None => TestResult::discard(),
                Some(sum) => TestResult::from_bool(sum == (a_big + b_big).to_u128().unwrap())
            }
        }
    }

    quickcheck! {
        fn quickcheck_add_primitive_signed(a: i128, b: i128) -> TestResult {
            let a = a;
            let b = b;
            let a_big = BigInt::from(a);
            let b_big = BigInt::from(b);
            //maximum value of u64 means no overflow
            match a.checked_add(b) {
                None => TestResult::discard(),
                Some(sum) => TestResult::from_bool(sum == (a_big + b_big).to_i128().unwrap())
            }
        }
    }

    quickcheck! {
        fn quickcheck_add_commutative(a: BigUint, b: BigUint) -> bool {
            a.clone() + b.clone() == b + a
        }
    }

    quickcheck! {
        fn quickcheck_add_commutative_signed(a: BigInt, b: BigInt) -> bool {
            a.clone() + b.clone() == b + a
        }
    }

    quickcheck! {
        fn quickcheck_add_zero(a: BigUint) -> bool {
            a.clone() == a + BigUint::from(0_u32)
        }
    }

    quickcheck! {
        fn quickcheck_add_associative(a: BigUint, b: BigUint, c: BigUint) -> bool {
            (a.clone() + b.clone()) + c.clone() == a + (b + c)
        }
    }

    quickcheck! {
        fn quickcheck_mul_primitive(a: u64, b: u64) -> bool {
            let a = a as u128;
            let b = b as u128;
            let a_big = BigUint::from(a);
            let b_big = BigUint::from(b);
            //maximum value of u64 means no overflow
            a * b == (a_big * b_big).to_u128().unwrap()
        }
    }

    quickcheck! {
        fn quickcheck_mul_commutative(a: BigUint, b: BigUint) -> bool {
            a.clone() * b.clone() == b * a
        }
    }

    quickcheck! {
        fn quickcheck_mul_associative(a: BigUint, b: BigUint, c: BigUint) -> bool {
            (a.clone() * b.clone()) * c.clone() == a * (b * c)
        }
    }

    quickcheck! {
        fn quickcheck_distributive(a: BigUint, b: BigUint, c: BigUint) -> bool {
            a.clone() * (b.clone() + c.clone()) == a.clone() * b + a * c
        }
    }

    quickcheck! {
        ///Tests that exactly one of a<b a>b a=b is true
        fn quickcheck_ge_le_eq_mut_exclusive(a: BigUint, b: BigUint) -> bool {
            let gt_lt_eq = vec![a > b, a < b, a == b];
            gt_lt_eq
                .iter()
                .fold(0, |acc, e| if *e { acc + 1 } else { acc })
                == 1
    }
    }

    quickcheck! {
        /// Tests correctness of subtraction assuming addition is correct
        fn quickcheck_sub(a: BigUint, b: BigUint) -> bool {
            if b < a {
                a.clone() - b.clone() + b == a
            } else {
                b.clone() - a.clone() + a == b
            }
        }
    }

    quickcheck! {
        fn quickcheck_pow_zero(a: BigUint) -> bool {
            a.pow(0_u32) == BigUint::one()
        }
    }

    quickcheck! {
        fn quickcheck_pow_one(a: BigUint) -> bool {
            a.pow(1_u32) == a
        }
    }

    quickcheck! {
        fn quickcheck_sqrt(a: BigUint) -> bool {
            (a.clone() * a.clone()).sqrt() == a
        }
    }

    quickcheck! {
        fn quickcheck_cbrt(a: BigUint) -> bool {
            (a.clone() * a.clone() * a.clone()).cbrt() == a
        }
    }

    quickcheck! {
        fn quickcheck_conversion(a:BigUint) -> bool {
            let mut success = true;
            for radix in 2..=36 {
                let string = a.to_str_radix(radix);
                success &= a == BigUint::from_str_radix(&string, radix).unwrap()
            }
            success
        }
    }

    //this test takes too long (no surprise)
    //quickcheck automaticly stops if a test takes more than 60s
    // fn quickcheck_pow_and_root(a: BigUint, n: u8) -> TestResult {
    //     match n {
    //         0 => TestResult::discard(),
    //         n => TestResult::from_bool(a.clone().pow(n).nth_root(n as u32) == a),
    //     }
    // }

    quickcheck! {
        fn quickcheck_modpow(a:BigUint, exp:u8, modulus:u64) -> TestResult {
            if modulus == 0 {
                return TestResult::discard()
            }
            let expected = a.clone().pow(exp) % modulus;
            let result = a.modpow(&BigUint::from(exp), &BigUint::from(modulus));
            TestResult::from_bool(expected == result)
        }
    }
}
