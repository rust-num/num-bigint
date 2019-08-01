#![cfg(feature = "quickcheck")]

extern crate num_bigint;
extern crate num_integer;
extern crate num_traits;

extern crate quickcheck;
#[macro_use]
extern crate quickcheck_macros;

use num_bigint::{BigInt, BigUint};
use num_traits::{Num, One, Pow, Zero};
use quickcheck::TestResult;

#[quickcheck]
fn quickcheck_unsigned_eq_reflexive(a: BigUint) -> bool {
    a == a
}

#[quickcheck]
fn quickcheck_signed_eq_reflexive(a: BigInt) -> bool {
    a == a
}

#[quickcheck]
fn quickcheck_unsigned_eq_symmetric(a: BigUint, b: BigUint) -> bool {
    if a == b {
        b == a
    } else {
        b != a
    }
}

#[quickcheck]
fn quickcheck_signed_eq_symmetric(a: BigInt, b: BigInt) -> bool {
    if a == b {
        b == a
    } else {
        b != a
    }
}

#[quickcheck]
fn quickcheck_unsigned_add_primitive(a: u128, b: u128) -> TestResult {
    let actual = BigUint::from(a) + BigUint::from(b);
    match a.checked_add(b) {
        None => TestResult::discard(),
        Some(expected) => TestResult::from_bool(BigUint::from(expected) == actual),
    }
}

#[quickcheck]
fn quickcheck_signed_add_primitive(a: i128, b: i128) -> TestResult {
    let actual = BigInt::from(a) + BigInt::from(b);
    match a.checked_add(b) {
        None => TestResult::discard(),
        Some(expected) => TestResult::from_bool(BigInt::from(expected) == actual),
    }
}

#[quickcheck]
fn quickcheck_unsigned_add_commutative(a: BigUint, b: BigUint) -> bool {
    a.clone() + b.clone() == b + a
}

#[quickcheck]
fn quickcheck_signed_add_commutative(a: BigInt, b: BigInt) -> bool {
    a.clone() + b.clone() == b + a
}

#[quickcheck]
fn quickcheck_unsigned_add_zero(a: BigUint) -> bool {
    a.clone() == a + BigUint::zero()
}

#[quickcheck]
fn quickcheck_signed_add_zero(a: BigInt) -> bool {
    a.clone() == a + BigInt::zero()
}

#[quickcheck]
fn quickcheck_unsigned_add_associative(a: BigUint, b: BigUint, c: BigUint) -> bool {
    (a.clone() + b.clone()) + c.clone() == a + (b + c)
}

#[quickcheck]
fn quickcheck_signed_add_associative(a: BigInt, b: BigInt, c: BigInt) -> bool {
    (a.clone() + b.clone()) + c.clone() == a + (b + c)
}

#[quickcheck]
fn quickcheck_unsigned_mul_primitive(a: u64, b: u64) -> bool {
    //maximum value of u64 means no overflow
    BigUint::from(a as u128 * b as u128) == BigUint::from(a) * BigUint::from(b)
}

#[quickcheck]
fn quickcheck_signed_mul_primitive(a: i64, b: i64) -> bool {
    //maximum value of i64 means no overflow
    BigInt::from(a as i128 * b as i128) == BigInt::from(a) * BigInt::from(b)
}

#[quickcheck]
fn quickcheck_unsigned_mul_zero(a: BigUint) -> bool {
    a * BigUint::zero() == BigUint::zero()
}

#[quickcheck]
fn quickcheck_signed_mul_zero(a: BigInt) -> bool {
    a * BigInt::zero() == BigInt::zero()
}

#[quickcheck]
fn quickcheck_unsigned_mul_one(a: BigUint) -> bool {
    a.clone() * BigUint::one() == a
}

#[quickcheck]
fn quickcheck_signed_mul_one(a: BigInt) -> bool {
    a.clone() * BigInt::one() == a
}

#[quickcheck]
fn quickcheck_unsigned_mul_commutative(a: BigUint, b: BigUint) -> bool {
    a.clone() * b.clone() == b * a
}

#[quickcheck]
fn quickcheck_signed_mul_commutative(a: BigInt, b: BigInt) -> bool {
    a.clone() * b.clone() == b * a
}

#[quickcheck]
fn quickcheck_unsigned_mul_associative(a: BigUint, b: BigUint, c: BigUint) -> bool {
    (a.clone() * b.clone()) * c.clone() == a * (b * c)
}

#[quickcheck]
fn quickcheck_signed_mul_associative(a: BigInt, b: BigInt, c: BigInt) -> bool {
    (a.clone() * b.clone()) * c.clone() == a * (b * c)
}

#[quickcheck]
fn quickcheck_unsigned_distributive(a: BigUint, b: BigUint, c: BigUint) -> bool {
    a.clone() * (b.clone() + c.clone()) == a.clone() * b + a * c
}

#[quickcheck]
fn quickcheck_signed_distributive(a: BigInt, b: BigInt, c: BigInt) -> bool {
    a.clone() * (b.clone() + c.clone()) == a.clone() * b + a * c
}

#[quickcheck]
///Tests that exactly one of a<b a>b a=b is true
fn quickcheck_unsigned_ge_le_eq_mut_exclusive(a: BigUint, b: BigUint) -> bool {
    let gt_lt_eq = vec![a > b, a < b, a == b];
    gt_lt_eq
        .iter()
        .fold(0, |acc, e| if *e { acc + 1 } else { acc })
        == 1
}

#[quickcheck]
///Tests that exactly one of a<b a>b a=b is true
fn quickcheck_signed_ge_le_eq_mut_exclusive(a: BigInt, b: BigInt) -> bool {
    let gt_lt_eq = vec![a > b, a < b, a == b];
    gt_lt_eq
        .iter()
        .fold(0, |acc, e| if *e { acc + 1 } else { acc })
        == 1
}

#[quickcheck]
fn quickcheck_unsigned_sub_primitive(a: u128, b: u128) -> bool {
    if b < a {
        BigUint::from(a - b) == BigUint::from(a) - BigUint::from(b)
    } else {
        BigUint::from(b - a) == BigUint::from(b) - BigUint::from(a)
    }
}

#[quickcheck]
fn quickcheck_signed_sub_primitive(a: i128, b: i128) -> bool {
    if b < a {
        BigInt::from(a - b) == BigInt::from(a) - BigInt::from(b)
    } else {
        BigInt::from(b - a) == BigInt::from(b) - BigInt::from(a)
    }
}

#[quickcheck]
/// Tests correctness of subtraction assuming addition is correct
fn quickcheck_unsigned_sub(a: BigUint, b: BigUint) -> bool {
    if b < a {
        a.clone() - b.clone() + b == a
    } else {
        b.clone() - a.clone() + a == b
    }
}

#[quickcheck]
/// Tests correctness of subtraction assuming addition is correct
fn quickcheck_signed_sub(a: BigInt, b: BigInt) -> bool {
    if b < a {
        a.clone() - b.clone() + b == a
    } else {
        b.clone() - a.clone() + a == b
    }
}

#[quickcheck]
fn quickcheck_unsigned_div_primitive(a: u128, b: u128) -> TestResult {
    if b == 0 {
        TestResult::discard()
    } else {
        TestResult::from_bool(BigUint::from(a / b) == BigUint::from(a) / BigUint::from(b))
    }
}

#[quickcheck]
fn quickcheck_signed_div_primitive(a: i128, b: i128) -> TestResult {
    if b == 0 {
        TestResult::discard()
    } else {
        TestResult::from_bool(BigInt::from(a / b) == BigInt::from(a) / BigInt::from(b))
    }
}

#[quickcheck]
fn quickcheck_unsigned_pow_zero(a: BigUint) -> bool {
    a.pow(0_u32) == BigUint::one()
}

#[quickcheck]
fn quickcheck_unsigned_pow_one(a: BigUint) -> bool {
    a.pow(1_u32) == a
}

#[quickcheck]
fn quickcheck_unsigned_sqrt(a: BigUint) -> bool {
    (a.clone() * a.clone()).sqrt() == a
}

#[quickcheck]
fn quickcheck_unsigned_cbrt(a: BigUint) -> bool {
    (a.clone() * a.clone() * a.clone()).cbrt() == a
}

#[quickcheck]
fn quickcheck_signed_cbrt(a: BigInt) -> bool {
    (a.clone() * a.clone() * a.clone()).cbrt() == a
}

#[quickcheck]
fn quickcheck_unsigned_conversion(a: BigUint) -> bool {
    let mut success = true;
    for radix in 2..=36 {
        let string = a.to_str_radix(radix);
        success &= a == BigUint::from_str_radix(&string, radix).unwrap()
    }
    success
}

#[quickcheck]
fn quickcheck_signed_conversion(a: BigInt) -> bool {
    let mut success = true;
    for radix in 2..=36 {
        let string = a.to_str_radix(radix);
        success &= a == BigInt::from_str_radix(&string, radix).unwrap()
    }
    success
}

#[quickcheck]
fn quickcheck_unsigned_modpow(a: BigUint, exp: u8, modulus: BigUint) -> TestResult {
    if modulus.clone().is_zero() {
        TestResult::discard()
    } else {
        let expected = a.clone().pow(exp) % modulus.clone();
        let result = a.modpow(&BigUint::from(exp), &BigUint::from(modulus));
        TestResult::from_bool(expected == result)
    }
}
