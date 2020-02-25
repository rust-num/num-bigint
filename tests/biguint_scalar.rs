use num_bigint::BigUint;
use num_traits::{ToPrimitive, Zero};

mod consts;
use crate::consts::*;

#[macro_use]
mod macros;

#[test]
fn test_scalar_add() {
    fn check(x: &BigUint, y: &BigUint, z: &BigUint) {
        let (x, y, z) = (x.clone(), y.clone(), z.clone());
        assert_unsigned_scalar_op!(x + y == z);
        assert_unsigned_scalar_assign_op!(x += y == z);
    }

    for elm in SUM_TRIPLES.iter() {
        let (a_vec, b_vec, c_vec) = *elm;
        let a = BigUint::from_slice(a_vec);
        let b = BigUint::from_slice(b_vec);
        let c = BigUint::from_slice(c_vec);

        check(&a, &b, &c);
        check(&b, &a, &c);
    }
}

#[test]
fn test_scalar_sub() {
    fn check(x: &BigUint, y: &BigUint, z: &BigUint) {
        let (x, y, z) = (x.clone(), y.clone(), z.clone());
        assert_unsigned_scalar_op!(x - y == z);
        assert_unsigned_scalar_assign_op!(x -= y == z);
    }

    for elm in SUM_TRIPLES.iter() {
        let (a_vec, b_vec, c_vec) = *elm;
        let a = BigUint::from_slice(a_vec);
        let b = BigUint::from_slice(b_vec);
        let c = BigUint::from_slice(c_vec);

        check(&c, &a, &b);
        check(&c, &b, &a);
    }
}

#[test]
fn test_scalar_mul() {
    fn check(x: &BigUint, y: &BigUint, z: &BigUint) {
        let (x, y, z) = (x.clone(), y.clone(), z.clone());
        assert_unsigned_scalar_op!(x * y == z);
        assert_unsigned_scalar_assign_op!(x *= y == z);
    }

    for elm in MUL_TRIPLES.iter() {
        let (a_vec, b_vec, c_vec) = *elm;
        let a = BigUint::from_slice(a_vec);
        let b = BigUint::from_slice(b_vec);
        let c = BigUint::from_slice(c_vec);

        check(&a, &b, &c);
        check(&b, &a, &c);
    }
}

#[test]
fn test_scalar_rem_noncommutative() {
    assert_eq!(5u8 % BigUint::from(7u8), 5u8);
    assert_eq!(BigUint::from(5u8) % 7u8, 5u8);
}

#[test]
fn test_scalar_div_rem() {
    fn check(x: &BigUint, y: &BigUint, z: &BigUint, r: &BigUint) {
        let (x, y, z, r) = (x.clone(), y.clone(), z.clone(), r.clone());
        assert_unsigned_scalar_op!(x / y == z);
        assert_unsigned_scalar_op!(x % y == r);
        assert_unsigned_scalar_assign_op!(x /= y == z);
        assert_unsigned_scalar_assign_op!(x %= y == r);
    }

    for elm in MUL_TRIPLES.iter() {
        let (a_vec, b_vec, c_vec) = *elm;
        let a = BigUint::from_slice(a_vec);
        let b = BigUint::from_slice(b_vec);
        let c = BigUint::from_slice(c_vec);

        if !a.is_zero() {
            check(&c, &a, &b, &Zero::zero());
        }

        if !b.is_zero() {
            check(&c, &b, &a, &Zero::zero());
        }
    }

    for elm in DIV_REM_QUADRUPLES.iter() {
        let (a_vec, b_vec, c_vec, d_vec) = *elm;
        let a = BigUint::from_slice(a_vec);
        let b = BigUint::from_slice(b_vec);
        let c = BigUint::from_slice(c_vec);
        let d = BigUint::from_slice(d_vec);

        if !b.is_zero() {
            check(&a, &b, &c, &d);
            assert_unsigned_scalar_op!(a / b == c);
            assert_unsigned_scalar_op!(a % b == d);
            assert_unsigned_scalar_assign_op!(a /= b == c);
            assert_unsigned_scalar_assign_op!(a %= b == d);
        }
    }
}

#[test]
fn test_biguint_scalar_cmp() {
    use std::cmp::Ordering::*;

    let zero = BigUint::from(0u32);
    let one = BigUint::from(1u32);
    let five = BigUint::from(5u32);

    fn scalar_cmp_asserts(num: &BigUint, scalar: i32) {
        assert_eq!(num.partial_cmp(&(scalar as u8)), Some(Equal));
        assert_eq!((scalar as u8).partial_cmp(num), Some(Equal));
        assert_eq!((scalar as u8 + 1).partial_cmp(num), Some(Greater));
        assert_eq!(num.partial_cmp(&(scalar as u8 + 1)), Some(Less));

        assert_eq!(num.partial_cmp(&(scalar as u16)), Some(Equal));
        assert_eq!((scalar as u16).partial_cmp(num), Some(Equal));
        assert_eq!((scalar as u16 + 1).partial_cmp(num), Some(Greater));
        assert_eq!(num.partial_cmp(&(scalar as u16 + 1)), Some(Less));

        assert_eq!(num.partial_cmp(&(scalar as u32)), Some(Equal));
        assert_eq!((scalar as u32).partial_cmp(num), Some(Equal));
        assert_eq!((scalar as u32 + 1).partial_cmp(num), Some(Greater));
        assert_eq!(num.partial_cmp(&(scalar as u32 + 1)), Some(Less));

        assert_eq!(num.partial_cmp(&(scalar as u64)), Some(Equal));
        assert_eq!((scalar as u64).partial_cmp(num), Some(Equal));
        assert_eq!((scalar as u64 + 1).partial_cmp(num), Some(Greater));
        assert_eq!(num.partial_cmp(&(scalar as u64 + 1)), Some(Less));

        assert_eq!(num.partial_cmp(&(scalar as u128)), Some(Equal));
        assert_eq!((scalar as u128).partial_cmp(num), Some(Equal));
        assert_eq!((scalar as u128 + 1).partial_cmp(num), Some(Greater));
        assert_eq!(num.partial_cmp(&(scalar as u128 + 1)), Some(Less));

        assert_eq!(num.partial_cmp(&(scalar as usize)), Some(Equal));
        assert_eq!((scalar as usize).partial_cmp(num), Some(Equal));
        assert_eq!((scalar as usize + 1).partial_cmp(num), Some(Greater));
        assert_eq!(num.partial_cmp(&(scalar as usize + 1)), Some(Less));
    }

    scalar_cmp_asserts(&zero, 0);
    scalar_cmp_asserts(&one, 1);
    scalar_cmp_asserts(&five, 5);

    let a = BigUint::from(10000000000u64);
    assert_eq!(a.partial_cmp(&10000000000u64), Some(Equal));
    assert_eq!(a.partial_cmp(&1000000000u64), Some(Greater));
    assert_eq!(a.partial_cmp(&100000000000u64), Some(Less));
}
