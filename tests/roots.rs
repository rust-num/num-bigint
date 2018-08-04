extern crate num_bigint;
extern crate num_integer;
extern crate num_traits;

mod biguint {
    use num_bigint::BigUint;
    use num_traits::Pow;
    use std::str::FromStr;

    fn check(x: u64, n: u32) {
        let big_x = BigUint::from(x);
        let res = big_x.nth_root(n);

        if n == 2 {
            assert_eq!(&res, &big_x.sqrt())
        } else if n == 3 {
            assert_eq!(&res, &big_x.cbrt())
        }

        assert!(res.pow(n) <= big_x);
        assert!((res + 1u32).pow(n) > big_x);
    }

    #[test]
    fn test_sqrt() {
        check(99, 2);
        check(100, 2);
        check(120, 2);
    }

    #[test]
    fn test_cbrt() {
        check(8, 3);
        check(26, 3);
    }

    #[test]
    fn test_nth_root() {
        check(0, 1);
        check(10, 1);
        check(100, 4);
    }

    #[test]
    #[should_panic]
    fn test_nth_root_n_is_zero() {
        check(4, 0);
    }

    #[test]
    fn test_nth_root_big() {
        let x = BigUint::from_str("123_456_789").unwrap();
        let expected = BigUint::from(6u32);

        assert_eq!(x.nth_root(10), expected);
    }
}

mod bigint {
    use num_bigint::BigInt;
    use num_traits::{Signed, Pow};

    fn check(x: i64, n: u32) {
        let big_x = BigInt::from(x);
        let res = big_x.nth_root(n);

        if n == 2 {
            assert_eq!(&res, &big_x.sqrt())
        } else if n == 3 {
            assert_eq!(&res, &big_x.cbrt())
        }

        if big_x.is_negative() {
            assert!(res.pow(n) >= big_x);
            assert!((res - 1u32).pow(n) < big_x);
        } else {
            assert!(res.pow(n) <= big_x);
            assert!((res + 1u32).pow(n) > big_x);
        }
    }

    #[test]
    fn test_nth_root() {
        check(-100, 3);
    }

    #[test]
    #[should_panic]
    fn test_nth_root_x_neg_n_even() {
        check(-100, 4);
    }

    #[test]
    #[should_panic]
    fn test_sqrt_x_neg() {
        check(-4, 2);
    }

    #[test]
    fn test_cbrt() {
        check(8, 3);
        check(-8, 3);
    }
}
