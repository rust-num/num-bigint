extern crate num_bigint;
extern crate num_integer;
extern crate num_traits;

mod biguint {
    use num_bigint::BigUint;
    use num_traits::FromPrimitive;
    use std::str::FromStr;

    fn check(x: i32, n: u32, expected: i32) {
        let big_x: BigUint = FromPrimitive::from_i32(x).unwrap();
        let big_expected: BigUint = FromPrimitive::from_i32(expected).unwrap();

        assert_eq!(big_x.nth_root(n), big_expected);
    }

    #[test]
    fn test_sqrt() {
        check(99, 2, 9);
        check(100, 2, 10);
        check(120, 2, 10);
    }

    #[test]
    fn test_cbrt() {
        check(8, 3, 2);
        check(26, 3, 2);
    }

    #[test]
    fn test_nth_root() {
        check(0, 1, 0);
        check(10, 1, 10);
        check(100, 4, 3);
    }

    #[test]
    #[should_panic]
    fn test_nth_root_n_is_zero() {
        check(4, 0, 0);
    }

    #[test]
    fn test_nth_root_big() {
        let x: BigUint = FromStr::from_str("123_456_789").unwrap();
        let expected : BigUint = FromPrimitive::from_i32(6).unwrap();

        assert_eq!(x.nth_root(10), expected);
    }
}

mod bigint {
    use num_bigint::BigInt;
    use num_traits::FromPrimitive;

    fn check(x: i32, n: u32, expected: i32) {
        let big_x: BigInt = FromPrimitive::from_i32(x).unwrap();
        let big_expected: BigInt = FromPrimitive::from_i32(expected).unwrap();

        assert_eq!(big_x.nth_root(n), big_expected);
    }

    #[test]
    fn test_nth_root() {
        check(-100, 3, -4);
    }

    #[test]
    #[should_panic]
    fn test_nth_root_x_neg_n_even() {
        check(-100, 4, 0);
    }

    #[test]
    #[should_panic]
    fn test_sqrt_x_neg() {
        check(-4, 2, -2);
    }

    #[test]
    fn test_cbrt() {
        check(-8, 3, -2);
    }
}
