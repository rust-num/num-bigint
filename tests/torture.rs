#![cfg(feature = "rand")]

extern crate num_bigint;
extern crate num_traits;
extern crate rand;

use num_bigint::RandBigInt;
use num_traits::Zero;
use rand::rngs::SmallRng;
use rand::prelude::*;

fn get_rng() -> SmallRng {
    let seed = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16];
    SmallRng::from_seed(seed)
}

fn test_mul_divide_torture_count(count: usize) {
    let bits_max = 1 << 12;
    let mut rng = get_rng();

    for _ in 0..count {
        // Test with numbers of random sizes:
        let xbits = rng.gen_range(0, bits_max);
        let ybits = rng.gen_range(0, bits_max);

        let x = rng.gen_biguint(xbits);
        let y = rng.gen_biguint(ybits);

        if x.is_zero() || y.is_zero() {
            continue;
        }

        let prod = &x * &y;
        assert_eq!(&prod / &x, y);
        assert_eq!(&prod / &y, x);
    }
}

#[test]
fn test_mul_divide_torture() {
    test_mul_divide_torture_count(1_000);
}

#[test]
#[ignore]
fn test_mul_divide_torture_long() {
    test_mul_divide_torture_count(1_000_000);
}

fn test_factored_mul_torture_count(count: usize) {
    let bits = 1 << 16;
    let mut rng = get_rng();

    for _ in 0..count {
        let w = rng.gen_biguint(bits);
        let x = rng.gen_biguint(bits);
        let y = rng.gen_biguint(bits);
        let z = rng.gen_biguint(bits);

        let prod1 = (&w * &x) * (&y * &z);
        let prod2 = (&w * &y) * (&x * &z);
        let prod3 = (&w * &z) * (&x * &y);
        assert_eq!(prod1, prod2);
        assert_eq!(prod2, prod3);
    }
}

#[test]
fn test_factored_mul_torture() {
    test_factored_mul_torture_count(50);
}

#[test]
#[ignore]
fn test_factored_mul_torture_long() {
    test_factored_mul_torture_count(1_000);
}
