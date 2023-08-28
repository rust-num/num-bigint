#![feature(test)]

extern crate test;

use num_bigint::BigUint;
use num_traits::One;
use std::ops::{Div, Mul};
use test::Bencher;

#[bench]
fn factorial_mul_biguint(b: &mut Bencher) {
    b.iter(|| {
        (1u32..1000)
            .map(BigUint::from)
            .fold(BigUint::one(), Mul::mul)
    });
}

fn factorial_product(l: usize, r: usize) -> BigUint {
    if l >= r {
        BigUint::from(l)
    } else {
        let m = (l+r)/2;
        factorial_product(l, m) * factorial_product(m+1, r)
    }
}

#[bench]
fn factorial_mul_biguint_dnc_10k(b: &mut Bencher) {
    b.iter(|| {
        factorial_product(1, 10_000);
    });
}

#[bench]
fn factorial_mul_biguint_dnc_100k(b: &mut Bencher) {
    b.iter(|| {
        factorial_product(1, 100_000);
    });
}

#[bench]
fn factorial_mul_biguint_dnc_300k(b: &mut Bencher) {
    b.iter(|| {
        factorial_product(1, 300_000);
    });
}

#[bench]
fn factorial_mul_u32(b: &mut Bencher) {
    b.iter(|| (1u32..1000).fold(BigUint::one(), Mul::mul));
}

// The division test is inspired by this blog comparison:
// <https://tiehuis.github.io/big-integers-in-zig#division-test-single-limb>

#[bench]
fn factorial_div_biguint(b: &mut Bencher) {
    let n: BigUint = (1u32..1000).fold(BigUint::one(), Mul::mul);
    b.iter(|| {
        (1u32..1000)
            .rev()
            .map(BigUint::from)
            .fold(n.clone(), Div::div)
    });
}

#[bench]
fn factorial_div_u32(b: &mut Bencher) {
    let n: BigUint = (1u32..1000).fold(BigUint::one(), Mul::mul);
    b.iter(|| (1u32..1000).rev().fold(n.clone(), Div::div));
}
