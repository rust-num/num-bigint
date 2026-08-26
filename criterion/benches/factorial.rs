use criterion::{criterion_group, criterion_main, Criterion};

use num_bigint::BigUint;
use num_traits::One;
use std::ops::{Div, Mul};

fn criterion_benchmark(c: &mut Criterion) {
    let mut group = c.benchmark_group("factorial");
    group.bench_function("mul_biguint", |b| {
        b.iter(|| {
            (1u32..1000)
                .map(BigUint::from)
                .fold(BigUint::one(), Mul::mul)
        });
    });

    group.bench_function("mul_u32", |b| {
        b.iter(|| (1u32..1000).fold(BigUint::one(), Mul::mul));
    });

    // The division test is inspired by this blog comparison:
    // <https://tiehuis.github.io/big-integers-in-zig#division-test-single-limb>

    group.bench_function("div_biguint", |b| {
        let n: BigUint = (1u32..1000).fold(BigUint::one(), Mul::mul);
        b.iter(|| {
            (1u32..1000)
                .rev()
                .map(BigUint::from)
                .fold(n.clone(), Div::div)
        });
    });

    group.bench_function("div_u32", |b| {
        let n: BigUint = (1u32..1000).fold(BigUint::one(), Mul::mul);
        b.iter(|| (1u32..1000).rev().fold(n.clone(), Div::div));
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
