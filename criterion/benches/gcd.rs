use criterion::{criterion_group, criterion_main, Bencher, Criterion};

use num_bigint::{BigUint, RandBigInt};
use num_integer::Integer;
use num_traits::Zero;

mod rng;
use rng::get_rng;

fn bench(b: &mut Bencher, bits: u64, gcd: fn(&BigUint, &BigUint) -> BigUint) {
    let mut rng = get_rng();
    let x = rng.gen_biguint(bits);
    let y = rng.gen_biguint(bits);

    assert_eq!(euclid(&x, &y), x.gcd(&y));

    b.iter(|| gcd(&x, &y));
}

fn euclid(x: &BigUint, y: &BigUint) -> BigUint {
    // Use Euclid's algorithm
    let mut m = x.clone();
    let mut n = y.clone();
    while !m.is_zero() {
        let temp = m;
        m = n % &temp;
        n = temp;
    }
    n
}

fn criterion_benchmark(c: &mut Criterion) {
    {
        let mut group = c.benchmark_group("gcd_euclid");
        group.bench_function("0064", |b| bench(b, 64, euclid));
        group.bench_function("0256", |b| bench(b, 256, euclid));
        group.bench_function("1024", |b| bench(b, 1024, euclid));
        group.bench_function("4096", |b| bench(b, 4096, euclid));
    }

    {
        let mut group = c.benchmark_group("gcd_stein");
        group.bench_function("0064", |b| bench(b, 64, BigUint::gcd));
        group.bench_function("0256", |b| bench(b, 256, BigUint::gcd));
        group.bench_function("1024", |b| bench(b, 1024, BigUint::gcd));
        group.bench_function("4096", |b| bench(b, 4096, BigUint::gcd));
    }
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
