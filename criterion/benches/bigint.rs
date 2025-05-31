use criterion::{criterion_group, criterion_main, Bencher, Criterion};
#[cfg(not(feature = "smallvec"))]
use num_bigint::{BigInt, BigUint, RandBigInt};
#[cfg(feature = "smallvec")]
use num_bigint_small::{BigInt, BigUint, RandBigInt};
use num_traits::{FromPrimitive, Num, One, Zero};
use std::mem::replace;

mod rng;
use rng::get_rng;

fn multiply_bench(b: &mut Bencher, xbits: u64, ybits: u64) {
    let mut rng = get_rng();
    let x = rng.gen_bigint(xbits);
    let y = rng.gen_bigint(ybits);

    b.iter(|| &x * &y);
}

fn divide_bench(b: &mut Bencher, xbits: u64, ybits: u64) {
    let mut rng = get_rng();
    let x = rng.gen_bigint(xbits);
    let y = rng.gen_bigint(ybits);

    b.iter(|| &x / &y);
}

fn remainder_bench(b: &mut Bencher, xbits: u64, ybits: u64) {
    let mut rng = get_rng();
    let x = rng.gen_bigint(xbits);
    let y = rng.gen_bigint(ybits);

    b.iter(|| &x % &y);
}

fn factorial(n: usize) -> BigUint {
    let mut f: BigUint = One::one();
    for i in 1..=n {
        let bu: BigUint = FromPrimitive::from_usize(i).unwrap();
        f *= bu;
    }
    f
}

/// Compute Fibonacci numbers
fn fib(n: usize) -> BigUint {
    let mut f0: BigUint = Zero::zero();
    let mut f1: BigUint = One::one();
    for _ in 0..n {
        let f2 = f0 + &f1;
        f0 = replace(&mut f1, f2);
    }
    f0
}

/// Compute Fibonacci numbers with two ops per iteration
/// (add and subtract, like issue #200)
fn fib2(n: usize) -> BigUint {
    let mut f0: BigUint = Zero::zero();
    let mut f1: BigUint = One::one();
    for _ in 0..n {
        f1 += &f0;
        f0 = &f1 - f0;
    }
    f0
}

fn fib_u128(n: usize) -> u128 {
    let mut f0: u128 = Zero::zero();
    let mut f1: u128 = One::one();
    for _ in 0..n {
        f1 = f1.checked_add(f0).unwrap();
        f0 = f1.checked_sub(f0).unwrap();
    }
    f0
}

fn to_str_radix_bench(b: &mut Bencher, radix: u32) {
    let mut rng = get_rng();
    let x = rng.gen_bigint(1009);
    b.iter(|| x.to_str_radix(radix));
}

fn from_str_radix_bench(b: &mut Bencher, radix: u32) {
    let mut rng = get_rng();
    let x = rng.gen_bigint(1009);
    let s = x.to_str_radix(radix);
    assert_eq!(x, BigInt::from_str_radix(&s, radix).unwrap());
    b.iter(|| BigInt::from_str_radix(&s, radix));
}

fn rand_bench(b: &mut Bencher, bits: u64) {
    let mut rng = get_rng();

    b.iter(|| rng.gen_bigint(bits));
}

/// This modulus is the prime from the 2048-bit MODP DH group:
/// https://tools.ietf.org/html/rfc3526#section-3
const RFC3526_2048BIT_MODP_GROUP: &str = "\
                                          FFFFFFFF_FFFFFFFF_C90FDAA2_2168C234_C4C6628B_80DC1CD1\
                                          29024E08_8A67CC74_020BBEA6_3B139B22_514A0879_8E3404DD\
                                          EF9519B3_CD3A431B_302B0A6D_F25F1437_4FE1356D_6D51C245\
                                          E485B576_625E7EC6_F44C42E9_A637ED6B_0BFF5CB6_F406B7ED\
                                          EE386BFB_5A899FA5_AE9F2411_7C4B1FE6_49286651_ECE45B3D\
                                          C2007CB8_A163BF05_98DA4836_1C55D39A_69163FA8_FD24CF5F\
                                          83655D23_DCA3AD96_1C62F356_208552BB_9ED52907_7096966D\
                                          670C354E_4ABC9804_F1746C08_CA18217C_32905E46_2E36CE3B\
                                          E39E772C_180E8603_9B2783A2_EC07A28F_B5C55DF0_6F4C52C9\
                                          DE2BCBF6_95581718_3995497C_EA956AE5_15D22618_98FA0510\
                                          15728E5A_8AACAA68_FFFFFFFF_FFFFFFFF";

fn criterion_benchmark(c: &mut Criterion) {
    {
        let mut group = c.benchmark_group("multiply");
        group.bench_function("0", |b| multiply_bench(b, 1 << 8, 1 << 8));
        group.bench_function("1", |b| multiply_bench(b, 1 << 8, 1 << 16));
        group.bench_function("2", |b| multiply_bench(b, 1 << 16, 1 << 16));
        group.bench_function("3", |b| multiply_bench(b, 1 << 16, 1 << 17));
    }

    {
        let mut group = c.benchmark_group("divide");
        group.bench_function("0", |b| divide_bench(b, 1 << 8, 1 << 6));
        group.bench_function("1", |b| divide_bench(b, 1 << 12, 1 << 8));
        group.bench_function("2", |b| divide_bench(b, 1 << 16, 1 << 12));
        group.bench_function("big_little", |b| divide_bench(b, 1 << 16, 1 << 4));
    }

    {
        let mut group = c.benchmark_group("remainder");
        group.bench_function("0", |b| remainder_bench(b, 1 << 8, 1 << 6));
        group.bench_function("1", |b| remainder_bench(b, 1 << 12, 1 << 8));
        group.bench_function("2", |b| remainder_bench(b, 1 << 16, 1 << 12));
        group.bench_function("big_little", |b| remainder_bench(b, 1 << 16, 1 << 4));
    }

    c.bench_function("factorial_100", |b| b.iter(|| factorial(100)));

    {
        let mut group = c.benchmark_group("fib");
        group.bench_function("100", |b| b.iter(|| fib(100)));
        group.bench_function("1000", |b| b.iter(|| fib(1000)));
        group.bench_function("10000", |b| b.iter(|| fib(10000)));
    }
    {
        let mut group = c.benchmark_group("fib2");
        group.bench_function("100", |b| b.iter(|| fib2(100)));
        group.bench_function("1000", |b| b.iter(|| fib2(1000)));
        group.bench_function("10000", |b| b.iter(|| fib2(10000)));
    }
    c.bench_function("fib_u128", |b| b.iter(|| fib_u128(100)));

    c.bench_function("fac_to_string", |b| {
        let fac = factorial(100);
        b.iter(|| fac.to_string())
    });

    c.bench_function("fib_to_string", |b| {
        let fib = fib(100);
        b.iter(|| fib.to_string());
    });

    {
        let mut group = c.benchmark_group("to_str_radix");
        group.bench_function("02", |b| to_str_radix_bench(b, 2));
        group.bench_function("08", |b| to_str_radix_bench(b, 8));
        group.bench_function("10", |b| to_str_radix_bench(b, 10));
        group.bench_function("16", |b| to_str_radix_bench(b, 16));
        group.bench_function("36", |b| to_str_radix_bench(b, 36));
    }

    {
        let mut group = c.benchmark_group("from_str_radix");
        group.bench_function("from_str_radix_02", |b| from_str_radix_bench(b, 2));
        group.bench_function("from_str_radix_08", |b| from_str_radix_bench(b, 8));
        group.bench_function("from_str_radix_10", |b| from_str_radix_bench(b, 10));
        group.bench_function("from_str_radix_16", |b| from_str_radix_bench(b, 16));
        group.bench_function("from_str_radix_36", |b| from_str_radix_bench(b, 36));
    }

    {
        let mut group = c.benchmark_group("rand");
        group.bench_function("64", |b| rand_bench(b, 1 << 6));
        group.bench_function("256", |b| rand_bench(b, 1 << 8));
        group.bench_function("1009", |b| rand_bench(b, 1009));
        group.bench_function("2048", |b| rand_bench(b, 1 << 11));
        group.bench_function("4096", |b| rand_bench(b, 1 << 12));
        group.bench_function("8192", |b| rand_bench(b, 1 << 13));
        group.bench_function("65536", |b| rand_bench(b, 1 << 16));
        group.bench_function("131072", |b| rand_bench(b, 1 << 17));
    }

    c.bench_function("shl", |b| {
        let n = BigUint::one() << 1000u32;
        let mut m = n.clone();
        b.iter(|| {
            m.clone_from(&n);
            for i in 0..50 {
                m <<= i;
            }
        })
    });

    c.bench_function("shr", |b| {
        let n = BigUint::one() << 2000u32;
        let mut m = n.clone();
        b.iter(|| {
            m.clone_from(&n);
            for i in 0..50 {
                m >>= i;
            }
        })
    });

    c.bench_function("hahs", |b| {
        use std::collections::HashSet;
        let mut rng = get_rng();
        let v: Vec<BigInt> = (1000..2000).map(|bits| rng.gen_bigint(bits)).collect();
        b.iter(|| {
            let h: HashSet<&BigInt> = v.iter().collect();
            assert_eq!(h.len(), v.len());
        });
    });

    c.bench_function("pow_bench", |b| {
        b.iter(|| {
            let upper = 100_u32;
            let mut i_big = BigUint::from(1u32);
            for _i in 2..=upper {
                i_big += 1u32;
                for j in 2..=upper {
                    i_big.pow(j);
                }
            }
        });
    });

    {
        let mut group = c.benchmark_group("pow");
        group.bench_function("pow_bench_bigexp", |b| {
            use num_traits::Pow;

            b.iter(|| {
                let upper = 100_u32;
                let mut i_big = BigUint::from(1u32);
                for _i in 2..=upper {
                    i_big += 1u32;
                    let mut j_big = BigUint::from(1u32);
                    for _j in 2..=upper {
                        j_big += 1u32;
                        Pow::pow(&i_big, &j_big);
                    }
                }
            });
        });

        group.bench_function("1e1000", |b| {
            b.iter(|| BigUint::from(10u32).pow(1_000));
        });

        group.bench_function("1e10000", |b| {
            b.iter(|| BigUint::from(10u32).pow(10_000));
        });

        group.bench_function("1e100000", |b| {
            b.iter(|| BigUint::from(10u32).pow(100_000));
        });

        group.bench_function("modpow", |b| {
            let mut rng = get_rng();
            let base = rng.gen_biguint(2048);
            let e = rng.gen_biguint(2048);
            let m = BigUint::from_str_radix(RFC3526_2048BIT_MODP_GROUP, 16).unwrap();

            b.iter(|| base.modpow(&e, &m));
        });

        group.bench_function("modpow_even", |b| {
            let mut rng = get_rng();
            let base = rng.gen_biguint(2048);
            let e = rng.gen_biguint(2048);
            // Make the modulus even, so monty (base-2^32) doesn't apply.
            let m = BigUint::from_str_radix(RFC3526_2048BIT_MODP_GROUP, 16).unwrap() - 1u32;

            b.iter(|| base.modpow(&e, &m));
        });
    }

    {
        let mut group = c.benchmark_group("iters");
        group.bench_function("to_u32_digits", |b| {
            let mut rng = get_rng();
            let n = rng.gen_biguint(2048);

            b.iter(|| n.to_u32_digits());
        });

        group.bench_function("iter_u32_digits", |b| {
            let mut rng = get_rng();
            let n = rng.gen_biguint(2048);

            b.iter(|| n.iter_u32_digits().max());
        });

        group.bench_function("to_u64_digits", |b| {
            let mut rng = get_rng();
            let n = rng.gen_biguint(2048);

            b.iter(|| n.to_u64_digits());
        });

        group.bench_function("iter_u64_digits", |b| {
            let mut rng = get_rng();
            let n = rng.gen_biguint(2048);

            b.iter(|| n.iter_u64_digits().max());
        });
    }
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
