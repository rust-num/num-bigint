#![feature(test)]

extern crate test;

// Run with: cargo bench --bench ntt_plan --features bench-internals

use num_bigint::__benchmark_ntt_plan_build;
use test::{black_box, Bencher};

macro_rules! ntt_plan_bench {
    ($name:ident, $prime:expr, $exponent:expr) => {
        #[bench]
        fn $name(bencher: &mut Bencher) {
            let min_len = 1usize << $exponent;
            bencher.iter(|| black_box(__benchmark_ntt_plan_build::<$prime>(black_box(min_len))));
        }
    };
}

ntt_plan_bench!(p1_2_09, 1, 9);
ntt_plan_bench!(p1_2_10, 1, 10);
ntt_plan_bench!(p1_2_11, 1, 11);
ntt_plan_bench!(p1_2_12, 1, 12);
ntt_plan_bench!(p1_2_13, 1, 13);
ntt_plan_bench!(p1_2_14, 1, 14);
ntt_plan_bench!(p1_2_15, 1, 15);
ntt_plan_bench!(p1_2_16, 1, 16);
ntt_plan_bench!(p1_2_17, 1, 17);
ntt_plan_bench!(p1_2_18, 1, 18);
ntt_plan_bench!(p1_2_19, 1, 19);
ntt_plan_bench!(p1_2_20, 1, 20);
ntt_plan_bench!(p2_2_09, 2, 9);
ntt_plan_bench!(p2_2_10, 2, 10);
ntt_plan_bench!(p2_2_11, 2, 11);
ntt_plan_bench!(p2_2_12, 2, 12);
ntt_plan_bench!(p2_2_13, 2, 13);
ntt_plan_bench!(p2_2_14, 2, 14);
ntt_plan_bench!(p2_2_15, 2, 15);
ntt_plan_bench!(p2_2_16, 2, 16);
ntt_plan_bench!(p2_2_17, 2, 17);
ntt_plan_bench!(p2_2_18, 2, 18);
ntt_plan_bench!(p2_2_19, 2, 19);
ntt_plan_bench!(p2_2_20, 2, 20);
ntt_plan_bench!(p3_2_09, 3, 9);
ntt_plan_bench!(p3_2_10, 3, 10);
ntt_plan_bench!(p3_2_11, 3, 11);
ntt_plan_bench!(p3_2_12, 3, 12);
ntt_plan_bench!(p3_2_13, 3, 13);
ntt_plan_bench!(p3_2_14, 3, 14);
ntt_plan_bench!(p3_2_15, 3, 15);
ntt_plan_bench!(p3_2_16, 3, 16);
ntt_plan_bench!(p3_2_17, 3, 17);
ntt_plan_bench!(p3_2_18, 3, 18);
ntt_plan_bench!(p3_2_19, 3, 19);
ntt_plan_bench!(p3_2_20, 3, 20);
