#![no_main]
use libfuzzer_sys::fuzz_target;
use num_bigint::BigUint;

fuzz_target!(|a: BigUint| {
    assert_eq!((&a * &a).sqrt(), a);
});
