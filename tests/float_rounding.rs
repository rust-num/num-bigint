use num_bigint::{BigInt, BigUint};
use num_traits::{One, ToPrimitive};

fn check_f64(n: BigUint, expected: u64) {
    assert_eq!(n.to_f64().unwrap().to_bits(), expected);
    let signed = BigInt::from(n);
    assert_eq!(signed.to_f64().unwrap().to_bits(), expected);
    assert_eq!((-signed).to_f64().unwrap().to_bits(), expected | (1 << 63));
}

fn check_f32(n: BigUint, expected: u32) {
    assert_eq!(n.to_f32().unwrap().to_bits(), expected);
    let signed = BigInt::from(n);
    assert_eq!(signed.to_f32().unwrap().to_bits(), expected);
    assert_eq!((-signed).to_f32().unwrap().to_bits(), expected | (1 << 31));
}

#[test]
fn test_f64_discarded_bits() {
    for exponent in [128usize, 192, 256] {
        let base = BigUint::one() << exponent;
        let half_ulp = BigUint::one() << (exponent - 53);
        let midpoint = &base + &half_ulp;
        let even = ((exponent as u64) + 1023) << 52;

        check_f64(base, even);
        check_f64(&midpoint - 1u8, even);
        check_f64(midpoint.clone(), even);
        // Each bit of a fully discarded 64-bit digit must still affect rounding.
        for bit in 0..64usize {
            check_f64(&midpoint + (BigUint::one() << bit), even + 1);
        }
        check_f64(&midpoint + &half_ulp, even + 1);
        check_f64(midpoint + (half_ulp << 1usize), even + 2);
    }
}

#[test]
fn test_f32_discarded_bits() {
    for exponent in [96usize, 127] {
        let base = BigUint::one() << exponent;
        let half_ulp = BigUint::one() << (exponent - 24);
        let midpoint = &base + &half_ulp;
        let even = ((exponent as u32) + 127) << 23;

        check_f32(base, even);
        check_f32(&midpoint - 1u8, even);
        check_f32(midpoint.clone(), even);
        // On 32-bit targets these bits are in a fully discarded digit.
        for bit in 0..32usize {
            check_f32(&midpoint + (BigUint::one() << bit), even + 1);
        }
        check_f32(&midpoint + &half_ulp, even + 1);
        check_f32(midpoint + (half_ulp << 1usize), even + 2);
    }
}
