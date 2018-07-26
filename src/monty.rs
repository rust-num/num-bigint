use std::ops::Shl;
use traits::{One, Zero};

use big_digit::{self, BigDigit, DoubleBigDigit, SignedDoubleBigDigit};
use biguint::BigUint;

struct MontyReducer<'a> {
    n: &'a BigUint,
    n0inv: BigDigit,
}

// k0 = -m**-1 mod 2**BITS. Algorithm from: Dumas, J.G. "On Newton–Raphson
// Iteration for Multiplicative Inverses Modulo Prime Powers".
fn inv_mod_alt(b: BigDigit) -> BigDigit {
    assert_ne!(b & 1, 0);

    let mut k0 = 2 - b as SignedDoubleBigDigit;
    let mut t = (b - 1) as SignedDoubleBigDigit;
    let mut i = 1;
    while i < big_digit::BITS {
        t = t.wrapping_mul(t);
        k0 = k0.wrapping_mul(t + 1);

        i <<= 1;
    }
    -k0 as BigDigit
}

// Calculate the modular inverse of `num`, using Extended GCD.
//
// Reference:
// Brent & Zimmermann, Modern Computer Arithmetic, v0.5.9, Algorithm 1.20
fn inv_mod(num: BigDigit) -> BigDigit {
    // num needs to be relatively prime to 2**{32|64} -- i.e. it must be odd.
    assert!(num % 2 != 0);

    let mut a: SignedDoubleBigDigit = num as SignedDoubleBigDigit;
    let mut b: SignedDoubleBigDigit = (BigDigit::max_value() as SignedDoubleBigDigit) + 1;

    // ExtendedGcd
    // Input: positive integers a and b
    // Output: integers (g, u, v) such that g = gcd(a, b) = ua + vb
    // As we don't need v for modular inverse, we don't calculate it.

    // 1: (u, w) <- (1, 0)
    let mut u = 1;
    let mut w = 0;
    // 3: while b != 0
    while b != 0 {
        // 4: (q, r) <- DivRem(a, b)
        let q = a / b;
        let r = a % b;
        // 5: (a, b) <- (b, r)
        a = b;
        b = r;
        // 6: (u, w) <- (w, u - qw)
        let m = u - w * q;
        u = w;
        w = m;
    }

    assert!(a == 1);
    // Downcasting acts like a mod 2^32 too.
    u as BigDigit
}

impl<'a> MontyReducer<'a> {
    fn new(n: &'a BigUint) -> Self {
        let n0inv = inv_mod_alt(n.data[0]);
        MontyReducer { n: n, n0inv: n0inv }
    }
}

// Montgomery Reduction
//
// Reference:
// Brent & Zimmermann, Modern Computer Arithmetic, v0.5.9, Algorithm 2.6
fn monty_redc(a: BigUint, mr: &MontyReducer) -> BigUint {
    // println!("redc: {}", &a);
    let mut c = a.data;
    let n = &mr.n.data;
    let n_size = n.len();

    // Allocate sufficient work space
    c.resize(2 * n_size + 2, 0);

    // β is the size of a word, in this case {32|64} bits. So "a mod β" is
    // equivalent to masking a to {32|64} bits.
    // mu <- -N^(-1) mod β
    let mu = mr.n0inv; //(0 as BigDigit).wrapping_sub(mr.n0inv);

    // 1: for i = 0 to (n-1)
    for i in 0..n_size {
        // 2: q_i <- mu*c_i mod β
        let q_i = c[i].wrapping_mul(mu as BigDigit);
        // println!("{} {}", q_i, c[i]);

        // 3: C <- C + q_i * N * β^i
        super::algorithms::mac_digit(&mut c[i..], n, q_i);
        // println!("updated: {:?}", c);
    }

    // println!("c: {:?}", c);
    // println!("n_size: {}\nmu: {:?}", n_size, mu);
    // println!("n: {:?}", n);

    // 4: R <- C * β^(-n)
    // This is an n-word bitshift, equivalent to skipping n words.
    let ret = BigUint::new_native(c[n_size..].to_vec());

    // 5: if R >= β^n then return R-N else return R.
    if &ret < mr.n {
        ret
    } else {
        ret - mr.n
    }
}

// Montgomery Multiplication
fn monty_mult(a: BigUint, b: &BigUint, mr: &MontyReducer) -> BigUint {
    // println!("mul: {} * {}", &a, &b);
    monty_redc(a * b, mr)
}

// Montgomery Squaring
fn monty_sqr(a: BigUint, mr: &MontyReducer) -> BigUint {
    // println!("sqr {} **2", &a);
    // TODO: Replace with an optimised squaring function
    monty_redc(&a * &a, mr)
}

/// Computes z mod m = x * y * 2 ** (-n*_W) mod m
/// assuming k = -1/m mod 2**_W
/// See Gueron, "Efficient Software Implementations of Modular Exponentiation".
/// https://eprint.iacr.org/2011/239.pdf
/// In the terminology of that paper, this is an "Almost Montgomery Multiplication":
/// x and y are required to satisfy 0 <= z < 2**(n*_W) and then the result
/// z is guaranteed to satisfy 0 <= z < 2**(n*_W), but it may not be < m.
fn montgomery(x: &BigUint, y: &BigUint, m: &BigUint, k: BigDigit, n: usize) -> BigUint {
    // println!("montgomery: {:?} {:?} {:?} {} {}", x, y, m, k, n);
    // This code assumes x, y, m are all the same length, n.
    // (required by addMulVVW and the for loop).
    // It also assumes that x, y are already reduced mod m,
    // or else the result will not be properly reduced.
    assert!(
        x.data.len() == n && y.data.len() == n && m.data.len() == n,
        "{:?} {:?} {:?} {}",
        x,
        y,
        m,
        n
    );

    let mut z = BigUint::zero();
    z.data.resize(n * 2, 0);

    let mut c: BigDigit = 0;
    for i in 0..n {
        let c2 = add_mul_vvw(&mut z.data[i..n + i], &x.data, y.data[i]);
        let t = z.data[i].wrapping_mul(k);
        let c3 = add_mul_vvw(&mut z.data[i..n + i], &m.data, t);
        let cx = c.wrapping_add(c2);
        let cy = cx.wrapping_add(c3);
        z.data[n + i] = cy;
        if cx < c2 || cy < c3 {
            c = 1;
        } else {
            c = 0;
        }
    }

    if c == 0 {
        z.data = z.data[n..].to_vec();
    } else {
        {
            let (mut first, second) = z.data.split_at_mut(n);
            sub_vv(&mut first, &second, &m.data);
        }
        z.data = z.data[..n].to_vec();
    }

    z
}

#[inline(always)]
fn add_mul_vvw(z: &mut [BigDigit], x: &[BigDigit], y: BigDigit) -> BigDigit {
    let mut c = 0;
    for (i, zi) in z.iter_mut().enumerate() {
        let (z1, z0) = mul_add_www(x[i], y, *zi);
        let (c_, zi_) = add_ww(z0, c, 0);
        *zi = zi_;
        c = c_ + z1;
    }

    c
}

/// The resulting carry c is either 0 or 1.
#[inline(always)]
fn sub_vv(z: &mut [BigDigit], x: &[BigDigit], y: &[BigDigit]) -> BigDigit {
    let mut c = 0;
    for (i, xi) in x.iter().enumerate().take(z.len()) {
        let yi = y[i];
        let zi = xi.wrapping_sub(yi).wrapping_sub(c);
        z[i] = zi;
        // see "Hacker's Delight", section 2-12 (overflow detection)
        c = ((yi & !xi) | ((yi | !xi) & zi)) >> (big_digit::BITS - 1)
    }

    c
}

/// z1<<_W + z0 = x+y+c, with c == 0 or 1
#[inline(always)]
fn add_ww(x: BigDigit, y: BigDigit, c: BigDigit) -> (BigDigit, BigDigit) {
    let yc = y.wrapping_add(c);
    let z0 = x.wrapping_add(yc);
    let mut z1 = 0;
    if z0 < x || yc < y {
        z1 = 1;
    }

    (z1, z0)
}

/// z1 << _W + z0 = x * y + c
#[inline(always)]
fn mul_add_www(x: BigDigit, y: BigDigit, c: BigDigit) -> (BigDigit, BigDigit) {
    let (z1, zz0) = mul_ww(x, y);
    let z0 = zz0.wrapping_add(c);
    if z0 < zz0 {
        (z1 + 1, z0)
    } else {
        (z1, z0)
    }
}

/// z1<<_W + z0 = x*y
/// Adapted from Warren, Hacker's Delight, p. 132.
#[inline(always)]
fn mul_ww(x: BigDigit, y: BigDigit) -> (BigDigit, BigDigit) {
    let x0 = x & big_digit::HALF_DIGIT_MASK as BigDigit;
    let x1 = x >> big_digit::HALF_WORD_SIZE;
    let y0 = y & big_digit::HALF_DIGIT_MASK as BigDigit;
    let y1 = y >> big_digit::HALF_WORD_SIZE;
    let w0 = x0.wrapping_mul(y0);
    let t = x1
        .wrapping_mul(y0)
        .wrapping_add(w0 >> big_digit::HALF_WORD_SIZE);

    let mut w1 = t & big_digit::HALF_DIGIT_MASK as BigDigit;
    let w2 = t >> big_digit::HALF_WORD_SIZE;
    w1 = w1.wrapping_add(x0.wrapping_mul(y1));
    let z1 = x1
        .wrapping_mul(y1)
        .wrapping_add(w2)
        .wrapping_add(w1 >> big_digit::HALF_WORD_SIZE);
    let z0 = x.wrapping_mul(y);

    (z1, z0)
}

/// Calculates x ** y mod m using a fixed, 4-bit window.
pub fn monty_modpow(x: &BigUint, y: &BigUint, m: &BigUint) -> BigUint {
    assert!(m.data[0] & 1 == 1);
    let mr = MontyReducer::new(m);
    let num_words = m.data.len();

    println!("modpow {:?} {:?} {:?}", x, y, m);
    println!("numWords {}", num_words);
    let mut x = x.clone();

    println!("inverse: {:?}", mr.n0inv);
    // We want the lengths of x and m to be equal.
    // It is OK if x >= m as long as len(x) == len(m).
    if x.data.len() > num_words {
        x %= m;
        // Note: now len(x) <= numWords, not guaranteed ==.
    }
    if x.data.len() < num_words {
        x.data.resize(num_words, 0);
    }

    // rr = 2**(2*_W*len(m)) mod m
    let mut rr = BigUint::one();
    rr = (rr.shl(2 * num_words * big_digit::BITS)) % m;
    if rr.data.len() < num_words {
        rr.data.resize(num_words, 0);
    }
    println!("rr: {:?}", rr);
    // one = 1, with equal length to that of m
    let mut one = BigUint::one();
    one.data.resize(num_words, 0);

    let n = 4;
    // powers[i] contains x^i
    let mut powers = Vec::with_capacity(1 << n);
    powers.push(montgomery(&one, &rr, m, mr.n0inv, num_words));
    powers.push(montgomery(&x, &rr, m, mr.n0inv, num_words));
    for i in 2..1 << n {
        let r = montgomery(&powers[i - 1], &powers[1], m, mr.n0inv, num_words);
        powers.push(r);
    }
    // println!("powers: {:?} {}", powers, 1 << n);
    // initialize z = 1 (Montgomery 1)
    let mut z = powers[0].clone();
    z.data.resize(num_words, 0);
    let mut zz = BigUint::zero();
    zz.data.resize(num_words, 0);

    println!("powers: {:?}", powers);
    // same windowed exponent, but with Montgomery multiplications
    for i in (0..y.data.len()).rev() {
        let mut yi = y.data[i];
        let mut j = 0;
        while j < big_digit::BITS {
            if i != y.data.len() - 1 || j != 0 {
                zz = montgomery(&z, &z, m, mr.n0inv, num_words);
                z = montgomery(&zz, &zz, m, mr.n0inv, num_words);
                zz = montgomery(&z, &z, m, mr.n0inv, num_words);
                z = montgomery(&zz, &zz, m, mr.n0inv, num_words);
            }
            zz = montgomery(
                &z,
                &powers[(yi >> (big_digit::BITS - n)) as usize],
                m,
                mr.n0inv,
                num_words,
            );
            ::std::mem::swap(&mut z, &mut zz);
            yi <<= n;
            j += n;
        }
    }

    // convert to regular number
    zz = montgomery(&z, &one, m, mr.n0inv, num_words);

    zz.normalize();
    // One last reduction, just in case.
    // See golang.org/issue/13907.
    if &zz >= m {
        // Common case is m has high bit set; in that case,
        // since zz is the same length as m, there can be just
        // one multiple of m to remove. Just subtract.
        // We think that the subtract should be sufficient in general,
        // so do that unconditionally, but double-check,
        // in case our beliefs are wrong.
        // The div is not expected to be reached.
        zz -= m;
        if &zz >= m {
            zz %= m;
        }
    }

    zz.normalize();
    zz
}

#[test]
fn test_inv_mod() {
    use big_digit::DoubleBigDigit;

    let modulus = (BigDigit::max_value() as DoubleBigDigit) + 1;

    for element in 0..100 {
        let element = element as BigDigit;
        if element % 2 == 0 {
            continue;
        }

        let inverse = inv_mod(element.clone());
        let cmp = (inverse as DoubleBigDigit * element as DoubleBigDigit) % modulus;

        assert_eq!(
            cmp as BigDigit, 1,
            "mod_inverse({}, {}) * {} % {} = {}, not 1",
            element, &modulus, element, &modulus, cmp,
        );
    }
}
