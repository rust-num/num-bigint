use integer::Integer;
use traits::{One, Zero};

use biguint::BigUint;
use big_digit::BITS;

struct MontyReducer<'a> {
    n: &'a BigUint,
    n0inv: u32
}

// Calculate the modular inverse of `num`, using Extended GCD.
//
// Reference:
// Brent & Zimmermann, Modern Computer Arithmetic, v0.5.9, Algorithm 1.20
fn inv_mod_u32(num: u32) -> u32 {
    // num needs to be relatively prime to 2**32 -- i.e. it must be odd.
    assert!(num % 2 != 0);

    let mut a: i64 = num as i64;
    let mut b: i64 = (u32::max_value() as i64) + 1;

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
        a = b; b = r;
        // 6: (u, w) <- (w, u - qw)
        let m = u - w*q;
        u = w; w = m;
    }

    assert!(a == 1);
    // Downcasting acts like a mod 2^32 too.
    u as u32
}

impl<'a> MontyReducer<'a> {
    fn new(n: &'a BigUint) -> Self {
        let n0inv = inv_mod_u32(n.data[0]);
        MontyReducer { n: n, n0inv: n0inv }
    }

    /// Map a number to the Montgomery domain
    fn map(&self, x: &BigUint) -> BigUint {
        let shift = self.n.data.len() * BITS;
        (x << shift) % self.n
    }

    /// Montgomery Reduction
    ///
    /// Reference:
    /// Brent & Zimmermann, Modern Computer Arithmetic, v0.5.9, Algorithm 2.6
    fn redc(&self, a: BigUint) -> BigUint {
        let mut c = a.data;
        let n = &self.n.data;
        let n_size = n.len();

        // Allocate sufficient work space
        c.resize(2 * n_size + 2, 0);

        // β is the size of a word, in this case 32 bits. So "a mod β" is
        // equivalent to masking a to 32 bits.
        // mu <- -N^(-1) mod β
        let mu = 0u32.wrapping_sub(self.n0inv);

        // 1: for i = 0 to (n-1)
        for i in 0..n_size {
            // 2: q_i <- mu*c_i mod β
            let q_i = c[i].wrapping_mul(mu);

            // 3: C <- C + q_i * N * β^i
            super::algorithms::mac_digit(&mut c[i..], n, q_i);
        }

        // 4: R <- C * β^(-n)
        // This is an n-word bitshift, equivalent to skipping n words.
        let ret = BigUint::new(c[n_size..].to_vec());

        // 5: if R >= β^n then return R-N else return R.
        if &ret < self.n {
            ret
        } else {
            ret - self.n
        }
    }

    /// Montgomery Multiplication
    fn mul(&self, a: BigUint, b: &BigUint) -> BigUint {
        self.redc(a * b)
    }

    /// Montgomery Squaring
    fn square(&self, a: BigUint) -> BigUint {
        // TODO: Replace with an optimised squaring function
        self.redc(&a * &a)
    }

    /// Montgomery Exponentiation
    fn pow(&self, mut base: BigUint, exp: &BigUint) -> BigUint {
        debug_assert!(!exp.is_zero());

        // Binary exponentiation
        let mut exp = exp.clone();
        while exp.is_even() {
            base = self.square(base);
            exp >>= 1;
        }
        if exp.is_one() { return base; }

        let mut acc = base.clone();
        while !exp.is_one() {
            exp >>= 1;
            base = self.square(base);
            if exp.is_odd() {
                acc = self.mul(acc, &base);
            }
        }

        acc
    }
}

pub fn monty_modpow(a: &BigUint, exp: &BigUint, modulus: &BigUint) -> BigUint{
    let mr = MontyReducer::new(modulus);

    // Map the base to the Montgomery domain
    let base = mr.map(a);

    // Do the computation
    let answer = mr.pow(base, exp);

    // Map the result back to the residues domain
    mr.redc(answer)
}
