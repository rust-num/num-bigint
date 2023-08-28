#![allow(clippy::cast_sign_loss)]
#![allow(clippy::cast_lossless)]
#![allow(clippy::cast_possible_truncation)]
#![allow(clippy::many_single_char_names)]
#![allow(clippy::needless_range_loop)]
#![allow(clippy::similar_names)]

mod arith {
    // Extended Euclid algorithm:
    //   (g, x, y) is a solution to ax + by = g, where g = gcd(a, b)
    pub const fn egcd(mut a: i128, mut b: i128) -> (i128, i128, i128) {
        if a < 0 { a = -a; }
        if b < 0 { b = -b; }
        assert!(a > 0 || b > 0);
        let mut c = [1, 0, 0, 1]; // treat as a row-major 2x2 matrix
        let (g, x, y) = loop {
            if a == 0 { break (b, 0, 1); }
            if b == 0 { break (a, 1, 0); }
            if a < b {
                let (q, r) = (b/a, b%a);
                b = r;
                c = [c[0], c[1] - q*c[0], c[2], c[3] - q*c[2]];
            } else {
                let (q, r) = (a/b, a%b);
                a = r;
                c = [c[0] - q*c[1], c[1], c[2] - q*c[3], c[3]];
            }
        };
        (g, c[0]*x + c[1]*y, c[2]*x + c[3]*y)
    }
    // Modular inverse: a^-1 mod modulus
    //   (m == 0 means m == 2^64)
    pub const fn invmod(a: u64, modulus: u64) -> u64 {
        let m = if modulus == 0 { 1i128 << 64 } else { modulus as i128 };
        let (g, mut x, _y) = egcd(a as i128, m);
        assert!(g == 1);
        if x < 0 { x += m; }
        assert!(x > 0 && x < 1i128 << 64);
        x as u64
    }
}

struct Arith<const P: u64> {}
impl<const P: u64> Arith<P> {
    pub const FACTOR_TWO: usize = (P-1).trailing_zeros() as usize;
    pub const FACTOR_THREE: usize = Self::factor_three();
    pub const FACTOR_FIVE: usize = Self::factor_five();
    pub const MAX_NTT_LEN: u64 = Self::max_ntt_len();
    pub const R: u64 = ((1u128 << 64) % P as u128) as u64; // 2^64 mod P
    pub const R2: u64 = ((Self::R as u128 * Self::R as u128) % P as u128) as u64; // R^2 mod P
    pub const R3: u64 = ((Self::R2 as u128 * Self::R as u128) % P as u128) as u64; // R^3 mod P
    pub const RNEG: u64 = P.wrapping_sub(Self::R); // -2^64 mod P
    pub const PINV: u64 = arith::invmod(P, 0); // P^-1 mod 2^64
    pub const ROOT: u64 = Self::ntt_root(); // MultiplicativeOrder[ROOT, P] == MAX_NTT_LEN
    pub const ROOTR: u64 = Self::mulmod_naive(Self::ROOT, Self::R); // ROOT * R mod P
    const fn factor_three() -> usize {
        let mut tmp = P-1;
        let mut ans = 0;
        while tmp % 3 == 0 { tmp /= 3; ans += 1; }
        ans
    }
    const fn factor_five() -> usize {
        let mut tmp = P-1;
        let mut ans = 0;
        while tmp % 5 == 0 { tmp /= 5; ans += 1; }
        ans
    }
    const fn max_ntt_len() -> u64 {
        let mut ans = 1u64 << Self::FACTOR_TWO;
        let mut i = 0;
        while i < Self::FACTOR_THREE { ans *= 3; i += 1; }
        let mut i = 0;
        while i < Self::FACTOR_FIVE { ans *= 5; i += 1; }
        assert!(ans % 4050 == 0);
        ans
    }
    const fn ntt_root() -> u64 {
        let mut p = 1;
        'outer: loop {
            /* ensure p is prime */
            p += 1;
            let mut i = 2;
            while i * i <= p {
                if p % i == 0 { continue 'outer; }
                i += 1;
            }
            let root = Self::powmod_naive(p, P/Self::MAX_NTT_LEN);
            let mut j = 0;
            while j <= Self::FACTOR_TWO {
                let mut k = 0;
                while k <= Self::FACTOR_THREE {
                    let mut l = 0;
                    while l <= Self::FACTOR_FIVE {
                        let p2 = Self::powmod_naive(2, j as u64);
                        let p3 = Self::powmod_naive(3, k as u64);
                        let p5 = Self::powmod_naive(5, l as u64);
                        let exponent = p2 * p3 * p5;
                        if exponent < Self::MAX_NTT_LEN && Self::powmod_naive(root, exponent) == 1 {
                            continue 'outer;
                        }
                        l += 1;
                    }
                    k += 1;
                }
                j += 1;
            }
            break root
        }
    }
    // Computes a * b mod P
    const fn mulmod_naive(a: u64, b: u64) -> u64 {
        ((a as u128 * b as u128) % P as u128) as u64
    }
    // Computes base^exponent mod P
    const fn powmod_naive(base: u64, exponent: u64) -> u64 {
        let mut cur = 1;
        let mut pow = base as u128;
        let mut p = exponent;
        while p > 0 {
            if p % 2 > 0 {
                cur = (cur * pow) % P as u128;
            }
            p /= 2;
            pow = (pow * pow) % P as u128;
        }
        cur as u64
    }
    // Multiplication with Montgomery reduction:
    //   a * b * R^-1 mod P
    pub const fn mmulmod(a: u64, b: u64) -> u64 {
        let x = a as u128 * b as u128;
        let m = (x as u64).wrapping_mul(Self::PINV);
        let y = ((m as u128 * P as u128) >> 64) as u64;
        let (out, overflow) = ((x >> 64) as u64).overflowing_sub(y);
        if overflow { out.wrapping_add(P) } else { out }
    }
    pub const fn mmulmod_cond<const INV: bool>(a: u64, b: u64) -> u64 {
        if INV { Self::mmulmod(a, b) } else { b }
    }
    // Fused-multiply-add with Montgomery reduction:
    //   a * b * R^-1 + c mod P
    pub const fn mmuladdmod(a: u64, b: u64, c: u64) -> u64 {
        let x = a as u128 * b as u128;
        let hi = Self::addmod((x >> 64) as u64, c);
        let m = (x as u64).wrapping_mul(Self::PINV);
        let y = ((m as u128 * P as u128) >> 64) as u64;
        let (out, overflow) = hi.overflowing_sub(y);
        if overflow { out.wrapping_add(P) } else { out }
    }
    // Fused-multiply-sub with Montgomery reduction:
    //   a * b * R^-1 - c mod P
    pub const fn mmulsubmod(a: u64, b: u64, c: u64) -> u64 {
        let x = a as u128 * b as u128;
        let hi = Self::submod((x >> 64) as u64, c);
        let m = (x as u64).wrapping_mul(Self::PINV);
        let y = ((m as u128 * P as u128) >> 64) as u64;
        let (out, overflow) = hi.overflowing_sub(y);
        if overflow { out.wrapping_add(P) } else { out }
    }
    // Computes base^exponent mod P with Montgomery reduction
    pub const fn mpowmod(base: u64, exponent: u64) -> u64 {
        let mut cur = Self::R;
        let mut pow = base;
        let mut p = exponent;
        while p > 0 {
            if p % 2 > 0 {
                cur = Self::mmulmod(cur, pow);
            }
            p /= 2;
            pow = Self::mmulmod(pow, pow);
        }
        cur
    }
    // Computes a + b mod P, output range [0, P)
    pub const fn addmod(a: u64, b: u64) -> u64 {
        Self::submod(a, P.wrapping_sub(b))
    }
    // Computes a + b mod P, output range [0, 2^64)
    pub const fn addmod64(a: u64, b: u64) -> u64 {
        let (out, overflow) = a.overflowing_add(b);
        if overflow { out.wrapping_sub(P) } else { out }
    }
    // Computes a + b mod P, selects addmod64 or addmod depending on INV
    pub const fn addmodopt<const INV: bool>(a: u64, b: u64) -> u64 {
        if INV { Self::addmod64(a, b) } else { Self::addmod(a, b) }
    }
    // Computes a - b mod P, output range [0, P)
    pub const fn submod(a: u64, b: u64) -> u64 {
        let (out, overflow) = a.overflowing_sub(b);
        if overflow { out.wrapping_add(P) } else { out }
    }
}

struct NttKernelImpl<const P: u64, const RADIX: usize, const INV: bool>;
impl<const P: u64, const RADIX: usize, const INV: bool> NttKernelImpl<P, RADIX, INV> {
    pub const ROOTR: u64 = Arith::<P>::mpowmod(Arith::<P>::ROOTR, if INV { Arith::<P>::MAX_NTT_LEN - 1 } else { 1 });
    pub const U3: u64 = Arith::<P>::mpowmod(Self::ROOTR, Arith::<P>::MAX_NTT_LEN/3);
    pub const U4: u64 = Arith::<P>::mpowmod(Self::ROOTR, Arith::<P>::MAX_NTT_LEN/4);
    pub const U5: u64 = Arith::<P>::mpowmod(Self::ROOTR, Arith::<P>::MAX_NTT_LEN/5);
    pub const U6: u64 = Arith::<P>::mpowmod(Self::ROOTR, Arith::<P>::MAX_NTT_LEN/6);
    pub const C51: u64 = Self::c5().0;
    pub const C52: u64 = Self::c5().1;
    pub const C53: u64 = Self::c5().2;
    pub const C54: u64 = Self::c5().3;
    pub const C55: u64 = Self::c5().4;
    const fn c5() -> (u64, u64, u64, u64, u64) {
        let w = Self::U5;
        let w2 = Arith::<P>::mpowmod(w, 2);
        let w4 = Arith::<P>::mpowmod(w, 4);
        let inv2 = Arith::<P>::mmulmod(Arith::<P>::R2, arith::invmod(2, P)); 
        let inv4 = Arith::<P>::mmulmod(Arith::<P>::R2, arith::invmod(4, P));
        let c51 = Arith::<P>::submod(Arith::<P>::submod(0, Arith::<P>::R), inv4); // (-1) + (-1) * 4^-1 mod P
        let c52 = Arith::<P>::addmod(Arith::<P>::mmulmod(inv2, Arith::<P>::addmod(w, w4)), inv4); // 4^-1 * (2*w + 2*w^4 + 1) mod P
        let c53 = Arith::<P>::mmulmod(inv2, Arith::<P>::submod(w, w4)); // 2^-1 * (w - w^4) mod P
        let c54 = Arith::<P>::addmod(Arith::<P>::addmod(w, w2), inv2); // 2^-1 * (2*w + 2*w^2 + 1) mod P
        let c55 = Arith::<P>::addmod(Arith::<P>::addmod(w2, w4), inv2); // 2^-1 * (2*w^2 + 2*w^4 + 1) mod P
        (c51, c52, c53, c54, c55)
    }
}

impl<const P: u64, const INV: bool> NttKernelImpl<P, 2, INV> {
    unsafe fn apply<'a>(n: usize, s: usize, eo: bool, x: &'a mut [u64], y: &'a mut [u64]) -> (usize, usize, bool, &'a mut [u64], &'a mut [u64]) {
        let mut src = x.as_ptr();
        let mut dst = y.as_mut_ptr();
        let omega1 = Arith::<P>::mpowmod(Self::ROOTR, Arith::<P>::MAX_NTT_LEN/n as u64);
        let (n1, n1s) = (n/2, n/2*s);
        let mut w1p = Arith::<P>::R;
        for _ in 0..n1 {
            for _ in 0..s {
                let a = *src.wrapping_add(0);
                let b = *src.wrapping_add(n1s);
                *dst.wrapping_add(0) = Arith::<P>::addmod(a, b);
                *dst.wrapping_add(s) = Arith::<P>::mmulmod(w1p, Arith::<P>::submod(a, b));
                src = src.wrapping_add(1);
                dst = dst.wrapping_add(1);
            }
            dst = dst.wrapping_add(s);
            w1p = Arith::<P>::mmulmod(w1p, omega1);
        }
        (n/2, s*2, !eo, y, x)
    }
    unsafe fn apply_last<'a>(n: usize, s: usize, eo: bool, x: &'a mut [u64], y: &'a mut [u64], mult: u64) -> (usize, usize, bool, &'a mut [u64], &'a mut [u64]) {
        assert_eq!(n, 2);
        let mut src = x.as_ptr();
        let mut dst = if eo { y.as_mut_ptr() } else { x.as_mut_ptr() };
        for _ in 0..s {
            let a = *src.wrapping_add(0);
            let b = *src.wrapping_add(s);
            *dst.wrapping_add(0) = Arith::<P>::mmulmod_cond::<INV>(mult, Arith::<P>::addmodopt::<INV>(a, b));
            *dst.wrapping_add(s) = Arith::<P>::mmulmod_cond::<INV>(mult, Arith::<P>::submod(a, b));
            src = src.wrapping_add(1);
            dst = dst.wrapping_add(1);
        }
        if eo { (n/2, s*2, !eo, y, x) } else { (n/2, s*2, eo, x, y) }
    }
}

impl<const P: u64, const INV: bool> NttKernelImpl<P, 3, INV> {
    unsafe fn apply<'a>(n: usize, s: usize, eo: bool, x: &'a mut [u64], y: &'a mut [u64]) -> (usize, usize, bool, &'a mut [u64], &'a mut [u64]) {
        let mut src = x.as_ptr();
        let mut dst = y.as_mut_ptr();
        let omega1 = Arith::<P>::mpowmod(Self::ROOTR, Arith::<P>::MAX_NTT_LEN/n as u64);
        let (n1, n1s) = (n/3, n/3*s);
        let (mut w1p, mut w2p) = (Arith::<P>::R, Arith::<P>::R);
        for _ in 0..n1 {
            for _ in 0..s {
                let a = *src.wrapping_add(0);
                let b = *src.wrapping_add(n1s);
                let c = *src.wrapping_add(2*n1s);
                let kbmc = Arith::<P>::mmulmod(Self::U3, Arith::<P>::submod(b, c));
                *dst.wrapping_add(0) = Arith::<P>::addmod(a, Arith::<P>::addmod(b, c));
                *dst.wrapping_add(s) = Arith::<P>::mmulmod(w1p, Arith::<P>::addmod64(Arith::<P>::submod(a, c), kbmc));
                *dst.wrapping_add(2*s) = Arith::<P>::mmulmod(w2p, Arith::<P>::submod(Arith::<P>::submod(a, b), kbmc));
                src = src.wrapping_add(1);
                dst = dst.wrapping_add(1);
            }
            dst = dst.wrapping_add(2*s);
            w1p = Arith::<P>::mmulmod(w1p, omega1);
            w2p = Arith::<P>::mmulmod(w1p, w1p);
        }
        (n/3, s*3, !eo, y, x)
    }
    unsafe fn apply_last<'a>(n: usize, s: usize, eo: bool, x: &'a mut [u64], y: &'a mut [u64], mult: u64) -> (usize, usize, bool, &'a mut [u64], &'a mut [u64]) {
        assert_eq!(n, 3);
        let mut src = x.as_ptr();
        let mut dst = if eo { y.as_mut_ptr() } else { x.as_mut_ptr() };
        for _ in 0..s {
            let a = *src.wrapping_add(0);
            let b = *src.wrapping_add(s);
            let c = *src.wrapping_add(2*s);
            let kbmc = Arith::<P>::mmulmod(Self::U3, Arith::<P>::submod(b, c));
            *dst.wrapping_add(0) = Arith::<P>::mmulmod_cond::<INV>(mult, Arith::<P>::addmodopt::<INV>(a, Arith::<P>::addmodopt::<INV>(b, c)));
            *dst.wrapping_add(s) = Arith::<P>::mmulmod_cond::<INV>(mult, Arith::<P>::addmodopt::<INV>(Arith::<P>::submod(a, c), kbmc));
            *dst.wrapping_add(2*s) = Arith::<P>::mmulmod_cond::<INV>(mult, Arith::<P>::submod(Arith::<P>::submod(a, b), kbmc));
            src = src.wrapping_add(1);
            dst = dst.wrapping_add(1);
        }
        if eo { (n/3, s*3, !eo, y, x) } else { (n/3, s*3, eo, x, y) }
    }
}

impl<const P: u64, const INV: bool> NttKernelImpl<P, 4, INV> {
    unsafe fn apply<'a>(n: usize, s: usize, eo: bool, x: &'a mut [u64], y: &'a mut [u64]) -> (usize, usize, bool, &'a mut [u64], &'a mut [u64]) {
        let mut src = x.as_ptr();
        let mut dst = y.as_mut_ptr();
        let omega1 = Arith::<P>::mpowmod(Self::ROOTR, Arith::<P>::MAX_NTT_LEN/n as u64);
        let (n1, n1s) = (n/4, n/4*s);
        let (mut w1p, mut w2p, mut w3p) = (Arith::<P>::R, Arith::<P>::R, P.wrapping_sub(Self::U4));
        for _ in 0..n1 {
            for _ in 0..s {
                let a = *src.wrapping_add(0);
                let b = *src.wrapping_add(n1s);
                let c = *src.wrapping_add(2*n1s);
                let d = *src.wrapping_add(3*n1s);
                let apc = Arith::<P>::addmod(a, c);
                let amc = Arith::<P>::mmulmod(w1p, Arith::<P>::submod(a, c));
                let bpd = Arith::<P>::addmod(b, d);
                let bmd = Arith::<P>::submod(b, d);
                let jbmd = Arith::<P>::mmulmod(w3p, bmd);
                *dst.wrapping_add(0) = Arith::<P>::addmod(apc, bpd);
                *dst.wrapping_add(s) = Arith::<P>::submod(amc, jbmd);
                *dst.wrapping_add(2*s) = Arith::<P>::mmulmod(w2p, Arith::<P>::submod(apc,  bpd));
                *dst.wrapping_add(3*s) = Arith::<P>::mmulmod(w2p, Arith::<P>::addmod64(amc, jbmd));
                src = src.wrapping_add(1);
                dst = dst.wrapping_add(1);
            }
            dst = dst.wrapping_add(3*s);
            w1p = Arith::<P>::mmulmod(w1p, omega1);
            w2p = Arith::<P>::mmulmod(w1p, w1p);
            w3p = Arith::<P>::mmulmod(w3p, omega1);
        }
        (n/4, s*4, !eo, y, x)
    }
    unsafe fn apply_last<'a>(n: usize, s: usize, eo: bool, x: &'a mut [u64], y: &'a mut [u64], mult: u64) -> (usize, usize, bool, &'a mut [u64], &'a mut [u64]) {
        assert_eq!(n, 4);
        let mut src = x.as_ptr();
        let mut dst = if eo { y.as_mut_ptr() } else { x.as_mut_ptr() };
        for _ in 0..s {
            let a = *src.wrapping_add(0);
            let b = *src.wrapping_add(s);
            let c = *src.wrapping_add(2*s);
            let d = *src.wrapping_add(3*s);
            let apc = Arith::<P>::addmod(a, c);
            let amc = Arith::<P>::submod(a, c);
            let bpd = Arith::<P>::addmod(b, d);
            let bmd = Arith::<P>::submod(b, d);
            let jbmd = Arith::<P>::mmulmod(bmd, P.wrapping_sub(Self::U4));
            *dst.wrapping_add(0) = Arith::<P>::mmulmod_cond::<INV>(mult, Arith::<P>::addmodopt::<INV>(apc, bpd));
            *dst.wrapping_add(s) = Arith::<P>::mmulmod_cond::<INV>(mult, Arith::<P>::submod(amc, jbmd));
            *dst.wrapping_add(2*s) = Arith::<P>::mmulmod_cond::<INV>(mult, Arith::<P>::submod(apc,  bpd));
            *dst.wrapping_add(3*s) = Arith::<P>::mmulmod_cond::<INV>(mult, Arith::<P>::addmodopt::<INV>(amc, jbmd));
            src = src.wrapping_add(1);
            dst = dst.wrapping_add(1);
        }
        if eo { (n/4, s*4, !eo, y, x) } else { (n/4, s*4, eo, x, y) }
    }
}

impl<const P: u64, const INV: bool> NttKernelImpl<P, 5, INV> {
    unsafe fn apply<'a>(n: usize, s: usize, eo: bool, x: &'a mut [u64], y: &'a mut [u64]) -> (usize, usize, bool, &'a mut [u64], &'a mut [u64]) {
        let mut src = x.as_ptr();
        let mut dst = y.as_mut_ptr();
        let omega1 = Arith::<P>::mpowmod(Self::ROOTR, Arith::<P>::MAX_NTT_LEN/n as u64);
        let (n1, n1s) = (n/5, n/5*s);
        let (mut w1p, mut w2p, mut w3p, mut w4p) = (Arith::<P>::R, Arith::<P>::RNEG, Arith::<P>::RNEG, Arith::<P>::R);
        for _ in 0..n1 {
            for _ in 0..s {
                let a = *src.wrapping_add(0);
                let b = *src.wrapping_add(n1s);
                let c = *src.wrapping_add(2*n1s);
                let d = *src.wrapping_add(3*n1s);
                let e = *src.wrapping_add(4*n1s);
                let t1 = Arith::<P>::addmod(b, e);
                let t2 = Arith::<P>::addmod(c, d);
                let t3 = Arith::<P>::submod(b, e);
                let t4= Arith::<P>::submod(d, c);
                let t5 = Arith::<P>::addmod(t1, t2);
                let t6 = Arith::<P>::submod(t1, t2);
                let t7 = Arith::<P>::addmod64(t3, t4);
                let m1 = Arith::<P>::addmod(a, t5);
                let m2 = Arith::<P>::mmulsubmod(P.wrapping_sub(Self::C51), t5, m1);
                let m3 = Arith::<P>::mmulmod(Self::C52, t6);
                let m4 = Arith::<P>::mmulmod(Self::C53, t7);
                let m5 = Arith::<P>::mmulsubmod(Self::C54, t4, m4);
                let m6 = Arith::<P>::mmulsubmod(P.wrapping_sub(Self::C55), t3, m4);
                let s2 = Arith::<P>::submod(m3, m2);
                let s4 = Arith::<P>::addmod64(m2, m3);
                *dst.wrapping_add(0) = m1;
                *dst.wrapping_add(s) = Arith::<P>::mmulmod(w1p, Arith::<P>::submod(s2, m5));
                *dst.wrapping_add(2*s) = Arith::<P>::mmulmod(w2p, Arith::<P>::addmod64(s4, m6));
                *dst.wrapping_add(3*s) = Arith::<P>::mmulmod(w3p, Arith::<P>::submod(s4, m6));
                *dst.wrapping_add(4*s) = Arith::<P>::mmulmod(w4p, Arith::<P>::addmod64(s2, m5));
                src = src.wrapping_add(1);
                dst = dst.wrapping_add(1);
            }
            dst = dst.wrapping_add(4*s);
            w1p = Arith::<P>::mmulmod(w1p, omega1);
            w2p = Arith::<P>::mmulmod(w1p, P.wrapping_sub(w1p));
            w3p = Arith::<P>::mmulmod(w1p, w2p);
            w4p = Arith::<P>::mmulmod(w2p, w2p);
        }
        (n/5, s*5, !eo, y, x)
    }
    unsafe fn apply_last<'a>(n: usize, s: usize, eo: bool, x: &'a mut [u64], y: &'a mut [u64], mult: u64) -> (usize, usize, bool, &'a mut [u64], &'a mut [u64]) {
        assert_eq!(n, 5);
        let mut src = x.as_ptr();
        let mut dst = if eo { y.as_mut_ptr() } else { x.as_mut_ptr() };
        for _ in 0..s {
            let a = *src.wrapping_add(0);
            let b = *src.wrapping_add(s);
            let c = *src.wrapping_add(2*s);
            let d = *src.wrapping_add(3*s);
            let e = *src.wrapping_add(4*s);
            let t1 = Arith::<P>::addmod(b, e);
            let t2 = Arith::<P>::addmod(c, d);
            let t3 = Arith::<P>::submod(b, e);
            let t4= Arith::<P>::submod(d, c);
            let t5 = Arith::<P>::addmod(t1, t2);
            let t6 = Arith::<P>::submod(t1, t2);
            let t7 = Arith::<P>::addmod64(t3, t4);
            let m1 = Arith::<P>::addmod(a, t5);
            let m2 = Arith::<P>::mmuladdmod(Self::C51, t5, m1);
            let m3 = Arith::<P>::mmulmod(Self::C52, t6);
            let m4 = Arith::<P>::mmulmod(Self::C53, t7);
            let m5 = Arith::<P>::mmulsubmod(Self::C54, t4, m4);
            let m6 = Arith::<P>::mmulsubmod(P.wrapping_sub(Self::C55), t3, m4);
            let s2 = Arith::<P>::addmod(m2, m3);
            let s4 = Arith::<P>::submod(m2, m3);
            *dst.wrapping_add(0) = Arith::<P>::mmulmod_cond::<INV>(mult, m1);
            *dst.wrapping_add(s) = Arith::<P>::mmulmod_cond::<INV>(mult, Arith::<P>::submod(s2, m5));
            *dst.wrapping_add(2*s) = Arith::<P>::mmulmod_cond::<INV>(mult, Arith::<P>::submod(s4, m6));
            *dst.wrapping_add(3*s) = Arith::<P>::mmulmod_cond::<INV>(mult, Arith::<P>::addmodopt::<INV>(s4, m6));
            *dst.wrapping_add(4*s) = Arith::<P>::mmulmod_cond::<INV>(mult, Arith::<P>::addmodopt::<INV>(s2, m5));
            src = src.wrapping_add(1);
            dst = dst.wrapping_add(1);
        }
        if eo { (n/5, s*5, !eo, y, x) } else { (n/5, s*5, eo, x, y) }
    }
}

impl<const P: u64, const INV: bool> NttKernelImpl<P, 6, INV> {
    unsafe fn apply<'a>(n: usize, s: usize, eo: bool, x: &'a mut [u64], y: &'a mut [u64]) -> (usize, usize, bool, &'a mut [u64], &'a mut [u64]) {
        let mut src = x.as_ptr();
        let mut dst = y.as_mut_ptr();
        let omega1 = Arith::<P>::mpowmod(Self::ROOTR, Arith::<P>::MAX_NTT_LEN/n as u64);
        let (n1, n1s) = (n/6, n/6*s);
        let (mut w1p, mut w2p, mut w3p, mut w4p, mut w5p) = (Arith::<P>::R, Arith::<P>::R, Arith::<P>::R, Arith::<P>::R, Arith::<P>::R);
        for _ in 0..n1 {
            for _ in 0..s {
                let mut a = *src.wrapping_add(0);
                let mut b = *src.wrapping_add(n1s);
                let mut c = *src.wrapping_add(2*n1s);
                let mut d = *src.wrapping_add(3*n1s);
                let mut e = *src.wrapping_add(4*n1s);
                let mut f = *src.wrapping_add(5*n1s);
                (a, d) = (Arith::<P>::addmod(a, d), Arith::<P>::submod(a, d));
                (b, e) = (Arith::<P>::addmod(b, e), Arith::<P>::submod(b, e));
                (c, f) = (Arith::<P>::addmod(c, f), Arith::<P>::submod(c, f));
                let lbmc = Arith::<P>::mmulmod(Self::U6, Arith::<P>::submod(b, c));
                *dst.wrapping_add(0) = Arith::<P>::addmod(a, Arith::<P>::addmod(b, c));
                *dst.wrapping_add(2*s) = Arith::<P>::mmulmod(w2p, Arith::<P>::addmod64(Arith::<P>::submod(a, b), lbmc));
                *dst.wrapping_add(4*s) = Arith::<P>::mmulmod(w4p, Arith::<P>::submod(Arith::<P>::submod(a, c), lbmc));
                let mlepf = Arith::<P>::mmulmod(P.wrapping_sub(Self::U6), Arith::<P>::addmod64(e, f));
                *dst.wrapping_add(s) = Arith::<P>::mmulmod(w1p, Arith::<P>::submod(Arith::<P>::submod(d, f), mlepf));
                *dst.wrapping_add(3*s) = Arith::<P>::mmulmod(w3p, Arith::<P>::submod(d, Arith::<P>::submod(e, f)));
                *dst.wrapping_add(5*s) = Arith::<P>::mmulmod(w5p, Arith::<P>::addmod64(Arith::<P>::addmod64(d, mlepf), e));
                src = src.wrapping_add(1);
                dst = dst.wrapping_add(1);
            }
            dst = dst.wrapping_add(5*s);
            w1p = Arith::<P>::mmulmod(w1p, omega1);
            w2p = Arith::<P>::mmulmod(w1p, w1p);
            w3p = Arith::<P>::mmulmod(w1p, w2p);
            w4p = Arith::<P>::mmulmod(w2p, w2p);
            w5p = Arith::<P>::mmulmod(w2p, w3p);
        }
        (n/6, s*6, !eo, y, x)
    }
    unsafe fn apply_last<'a>(n: usize, s: usize, eo: bool, x: &'a mut [u64], y: &'a mut [u64], mult: u64) -> (usize, usize, bool, &'a mut [u64], &'a mut [u64]) {
        assert_eq!(n, 6);
        let mut src = x.as_ptr();
        let mut dst = if eo { y.as_mut_ptr() } else { x.as_mut_ptr() };
        for _ in 0..s {
            let mut a = *src.wrapping_add(0);
            let mut b = *src.wrapping_add(s);
            let mut c = *src.wrapping_add(2*s);
            let mut d = *src.wrapping_add(3*s);
            let mut e = *src.wrapping_add(4*s);
            let mut f = *src.wrapping_add(5*s);
            (a, d) = (Arith::<P>::addmod(a, d), Arith::<P>::submod(a, d));
            (b, e) = (Arith::<P>::addmod(b, e), Arith::<P>::submod(b, e));
            (c, f) = (Arith::<P>::addmod(c, f), Arith::<P>::submod(c, f));
            let lbmc = Arith::<P>::mmulmod(Self::U6, Arith::<P>::submod(b, c));
            *dst.wrapping_add(0) = Arith::<P>::mmulmod_cond::<INV>(mult, Arith::<P>::addmodopt::<INV>(a, Arith::<P>::addmodopt::<INV>(b, c)));
            *dst.wrapping_add(2*s) = Arith::<P>::mmulmod_cond::<INV>(mult, Arith::<P>::addmodopt::<INV>(Arith::<P>::submod(a, b), lbmc));
            *dst.wrapping_add(4*s) = Arith::<P>::mmulmod_cond::<INV>(mult, Arith::<P>::submod(Arith::<P>::submod(a, c), lbmc));
            let lepf = Arith::<P>::mmulmod(Self::U6, Arith::<P>::addmod64(e, f));
            *dst.wrapping_add(s) = Arith::<P>::mmulmod_cond::<INV>(mult, Arith::<P>::addmodopt::<INV>(Arith::<P>::submod(d, f), lepf));
            *dst.wrapping_add(3*s) = Arith::<P>::mmulmod_cond::<INV>(mult, Arith::<P>::submod(d, Arith::<P>::submod(e, f)));
            *dst.wrapping_add(5*s) = Arith::<P>::mmulmod_cond::<INV>(mult, Arith::<P>::addmodopt::<INV>(Arith::<P>::submod(d, lepf), e));
            src = src.wrapping_add(1);
            dst = dst.wrapping_add(1);
        }
        if eo { (n/6, s*6, !eo, y, x) } else { (n/6, s*6, eo, x, y) }
    }
}

fn ntt_stockham<const P: u64, const INV: bool>(input: &mut [u64], buf: &mut [u64]) {
    let (mut n, mut s, mut eo, mut x, mut y) = (input.len(), 1, false, input, buf);
    assert!(Arith::<P>::MAX_NTT_LEN % n as u64 == 0);
    let inv_p2 = Arith::<P>::mmulmod(Arith::<P>::R3, Arith::<P>::submod(0, (P-1)/n as u64));
    if n == 1 {
        x[0] = Arith::<P>::mmulmod_cond::<INV>(inv_p2, x[0]);
        return;
    }
    let (mut cnt6, mut cnt5, mut cnt4, mut cnt3, mut cnt2) = (0, 0, 0, 0, 0);
    let mut tmp = n;
    while tmp % 6 == 0 { tmp /= 6; cnt6 += 1; }
    while tmp % 5 == 0 { tmp /= 5; cnt5 += 1; }
    while tmp % 4 == 0 { tmp /= 4; cnt4 += 1; }
    while tmp % 3 == 0 { tmp /= 3; cnt3 += 1; }
    while tmp % 2 == 0 { tmp /= 2; cnt2 += 1; }
    while cnt6 > 0 && cnt2 > 0 { cnt6 -= 1; cnt2 -= 1; cnt4 += 1; cnt3 += 1; }
    unsafe {
        while cnt2 > 0 {
            (n, s, eo, x, y) = if n > 2 {
                NttKernelImpl::<P, 2, INV>::apply(n, s, eo, x, y)
            } else {
                NttKernelImpl::<P, 2, INV>::apply_last(n, s, eo, x, y, inv_p2)
            };
            cnt2 -= 1;
        }
        while cnt3 > 0 {
            (n, s, eo, x, y) = if n > 3 {
                NttKernelImpl::<P, 3, INV>::apply(n, s, eo, x, y)
            } else {
                NttKernelImpl::<P, 3, INV>::apply_last(n, s, eo, x, y, inv_p2)
            };
            cnt3 -= 1;
        }
        while cnt4 > 0 {
            (n, s, eo, x, y) = if n > 4 {
                NttKernelImpl::<P, 4, INV>::apply(n, s, eo, x, y)
            } else {
                NttKernelImpl::<P, 4, INV>::apply_last(n, s, eo, x, y, inv_p2)
            };
            cnt4 -= 1;
        }
        while cnt5 > 0 {
            (n, s, eo, x, y) = if n > 5 {
                NttKernelImpl::<P, 5, INV>::apply(n, s, eo, x, y)
            } else {
                NttKernelImpl::<P, 5, INV>::apply_last(n, s, eo, x, y, inv_p2)
            };
            cnt5 -= 1;
        }
        while cnt6 > 0 {
            (n, s, eo, x, y) = if n > 6 {
                NttKernelImpl::<P, 6, INV>::apply(n, s, eo, x, y)
            } else {
                NttKernelImpl::<P, 6, INV>::apply_last(n, s, eo, x, y, inv_p2)
            };
            cnt6 -= 1;
        }
    }
}

fn plan_ntt<const P: u64>(min_len: usize) -> (usize, usize) {
    assert!(min_len as u64 <= Arith::<P>::MAX_NTT_LEN);
    let (mut len_max, mut len_max_cost) = (0usize, usize::MAX);
    let mut len5 = 1;
    for _ in 0..=Arith::<P>::FACTOR_FIVE {
        let mut len35 = len5;
        for _ in 0..=Arith::<P>::FACTOR_THREE {
            let mut len = len35;
            let mut i = 0;
            while len < min_len && i < Arith::<P>::FACTOR_TWO { len *= 2; i += 1; }
            if len >= min_len && len < len_max_cost {
                let (mut tmp, mut cost) = (len, 0);
                while tmp % 6 == 0 { (tmp, cost) = (tmp/6, cost + len); }
                while tmp % 5 == 0 { (tmp, cost) = (tmp/5, cost + len + len/5); }
                while tmp % 4 == 0 { (tmp, cost) = (tmp/4, cost + len); }
                while tmp % 3 == 0 { (tmp, cost) = (tmp/3, cost + len); }
                while tmp % 2 == 0 { (tmp, cost) = (tmp/2, cost + len); }
                if cost < len_max_cost { (len_max, len_max_cost) = (len, cost); }
            }
            len35 *= 3;
        }
        len5 *= 5;
    }
    (len_max, len_max_cost)
}

// Performs (cyclic) integer convolution modulo P using NTT.
// Modifies the three buffers in-place.
// The output is saved in the slice `x`.
// The three slices must have the same length which divides `Arith::<P>::MAX_NTT_LEN`.
fn conv<const P: u64>(x: &mut [u64], y: &mut [u64], buf: &mut [u64]) {
    assert!(!x.is_empty() && x.len() == y.len() && y.len() == buf.len());
    ntt_stockham::<P, false>(x, buf);
    ntt_stockham::<P, false>(y, buf);
    for i in 0..x.len() { x[i] = Arith::<P>::mmulmod(x[i], y[i]); }
    ntt_stockham::<P, true>(x, buf);
}

////////////////////////////////////////////////////////////////////////////////

use core::cmp::{min, max};
use crate::big_digit::BigDigit;

const P1: u64 = 10_237_243_632_176_332_801; // Max NTT length = 2^24 * 3^20 * 5^2 = 1_462_463_376_025_190_400
const P2: u64 = 13_649_658_176_235_110_401; // Max NTT length = 2^26 * 3^19 * 5^2 = 1_949_951_168_033_587_200
const P3: u64 = 14_259_017_916_245_606_401; // Max NTT length = 2^22 * 3^21 * 5^2 = 1_096_847_532_018_892_800

const P1INV_R_MOD_P2: u64 = Arith::<P2>::mmulmod(Arith::<P2>::R2, arith::invmod(P1, P2));
const P1P2INV_R_MOD_P3: u64 = Arith::<P3>::mmulmod(
    Arith::<P3>::R3,
    Arith::<P3>::mmulmod(
        arith::invmod(P1, P3),
        arith::invmod(P2, P3)
    )
);
const P1_R_MOD_P3: u64 = Arith::<P3>::mmulmod(Arith::<P3>::R2, P1);
const P1P2_LO: u64 = (P1 as u128 * P2 as u128) as u64;
const P1P2_HI: u64 = ((P1 as u128 * P2 as u128) >> 64) as u64;

fn mac3_two_primes(acc: &mut [u64], b: &[u64], c: &[u64], bits: u64) {
    let min_len = b.len() + c.len();
    let len_max_1 = plan_ntt::<P1>(min_len).0;
    let len_max_2 = plan_ntt::<P2>(min_len).0;
    let len_max = max(len_max_1, len_max_2);
    let mut x = vec![0u64; len_max_1];
    let mut y = vec![0u64; len_max_2];
    let mut r = vec![0u64; len_max];
    let mut s = vec![0u64; len_max];

    /* convolution with modulo P1 */
    for i in 0..b.len() { x[i] = if b[i] >= P1 { b[i] - P1 } else { b[i] }; }
    for i in 0..c.len() { r[i] = if c[i] >= P1 { c[i] - P1 } else { c[i] }; }
    r[c.len()..len_max_1].fill(0u64);
    conv::<P1>(&mut x, &mut r[..len_max_1], &mut s[..len_max_1]);

    /* convolution with modulo P2 */
    for i in 0..b.len() { y[i] = if b[i] >= P2 { b[i] - P2 } else { b[i] }; }
    for i in 0..c.len() { r[i] = if c[i] >= P2 { c[i] - P2 } else { c[i] }; }
    r[c.len()..len_max_2].fill(0u64);
    conv::<P2>(&mut y, &mut r[..len_max_2], &mut s[..len_max_2]);

    /* merge the results in {x, y} into r (process carry along the way) */
    let mask = (1u64 << bits) - 1;
    let mut carry: u128 = 0;
    let (mut j, mut p) = (0usize, 0u64);
    for i in 0..min_len {
        /* extract the convolution result */
        let (a, b) = (x[i], y[i]);
        let bma = Arith::<P2>::submod(b, a);
        let u = Arith::<P2>::mmulmod(bma, P1INV_R_MOD_P2);
        let v = a as u128 + P1 as u128 * u as u128 + carry;
        carry = v >> bits;

        /* write to r */
        let out = (v as u64) & mask;
        r[j] = (r[j] & ((1u64 << p) - 1)) | (out << p);
        p += bits;
        if p >= 64 {
            (j, p) = (j+1, p-64);
            r[j] = out >> (bits - p);
        }
    }

    /* add r to acc */
    let mut carry: u64 = 0;
    for i in 0..min(acc.len(), j+1) {
        let w = r[i];
        let (v, overflow1) = acc[i].overflowing_add(w);
        let (v, overflow2) = v.overflowing_add(carry);
        acc[i] = v;
        carry = u64::from(overflow1 || overflow2);
    }
}

fn mac3_three_primes(acc: &mut [u64], b: &[u64], c: &[u64]) {
    let min_len = b.len() + c.len();
    let len_max_1 = plan_ntt::<P1>(min_len).0;
    let len_max_2 = plan_ntt::<P2>(min_len).0;
    let len_max_3 = plan_ntt::<P3>(min_len).0;
    let len_max = max(len_max_1, max(len_max_2, len_max_3));
    let mut x = vec![0u64; len_max_1];
    let mut y = vec![0u64; len_max_2];
    let mut z = vec![0u64; len_max_3];
    let mut r = vec![0u64; len_max];
    let mut s = vec![0u64; len_max];

    /* convolution with modulo P1 */
    for i in 0..b.len() { x[i] = if b[i] >= P1 { b[i] - P1 } else { b[i] }; }
    for i in 0..c.len() { r[i] = if c[i] >= P1 { c[i] - P1 } else { c[i] }; }
    r[c.len()..len_max_1].fill(0u64);
    conv::<P1>(&mut x, &mut r[..len_max_1], &mut s[..len_max_1]);

    /* convolution with modulo P2 */
    for i in 0..b.len() { y[i] = if b[i] >= P2 { b[i] - P2 } else { b[i] }; }
    for i in 0..c.len() { r[i] = if c[i] >= P2 { c[i] - P2 } else { c[i] }; }
    r[c.len()..len_max_2].fill(0u64);
    conv::<P2>(&mut y, &mut r[..len_max_2], &mut s[..len_max_2]);

    /* convolution with modulo P3 */
    for i in 0..b.len() { z[i] = if b[i] >= P3 { b[i] - P3 } else { b[i] }; }
    for i in 0..c.len() { r[i] = if c[i] >= P3 { c[i] - P3 } else { c[i] }; }
    r[c.len()..len_max_3].fill(0u64);
    conv::<P3>(&mut z, &mut r[..len_max_3], &mut s[..len_max_3]);

    /* merge the results in {x, y, z} into acc (process carry along the way) */
    let mut carry: u128 = 0;
    for i in 0..min_len {
        let (a, b, c) = (x[i], y[i], z[i]);
        // We need to solve the following system of linear congruences:
        //     x === a mod P1,
        //     x === b mod P2,
        //     x === c mod P3.
        // The first two equations are equivalent to
        //     x === a + P1 * (U * (b-a) mod P2) mod P1P2,
        // where U is the solution to
        //     P1 * U === 1 mod P2.
        let bma = Arith::<P2>::submod(b, a);
        let u = Arith::<P2>::mmulmod(bma, P1INV_R_MOD_P2);
        let v = a as u128 + P1 as u128 * u as u128;
        let v_mod_p3 = Arith::<P3>::addmod(a, Arith::<P3>::mmulmod(P1_R_MOD_P3, u));
        // Now we have reduced the congruences into two:
        //     x === v mod P1P2,
        //     x === c mod P3.
        // The solution is
        //     x === v + P1P2 * (V * (c-v) mod P3) mod P1P2P3,
        // where V is the solution to
        //     P1P2 * V === 1 mod P3.
        let cmv = Arith::<P3>::submod(c, v_mod_p3);
        let vcmv = Arith::<P3>::mmulmod(cmv, P1P2INV_R_MOD_P3);
        let (out_01, overflow) = carry.overflowing_add(v + P1P2_LO as u128 * vcmv as u128);
        let out_0 = out_01 as u64;
        let out_12 = P1P2_HI as u128 * vcmv as u128 + (out_01 >> 64) +
            if overflow { 1u128 << 64 } else { 0 };
        let out_1 = out_12 as u64;
        let out_2 = (out_12 >> 64) as u64;

        let (v, overflow) = acc[i].overflowing_add(out_0);
        acc[i] = v;
        carry = out_1 as u128 + ((out_2 as u128) << 64) + u128::from(overflow);
    }
    let mut carry = carry as u64;
    for i in min_len..acc.len() {
        let (v, overflow) = acc[i].overflowing_add(carry);
        acc[i] = v;
        carry = u64::from(overflow);
    }
}

fn mac3_u64(acc: &mut [u64], bb: &[u64], cc: &[u64], split_unbalanced: bool) {
    let (b, c) = if bb.len() < cc.len() { (bb, cc) } else { (cc, bb) };
    if split_unbalanced && b.len() * 2 <= c.len() {
        /* special handling for unbalanced multiplication:
           we reduce it to about `c.len()/b.len()` balanced multiplications */
        let mut i = 0usize;
        let mut carry = 0u64;
        while i < c.len() {
            let j = min(i + b.len(), c.len());
            let k = j + b.len();
            let tmp = acc[k];
            acc[k] = 0;
            mac3_u64(&mut acc[i..=k], b, &c[i..j], false);
            let mut l = j;
            while carry > 0 && l < k {
                let (v, overflow) = acc[l].overflowing_add(carry);
                acc[l] = v;
                carry = u64::from(overflow);
                l += 1;
            }
            i = j;
            carry += tmp;
        }
        i += b.len();
        while i < acc.len() {
            let (v, overflow) = acc[i].overflowing_add(carry);
            acc[i] = v;
            carry = u64::from(overflow);
            i += 1;
        }
        return;
    }

    let max_cnt = max(b.len(), c.len()) as u64;
    let mut bits = 0u64;
    while 1u64 << (2*bits) < max_cnt { bits += 1; }
    bits = 63 - bits;
    if bits >= 44 {
        /* can pack more effective bits per u64 with two primes than with three primes */
        fn pack_into(src: &[u64], dst: &mut [u64], bits: u64) -> usize {
            let (mut j, mut p) = (0usize, 0u64);
            for i in 0..src.len() {
                let mut k = 0;
                while k < 64 {
                    let bits_this_time = min(64 - k, bits - p);
                    dst[j] = (dst[j] & ((1u64 << p) - 1)) | (((src[i] >> k) & ((1u64 << bits_this_time) - 1)) << p);
                    k += bits_this_time;
                    p += bits_this_time;
                    if p == bits { (j, p) = (j+1, 0); }
                }
            }
            if p == 0 { j } else { j+1 }
        }
        let mut b2 = vec![0u64; ((64 * b.len() as u64 + bits - 1) / bits) as usize];
        let mut c2 = vec![0u64; ((64 * c.len() as u64 + bits - 1) / bits) as usize];
        let b2_len = pack_into(b, &mut b2, bits);
        let c2_len = pack_into(c, &mut c2, bits);
        mac3_two_primes(acc, &b2[..b2_len], &c2[..c2_len], bits);
    } else {
        /* can pack at most 21 effective bits per u64, which is worse than
           64/3 = 21.3333.. effective bits per u64 achieved with three primes */
        mac3_three_primes(acc, b, c);
    }
}

////////////////////////////////////////////////////////////////////////////////

#[cfg(u64_digit)]
pub fn mac3(acc: &mut [BigDigit], b: &[BigDigit], c: &[BigDigit]) {
    mac3_u64(acc, b, c, true);
}

#[cfg(not(u64_digit))]
pub fn mac3(acc: &mut [BigDigit], b: &[BigDigit], c: &[BigDigit]) {
    fn bigdigit_to_u64(src: &[BigDigit]) -> crate::biguint::Vec::<u64> {
        let mut out = vec![0u64; (src.len() + 1) / 2];
        for i in 0..src.len()/2 {
            out[i] = (src[2*i] as u64) | ((src[2*i+1] as u64) << 32);
        }
        if src.len() % 2 == 1 {
            out[src.len()/2] = src[src.len()-1] as u64;
        }
        out
    }
    fn u64_to_bigdigit(src: &[u64], dst: &mut [BigDigit]) {
        for i in 0..dst.len()/2 {
            dst[2*i] = src[i] as BigDigit;
            dst[2*i+1] = (src[i] >> 32) as BigDigit;
        }
        if dst.len() % 2 == 1 {
            dst[dst.len()-1] = src[src.len()-1] as BigDigit;
        }
    }

    /* convert to u64 => process => convert back to BigDigit (u32) */
    let mut acc_u64 = bigdigit_to_u64(acc);
    let b_u64 = bigdigit_to_u64(b);
    let c_u64 = bigdigit_to_u64(c);
    mac3_u64(&mut acc_u64, &b_u64, &c_u64, true);
    u64_to_bigdigit(&acc_u64, acc);
}