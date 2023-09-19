#![allow(clippy::cast_sign_loss)]
#![allow(clippy::cast_lossless)]
#![allow(clippy::cast_possible_truncation)]
#![allow(clippy::many_single_char_names)]
#![allow(clippy::needless_range_loop)]
#![allow(clippy::too_many_arguments)]
#![allow(clippy::similar_names)]

use crate::biguint::Vec;

mod arith {
    // Extended Euclid algorithm:
    //   (g, x, y) is a solution to ax + by = g, where g = gcd(a, b)
    const fn egcd(mut a: i128, mut b: i128) -> (i128, i128, i128) {
        assert!(a > 0 && b > 0);
        let mut c = if a > b { (a, b) = (b, a); [0, 1, 1, 0] } else { [1, 0, 0, 1] }; // treat as a row-major 2x2 matrix
        loop {
            if a == 0 { break (b, c[1], c[3]) }
            let (q, r) = (b/a, b%a);
            (a, b) = (r, a);
            c = [c[1] - q*c[0], c[0], c[3] - q*c[2], c[2]];
        }
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
    pub const FACTOR_TWO: u32 = (P-1).trailing_zeros();
    pub const FACTOR_THREE: u32 = Self::factor_three();
    pub const FACTOR_FIVE: u32 = Self::factor_five();
    pub const MAX_NTT_LEN: u64 = Self::max_ntt_len();
    pub const R: u64 = ((1u128 << 64) % P as u128) as u64; // 2^64 mod P
    pub const R2: u64 = (Self::R as u128 * Self::R as u128 % P as u128) as u64; // R^2 mod P
    pub const R3: u64 = (Self::R2 as u128 * Self::R as u128 % P as u128) as u64; // R^3 mod P
    pub const PINV: u64 = arith::invmod(P, 0); // P^-1 mod 2^64
    pub const ROOT: u64 = Self::ntt_root(); // MultiplicativeOrder[ROOT, P] == MAX_NTT_LEN
    pub const ROOTR: u64 = (Self::ROOT as u128 * Self::R as u128 % P as u128) as u64; // ROOT * R mod P
    const fn factor_three() -> u32 {
        let mut tmp = P-1;
        let mut ans = 0;
        while tmp % 3 == 0 { tmp /= 3; ans += 1; }
        ans
    }
    const fn factor_five() -> u32 {
        let mut tmp = P-1;
        let mut ans = 0;
        while tmp % 5 == 0 { tmp /= 5; ans += 1; }
        ans
    }
    const fn max_ntt_len() -> u64 {
        let ans = 2u64.pow(Self::FACTOR_TWO) * 3u64.pow(Self::FACTOR_THREE) * 5u64.pow(Self::FACTOR_FIVE);
        assert!(ans % 4050 == 0);
        ans
    }
    const fn ntt_root() -> u64 {
        let mut p = 2;
        'outer: loop {
            let root = Self::powmod_naive(p, P/Self::MAX_NTT_LEN);
            let mut j = 0;
            while j <= Self::FACTOR_TWO {
                let mut k = 0;
                while k <= Self::FACTOR_THREE {
                    let mut l = 0;
                    while l <= Self::FACTOR_FIVE {
                        let exponent = 2u64.pow(j) * 3u64.pow(k) * 5u64.pow(l);
                        if exponent < Self::MAX_NTT_LEN && Self::powmod_naive(root, exponent) == 1 {
                            p += 1;
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
    // Montgomery reduction:
    //   x * R^-1 mod P
    pub const fn mreduce(x: u128) -> u64 {
        let m = (x as u64).wrapping_mul(Self::PINV);
        let y = ((m as u128 * P as u128) >> 64) as u64;
        let (out, overflow) = ((x >> 64) as u64).overflowing_sub(y);
        if overflow { out.wrapping_add(P) } else { out }
    }
    // Multiplication with Montgomery reduction:
    //   a * b * R^-1 mod P
    pub const fn mmulmod(a: u64, b: u64) -> u64 {
        Self::mreduce(a as u128 * b as u128)
    }
    // Multiplication with Montgomery reduction:
    //   a * b * R^-1 mod P
    // This function only applies the multiplication when INV && TWIDDLE,
    //   otherwise it just returns b.
    pub const fn mmulmod_invtw<const INV: bool, const TWIDDLE: bool>(a: u64, b: u64) -> u64 {
        if INV && TWIDDLE { Self::mmulmod(a, b) } else { b }
    }
    // Fused-multiply-sub with Montgomery reduction:
    //   a * b * R^-1 - c mod P
    pub const fn mmulsubmod(a: u64, b: u64, c: u64) -> u64 {
        let x = a as u128 * b as u128;
        let lo = x as u64;
        let hi = Self::submod((x >> 64) as u64, c);
        Self::mreduce(lo as u128 | ((hi as u128) << 64))
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
    // Computes c as u128 * mreduce(v) as u128,
    //   using d: u64 = mmulmod(P-1, c).
    // It is caller's responsibility to ensure that d is correct.
    // Note that d can be computed by calling mreducelo(c).
    pub const fn mmulmod_noreduce(v: u128, c: u64, d: u64) -> u128 {
        let a: u128 = c as u128 * (v >> 64);
        let b: u128 = d as u128 * (v as u64 as u128);
        let (w, overflow) = a.overflowing_sub(b);
        if overflow { w.wrapping_add((P as u128) << 64) } else { w }
    }
    // Computes submod(0, mreduce(x as u128)) fast.
    pub const fn mreducelo(x: u64) -> u64 {
        let m = x.wrapping_mul(Self::PINV);
        ((m as u128 * P as u128) >> 64) as u64
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
    // Computes a + b mod P, selects addmod64 or addmod depending on INV && TWIDDLE
    pub const fn addmodopt_invtw<const INV: bool, const TWIDDLE: bool>(a: u64, b: u64) -> u64 {
        if INV && TWIDDLE { Self::addmod64(a, b) } else { Self::addmod(a, b) }
    }
    // Computes a - b mod P, output range [0, P)
    pub const fn submod(a: u64, b: u64) -> u64 {
        let (out, overflow) = a.overflowing_sub(b);
        if overflow { out.wrapping_add(P) } else { out }
    }
}

struct NttPlan {
    pub n: usize,   // n == g*m
    pub g: usize,   // g <= NttPlan::GMAX
    pub m: usize,   // m divides Arith::<P>::MAX_NTT_LEN
    pub cost: usize,
    pub last_radix: usize,
    pub s_list: Vec<(usize, usize)>,
}
impl NttPlan {
    pub fn build<const P: u64>(min_len: usize) -> Self {
        assert!(min_len as u64 <= Arith::<P>::MAX_NTT_LEN);
        let (mut len_max, mut len_max_cost, mut g) = (usize::MAX, usize::MAX, 1);
        for m7 in 0..=1 {
            for m5 in 0..=Arith::<P>::FACTOR_FIVE {
                for m3 in 0..=Arith::<P>::FACTOR_THREE {
                    let len = 7u64.pow(m7) * 5u64.pow(m5) * 3u64.pow(m3);
                    if len >= 2 * min_len as u64 { continue; }
                    let (mut len, mut m2) = (len as usize, 0);
                    while len < min_len && m2 < Arith::<P>::FACTOR_TWO { len *= 2; m2 += 1; }
                    if len >= min_len && len < len_max_cost {
                        let (mut tmp, mut cost) = (len, 0);
                        let mut g_new = 1;
                        if len % 7 == 0 {
                            (g_new, tmp, cost) = (7, tmp/7, cost + len*115/100);
                        } else if len % 5 == 0 {
                            (g_new, tmp, cost) = (5, tmp/5, cost + len*89/100);
                        } else if m2 >= m3 + 3 {
                            (g_new, tmp, cost) = (8, tmp/8, cost + len*130/100);
                        } else if m2 >= m3 + 2 {
                            (g_new, tmp, cost) = (4, tmp/4, cost + len*87/100);
                        } else if m2 == 0 && m3 >= 1 {
                            (g_new, tmp, cost) = (3, tmp/3, cost + len*86/100);
                        } else if m3 == 0 && m2 >= 1 {
                            (g_new, tmp, cost) = (2, tmp/2, cost + len*86/100);
                        } else if len % 6 == 0 {
                            (g_new, tmp, cost) = (6, tmp/6, cost + len*91/100);
                        }
                        let (mut b6, mut b2) = (false, false);
                        while tmp % 6 == 0 { (tmp, cost) = (tmp/6, cost + len*106/100); b6 = true; }
                        while tmp % 5 == 0 { (tmp, cost) = (tmp/5, cost + len*131/100); }
                        while tmp % 4 == 0 { (tmp, cost) = (tmp/4, cost + len); }
                        while tmp % 3 == 0 { (tmp, cost) = (tmp/3, cost + len); }
                        while tmp % 2 == 0 { (tmp, cost) = (tmp/2, cost + len); b2 = true; }
                        if b6 && b2 { cost -= len*6/100; }
                        if cost < len_max_cost { (len_max, len_max_cost, g) = (len, cost, g_new); }
                    }
                }
            }
        }
        let (mut cnt6, mut cnt5, mut cnt4, mut cnt3, mut cnt2) = (0, 0, 0, 0, 0);
        let mut tmp = len_max/g;
        while tmp % 6 == 0 { tmp /= 6; cnt6 += 1; }
        while tmp % 5 == 0 { tmp /= 5; cnt5 += 1; }
        while tmp % 4 == 0 { tmp /= 4; cnt4 += 1; }
        while tmp % 3 == 0 { tmp /= 3; cnt3 += 1; }
        while tmp % 2 == 0 { tmp /= 2; cnt2 += 1; }
        while cnt6 > 0 && cnt2 > 0 { cnt6 -= 1; cnt2 -= 1; cnt4 += 1; cnt3 += 1; }
        let s_list = {
            let mut out = vec![];
            let mut tmp = len_max;
            for _ in 0..cnt2 { out.push((tmp, 2)); tmp /= 2; }
            for _ in 0..cnt3 { out.push((tmp, 3)); tmp /= 3; }
            for _ in 0..cnt4 { out.push((tmp, 4)); tmp /= 4; }
            for _ in 0..cnt5 { out.push((tmp, 5)); tmp /= 5; }
            for _ in 0..cnt6 { out.push((tmp, 6)); tmp /= 6; }
            out
        };
        Self {
            n: len_max,
            g,
            m: len_max / g,
            cost: len_max_cost,
            last_radix: s_list.last().unwrap_or(&(1, 1)).1,
            s_list,
        }
    }
}
fn conv_base<const P: u64>(n: usize, x: *mut u64, y: *mut u64, c: u64) {
    unsafe {
        let c2 = Arith::<P>::mreducelo(c);
        let out = x.wrapping_sub(n);
        for i in 0..n {
            let mut v: u128 = 0;
            for j in i+1..n {
                let (w, overflow) = v.overflowing_sub(*x.wrapping_add(j) as u128 * *y.wrapping_add(i+n-j) as u128);
                v = if overflow { w.wrapping_add((P as u128) << 64) } else { w };
            }
            v = Arith::<P>::mmulmod_noreduce(v, c, c2);
            for j in 0..=i {
                let (w, overflow) = v.overflowing_sub(*x.wrapping_add(j) as u128 * *y.wrapping_add(i-j) as u128);
                v = if overflow { w.wrapping_add((P as u128) << 64) } else { w };
            }
            *out.wrapping_add(i) = Arith::<P>::mreduce(v);
        }
    }
}

struct NttKernelImpl<const P: u64, const INV: bool>;
impl<const P: u64, const INV: bool> NttKernelImpl<P, INV> {
    pub const ROOTR: u64 = Arith::<P>::mpowmod(Arith::<P>::ROOTR, if INV { Arith::<P>::MAX_NTT_LEN - 1 } else { 1 });
    pub const U2: u64 = Arith::<P>::mpowmod(Self::ROOTR, Arith::<P>::MAX_NTT_LEN/2); // U2 == P - Arith::<P>::R
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
const fn ntt2_kernel<const P: u64, const INV: bool, const TWIDDLE: bool>(
    w1: u64,
    a: u64, mut b: u64) -> (u64, u64) {
    if !INV && TWIDDLE {
        b = Arith::<P>::mmulmod(w1, b);
    }
    let out0 = Arith::<P>::addmod(a, b);
    let out1 = Arith::<P>::mmulmod_invtw::<INV, TWIDDLE>(w1, Arith::<P>::submod(a, b));
    (out0, out1)
}
fn ntt2_single_block<const P: u64, const INV: bool, const TWIDDLE: bool>(
    s1: usize, mut px: *mut u64, ptf: *const u64) -> (*mut u64, *const u64) {
    unsafe {
        let w1 = if TWIDDLE { *ptf } else { 0 };
        for _ in 0..s1 {
            (*px, *px.wrapping_add(s1)) =
                ntt2_kernel::<P, INV, TWIDDLE>(w1,
                    *px, *px.wrapping_add(s1));
            px = px.wrapping_add(1);
        }
    }
    (px.wrapping_add(s1), ptf.wrapping_add(1))
}
const fn ntt3_kernel<const P: u64, const INV: bool, const TWIDDLE: bool>(
    w1: u64, w2: u64,
    a: u64, mut b: u64, mut c: u64) -> (u64, u64, u64) {
    if !INV && TWIDDLE {
        b = Arith::<P>::mmulmod(w1, b);
        c = Arith::<P>::mmulmod(w2, c);
    }
    let kbmc = Arith::<P>::mmulmod(NttKernelImpl::<P, INV>::U3, Arith::<P>::submod(b, c));
    let out0 = Arith::<P>::addmod(a, Arith::<P>::addmod(b, c));
    let out1 = Arith::<P>::mmulmod_invtw::<INV, TWIDDLE>(w1, Arith::<P>::submod(a, Arith::<P>::submod(c, kbmc)));
    let out2 = Arith::<P>::mmulmod_invtw::<INV, TWIDDLE>(w2, Arith::<P>::submod(Arith::<P>::submod(a, b), kbmc));
    (out0, out1, out2)
}
fn ntt3_single_block<const P: u64, const INV: bool, const TWIDDLE: bool>(
    s1: usize, mut px: *mut u64, ptf: *const u64) -> (*mut u64, *const u64) {
    unsafe {
        let (w1, w2) = if TWIDDLE {
            let w1 = *ptf;
            let w2 = Arith::<P>::mmulmod(w1, w1);
            (w1, w2)
        } else {
            (0, 0)
        };
        for _ in 0..s1 {
            (*px, *px.wrapping_add(s1), *px.wrapping_add(2*s1)) =
                ntt3_kernel::<P, INV, TWIDDLE>(w1, w2,
                    *px, *px.wrapping_add(s1), *px.wrapping_add(2*s1));
            px = px.wrapping_add(1);
        }
    }
    (px.wrapping_add(2*s1), ptf.wrapping_add(1))
}
const fn ntt4_kernel<const P: u64, const INV: bool, const TWIDDLE: bool>(
    w1: u64, w2: u64, w3: u64,
    a: u64, mut b: u64, mut c: u64, mut d: u64) -> (u64, u64, u64, u64) {
    if !INV && TWIDDLE {
        b = Arith::<P>::mmulmod(w1, b);
        c = Arith::<P>::mmulmod(w2, c);
        d = Arith::<P>::mmulmod(w3, d);
    }
    let apc = Arith::<P>::addmod(a, c);
    let amc = Arith::<P>::submod(a, c);
    let bpd = Arith::<P>::addmod(b, d);
    let bmd = Arith::<P>::submod(b, d);
    let jbmd = Arith::<P>::mmulmod(NttKernelImpl::<P, INV>::U4, bmd);
    let out0 = Arith::<P>::addmod(apc, bpd);
    let out1 = Arith::<P>::mmulmod_invtw::<INV, TWIDDLE>(w1, Arith::<P>::addmodopt_invtw::<INV, TWIDDLE>(amc, jbmd));
    let out2 = Arith::<P>::mmulmod_invtw::<INV, TWIDDLE>(w2, Arith::<P>::submod(apc,  bpd));
    let out3 = Arith::<P>::mmulmod_invtw::<INV, TWIDDLE>(w3, Arith::<P>::submod(amc, jbmd));
    (out0, out1, out2, out3)
}
fn ntt4_single_block<const P: u64, const INV: bool, const TWIDDLE: bool>(
    s1: usize, mut px: *mut u64, ptf: *const u64) -> (*mut u64, *const u64) {
    unsafe {
        let (w1, w2, w3) = if TWIDDLE {
            let w1 = *ptf;
            let w2 = Arith::<P>::mmulmod(w1, w1);
            let w3 = Arith::<P>::mmulmod(w1, w2);
            (w1, w2, w3)
        } else {
            (0, 0, 0)
        };
        for _ in 0..s1 {
            (*px, *px.wrapping_add(s1), *px.wrapping_add(2*s1),
            *px.wrapping_add(3*s1)) =
                ntt4_kernel::<P, INV, TWIDDLE>(w1, w2, w3,
                    *px, *px.wrapping_add(s1), *px.wrapping_add(2*s1),
                    *px.wrapping_add(3*s1));
            px = px.wrapping_add(1);
        }
    }
    (px.wrapping_add(3*s1), ptf.wrapping_add(1))
}
const fn ntt5_kernel<const P: u64, const INV: bool, const TWIDDLE: bool>(
    w1: u64, w2: u64, w3: u64, w4: u64,
    a: u64, mut b: u64, mut c: u64, mut d: u64, mut e: u64) -> (u64, u64, u64, u64, u64) {
    if !INV && TWIDDLE {
        b = Arith::<P>::mmulmod(w1, b);
        c = Arith::<P>::mmulmod(w2, c);
        d = Arith::<P>::mmulmod(w3, d);
        e = Arith::<P>::mmulmod(w4, e);
    }
    let t1 = Arith::<P>::addmod(b, e);
    let t2 = Arith::<P>::addmod(c, d);
    let t3 = Arith::<P>::submod(b, e);
    let t4 = Arith::<P>::submod(d, c);
    let t5 = Arith::<P>::addmod(t1, t2);
    let t6 = Arith::<P>::submod(t1, t2);
    let t7 = Arith::<P>::addmod64(t3, t4);
    let m1 = Arith::<P>::addmod(a, t5);
    let m2 = Arith::<P>::mmulsubmod(P.wrapping_sub(NttKernelImpl::<P, INV>::C51), t5, m1);
    let m3 = Arith::<P>::mmulmod(NttKernelImpl::<P, INV>::C52, t6);
    let m4 = Arith::<P>::mmulmod(NttKernelImpl::<P, INV>::C53, t7);
    let m5 = Arith::<P>::mmulsubmod(NttKernelImpl::<P, INV>::C54, t4, m4);
    let m6 = Arith::<P>::mmulsubmod(P.wrapping_sub(NttKernelImpl::<P, INV>::C55), t3, m4);
    let s1 = Arith::<P>::submod(m3, m2);
    let s2 = Arith::<P>::addmod(m2, m3);
    let out0 = m1;
    let out1 = Arith::<P>::mmulmod_invtw::<INV, TWIDDLE>(w1, Arith::<P>::submod(s1, m5));
    let out2 = Arith::<P>::mmulmod_invtw::<INV, TWIDDLE>(w2, Arith::<P>::submod(0, Arith::<P>::addmod(s2, m6)));
    let out3 = Arith::<P>::mmulmod_invtw::<INV, TWIDDLE>(w3, Arith::<P>::submod(m6, s2));
    let out4 = Arith::<P>::mmulmod_invtw::<INV, TWIDDLE>(w4, Arith::<P>::addmodopt_invtw::<INV, TWIDDLE>(s1, m5));
    (out0, out1, out2, out3, out4)
}
fn ntt5_single_block<const P: u64, const INV: bool, const TWIDDLE: bool>(
    s1: usize, mut px: *mut u64, ptf: *const u64) -> (*mut u64, *const u64) {
    unsafe {
        let (w1, w2, w3, w4) = if TWIDDLE {
            let w1 = *ptf;
            let w2 = Arith::<P>::mmulmod(w1, w1);
            let w3 = Arith::<P>::mmulmod(w1, w2);
            let w4 = Arith::<P>::mmulmod(w2, w2);
            (w1, w2, w3, w4)
        } else {
            (0, 0, 0, 0)
        };
        for _ in 0..s1 {
            (*px, *px.wrapping_add(s1), *px.wrapping_add(2*s1),
            *px.wrapping_add(3*s1), *px.wrapping_add(4*s1)) =
                ntt5_kernel::<P, INV, TWIDDLE>(w1, w2, w3, w4,
                    *px, *px.wrapping_add(s1), *px.wrapping_add(2*s1),
                    *px.wrapping_add(3*s1), *px.wrapping_add(4*s1));
            px = px.wrapping_add(1);
        }
    }
    (px.wrapping_add(4*s1), ptf.wrapping_add(1))
}
const fn ntt6_kernel<const P: u64, const INV: bool, const TWIDDLE: bool>(
    w1: u64, w2: u64, w3: u64, w4: u64, w5: u64,
    mut a: u64, mut b: u64, mut c: u64, mut d: u64, mut e: u64, mut f: u64) -> (u64, u64, u64, u64, u64, u64) {
    if !INV && TWIDDLE {
        b = Arith::<P>::mmulmod(w1, b);
        c = Arith::<P>::mmulmod(w2, c);
        d = Arith::<P>::mmulmod(w3, d);
        e = Arith::<P>::mmulmod(w4, e);
        f = Arith::<P>::mmulmod(w5, f);
    }
    (a, d) = (Arith::<P>::addmod(a, d), Arith::<P>::submod(a, d));
    (b, e) = (Arith::<P>::addmod(b, e), Arith::<P>::submod(b, e));
    (c, f) = (Arith::<P>::addmod(c, f), Arith::<P>::submod(c, f));
    let lbmc = Arith::<P>::mmulmod(NttKernelImpl::<P, INV>::U6, Arith::<P>::submod(b, c));
    let out0 = Arith::<P>::addmod(a, Arith::<P>::addmod(b, c));
    let out2 = Arith::<P>::mmulmod_invtw::<INV, TWIDDLE>(w2, Arith::<P>::submod(a, Arith::<P>::submod(b, lbmc)));
    let out4 = Arith::<P>::mmulmod_invtw::<INV, TWIDDLE>(w4, Arith::<P>::submod(Arith::<P>::submod(a, c), lbmc));
    let lepf = Arith::<P>::mmulmod(NttKernelImpl::<P, INV>::U6, Arith::<P>::addmod64(e, f));
    let out1 = Arith::<P>::mmulmod_invtw::<INV, TWIDDLE>(w1, Arith::<P>::submod(d, Arith::<P>::submod(f, lepf)));
    let out3 = Arith::<P>::mmulmod_invtw::<INV, TWIDDLE>(w3, Arith::<P>::submod(d, Arith::<P>::submod(e, f)));
    let out5 = Arith::<P>::mmulmod_invtw::<INV, TWIDDLE>(w5, Arith::<P>::submod(d, Arith::<P>::submod(lepf, e)));
    (out0, out1, out2, out3, out4, out5)
}
fn ntt6_single_block<const P: u64, const INV: bool, const TWIDDLE: bool>(
    s1: usize, mut px: *mut u64, ptf: *const u64) -> (*mut u64, *const u64) {
    unsafe {
        let (w1, w2, w3, w4, w5) = if TWIDDLE {
            let w1 = *ptf;
            let w2 = Arith::<P>::mmulmod(w1, w1);
            let w3 = Arith::<P>::mmulmod(w1, w2);
            let w4 = Arith::<P>::mmulmod(w2, w2);
            let w5 = Arith::<P>::mmulmod(w2, w3);
            (w1, w2, w3, w4, w5)
        } else {
            (0, 0, 0, 0, 0)
        };
        for _ in 0..s1 {
            (*px, *px.wrapping_add(s1), *px.wrapping_add(2*s1),
            *px.wrapping_add(3*s1), *px.wrapping_add(4*s1), *px.wrapping_add(5*s1)) =
                ntt6_kernel::<P, INV, TWIDDLE>(w1, w2, w3, w4, w5,
                    *px, *px.wrapping_add(s1), *px.wrapping_add(2*s1),
                    *px.wrapping_add(3*s1), *px.wrapping_add(4*s1), *px.wrapping_add(5*s1));
            px = px.wrapping_add(1);
        }
    }
    (px.wrapping_add(5*s1), ptf.wrapping_add(1))
}

fn ntt_dif_dit<const P: u64, const INV: bool>(plan: &NttPlan, x: &mut [u64], tf_list: &[u64]) {
    let mut i_list: Vec<_> = (0..plan.s_list.len()).collect();
    if INV { i_list.reverse(); }
    let mut ptf = tf_list.as_ptr();
    for i in i_list {
        let (s, radix) = plan.s_list[i];
        let s1 = s/radix;
        let mut px = x.as_mut_ptr();
        let px_end = x.as_mut_ptr().wrapping_add(plan.n);
        match radix {
            2 => {
                (px, ptf) = ntt2_single_block::<P, INV, false>(s1, px, ptf);
                while px < px_end {
                    (px, ptf) = ntt2_single_block::<P, INV, true>(s1, px, ptf);
                }
            },
            3 => {
                (px, ptf) = ntt3_single_block::<P, INV, false>(s1, px, ptf);
                while px < px_end {
                    (px, ptf) = ntt3_single_block::<P, INV, true>(s1, px, ptf);
                }
            },
            4 => {
                (px, ptf) = ntt4_single_block::<P, INV, false>(s1, px, ptf);
                while px < px_end {
                    (px, ptf) = ntt4_single_block::<P, INV, true>(s1, px, ptf);
                }
            },
            5 => {
                (px, ptf) = ntt5_single_block::<P, INV, false>(s1, px, ptf);
                while px < px_end {
                    (px, ptf) = ntt5_single_block::<P, INV, true>(s1, px, ptf);
                }
            },
            6 => {
                (px, ptf) = ntt6_single_block::<P, INV, false>(s1, px, ptf);
                while px < px_end {
                    (px, ptf) = ntt6_single_block::<P, INV, true>(s1, px, ptf);
                }
            },
            _ => { unreachable!() }
        }
    }
}

fn compute_twiddle_factors<const P: u64, const INV: bool>(s_list: &[(usize, usize)], out: &mut [u64]) -> usize {
    let r = s_list.last().unwrap_or(&(1, 1)).1;
    let mut p = 1;
    out[0] = Arith::<P>::R;
    for i in (1..s_list.len()).rev() {
        let radix = s_list[i-1].1;
        let w = Arith::<P>::mpowmod(NttKernelImpl::<P, INV>::ROOTR, Arith::<P>::MAX_NTT_LEN/(p as u64 * radix as u64 * r as u64));
        for j in p..radix*p {
            out[j] = Arith::<P>::mmulmod(w, out[j - p]);
        }
        p *= radix;
    }
    p
}

// Performs (cyclic) integer convolution modulo P using NTT.
// Modifies the input buffers in-place.
// The output is saved in the slice `x`.
// The input slices must have the same length.
fn conv<const P: u64>(plan: &NttPlan, x: &mut [u64], xlen: usize, y: &mut [u64], ylen: usize, mut mult: u64) {
    assert!(!x.is_empty() && x.len() == y.len());

    let (_n, g, m) = (plan.n, plan.g, plan.m);
    let last_radix = plan.last_radix;

    /* multiply by a constant in advance */
    let len_inv = Arith::<P>::mmulmod(Arith::<P>::R3, (P-1)/m as u64);
    mult = Arith::<P>::mmulmod(Arith::<P>::mmulmod(Arith::<P>::R2, mult), len_inv);
    for v in if xlen < ylen { &mut x[g..g+xlen] } else { &mut y[g..g+ylen] } {
        *v = Arith::<P>::mmulmod(*v, mult);
    }

    /* compute the total space needed for twiddle factors */
    let tf_all_count = {
        let (mut radix_cumul, mut out) = (1, 0);
        for &(_, radix) in &plan.s_list {
            out += radix_cumul;
            radix_cumul *= radix;
        }
        core::cmp::max(out, 1)
    };

    /* build twiddle factors */
    let mut tf_list = vec![0u64; tf_all_count];
    tf_list[0] = Arith::<P>::R;
    let mut tf_last_start = core::cmp::min(tf_all_count - 1, 1);
    for i in 1..plan.s_list.len() {
        let x = compute_twiddle_factors::<P, false>(&plan.s_list[0..=i], &mut tf_list[tf_last_start..]);
        if i + 1 < plan.s_list.len() { tf_last_start += x; }
    }

    /* dif fft */
    ntt_dif_dit::<P, false>(plan, &mut x[g..], &tf_list);
    ntt_dif_dit::<P, false>(plan, &mut y[g..], &tf_list);

    /* naive multiplication */
    let mut i = g;
    let (mut ii, mut ii_mod_last_radix) = (0, 0);
    let tf = &tf_list[tf_last_start..];
    let mut tf_current = tf[0];
    let tf_mult = match plan.last_radix {
        2 => NttKernelImpl::<P, false>::U2,
        3 => NttKernelImpl::<P, false>::U3,
        4 => NttKernelImpl::<P, false>::U4,
        5 => NttKernelImpl::<P, false>::U5,
        6 => NttKernelImpl::<P, false>::U6,
        _ => Arith::<P>::R
    };
    while i < g + plan.n {
        if ii_mod_last_radix == 0 {
            tf_current = tf[ii];
        } else {
            tf_current = Arith::<P>::mmulmod(tf_current, tf_mult);
        }
        conv_base::<P>(g, x.as_mut_ptr().wrapping_add(i), y.as_mut_ptr().wrapping_add(i),
            tf_current);
        i += g;
        ii_mod_last_radix += 1;
        if ii_mod_last_radix == last_radix {
            ii += 1;
            ii_mod_last_radix = 0;
        }
    }

    /* dit fft */
    let mut tf_last_start = 0;
    for i in (1..plan.s_list.len()).rev() {
        tf_last_start += compute_twiddle_factors::<P, true>(&plan.s_list[0..=i], &mut tf_list[tf_last_start..]);
    }
    tf_list[tf_last_start] = Arith::<P>::R;
    ntt_dif_dit::<P, true>(plan, x, &tf_list);
}

////////////////////////////////////////////////////////////////////////////////

use core::cmp::{min, max};
use crate::big_digit::BigDigit;

const P1: u64 = 14_259_017_916_245_606_401; // Max NTT length = 2^22 * 3^21 * 5^2 = 1_096_847_532_018_892_800
const P2: u64 = 17_984_575_660_032_000_001; // Max NTT length = 2^19 * 3^17 * 5^6 = 1_057_916_215_296_000_000
const P3: u64 = 17_995_154_822_184_960_001; // Max NTT length = 2^17 * 3^22 * 5^4 = 2_570_736_403_169_280_000

const P2P3: u128 = P2 as u128 * P3 as u128;
const P1INV_R_MOD_P2: u64 = Arith::<P2>::mmulmod(Arith::<P2>::R2, arith::invmod(P1, P2));
const P1P2INV_R_MOD_P3: u64 = Arith::<P3>::mmulmod(Arith::<P3>::R2, arith::invmod((P1 as u128 * P2 as u128 % P3 as u128) as u64, P3));
const P1_R_MOD_P3: u64 = Arith::<P3>::mmulmod(Arith::<P3>::R2, P1);
const P1P2_LO: u64 = (P1 as u128 * P2 as u128) as u64;
const P1P2_HI: u64 = ((P1 as u128 * P2 as u128) >> 64) as u64;

fn mac3_two_primes(acc: &mut [u64], b: &[u64], c: &[u64], bits: u64) {
    fn pack_into(src: &[u64], dst1: &mut [u64], dst2: &mut [u64], bits: u64) {
        let mut p = 0u64;
        let mut pdst1 = dst1.as_mut_ptr();
        let mut pdst2 = dst2.as_mut_ptr();
        let mut x = 0u64;
        let mask = (1u64 << bits).wrapping_sub(1);
        for v in src {
            let mut k = 0;
            while k < 64 {
                x |= (v >> k) << p;
                let q = 64 - k;
                if p + q >= bits {
                    unsafe { let out = x & mask; (*pdst1, *pdst2) = (out, out); }
                    x = 0;
                    (pdst1, pdst2, k, p) = (pdst1.wrapping_add(1), pdst2.wrapping_add(1), k + bits - p, 0);
                } else {
                    p += q;
                    break;
                }
            }
        }
        unsafe {
            if p > 0 { let out = x & mask; (*pdst1, *pdst2) = (out, out); }
        }
    }

    assert!(bits < 63);
    let b_len = ((64 * b.len() as u64 + bits - 1) / bits) as usize;
    let c_len = ((64 * c.len() as u64 + bits - 1) / bits) as usize;
    let min_len = b_len + c_len;
    let plan_x = NttPlan::build::<P2>(min_len);
    let plan_y = NttPlan::build::<P3>(min_len);

    let mut x = vec![0u64; plan_x.g + plan_x.n];
    let mut y = vec![0u64; plan_y.g + plan_y.n];
    let mut r = vec![0u64; plan_x.g + plan_x.n];
    let mut s = vec![0u64; plan_y.g + plan_y.n];
    pack_into(b, &mut x[plan_x.g..], &mut y[plan_y.g..], bits);
    pack_into(c, &mut r[plan_x.g..], &mut s[plan_y.g..], bits);
    conv::<P2>(&plan_x, &mut x, b_len, &mut r[..plan_x.g+plan_x.n], c_len, arith::invmod(P3, P2));
    conv::<P3>(&plan_y, &mut y, b_len, &mut s[..plan_y.g+plan_y.n], c_len, Arith::<P3>::submod(0, arith::invmod(P2, P3)));

    /* merge the results in {x, y} into r (process carry along the way) */
    let mask = (1u64 << bits) - 1;
    let mut carry: u128 = 0;
    let (mut j, mut p) = (0usize, 0u64);
    let mut s: u64 = 0;
    let mut carry_acc: u64 = 0;
    for i in 0..min_len {
        /* extract the convolution result */
        let (a, b) = (x[i], y[i]);
        let (mut v, overflow) = (a as u128 * P3 as u128 + carry).overflowing_sub(b as u128 * P2 as u128);
        if overflow { v = v.wrapping_add(P2P3); }
        carry = v >> bits;

        /* write to s */
        let out = (v as u64) & mask;
        s |= out << p;
        p += bits;
        if p >= 64 {
            /* flush s to the output buffer */
            let (w, overflow1) = s.overflowing_add(carry_acc);
            let (w, overflow2) = acc[j].overflowing_add(w);
            acc[j] = w;
            carry_acc = u64::from(overflow1 || overflow2);

            /* roll-over */
            (j, p) = (j+1, p-64);
            s = out >> (bits - p);
        }
    }

    /* process remaining carries */
    carry_acc += s;
    while j < acc.len() {
        let (w, overflow) = acc[j].overflowing_add(carry_acc);
        acc[j] = w;
        carry_acc = u64::from(overflow);
        j += 1;
    }
}

fn mac3_three_primes(acc: &mut [u64], b: &[u64], c: &[u64]) {
    let min_len = b.len() + c.len();
    let plan_x = NttPlan::build::<P1>(min_len);
    let plan_y = NttPlan::build::<P2>(min_len);
    let plan_z = NttPlan::build::<P3>(min_len);
    let mut x = vec![0u64; plan_x.g + plan_x.n];
    let mut y = vec![0u64; plan_y.g + plan_y.n];
    let mut z = vec![0u64; plan_z.g + plan_z.n];
    let mut r = vec![0u64; max(x.len(), max(y.len(), z.len()))];

    /* convolution with modulo P1 */
    for i in 0..b.len() { x[plan_x.g + i] = if b[i] >= P1 { b[i] - P1 } else { b[i] }; }
    for i in 0..c.len() { r[plan_x.g + i] = if c[i] >= P1 { c[i] - P1 } else { c[i] }; }
    conv::<P1>(&plan_x, &mut x, b.len(), &mut r[..plan_x.g+plan_x.n], c.len(), 1);

    /* convolution with modulo P2 */
    for i in 0..b.len() { y[plan_y.g + i] = if b[i] >= P2 { b[i] - P2 } else { b[i] }; }
    for i in 0..c.len() { r[plan_y.g + i] = if c[i] >= P2 { c[i] - P2 } else { c[i] }; }
    (&mut r[plan_y.g..])[c.len()..plan_y.n].fill(0u64);
    conv::<P2>(&plan_y, &mut y, b.len(), &mut r[..plan_y.g+plan_y.n], c.len(), 1);

    /* convolution with modulo P3 */
    for i in 0..b.len() { z[plan_z.g + i] = if b[i] >= P3 { b[i] - P3 } else { b[i] }; }
    for i in 0..c.len() { r[plan_z.g + i] = if c[i] >= P3 { c[i] - P3 } else { c[i] }; }
    (&mut r[plan_z.g..])[c.len()..plan_z.n].fill(0u64);
    conv::<P3>(&plan_z, &mut z, b.len(), &mut r[..plan_z.g+plan_z.n], c.len(), 1);

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

fn mac3_u64(acc: &mut [u64], b: &[u64], c: &[u64]) {
    const fn compute_bits(l: u64) -> u64 {
        let total_bits = l * 64;
        let (mut lo, mut hi) = (42, 62);
        while lo < hi {
            let mid = (lo + hi + 1) / 2;
            let single_digit_max_val = (1u64 << mid) - 1;
            let l_corrected = (total_bits + mid - 1) / mid;
            let (lhs, overflow) = (single_digit_max_val as u128 * single_digit_max_val as u128).overflowing_mul(l_corrected as u128);
            if !overflow && lhs < P2 as u128 * P3 as u128 { lo = mid; }
            else { hi = mid - 1; }
        }
        lo
    }

    let (b, c) = if b.len() < c.len() { (b, c) } else { (c, b) };
    let naive_cost = NttPlan::build::<P1>(b.len() + c.len()).cost;
    let split_cost = NttPlan::build::<P1>(b.len() + b.len()).cost * (c.len() / b.len()) +
        if c.len() % b.len() > 0 { NttPlan::build::<P1>(b.len() + (c.len() % b.len())).cost } else { 0 };
    if b.len() >= 128 && split_cost < naive_cost {
        /* special handling for unbalanced multiplication:
           we reduce it to about `c.len()/b.len()` balanced multiplications */
        let mut i = 0usize;
        let mut carry = 0u64;
        while i < c.len() {
            let j = min(i + b.len(), c.len());
            let k = j + b.len();
            let tmp = acc[k];
            acc[k] = 0;
            mac3_u64(&mut acc[i..=k], b, &c[i..j]);
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
        while carry > 0 && i < acc.len() {
            let (v, overflow) = acc[i].overflowing_add(carry);
            acc[i] = v;
            carry = u64::from(overflow);
            i += 1;
        }
        return;
    }

    // We have two choices:
    //     1. NTT with two primes.
    //     2. NTT with three primes.
    // Obviously we want to do only two passes for efficiency, not three.
    // However, the number of bits per u64 we can pack for NTT
    // depends on the length of the arrays being multiplied (convolved).
    // If the arrays are too long, the resulting values may exceed the
    // modulus range P2 * P3, which leads to incorrect results.
    // Hence, we compute the number of bits required by the length of NTT,
    // and use it to determine whether to use two-prime or three-prime.
    // Since we can pack 64 bits per u64 in three-prime NTT, the effective
    // number of bits in three-prime NTT is 64/3 = 21.3333..., which means
    // two-prime NTT can only do better when at least 43 bits per u64 can
    // be packed into each u64.
    let max_cnt = max(b.len(), c.len()) as u64;
    let bits = compute_bits(max_cnt);
    if bits >= 43 {
        /* can pack more effective bits per u64 with two primes than with three primes */
        mac3_two_primes(acc, b, c, bits);
    } else {
        /* can pack at most 21 effective bits per u64, which is worse than
           64/3 = 21.3333.. effective bits per u64 achieved with three primes */
        mac3_three_primes(acc, b, c);
    }
}

////////////////////////////////////////////////////////////////////////////////

#[cfg(u64_digit)]
pub fn mac3(acc: &mut [BigDigit], b: &[BigDigit], c: &[BigDigit]) {
    mac3_u64(acc, b, c);
}

#[cfg(not(u64_digit))]
pub fn mac3(acc: &mut [BigDigit], b: &[BigDigit], c: &[BigDigit]) {
    fn bigdigit_to_u64(src: &[BigDigit]) -> Vec<u64> {
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
    mac3_u64(&mut acc_u64, &b_u64, &c_u64);
    u64_to_bigdigit(&acc_u64, acc);
}