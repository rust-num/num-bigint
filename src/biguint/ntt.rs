#![forbid(unsafe_code)]

use super::{BigDigit, Vec};
use core::cmp::max;
use core::iter::zip;

mod arith {
    // Extended Euclidean algorithm:
    //   (g, x, y) is a solution to ax + by = g, where g = gcd(a, b)
    const fn egcd(mut a: i128, mut b: i128) -> (i128, i128, i128) {
        assert!(a > 0 && b > 0);

        // treat as a row-major 2x2 matrix
        let mut c = if a > b {
            (a, b) = (b, a);
            [0, 1, 1, 0]
        } else {
            [1, 0, 0, 1]
        };
        loop {
            if a == 0 {
                break (b, c[1], c[3]);
            }
            let (q, r) = (b / a, b % a);
            (a, b) = (r, a);
            c = [c[1] - q * c[0], c[0], c[3] - q * c[2], c[2]];
        }
    }

    // Modular inverse: a^-1 mod modulus
    //   (m == 0 means m == 2^64)
    const fn invmod_inner(a: u64, modulus: u64) -> u64 {
        let m = if modulus == 0 {
            1i128 << 64
        } else {
            modulus as i128
        };
        let (g, mut x, _y) = egcd(a as i128, m);
        assert!(g == 1);
        if x < 0 {
            x += m;
        }
        assert!(x > 0 && x < 1i128 << 64);
        x as u64
    }

    // Modular inverse: a^-1 mod modulus (modulus != 0)
    pub const fn invmod(a: u64, modulus: u64) -> u64 {
        assert!(modulus != 0);
        invmod_inner(a, modulus)
    }

    // Modular inverse: a^-1 mod 2^64
    pub const fn invmod_2pow64(a: u64) -> u64 {
        invmod_inner(a, 0)
    }
}

// NTT-local fixed-width modular arithmetic.
//
// Although this uses Montgomery reduction, it intentionally does not use
// `biguint::monty`: P is a compile-time u64 modulus and R is 2^64. The
// subtractive reduction remains const-evaluable, handles primes close to
// 2^64 without a 129-bit addition, and supports the fused and delayed
// reductions required by the NTT kernels.
struct Arith<const P: u64> {}

impl<const P: u64> Arith<P> {
    const R: u64 = ((1u128 << 64) % P as u128) as u64; // 2^64 mod P
    const R2: u64 = (Self::R as u128 * Self::R as u128 % P as u128) as u64; // R^2 mod P
    const PINV: u64 = arith::invmod_2pow64(P); // P^-1 mod 2^64

    const MAX_NTT_LEN: u64 =
        2u64.pow(Self::FACTORS_2) * 3u64.pow(Self::FACTORS_3) * 5u64.pow(Self::FACTORS_5);

    const ROOTR: u64 = {
        // ROOT * R mod P (ROOT: MAX_NTT_LEN divides MultiplicativeOrder[ROOT, P])
        //
        // Find an element p whose multiplicative order is a multiple of
        // MAX_NTT_LEN, then compute ROOTR = p^((P-1)/MAX_NTT_LEN) to obtain
        // an element of exact order MAX_NTT_LEN.
        //
        // We only test p^((P-1)/q) != 1 for q in {2, 3, 5} — the prime
        // factors of MAX_NTT_LEN. This suffices because:
        //   - FACTORS_2/3/5 extract all factors of 2, 3, 5 from P-1, so
        //     Q = (P-1)/MAX_NTT_LEN is coprime to MAX_NTT_LEN by construction.
        //   - Any remaining prime factors of P-1 (e.g. 7, 13, 17 for the
        //     current primes) are absorbed by the p^Q exponentiation and
        //     do not affect the order of ROOTR.
        //
        // P must be prime for Z/PZ to be a field; this is not checked here.
        assert!(Self::MAX_NTT_LEN % 4050 == 0);
        let mut p = Self::R;
        loop {
            if Self::mpowmod(p, P / 2) != Self::R
                && Self::mpowmod(p, P / 3) != Self::R
                && Self::mpowmod(p, P / 5) != Self::R
            {
                break Self::mpowmod(p, P / Self::MAX_NTT_LEN);
            }
            p = Self::addmod(p, Self::R);
        }
    };

    // Explicit instantiations to ensure const evaluation at compile time.
    const FACTORS_2: u32 = Self::factors(2);
    const FACTORS_3: u32 = Self::factors(3);
    const FACTORS_5: u32 = Self::factors(5);

    // Counts the number of `divisor` factors in P-1.
    const fn factors(divisor: u64) -> u32 {
        let (mut tmp, mut ans) = (P - 1, 0);
        while tmp % divisor == 0 {
            tmp /= divisor;
            ans += 1;
        }
        ans
    }

    // Montgomery reduction:
    //   x * R^-1 mod P
    const fn mreduce(x: u128) -> u64 {
        let m = (x as u64).wrapping_mul(Self::PINV);
        let y = ((m as u128 * P as u128) >> 64) as u64;
        let (out, overflow) = ((x >> 64) as u64).overflowing_sub(y);
        if overflow {
            out.wrapping_add(P)
        } else {
            out
        }
    }

    // Multiplication with Montgomery reduction:
    //   a * b * R^-1 mod P
    const fn mmulmod(a: u64, b: u64) -> u64 {
        Self::mreduce(a as u128 * b as u128)
    }

    // Fused-multiply-sub with Montgomery reduction:
    //   a * b * R^-1 - c mod P
    const fn mmulsubmod(a: u64, b: u64, c: u64) -> u64 {
        let x = a as u128 * b as u128;
        let lo = x as u64;
        let hi = Self::submod((x >> 64) as u64, c);
        Self::mreduce(lo as u128 | ((hi as u128) << 64))
    }

    // Computes base^exponent mod P with Montgomery reduction.
    const fn mpowmod(mut base: u64, mut exponent: u64) -> u64 {
        let mut cur = Self::R;
        while exponent > 0 {
            if exponent % 2 > 0 {
                cur = Self::mmulmod(cur, base);
            }
            exponent /= 2;
            base = Self::mmulmod(base, base);
        }
        cur
    }

    // Computes c as u128 * mreduce(v) as u128, using d: u64 = mmulmod(P-1, c).
    //
    // It is the caller's responsibility to ensure that d is correct.
    // Note that d can be computed by calling mreducelo(c).
    const fn mmulmod_noreduce(v: u128, c: u64, d: u64) -> u128 {
        let a: u128 = c as u128 * (v >> 64);
        let b: u128 = d as u128 * (v as u64 as u128);
        let (w, overflow) = a.overflowing_sub(b);
        if overflow {
            w.wrapping_add((P as u128) << 64)
        } else {
            w
        }
    }

    // Computes submod(0, mreduce(x as u128)) fast.
    const fn mreducelo(x: u64) -> u64 {
        let m = x.wrapping_mul(Self::PINV);
        ((m as u128 * P as u128) >> 64) as u64
    }

    // Computes a + b mod P, output range [0, P).
    const fn addmod(a: u64, b: u64) -> u64 {
        Self::submod(a, P.wrapping_sub(b))
    }

    // Computes a + b mod P, output range [0, 2^64).
    const fn addmod64(a: u64, b: u64) -> u64 {
        let (out, overflow) = a.overflowing_add(b);
        if overflow {
            out.wrapping_sub(P)
        } else {
            out
        }
    }

    // Computes a + b mod P, selects addmod64 or addmod depending on INV && TWIDDLE.
    const fn addmodopt_invtw<const INV: bool, const TWIDDLE: bool>(a: u64, b: u64) -> u64 {
        if INV && TWIDDLE {
            Self::addmod64(a, b)
        } else {
            Self::addmod(a, b)
        }
    }

    // Computes a - b mod P, output range [0, P).
    const fn submod(a: u64, b: u64) -> u64 {
        let (out, overflow) = a.overflowing_sub(b);
        if overflow {
            out.wrapping_add(P)
        } else {
            out
        }
    }
}

struct NttPlan {
    // The total convolution length, which equals g * m.
    n: usize,
    // The length of each base case handled by naive multiplication.
    g: usize,
    // The product of radices processed by NTT. Should divide `Arith::<P>::MAX_NTT_LEN`.
    m: usize,
    // The NTT radix scheduled closest to the naive multiplication.
    last_radix: usize,
    // The list of tuples (current block size, radix) in DIF (forward) order.
    s_list: Vec<(usize, usize)>,
}

// Planner costs normally fit in usize. The wide variant keeps comparisons exact
// for theoretical lengths where multiplying by the unit cost would overflow.
// Its declaration order also matches its numeric order: every Large value is
// greater than every Small value.
#[derive(Clone, Copy, Eq, Ord, PartialEq, PartialOrd)]
enum PlanCost {
    Small(usize),
    Large(u128),
}

impl PlanCost {
    fn new(len: usize, unit_cost: u32) -> Self {
        if let Ok(unit_cost_usize) = usize::try_from(unit_cost) {
            if let Some(cost) = len.checked_mul(unit_cost_usize) {
                return Self::Small(cost);
            }
        }
        Self::Large(len as u128 * unit_cost as u128)
    }
}

impl NttPlan {
    fn ntt_len_to_radices(mut ntt_len: usize) -> [u32; 5] {
        let m2 = ntt_len.trailing_zeros();
        ntt_len >>= m2;
        let mut m3 = 0;
        while ntt_len % 3 == 0 {
            ntt_len /= 3;
            m3 += 1;
        }
        let mut m5 = 0;
        while ntt_len % 5 == 0 {
            ntt_len /= 5;
            m5 += 1;
        }

        // Make sure `ntt_len` has no other factors besides 2, 3, and 5
        assert_eq!(ntt_len, 1);

        Self::ntt_len_factors_to_radices(m2, m3, m5)
    }

    fn ntt_len_factors_to_radices(mut m2: u32, mut m3: u32, m5: u32) -> [u32; 5] {
        // radix-5 does not compete with other radices (2, 3, 4, 6).
        let cnt5 = m5;

        // Pick radix-6 and radix-4 greedily and assign the remaining factors to radix-2 and radix-3.
        let cnt6 = m2.min(m3);
        m2 -= cnt6;
        m3 -= cnt6;
        let cnt4 = m2 / 2;
        m2 -= cnt4 * 2;
        let [cnt2, cnt3] = [m2, m3];

        [cnt2, cnt3, cnt4, cnt5, cnt6]
    }

    fn build<const P: u64>(min_len: usize) -> Self {
        // Allowed values for naive multiplication chunk length (base case).
        //
        // Under the cost model, g = 2 and g = 3 cannot be optimal. If the NTT radix
        // decomposition contains one of {4, 5, 6}, it can be exchanged with g to
        // improve the cost. The combination g = 3 and NTT radix 4 is an exception,
        // but g = 6 with NTT radix 2 is more favorable. If g is 2 or 3 and the
        // decomposition contains only {2, 3}, those radices can instead be absorbed
        // into g.
        //
        // Hence, we only consider g >= 4.
        const G_MIN: usize = 4;
        const G_MAX: usize = 9;

        // Cost model for NTT length selection.
        //
        // For a candidate convolution length `len = g * m`, the estimated cost is
        // `len * (G_COSTS[g - G_MIN] + sum(radix_counts[i] * RADIX_COSTS[i]))`.
        // Each radix-r stage processes all `len` elements in `len / r` groups of
        // r-point butterflies, so the cost of a stage is modeled as linear in
        // `len`. The constants below represent the relative per-element costs of
        // the naive base cases and radix butterflies and were calibrated
        // empirically. `G_COSTS` covers g = G_MIN..=G_MAX, while `RADIX_COSTS`
        // covers radix-2 through radix-6.
        const G_COSTS: [u32; G_MAX - G_MIN + 1] = [218, 287, 337, 433, 564, 663];
        const RADIX_COSTS: [u32; 5] = [458, 791, 756, 1230, 915];

        // Short cases do not require an NTT.
        // The threshold could be larger, but we keep it minimal for simplicity.
        // (Short cases will not invoke the NTT anyway.)
        if min_len <= G_MAX {
            // Special case for short `min_len`.
            let g = min_len.max(1);
            return Self::from_params::<P>(g, g);
        }

        let max_ntt_len = usize::try_from(Arith::<P>::MAX_NTT_LEN).unwrap_or(usize::MAX);
        assert!(min_len <= max_ntt_len.saturating_mul(G_MAX));

        // The NTT length `m` contains only factors of 2, 3, and 5, but the single
        // base-case length `g` may have any factorization: it is handled directly
        // by the naive quadratic kernel rather than by radix-specific NTT stages.
        // Thus `g = 7` can contribute one factor of 7 to `len = g * m`, reducing
        // padding for some input lengths without requiring a radix-7 NTT kernel.
        let mut best: Option<(usize, usize, PlanCost)> = None;
        for g in G_MIN..=G_MAX {
            let ntt_min_len = min_len / g + usize::from(min_len % g != 0);
            for m5 in 0..=Arith::<P>::FACTORS_5 {
                let pow5 = match 5usize.checked_pow(m5) {
                    Some(pow5) => pow5,
                    None => break,
                };
                for m3 in 0..=Arith::<P>::FACTORS_3 {
                    // Compute the length of the transform
                    let pow3 = match 3usize.checked_pow(m3) {
                        Some(pow3) => pow3,
                        None => break,
                    };
                    let mut ntt_len = match pow5.checked_mul(pow3) {
                        Some(ntt_len) => ntt_len,
                        None => break,
                    };
                    if ntt_len / 2 >= ntt_min_len {
                        // Too much padding will not be computationally optimal
                        break;
                    }
                    let mut m2 = 0;
                    while ntt_len < ntt_min_len && m2 < Arith::<P>::FACTORS_2 {
                        ntt_len = match ntt_len.checked_mul(2) {
                            Some(ntt_len) => ntt_len,
                            None => break,
                        };
                        m2 += 1;
                    }
                    if ntt_len < ntt_min_len {
                        continue;
                    }
                    let len = match g.checked_mul(ntt_len) {
                        Some(len) => len,
                        None => continue,
                    };

                    // Compute the cost of the transform
                    let radix_counts = Self::ntt_len_factors_to_radices(m2, m3, m5);
                    let unit_cost = G_COSTS[g - G_MIN]
                        + radix_counts
                            .into_iter()
                            .zip(RADIX_COSTS)
                            .map(|(radix_count, radix_cost)| radix_count * radix_cost)
                            .sum::<u32>();
                    let cost = PlanCost::new(len, unit_cost);

                    // Record the transform with the minimum cost
                    if best.map_or(true, |(_, _, best_cost)| cost < best_cost) {
                        best = Some((len, g, cost));
                    }
                }
            }
        }
        let (best_len, best_g, _) = best.expect("NTT length is too large");
        Self::from_params::<P>(best_len, best_g)
    }

    fn from_params<const P: u64>(len: usize, g: usize) -> Self {
        assert!(len % g == 0);
        let ntt_len = len / g;
        let [cnt2, cnt3, cnt4, cnt5, cnt6] = Self::ntt_len_to_radices(ntt_len);
        let s_list = {
            let capacity = (cnt2 + cnt3 + cnt4 + cnt5 + cnt6) as usize;
            let mut out = Vec::with_capacity(capacity);
            let mut tmp = len;
            for _ in 0..cnt2 {
                out.push((tmp, 2));
                tmp /= 2;
            }
            for _ in 0..cnt3 {
                out.push((tmp, 3));
                tmp /= 3;
            }
            for _ in 0..cnt4 {
                out.push((tmp, 4));
                tmp /= 4;
            }
            for _ in 0..cnt5 {
                out.push((tmp, 5));
                tmp /= 5;
            }
            for _ in 0..cnt6 {
                out.push((tmp, 6));
                tmp /= 6;
            }
            out
        };

        Self {
            n: len,
            g,
            m: ntt_len,
            last_radix: s_list.last().unwrap_or(&(1, 1)).1,
            s_list,
        }
    }
}

fn conv_base<const P: u64>(out: &mut [u64], x: &[u64], y: &[u64], c: u64) {
    assert_eq!(out.len(), x.len());
    assert_eq!(x.len(), y.len());

    let n = out.len();
    let c2 = Arith::<P>::mreducelo(c);
    for i in 0..n {
        let mut v: u128 = 0;
        for j in i + 1..n {
            let (w, overflow) = v.overflowing_sub(x[j] as u128 * y[i + n - j] as u128);
            v = if overflow {
                w.wrapping_add((P as u128) << 64)
            } else {
                w
            };
        }

        v = Arith::<P>::mmulmod_noreduce(v, c, c2);
        for j in 0..=i {
            let (w, overflow) = v.overflowing_sub(x[j] as u128 * y[i - j] as u128);
            v = if overflow {
                w.wrapping_add((P as u128) << 64)
            } else {
                w
            };
        }

        out[i] = Arith::<P>::mreduce(v);
    }
}

struct NttKernelImpl<const P: u64, const INV: bool>;

impl<const P: u64, const INV: bool> NttKernelImpl<P, INV> {
    const ROOTR: u64 = Arith::<P>::mpowmod(
        Arith::<P>::ROOTR,
        if INV { Arith::<P>::MAX_NTT_LEN - 1 } else { 1 },
    );

    const U3: u64 = Arith::<P>::mpowmod(Self::ROOTR, Arith::<P>::MAX_NTT_LEN / 3);
    const U4: u64 = Arith::<P>::mpowmod(Self::ROOTR, Arith::<P>::MAX_NTT_LEN / 4);
    const U5: u64 = Arith::<P>::mpowmod(Self::ROOTR, Arith::<P>::MAX_NTT_LEN / 5);
    const U6: u64 = Arith::<P>::mpowmod(Self::ROOTR, Arith::<P>::MAX_NTT_LEN / 6);

    const C5: (u64, u64, u64, u64, u64, u64) = {
        let w = Self::U5;
        let w2 = Arith::<P>::mpowmod(w, 2);
        let w4 = Arith::<P>::mpowmod(w, 4);
        let inv2 = Arith::<P>::mmulmod(Arith::<P>::R2, arith::invmod(2, P));
        let inv4 = Arith::<P>::mmulmod(Arith::<P>::R2, arith::invmod(4, P));
        let c51 = Arith::<P>::addmod(Arith::<P>::R, inv4); // 1 + 4^-1 mod P
        let c52 = Arith::<P>::addmod(Arith::<P>::mmulmod(inv2, Arith::<P>::addmod(w, w4)), inv4); // 4^-1 * (2*w + 2*w^4 + 1) mod P
        let c53 = Arith::<P>::mmulmod(inv2, Arith::<P>::submod(w, w4)); // 2^-1 * (w - w^4) mod P
        let c54 = Arith::<P>::addmod(Arith::<P>::addmod(w, w2), inv2); // 2^-1 * (2*w + 2*w^2 + 1) mod P
        let c55 = Arith::<P>::addmod(Arith::<P>::addmod(w2, w4), inv2); // 2^-1 * (2*w^2 + 2*w^4 + 1) mod P
        (0, c51, c52, c53, c54, c55)
    };
}

const fn ntt2_kernel<const P: u64, const INV: bool, const TWIDDLE: bool>(
    a: u64,
    mut b: u64,
    w1: u64,
) -> [u64; 2] {
    if !INV && TWIDDLE {
        b = Arith::<P>::mmulmod(w1, b);
    }
    let out0 = Arith::<P>::addmod(a, b);
    let mut out1 = Arith::<P>::submod(a, b);
    if INV && TWIDDLE {
        out1 = Arith::<P>::mmulmod(w1, out1);
    }
    [out0, out1]
}

fn ntt2_single_block<const P: u64, const INV: bool, const TWIDDLE: bool>(px: &mut [u64], w1: u64) {
    let w1 = if TWIDDLE { w1 } else { 0 };
    let s1 = px.len() / 2;
    let (a, b) = px.split_at_mut(s1);
    zip(a, b).for_each(|(a, b)| {
        [*a, *b] = ntt2_kernel::<P, INV, TWIDDLE>(*a, *b, w1);
    });
}

const fn ntt3_kernel<const P: u64, const INV: bool, const TWIDDLE: bool>(
    a: u64,
    mut b: u64,
    mut c: u64,
    [w1, w2]: [u64; 2],
) -> [u64; 3] {
    if !INV && TWIDDLE {
        b = Arith::<P>::mmulmod(w1, b);
        c = Arith::<P>::mmulmod(w2, c);
    }
    let kbmc = Arith::<P>::mmulmod(NttKernelImpl::<P, INV>::U3, Arith::<P>::submod(b, c));
    let out0 = Arith::<P>::addmod(a, Arith::<P>::addmod(b, c));
    let mut out1 = Arith::<P>::submod(a, Arith::<P>::submod(c, kbmc));
    let mut out2 = Arith::<P>::submod(Arith::<P>::submod(a, b), kbmc);
    if INV && TWIDDLE {
        out1 = Arith::<P>::mmulmod(w1, out1);
        out2 = Arith::<P>::mmulmod(w2, out2);
    }
    [out0, out1, out2]
}

fn ntt3_single_block<const P: u64, const INV: bool, const TWIDDLE: bool>(px: &mut [u64], w1: u64) {
    let w = if TWIDDLE {
        let w2 = Arith::<P>::mmulmod(w1, w1);
        [w1, w2]
    } else {
        [0; 2]
    };
    let s1 = px.len() / 3;
    let (a, rest) = px.split_at_mut(s1);
    let (b, c) = rest.split_at_mut(s1);
    zip(zip(a, b), c).for_each(|((a, b), c)| {
        [*a, *b, *c] = ntt3_kernel::<P, INV, TWIDDLE>(*a, *b, *c, w);
    });
}

const fn ntt4_kernel<const P: u64, const INV: bool, const TWIDDLE: bool>(
    a: u64,
    mut b: u64,
    mut c: u64,
    mut d: u64,
    [w1, w2, w3]: [u64; 3],
) -> [u64; 4] {
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
    let mut out1 = Arith::<P>::addmodopt_invtw::<INV, TWIDDLE>(amc, jbmd);
    let mut out2 = Arith::<P>::submod(apc, bpd);
    let mut out3 = Arith::<P>::submod(amc, jbmd);
    if INV && TWIDDLE {
        out1 = Arith::<P>::mmulmod(w1, out1);
        out2 = Arith::<P>::mmulmod(w2, out2);
        out3 = Arith::<P>::mmulmod(w3, out3);
    }
    [out0, out1, out2, out3]
}

fn ntt4_single_block<const P: u64, const INV: bool, const TWIDDLE: bool>(px: &mut [u64], w1: u64) {
    let w = if TWIDDLE {
        let w2 = Arith::<P>::mmulmod(w1, w1);
        let w3 = Arith::<P>::mmulmod(w1, w2);
        [w1, w2, w3]
    } else {
        [0; 3]
    };
    let s1 = px.len() / 4;
    let (a, rest) = px.split_at_mut(s1);
    let (b, rest) = rest.split_at_mut(s1);
    let (c, d) = rest.split_at_mut(s1);
    zip(zip(zip(a, b), c), d).for_each(|(((a, b), c), d)| {
        [*a, *b, *c, *d] = ntt4_kernel::<P, INV, TWIDDLE>(*a, *b, *c, *d, w);
    });
}

const fn ntt5_kernel<const P: u64, const INV: bool, const TWIDDLE: bool>(
    a: u64,
    mut b: u64,
    mut c: u64,
    mut d: u64,
    mut e: u64,
    [w1, w2, w3, w4]: [u64; 4],
) -> [u64; 5] {
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
    let m2 = Arith::<P>::mmulsubmod(NttKernelImpl::<P, INV>::C5.1, t5, m1);
    let m3 = Arith::<P>::mmulmod(NttKernelImpl::<P, INV>::C5.2, t6);
    let m4 = Arith::<P>::mmulmod(NttKernelImpl::<P, INV>::C5.3, t7);
    let m5 = Arith::<P>::mmulsubmod(NttKernelImpl::<P, INV>::C5.4, t4, m4);
    let m6 = Arith::<P>::mmulsubmod(P.wrapping_sub(NttKernelImpl::<P, INV>::C5.5), t3, m4);
    let s1 = Arith::<P>::submod(m3, m2);
    let s2 = Arith::<P>::addmod(m2, m3);
    let out0 = m1;
    let mut out1 = Arith::<P>::submod(s1, m5);
    let mut out2 = Arith::<P>::submod(Arith::<P>::submod(0, s2), m6);
    let mut out3 = Arith::<P>::submod(m6, s2);
    let mut out4 = Arith::<P>::addmodopt_invtw::<INV, TWIDDLE>(s1, m5);
    if INV && TWIDDLE {
        out1 = Arith::<P>::mmulmod(w1, out1);
        out2 = Arith::<P>::mmulmod(w2, out2);
        out3 = Arith::<P>::mmulmod(w3, out3);
        out4 = Arith::<P>::mmulmod(w4, out4);
    }
    [out0, out1, out2, out3, out4]
}

fn ntt5_single_block<const P: u64, const INV: bool, const TWIDDLE: bool>(px: &mut [u64], w1: u64) {
    let w = if TWIDDLE {
        let w2 = Arith::<P>::mmulmod(w1, w1);
        let w3 = Arith::<P>::mmulmod(w1, w2);
        let w4 = Arith::<P>::mmulmod(w2, w2);
        [w1, w2, w3, w4]
    } else {
        [0; 4]
    };
    let s1 = px.len() / 5;
    let (a, rest) = px.split_at_mut(s1);
    let (b, rest) = rest.split_at_mut(s1);
    let (c, rest) = rest.split_at_mut(s1);
    let (d, e) = rest.split_at_mut(s1);
    zip(zip(zip(zip(a, b), c), d), e).for_each(|((((a, b), c), d), e)| {
        [*a, *b, *c, *d, *e] = ntt5_kernel::<P, INV, TWIDDLE>(*a, *b, *c, *d, *e, w);
    });
}

const fn ntt6_kernel<const P: u64, const INV: bool, const TWIDDLE: bool>(
    mut a: u64,
    mut b: u64,
    mut c: u64,
    mut d: u64,
    mut e: u64,
    mut f: u64,
    [w1, w2, w3, w4, w5]: [u64; 5],
) -> [u64; 6] {
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
    let mut out2 = Arith::<P>::submod(a, Arith::<P>::submod(b, lbmc));
    let mut out4 = Arith::<P>::submod(Arith::<P>::submod(a, c), lbmc);
    let lepf = Arith::<P>::mmulmod(NttKernelImpl::<P, INV>::U6, Arith::<P>::addmod64(e, f));
    let mut out1 = Arith::<P>::submod(d, Arith::<P>::submod(f, lepf));
    let mut out3 = Arith::<P>::submod(d, Arith::<P>::submod(e, f));
    let mut out5 = Arith::<P>::submod(d, Arith::<P>::submod(lepf, e));
    if INV && TWIDDLE {
        out1 = Arith::<P>::mmulmod(w1, out1);
        out2 = Arith::<P>::mmulmod(w2, out2);
        out3 = Arith::<P>::mmulmod(w3, out3);
        out4 = Arith::<P>::mmulmod(w4, out4);
        out5 = Arith::<P>::mmulmod(w5, out5);
    }
    [out0, out1, out2, out3, out4, out5]
}

fn ntt6_single_block<const P: u64, const INV: bool, const TWIDDLE: bool>(px: &mut [u64], w1: u64) {
    let w = if TWIDDLE {
        let w2 = Arith::<P>::mmulmod(w1, w1);
        let w3 = Arith::<P>::mmulmod(w1, w2);
        let w4 = Arith::<P>::mmulmod(w2, w2);
        let w5 = Arith::<P>::mmulmod(w2, w3);
        [w1, w2, w3, w4, w5]
    } else {
        [0; 5]
    };
    let s1 = px.len() / 6;
    let (a, rest) = px.split_at_mut(s1);
    let (b, rest) = rest.split_at_mut(s1);
    let (c, rest) = rest.split_at_mut(s1);
    let (d, rest) = rest.split_at_mut(s1);
    let (e, f) = rest.split_at_mut(s1);
    zip(zip(zip(zip(zip(a, b), c), d), e), f).for_each(|(((((a, b), c), d), e), f)| {
        [*a, *b, *c, *d, *e, *f] = ntt6_kernel::<P, INV, TWIDDLE>(*a, *b, *c, *d, *e, *f, w);
    });
}

fn ntt_dif_dit<const P: u64, const INV: bool>(plan: &NttPlan, x: &mut [u64], twiddles: &[u64]) {
    let mut twiddle_index = 0;
    for step in 0..plan.s_list.len() {
        let i = if INV {
            plan.s_list.len() - 1 - step
        } else {
            step
        };
        let (s, radix) = plan.s_list[i];
        debug_assert_eq!(plan.n % s, 0);
        debug_assert_eq!(s % radix, 0);

        let mut x_iter = x[..plan.n].chunks_exact_mut(s);
        match radix {
            2 => {
                ntt2_single_block::<P, INV, false>(x_iter.next().unwrap(), 0);
                twiddle_index += 1;
                for px in x_iter {
                    ntt2_single_block::<P, INV, true>(px, twiddles[twiddle_index]);
                    twiddle_index += 1;
                }
            }
            3 => {
                ntt3_single_block::<P, INV, false>(x_iter.next().unwrap(), 0);
                twiddle_index += 1;
                for px in x_iter {
                    ntt3_single_block::<P, INV, true>(px, twiddles[twiddle_index]);
                    twiddle_index += 1;
                }
            }
            4 => {
                ntt4_single_block::<P, INV, false>(x_iter.next().unwrap(), 0);
                twiddle_index += 1;
                for px in x_iter {
                    ntt4_single_block::<P, INV, true>(px, twiddles[twiddle_index]);
                    twiddle_index += 1;
                }
            }
            5 => {
                ntt5_single_block::<P, INV, false>(x_iter.next().unwrap(), 0);
                twiddle_index += 1;
                for px in x_iter {
                    ntt5_single_block::<P, INV, true>(px, twiddles[twiddle_index]);
                    twiddle_index += 1;
                }
            }
            6 => {
                ntt6_single_block::<P, INV, false>(x_iter.next().unwrap(), 0);
                twiddle_index += 1;
                for px in x_iter {
                    ntt6_single_block::<P, INV, true>(px, twiddles[twiddle_index]);
                    twiddle_index += 1;
                }
            }
            _ => {
                unreachable!()
            }
        }
    }
}

fn calc_twiddle_factors<const P: u64, const INV: bool>(
    s_list: &[(usize, usize)],
    out: &mut [u64],
) -> usize {
    let mut p = 1;
    out[0] = Arith::<P>::R;
    for i in (1..s_list.len()).rev() {
        let radix = s_list[i - 1].1;
        let w = Arith::<P>::mpowmod(
            NttKernelImpl::<P, INV>::ROOTR,
            Arith::<P>::MAX_NTT_LEN / (p * radix * s_list.last().unwrap().1) as u64,
        );
        for j in p..radix * p {
            out[j] = Arith::<P>::mmulmod(w, out[j - p]);
        }
        p *= radix;
    }
    p
}

// Performs (cyclic) integer convolution modulo P using NTT.
// Modifies the input buffers in place.
// The output is saved in the slice `x`.
// The input slices must have the same length.
fn conv<const P: u64>(
    plan: &NttPlan,
    x: &mut [u64],
    xlen: usize,
    y: &mut [u64],
    ylen: usize,
    mut mult: u64,
) {
    assert!(!x.is_empty() && x.len() == y.len());
    let (_n, g, m, last_radix) = (plan.n, plan.g, plan.m, plan.last_radix as u64);

    // multiply by a constant in advance
    mult = Arith::<P>::mmulmod(
        Arith::<P>::mpowmod(Arith::<P>::R2, 3),
        Arith::<P>::mmulmod(mult, (P - 1) / m as u64),
    );
    for v in if xlen < ylen {
        &mut x[g..g + xlen]
    } else {
        &mut y[g..g + ylen]
    } {
        *v = Arith::<P>::mmulmod(*v, mult);
    }

    // compute the total space needed for twiddle factors
    let (mut radix_prefix_product, mut twiddle_count) = (1, 2); // 2 extra slots
    for &(_, radix) in &plan.s_list {
        twiddle_count += radix_prefix_product;
        radix_prefix_product *= radix;
    }

    // build twiddle factors
    let mut twiddles = vec![0u64; twiddle_count];
    let mut last_stage_twiddle_start = 0;
    for i in 0..plan.s_list.len() {
        let stage_twiddle_count = calc_twiddle_factors::<P, false>(
            &plan.s_list[0..=i],
            &mut twiddles[last_stage_twiddle_start..],
        );
        if i + 1 < plan.s_list.len() {
            last_stage_twiddle_start += stage_twiddle_count;
        }
    }

    // dif fft
    ntt_dif_dit::<P, false>(plan, &mut x[g..], &twiddles);
    ntt_dif_dit::<P, false>(plan, &mut y[g..], &twiddles);

    // naive multiplication
    let (mut i, mut last_stage_twiddle_index, mut position_in_last_radix) =
        (0, last_stage_twiddle_start, 0);
    let mut current_twiddle = Arith::<P>::R;
    let twiddle_step = Arith::<P>::mpowmod(
        NttKernelImpl::<P, false>::ROOTR,
        Arith::<P>::MAX_NTT_LEN / last_radix,
    );

    while i < plan.n {
        // We accumulate the results just before `x_block` for cache friendliness.
        let (out_block, x_block) = x[i..i + 2 * g].split_at_mut(g);
        let y_block = &y[i + g..i + 2 * g];
        conv_base::<P>(out_block, x_block, y_block, current_twiddle);
        i += g;
        position_in_last_radix += 1;
        if position_in_last_radix == last_radix {
            last_stage_twiddle_index += 1;
            position_in_last_radix = 0;
            current_twiddle = twiddles[last_stage_twiddle_index];
        } else {
            current_twiddle = Arith::<P>::mmulmod(current_twiddle, twiddle_step);
        }
    }

    // dit fft
    let mut twiddle_offset = 0;
    for i in (0..plan.s_list.len()).rev() {
        twiddle_offset +=
            calc_twiddle_factors::<P, true>(&plan.s_list[0..=i], &mut twiddles[twiddle_offset..]);
    }
    ntt_dif_dit::<P, true>(plan, x, &twiddles);
}

////////////////////////////////////////////////////////////////////////////////

// Prime selection balances two goals: a modulus close to 2^64 increases the
// packing density, while large powers of 2, 3, and 5 dividing P - 1 provide
// more supported NTT lengths. The NTT primes used by this module were found
// by exhaustive search subject to both requirements.
const P1: u64 = 14_259_017_916_245_606_401; // Max NTT length = 2^22 * 3^21 * 5^2 = 1_096_847_532_018_892_800
const P2: u64 = 17_984_575_660_032_000_001; // Max NTT length = 2^19 * 3^17 * 5^6 = 1_057_916_215_296_000_000
const P3: u64 = 17_995_154_822_184_960_001; // Max NTT length = 2^17 * 3^22 * 5^4 = 2_570_736_403_169_280_000

const P1INV_R_MOD_P2: u64 = Arith::<P2>::mmulmod(Arith::<P2>::R2, arith::invmod(P1, P2));

const P1P2INV_R_MOD_P3: u64 = Arith::<P3>::mmulmod(
    Arith::<P3>::R2,
    arith::invmod((P1 as u128 * P2 as u128 % P3 as u128) as u64, P3),
);

const P1_R_MOD_P3: u64 = Arith::<P3>::mmulmod(Arith::<P3>::R2, P1);

const P1P2_LO: u64 = (P1 as u128 * P2 as u128) as u64;
const P1P2_HI: u64 = ((P1 as u128 * P2 as u128) >> 64) as u64;

// Propagates carry from the beginning to the end of acc,
// and returns the resulting carry if it is nonzero.
fn propagate_carry(acc: &mut [u64], mut carry: u64) -> u64 {
    for x in acc {
        let (v, overflow) = x.overflowing_add(carry);
        (*x, carry) = (v, u64::from(overflow));
        if !overflow {
            break;
        }
    }
    carry
}

fn mac3_two_primes(acc: &mut [u64], b: &[u64], c: &[u64], bits: u64) {
    fn pack_into(src: &[u64], dst1: &mut [u64], dst2: &mut [u64], bits: u64) {
        let mut p = 0u64;
        let mut dst = zip(dst1, dst2);
        let mut x = 0u64;
        let mask = (1u64 << bits) - 1;
        for v in src {
            let mut k = 0;
            while k < 64 {
                x |= (v >> k) << p;
                let q = 64 - k;
                if p + q >= bits {
                    let out = x & mask;
                    let (pdst1, pdst2) = dst.next().unwrap();
                    *pdst1 = out;
                    *pdst2 = out;
                    x = 0;
                    k = k + bits - p;
                    p = 0;
                } else {
                    p += q;
                    break;
                }
            }
        }
        if p > 0 {
            let out = x & mask;
            let (pdst1, pdst2) = dst.next().unwrap();
            *pdst1 = out;
            *pdst2 = out;
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

    conv::<P2>(
        &plan_x,
        &mut x,
        b_len,
        &mut r[..plan_x.g + plan_x.n],
        c_len,
        arith::invmod(P3, P2),
    );

    conv::<P3>(
        &plan_y,
        &mut y,
        b_len,
        &mut s[..plan_y.g + plan_y.n],
        c_len,
        Arith::<P3>::submod(0, arith::invmod(P2, P3)),
    );

    // merge the results in {x, y} into r (process carry along the way)
    let mask = (1u64 << bits) - 1;
    let mut carry: u128 = 0;
    let (mut j, mut p) = (0usize, 0u64);
    let mut bitbuf: u64 = 0;
    let mut carry_acc: u64 = 0;
    for i in 0..min_len {
        // extract the convolution result
        let (a, b) = (x[i], y[i]);
        let (mut v, overflow) =
            (a as u128 * P3 as u128 + carry).overflowing_sub(b as u128 * P2 as u128);
        if overflow {
            v = v.wrapping_add(P2 as u128 * P3 as u128);
        }
        carry = v >> bits;

        // write to bitbuf
        let out = (v as u64) & mask;
        bitbuf |= out << p;
        p += bits;
        if p >= 64 {
            // flush bitbuf to the output buffer
            let (w, overflow1) = bitbuf.overflowing_add(carry_acc);
            let (w, overflow2) = acc[j].overflowing_add(w);
            acc[j] = w;
            carry_acc = u64::from(overflow1 || overflow2);

            // roll over
            (j, p) = (j + 1, p - 64);
            bitbuf = out >> (bits - p);
        }
    }

    // Process remaining carries. The addition carry_acc + bitbuf should not overflow
    // since bitbuf is underfilled and carry_acc is always 0 or 1.
    propagate_carry(&mut acc[j..], carry_acc + bitbuf);
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

    // convolution with modulo P1
    for i in 0..b.len() {
        x[plan_x.g + i] = if b[i] >= P1 { b[i] - P1 } else { b[i] };
    }
    for i in 0..c.len() {
        r[plan_x.g + i] = if c[i] >= P1 { c[i] - P1 } else { c[i] };
    }
    conv::<P1>(
        &plan_x,
        &mut x,
        b.len(),
        &mut r[..plan_x.g + plan_x.n],
        c.len(),
        1,
    );

    // convolution with modulo P2
    for i in 0..b.len() {
        y[plan_y.g + i] = if b[i] >= P2 { b[i] - P2 } else { b[i] };
    }
    for i in 0..c.len() {
        r[plan_y.g + i] = if c[i] >= P2 { c[i] - P2 } else { c[i] };
    }
    (&mut r[plan_y.g..])[c.len()..plan_y.n].fill(0u64);
    conv::<P2>(
        &plan_y,
        &mut y,
        b.len(),
        &mut r[..plan_y.g + plan_y.n],
        c.len(),
        1,
    );

    // convolution with modulo P3
    for i in 0..b.len() {
        z[plan_z.g + i] = if b[i] >= P3 { b[i] - P3 } else { b[i] };
    }
    for i in 0..c.len() {
        r[plan_z.g + i] = if c[i] >= P3 { c[i] - P3 } else { c[i] };
    }
    (&mut r[plan_z.g..])[c.len()..plan_z.n].fill(0u64);
    conv::<P3>(
        &plan_z,
        &mut z,
        b.len(),
        &mut r[..plan_z.g + plan_z.n],
        c.len(),
        1,
    );

    // merge the results in {x, y, z} into acc (process carry along the way)
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
        let out_12 = P1P2_HI as u128 * vcmv as u128
            + (out_01 >> 64)
            + if overflow { 1u128 << 64 } else { 0 };

        let (v, overflow) = acc[i].overflowing_add(out_0);
        acc[i] = v;
        carry = out_12 + u128::from(overflow);
    }
    propagate_carry(&mut acc[min_len..], carry as u64);
}

fn mac3_u64(acc: &mut [u64], b: &[u64], c: &[u64]) {
    // Computes the largest safe packed-digit width from the maximum operand length,
    // ensuring the convolution coefficients fit below P2 * P3.
    const fn compute_bits(l: u64) -> u64 {
        let total_bits = l * 64;
        let (mut lo, mut hi) = (42, 62);
        while lo < hi {
            let mid = (lo + hi + 1) / 2;
            let single_digit_max_val = (1u64 << mid) - 1;
            let l_corrected = (total_bits + mid - 1) / mid;
            let (lhs, overflow) = (single_digit_max_val as u128)
                .pow(2)
                .overflowing_mul(l_corrected as u128);
            if !overflow && lhs < P2 as u128 * P3 as u128 {
                lo = mid;
            } else {
                hi = mid - 1;
            }
        }
        lo
    }

    // We have two choices:
    //     1. NTT with two primes.
    //     2. NTT with three primes.
    // Obviously we want to do only two passes for efficiency, not three.
    // However, the number of bits per u64 we can pack for NTT
    // depends on the length of the arrays being multiplied (convolved).
    // If the arrays are too long, the resulting values may exceed the
    // modulus range P2 * P3, which leads to incorrect results.
    // Hence, we compute the maximum number of bits that can be packed in one
    // u64, and use it to determine whether to use two-prime or three-prime.
    // Since we can pack 64 bits per u64 in three-prime NTT, the effective
    // number of bits in three-prime NTT is 64/3 = 21.3333..., which means
    // two-prime NTT can only do better when at least 43 bits per u64 can
    // be packed into each u64.
    let max_cnt = max(b.len(), c.len()) as u64;
    let bits = compute_bits(max_cnt);
    if bits >= 43 {
        // can pack more effective bits per u64 with two primes than with three primes
        mac3_two_primes(acc, b, c, bits);
    } else {
        // can pack at most 21 effective bits per u64, which is worse than
        // 64/3 = 21.3333.. effective bits per u64 achieved with three primes
        mac3_three_primes(acc, b, c);
    }
}

////////////////////////////////////////////////////////////////////////////////

cfg_digit! {
    pub fn mac3(acc: &mut [BigDigit], b: &[BigDigit], c: &[BigDigit]) {
        fn bigdigit_to_u64(src: &[BigDigit], is_acc: bool) -> Vec<u64> {
            let mut out = vec![0u64; (src.len() + 1) / 2 + is_acc as usize];
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

        // convert to u64 => process => convert back to BigDigit (u32)
        let mut acc_u64 = bigdigit_to_u64(acc, true);
        mac3_u64(&mut acc_u64, &bigdigit_to_u64(b, false), &bigdigit_to_u64(c, false));
        u64_to_bigdigit(&acc_u64, acc);
    }
    pub fn mac3(acc: &mut [BigDigit], b: &[BigDigit], c: &[BigDigit]) {
       mac3_u64(acc, b, c);
    }
}

////////////////////////////////////////////////////////////////////////////////

#[cfg(feature = "bench-internals")]
pub(super) fn benchmark_plan_build<const PRIME_INDEX: usize>(
    min_len: usize,
) -> (usize, usize, usize) {
    let plan = match PRIME_INDEX {
        1 => NttPlan::build::<P1>(min_len),
        2 => NttPlan::build::<P2>(min_len),
        3 => NttPlan::build::<P3>(min_len),
        _ => panic!("PRIME_INDEX must be 1, 2, or 3"),
    };
    (plan.n, plan.g, plan.s_list.len())
}
