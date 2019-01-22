use std::cmp;
use std::iter::repeat;

use crate::algorithms::{adc, add2, sub2, sub_sign};
use crate::big_digit::{BigDigit, DoubleBigDigit, BITS};
use crate::bigint::Sign::{Minus, NoSign, Plus};
use crate::{BigInt, BigUint};

#[inline]
pub fn mac_with_carry(a: BigDigit, b: BigDigit, c: BigDigit, acc: &mut DoubleBigDigit) -> BigDigit {
    *acc += a as DoubleBigDigit;
    *acc += (b as DoubleBigDigit) * (c as DoubleBigDigit);
    let lo = *acc as BigDigit;
    *acc >>= BITS;
    lo
}

/// Three argument multiply accumulate:
/// acc += b * c
pub fn mac_digit(acc: &mut [BigDigit], b: &[BigDigit], c: BigDigit) {
    if c == 0 {
        return;
    }

    let mut carry = 0;
    let (a_lo, a_hi) = acc.split_at_mut(b.len());

    for (a, &b) in a_lo.iter_mut().zip(b) {
        *a = mac_with_carry(*a, b, c, &mut carry);
    }

    let mut a = a_hi.iter_mut();
    while carry != 0 {
        let a = a.next().expect("carry overflow during multiplication!");
        *a = adc(*a, 0, &mut carry);
    }
}

/// Three argument multiply accumulate:
/// acc += b * c
pub fn mac3(acc: &mut [BigDigit], b: &[BigDigit], c: &[BigDigit]) {
    let (x, y) = if b.len() < c.len() { (b, c) } else { (c, b) };

    // We use three algorithms for different input sizes.
    //
    // - For small inputs, long multiplication is fastest.
    // - Next we use Karatsuba multiplication (Toom-2), which we have optimized
    //   to avoid unnecessary allocations for intermediate values.
    // - For the largest inputs we use Toom-3, which better optimizes the
    //   number of operations, but uses more temporary allocations.
    //
    // The thresholds are somewhat arbitrary, chosen by evaluating the results
    // of `cargo bench --bench bigint multiply`.

    if x.len() <= 32 {
        // Long multiplication:
        for (i, xi) in x.iter().enumerate() {
            mac_digit(&mut acc[i..], y, *xi);
        }
    } else if x.len() <= 256 {
        /*
         * Karatsuba multiplication:
         *
         * The idea is that we break x and y up into two smaller numbers that each have about half
         * as many digits, like so (note that multiplying by b is just a shift):
         *
         * x = x0 + x1 * b
         * y = y0 + y1 * b
         *
         * With some algebra, we can compute x * y with three smaller products, where the inputs to
         * each of the smaller products have only about half as many digits as x and y:
         *
         * x * y = (x0 + x1 * b) * (y0 + y1 * b)
         *
         * x * y = x0 * y0
         *       + x0 * y1 * b
         *       + x1 * y0 * b
         *       + x1 * y1 * b^2
         *
         * Let p0 = x0 * y0 and p2 = x1 * y1:
         *
         * x * y = p0
         *       + (x0 * y1 + x1 * y0) * b
         *       + p2 * b^2
         *
         * The real trick is that middle term:
         *
         *         x0 * y1 + x1 * y0
         *
         *       = x0 * y1 + x1 * y0 - p0 + p0 - p2 + p2
         *
         *       = x0 * y1 + x1 * y0 - x0 * y0 - x1 * y1 + p0 + p2
         *
         * Now we complete the square:
         *
         *       = -(x0 * y0 - x0 * y1 - x1 * y0 + x1 * y1) + p0 + p2
         *
         *       = -((x1 - x0) * (y1 - y0)) + p0 + p2
         *
         * Let p1 = (x1 - x0) * (y1 - y0), and substitute back into our original formula:
         *
         * x * y = p0
         *       + (p0 + p2 - p1) * b
         *       + p2 * b^2
         *
         * Where the three intermediate products are:
         *
         * p0 = x0 * y0
         * p1 = (x1 - x0) * (y1 - y0)
         * p2 = x1 * y1
         *
         * In doing the computation, we take great care to avoid unnecessary temporary variables
         * (since creating a BigUint requires a heap allocation): thus, we rearrange the formula a
         * bit so we can use the same temporary variable for all the intermediate products:
         *
         * x * y = p2 * b^2 + p2 * b
         *       + p0 * b + p0
         *       - p1 * b
         *
         * The other trick we use is instead of doing explicit shifts, we slice acc at the
         * appropriate offset when doing the add.
         */

        /*
         * When x is smaller than y, it's significantly faster to pick b such that x is split in
         * half, not y:
         */
        let b = x.len() / 2;
        let (x0, x1) = x.split_at(b);
        let (y0, y1) = y.split_at(b);

        /*
         * We reuse the same BigUint for all the intermediate multiplies and have to size p
         * appropriately here: x1.len() >= x0.len and y1.len() >= y0.len():
         */
        let len = x1.len() + y1.len() + 1;
        let mut p = BigUint {
            data: smallvec![0; len],
        };

        // p2 = x1 * y1
        mac3(&mut p.data[..], x1, y1);

        // Not required, but the adds go faster if we drop any unneeded 0s from the end:
        p.normalize();

        add2(&mut acc[b..], &p.data[..]);
        add2(&mut acc[b * 2..], &p.data[..]);

        // Zero out p before the next multiply:
        p.data.truncate(0);
        p.data.extend(repeat(0).take(len));

        // p0 = x0 * y0
        mac3(&mut p.data[..], x0, y0);
        p.normalize();

        add2(&mut acc[..], &p.data[..]);
        add2(&mut acc[b..], &p.data[..]);

        // p1 = (x1 - x0) * (y1 - y0)
        // We do this one last, since it may be negative and acc can't ever be negative:
        let (j0_sign, j0) = sub_sign(x1, x0);
        let (j1_sign, j1) = sub_sign(y1, y0);

        match j0_sign * j1_sign {
            Plus => {
                p.data.truncate(0);
                p.data.extend(repeat(0).take(len));

                mac3(&mut p.data[..], &j0.data[..], &j1.data[..]);
                p.normalize();

                sub2(&mut acc[b..], &p.data[..]);
            }
            Minus => {
                mac3(&mut acc[b..], &j0.data[..], &j1.data[..]);
            }
            NoSign => (),
        }
    } else {
        // Toom-3 multiplication:
        //
        // Toom-3 is like Karatsuba above, but dividing the inputs into three parts.
        // Both are instances of Toom-Cook, using `k=3` and `k=2` respectively.
        //
        // The general idea is to treat the large integers digits as
        // polynomials of a certain degree and determine the coefficients/digits
        // of the product of the two via interpolation of the polynomial product.
        let i = y.len() / 3 + 1;

        let x0_len = cmp::min(x.len(), i);
        let x1_len = cmp::min(x.len() - x0_len, i);

        let y0_len = i;
        let y1_len = cmp::min(y.len() - y0_len, i);

        // Break x and y into three parts, representating an order two polynomial.
        // t is chosen to be the size of a digit so we can use faster shifts
        // in place of multiplications.
        //
        // x(t) = x2*t^2 + x1*t + x0
        let x0 = BigInt::from_slice_native(Plus, &x[..x0_len]);
        let x1 = BigInt::from_slice_native(Plus, &x[x0_len..x0_len + x1_len]);
        let x2 = BigInt::from_slice_native(Plus, &x[x0_len + x1_len..]);

        // y(t) = y2*t^2 + y1*t + y0
        let y0 = BigInt::from_slice_native(Plus, &y[..y0_len]);
        let y1 = BigInt::from_slice_native(Plus, &y[y0_len..y0_len + y1_len]);
        let y2 = BigInt::from_slice_native(Plus, &y[y0_len + y1_len..]);

        // Let w(t) = x(t) * y(t)
        //
        // This gives us the following order-4 polynomial.
        //
        // w(t) = w4*t^4 + w3*t^3 + w2*t^2 + w1*t + w0
        //
        // We need to find the coefficients w4, w3, w2, w1 and w0. Instead
        // of simply multiplying the x and y in total, we can evaluate w
        // at 5 points. An n-degree polynomial is uniquely identified by (n + 1)
        // points.
        //
        // It is arbitrary as to what points we evaluate w at but we use the
        // following.
        //
        // w(t) at t = 0, 1, -1, -2 and inf
        //
        // The values for w(t) in terms of x(t)*y(t) at these points are:
        //
        // let a = w(0)   = x0 * y0
        // let b = w(1)   = (x2 + x1 + x0) * (y2 + y1 + y0)
        // let c = w(-1)  = (x2 - x1 + x0) * (y2 - y1 + y0)
        // let d = w(-2)  = (4*x2 - 2*x1 + x0) * (4*y2 - 2*y1 + y0)
        // let e = w(inf) = x2 * y2 as t -> inf

        // x0 + x2, avoiding temporaries
        let p = &x0 + &x2;

        // y0 + y2, avoiding temporaries
        let q = &y0 + &y2;

        // x2 - x1 + x0, avoiding temporaries
        let p2 = &p - &x1;

        // y2 - y1 + y0, avoiding temporaries
        let q2 = &q - &y1;

        // w(0)
        let r0 = &x0 * &y0;

        // w(inf)
        let r4 = &x2 * &y2;

        // w(1)
        let r1 = (p + x1) * (q + y1);

        // w(-1)
        let r2 = &p2 * &q2;

        // w(-2)
        let r3 = ((p2 + x2) * 2 - x0) * ((q2 + y2) * 2 - y0);

        // Evaluating these points gives us the following system of linear equations.
        //
        //  0  0  0  0  1 | a
        //  1  1  1  1  1 | b
        //  1 -1  1 -1  1 | c
        // 16 -8  4 -2  1 | d
        //  1  0  0  0  0 | e
        //
        // The solved equation (after gaussian elimination or similar)
        // in terms of its coefficients:
        //
        // w0 = w(0)
        // w1 = w(0)/2 + w(1)/3 - w(-1) + w(2)/6 - 2*w(inf)
        // w2 = -w(0) + w(1)/2 + w(-1)/2 - w(inf)
        // w3 = -w(0)/2 + w(1)/6 + w(-1)/2 - w(1)/6
        // w4 = w(inf)
        //
        // This particular sequence is given by Bodrato and is an interpolation
        // of the above equations.
        let mut comp3: BigInt = (r3 - &r1) / 3;
        let mut comp1: BigInt = (r1 - &r2) / 2;
        let mut comp2: BigInt = r2 - &r0;
        comp3 = (&comp2 - comp3) / 2 + &r4 * 2;
        comp2 = comp2 + &comp1 - &r4;
        comp1 = comp1 - &comp3;

        // Recomposition. The coefficients of the polynomial are now known.
        //
        // Evaluate at w(t) where t is our given base to get the result.
        let result = r0
            + (comp1 << 32 * i)
            + (comp2 << 2 * 32 * i)
            + (comp3 << 3 * 32 * i)
            + (r4 << 4 * 32 * i);
        let result_pos = result.to_biguint().unwrap();
        add2(&mut acc[..], &result_pos.data);
    }
}
