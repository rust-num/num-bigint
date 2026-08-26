use crate::BigUint;
use num_traits::One;

/// All factorials that fit into a u64
const SMALL_FACTORIALS: [u64; 21] = [
    1,                   // 0! = 1
    1,                   // 1! = 1
    2,                   // 2! = 2
    6,                   // 3! = 6
    24,                  // 4! = 24
    120,                 // 5! = 120
    720,                 // 6! = 720
    5040,                // 7! = 5040
    40320,               // 8! = 40320
    362880,              // 9! = 362880
    3628800,             // 10! = 3628800
    39916800,            // 11! = 39916800
    479001600,           // 12! = 479001600
    6227020800,          // 13! = 6227020800
    87178291200,         // 14! = 87178291200
    1307674368000,       // 15! = 1307674368000
    20922789888000,      // 16! = 20922789888000
    355687428096000,     // 17! = 355687428096000
    6402373705728000,    // 18! = 6402373705728000
    121645100408832000,  // 19! = 121645100408832000
    2432902008176640000, // 20! = 2432902008176640000
];

/// The odd factors of each factorial: Each number in this table is odd, and
/// self[i] << k == i! for some k
const ODD_FACTORIALS: [u64; 25] = [
    1,                  // 0! = 1 << 0
    1,                  // 1! = 1 << 0
    1,                  // 2! = 1 << 1
    3,                  // 3! = 3 << 1
    3,                  // 4! = 3 << 3
    15,                 // 5! = 15 << 3
    45,                 // 6! = 45 << 4
    315,                // 7! = 315 << 4
    315,                // 8! = 315 << 7
    2835,               // 9! = 2835 << 7
    14175,              // 10! = 14175 << 8
    155925,             // 11! = 155925 << 8
    467775,             // 12! = 467775 << 10
    6081075,            // 13! = 6081075 << 10
    42567525,           // 14! = 42567525 << 11
    638512875,          // 15! = 638512875 << 11
    638512875,          // 16! = 638512875 << 15
    10854718875,        // 17! = 10854718875 << 15
    97692469875,        // 18! = 97692469875 << 16
    1856156927625,      // 19! = 1856156927625 << 16
    9280784638125,      // 20! = 9280784638125 << 18
    194896477400625,    // 21! = 194896477400625 << 18
    2143861251406875,   // 22! = 2143861251406875 << 19
    49308808782358125,  // 23! = 49308808782358125 << 19
    147926426347074375, // 24! = 147926426347074375 << 22
];

pub fn factorial(x: usize) -> BigUint {
    if let Some(&x) = SMALL_FACTORIALS.get(x) {
        x.into()
    } else {
        // This uses a neat trick: The number of 2s that are a factor of x! is
        // equal to x - x.count_ones()
        let result = odd_factorial(x);
        let power_two_count = x - x.count_ones() as usize;
        result << power_two_count
    }
}

#[inline]
fn odd_factorial(mut x: usize) -> BigUint {
    if let Option::Some(&x) = ODD_FACTORIALS.get(x) {
        return x.into();
    }
    // Algorithm: `odd_factorial(n) = odds(1..=n) * odd_factorial(n / 2)
    // Prime sieving could be used to improve this in the future, see
    // https://gmplib.org/manual/Factorial-Algorithm
    // This uses a loop instead of recursion, to improve speed & reduce stack
    // space
    let mut result = BigUint::one();
    while x as usize >= ODD_FACTORIALS.len() {
        let mut i = 3;
        while i <= x {
            result *= i;
            i += 2;
        }
        x /= 2;
    }
    result * ODD_FACTORIALS[x as usize]
}
