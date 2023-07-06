use crate::{BigInt, Sign::*};

impl BigInt {
    pub fn true_div(&self, other: &BigInt) -> f64 {
        let d = self.data.true_div(&other.data);
        match (self.sign, other.sign) {
            (Plus, Plus) | (NoSign, Plus) | (Minus, Minus) => d,
            (Plus, Minus) | (NoSign, Minus) | (Minus, Plus) => -d,
            (_, NoSign) => unreachable!(),
        }
    }
}
