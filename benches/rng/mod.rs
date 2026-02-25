use rand_0_10::rngs::Xoshiro128PlusPlus;
use rand_0_10::{Rng, SeedableRng};

pub(crate) fn get_rng() -> impl Rng {
    Xoshiro128PlusPlus::from_seed(*b"benching bigints")
}
