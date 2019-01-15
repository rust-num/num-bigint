/// Generic trait to implement modular inverse
pub trait ModInverse<R: Sized>: Sized {
    /// Function to calculate the [modular multiplicative
    /// inverse](https://en.wikipedia.org/wiki/Modular_multiplicative_inverse) of an integer *a* modulo *m*.
    ///
    /// TODO: references
    /// Returns the modular inverse of `self`.
    /// If none exists it returns `None`.
    fn mod_inverse(self, m: R) -> Option<Self>;
}
