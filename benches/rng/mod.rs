use rand_core::{Rng, TryRng};

pub(crate) fn get_rng() -> impl Rng {
    XorShiftStar {
        a: 0x0123_4567_89AB_CDEF,
    }
}

/// Simple `Rng` for benchmarking without additional dependencies
struct XorShiftStar {
    a: u64,
}

impl TryRng for XorShiftStar {
    type Error = std::convert::Infallible;
    fn try_next_u32(&mut self) -> Result<u32, Self::Error> {
        Ok(self.try_next_u64()? as u32)
    }

    fn try_next_u64(&mut self) -> Result<u64, Self::Error> {
        // https://en.wikipedia.org/wiki/Xorshift#xorshift*
        self.a ^= self.a >> 12;
        self.a ^= self.a << 25;
        self.a ^= self.a >> 27;
        Ok(self.a.wrapping_mul(0x2545_F491_4F6C_DD1D))
    }

    fn try_fill_bytes(&mut self, dest: &mut [u8]) -> Result<(), Self::Error> {
        for chunk in dest.chunks_mut(8) {
            let bytes = self.try_next_u64()?.to_le_bytes();
            let slice = &bytes[..chunk.len()];
            chunk.copy_from_slice(slice);
        }
        Ok(())
    }

    // fn try_fill_bytes(&mut self, dest: &mut [u8]) -> Result<(), rand::Error> {
    //     Ok(self.fill_bytes(dest))
    // }
}
