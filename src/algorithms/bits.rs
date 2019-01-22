use std::mem;

/// Find last set bit
/// fls(0) == 0, fls(u32::MAX) == 32
pub fn fls<T: num_traits::PrimInt>(v: T) -> usize {
    mem::size_of::<T>() * 8 - v.leading_zeros() as usize
}

pub fn ilog2<T: num_traits::PrimInt>(v: T) -> usize {
    fls(v) - 1
}
