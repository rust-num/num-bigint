use num_bigint::BigInt;

#[global_allocator]
static ALLOC: dhat::Alloc = dhat::Alloc;

#[test]
fn test_biguint_allocation_size() {
    let _profiler = dhat::Profiler::builder().testing().build();
    let big: BigInt = "-1234567898765432123456789876543212345678987654321".parse().unwrap();
    let stats = dhat::HeapStats::get();
    assert_eq!(stats.curr_bytes, big.allocation_size());
}
