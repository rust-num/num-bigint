use std::env;

fn main() {
    let ptr_width = env::var("CARGO_CFG_TARGET_POINTER_WIDTH");
    let u64_digit = ptr_width
        .as_ref()
        .map(|x| x == "64" || x == "128")
        .unwrap_or(false);

    if u64_digit {
        println!("cargo:rustc-cfg=u64_digit");
    }

    println!("cargo:rerun-if-changed=build.rs");
}
