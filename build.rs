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

    if let Ok(arch) = env::var("CARGO_CFG_TARGET_ARCH") {
        if arch == "x86_64" || arch == "x86" {
            println!("cargo:rustc-cfg=use_addcarry");
        }
    }

    println!("cargo:rerun-if-changed=build.rs");
}
