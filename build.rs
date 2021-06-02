extern crate autocfg;
use std::env;

fn main() {
    let ac = autocfg::new();

    if ac.probe_type("i128") || env::var_os("CARGO_FEATURE_I128").is_some() {
        println!("cargo:rustc-cfg=has_i128");
    }
}
