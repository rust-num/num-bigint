extern crate autocfg;

fn main() {
    autocfg::rerun_path("build.rs");
    autocfg::emit_possibility(HAS_DERIVE);
    if std::env::var_os("CARGO_FEATURE_RUSTC_SERIALIZE").is_some() {
        let ac = autocfg::new();

        // These built-in derives are being removed! (rust-lang/rust#134272)
        //
        // It's hard to directly probe for `derive(RustcDecodable, RustcEncodable)`, because that
        // depends on the external `rustc-serialize` dependency. They're in `prelude::v1` where we
        // can probe by path, but ironically only on relatively new versions, so we're also using
        // *inaccessible* `rust_2024` as a proxy for older versions.
        if ac.probe_raw(PRELUDE_DERIVE).is_ok() || !ac.probe_path(RUST_2024) {
            autocfg::emit(HAS_DERIVE);
        } else {
            println!("cargo:warning=rustc-serialize is not supported by the current compiler");
        }
    }
}

const HAS_DERIVE: &str = "has_derive_rustc_serialize";

const PRELUDE_DERIVE: &str = "
#[allow(soft_unstable, deprecated)]
pub use std::prelude::v1::{RustcDecodable, RustcEncodable};
";

const RUST_2024: &str = "std::prelude::rust_2024";
