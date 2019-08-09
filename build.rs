extern crate autocfg;

use std::env;
use std::error::Error;
use std::fs::File;
use std::io::Write;
use std::path::Path;

fn main() {
    let ac = autocfg::new();

    if ac.probe_type("i128") {
        autocfg::emit("has_i128");

        let pointer_width = env::var("CARGO_CFG_TARGET_POINTER_WIDTH");
        if pointer_width.as_ref().map(String::as_str) == Ok("64") {
            autocfg::emit("u64_digit");
        }
    } else if env::var_os("CARGO_FEATURE_I128").is_some() {
        panic!("i128 support was not detected!");
    }

    autocfg::rerun_path(file!());

    write_radix_bases().unwrap();
}

#[allow(unknown_lints, bare_trait_objects)]
fn write_radix_bases() -> Result<(), Box<Error>> {
    let out_dir = env::var("OUT_DIR")?;
    let dest_path = Path::new(&out_dir).join("radix_bases.rs");
    let mut f = File::create(&dest_path)?;

    for &bits in &[16, 32, 64] {
        let max = if bits < 64 {
            (1 << bits) - 1
        } else {
            std::u64::MAX
        };

        writeln!(f, "#[deny(overflowing_literals)]")?;
        writeln!(
            f,
            "pub static BASES_{bits}: [(u{bits}, usize); 257] = [",
            bits = bits
        )?;
        for radix in 0u64..257 {
            let (base, power) = if radix == 0 || radix.is_power_of_two() {
                (0, 0)
            } else {
                let mut power = 1;
                let mut base = radix;

                while let Some(b) = base.checked_mul(radix) {
                    if b > max {
                        break;
                    }
                    base = b;
                    power += 1;
                }
                (base, power)
            };
            writeln!(f, "    ({}, {}), // {}", base, power, radix)?;
        }
        writeln!(f, "];")?;
    }

    Ok(())
}
