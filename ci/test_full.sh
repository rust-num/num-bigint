#!/bin/bash

set -ex

echo Testing num-bigint on rustc ${TRAVIS_RUST_VERSION}

FEATURES="serde i128 u64_digit"

# num-bigint should build and test everywhere.
cargo build --verbose
cargo test --verbose

# It should build with minimal features too.
cargo build --no-default-features --features="std"
cargo test --no-default-features --features="std"

# Each isolated feature should also work everywhere.
for feature in $FEATURES; do
  cargo build --verbose --no-default-features --features="std $feature"
  cargo test --verbose --no-default-features --features="std $feature"
done

# test all supported features together
cargo build --features="std $FEATURES"
cargo test --features="std $FEATURES"

# make sure benchmarks can be built
if [[ "$TRAVIS_RUST_VERSION" == "nightly" ]]; then
  cargo bench --all-features --no-run
fi
