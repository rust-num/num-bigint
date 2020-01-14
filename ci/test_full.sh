#!/bin/bash

set -ex

echo Testing num-bigint on rustc ${TRAVIS_RUST_VERSION}

case "$TRAVIS_RUST_VERSION" in
  1.31.*) STD_FEATURES="serde" ;;
  1.3[23].*) STD_FEATURES="serde rand" ;;
  *) STD_FEATURES="serde rand quickcheck" ;;
esac

case "$TRAVIS_RUST_VERSION" in
  1.3[1-5].*) ;;
  *) NO_STD_FEATURES="serde rand" ;;
esac

# num-bigint should build and test everywhere.
cargo build --verbose
cargo test --verbose

# It should build with minimal features too.
cargo build --no-default-features --features="std"
cargo test --no-default-features --features="std"

# Each isolated feature should also work everywhere.
for feature in $STD_FEATURES; do
  cargo build --verbose --no-default-features --features="std $feature"
  cargo test --verbose --no-default-features --features="std $feature"
done

# test all supported features together
cargo build --features="std $STD_FEATURES"
cargo test --features="std $STD_FEATURES"

if test -n "${NO_STD_FEATURES:+true}"; then
  # It should build with minimal features too.
  cargo build --no-default-features
  cargo test --no-default-features

  # Each isolated feature should also work everywhere.
  for feature in $NO_STD_FEATURES; do
    cargo build --verbose --no-default-features --features="$feature"
    cargo test --verbose --no-default-features --features="$feature"
  done

  # test all supported features together
  cargo build --no-default-features --features="$NO_STD_FEATURES"
  cargo test --no-default-features --features="$NO_STD_FEATURES"
fi

# make sure benchmarks can be built
if [[ "$TRAVIS_RUST_VERSION" == "nightly" ]]; then
  cargo bench --all-features --no-run
fi

case "$STD_FEATURES" in
  *serde*) cargo test --manifest-path ci/big_serde/Cargo.toml ;;&
  *rand*) cargo test --manifest-path ci/big_rand/Cargo.toml ;;&
  *quickcheck*) cargo test --manifest-path ci/big_quickcheck/Cargo.toml ;;&
esac
