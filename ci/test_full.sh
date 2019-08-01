#!/bin/bash

set -ex

echo Testing num-bigint on rustc ${TRAVIS_RUST_VERSION}

case "$TRAVIS_RUST_VERSION" in
  1.1[5-9].* | 1.2[0-1].*) STD_FEATURES="serde" ;;
  1.2[2-5].*) STD_FEATURES="serde rand" ;;
  1.2[6-9].* | 1.30.*) STD_FEATURES="serde rand i128" ;;
  *) STD_FEATURES="serde rand i128 quickcheck quickcheck_macros" ;;
esac

case "$TRAVIS_RUST_VERSION" in
  1.1[5-9].* | 1.2[0-9].* | 1.3[0-5].*) ;;
  *) NO_STD_FEATURES="i128" ;;
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
