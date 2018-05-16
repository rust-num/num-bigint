#!/bin/bash

set -ex

echo Testing num-bigint on rustc ${TRAVIS_RUST_VERSION}

# num-bigint should build and test everywhere.
cargo build --verbose
cargo test --verbose

# It should build with minimal features too.
cargo build --no-default-features
cargo test --no-default-features

# Each isolated feature should also work everywhere.
for feature in rand serde; do
  cargo build --verbose --no-default-features --features="$feature"
  cargo test --verbose --no-default-features --features="$feature"
done

# test `i128` and all features together
if [[ "$TRAVIS_RUST_VERSION" =~ ^(nightly|beta|stable)$ ]]; then
  cargo build --verbose --features=i128
  cargo test --verbose --features=i128

  cargo build --all-features
  cargo test --all-features
else
  # all except `i128`
  cargo build --features="rand serde"
  cargo test --features="rand serde"
fi
