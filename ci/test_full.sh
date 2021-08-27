#!/bin/bash

set -e

CRATE=num-bigint
MSRV=1.31

get_rust_version() {
  local array=($(rustc --version));
  echo "${array[1]}";
  return 0;
}
RUST_VERSION=$(get_rust_version)

check_version() {
  IFS=. read -ra rust <<< "$RUST_VERSION"
  IFS=. read -ra want <<< "$1"
  [[ "${rust[0]}" -gt "${want[0]}" ||
   ( "${rust[0]}" -eq "${want[0]}" &&
     "${rust[1]}" -ge "${want[1]}" )
  ]]
}

echo "Testing $CRATE on rustc $RUST_VERSION"
if ! check_version $MSRV ; then
  echo "The minimum for $CRATE is rustc $MSRV"
  exit 1
fi

STD_FEATURES=(serde)
check_version 1.36 && STD_FEATURES+=(rand)
check_version 1.36 && NO_STD_FEATURES=(serde rand)
check_version 1.40 && STD_FEATURES+=(arbitrary)
check_version 1.46 && STD_FEATURES+=(quickcheck)
echo "Testing supported features: ${STD_FEATURES[*]}"
if [ -n "${NO_STD_FEATURES[*]}" ]; then
  echo " no_std supported features: ${NO_STD_FEATURES[*]}"
fi

# arbitrary 1.0.1 added const-generic arrays, which requires Rust 1.51
check_version 1.51.0 || cargo update -p arbitrary --precise 1.0.0

set -x

# test the default with std
cargo build
cargo test

# test each isolated feature with std
for feature in ${STD_FEATURES[*]}; do
  cargo build --no-default-features --features="std $feature"
  cargo test --no-default-features --features="std $feature"
done

# test all supported features with std
cargo build --no-default-features --features="std ${STD_FEATURES[*]}"
cargo test --no-default-features --features="std ${STD_FEATURES[*]}"


if [ -n "${NO_STD_FEATURES[*]}" ]; then
  # test minimal `no_std`
  cargo build --no-default-features
  cargo test --no-default-features

  # test each isolated feature without std
  for feature in ${NO_STD_FEATURES[*]}; do
    cargo build --no-default-features --features="$feature"
    cargo test --no-default-features --features="$feature"
  done

  # test all supported features without std
  cargo build --no-default-features --features="${NO_STD_FEATURES[*]}"
  cargo test --no-default-features --features="${NO_STD_FEATURES[*]}"
fi


# make sure benchmarks can be built and sanity-tested
if rustc --version | grep -q nightly; then
    cargo test --all-features --benches
fi

case "${STD_FEATURES[*]}" in
  *serde*) cargo test --manifest-path ci/big_serde/Cargo.toml ;;&
  *rand*) cargo test --manifest-path ci/big_rand/Cargo.toml ;;&
  *quickcheck*) cargo test --manifest-path ci/big_quickcheck/Cargo.toml ;;&
esac
