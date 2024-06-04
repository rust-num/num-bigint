#!/bin/bash

set -e

CRATE=num-bigint
MSRV=1.60

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

STD_FEATURES=(inline arbitrary quickcheck rand serde)
NO_STD_FEATURES=(inline serde rand)
echo "Testing supported features: ${STD_FEATURES[*]}"
if [ -n "${NO_STD_FEATURES[*]}" ]; then
  echo " no_std supported features: ${NO_STD_FEATURES[*]}"
fi

# arbitrary 1.1.4 started using array::from_fn
check_version 1.63.0 || cargo update -p arbitrary --precise 1.1.3

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
  *serde*) (
      cd ci/big_serde
      cargo test
    ) ;;&
  *rand*) cargo test --manifest-path ci/big_rand/Cargo.toml ;;&
  *quickcheck*) (
      cd ci/big_quickcheck
      cargo test
    ) ;;&
esac
