#!/bin/sh
# Use rustup to locally run the same suite of tests as .github/workflows/
# (You should first install/update all of the versions below.)

set -ex

ci=$(dirname $0)
for version in 1.31.0 1.32.0 1.34.0 1.36.0 stable beta nightly; do
    rustup run "$version" "$ci/test_full.sh"
done
