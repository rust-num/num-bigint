# Release 0.1.42

- [num-bigint now has its own source repository][num-356] at [rust-num/num-bigint][home].
- [`lcm` now avoids creating a large intermediate product][num-350].
- [`gcd` now uses Stein's algorithm][15] with faster shifts instead of division.

**Contributors**: @cuviper, @Emerentius, @mhogrefe

[home]: https://github.com/rust-num/num-bigint
[num-350]: https://github.com/rust-num/num/pull/350
[num-356]: https://github.com/rust-num/num/pull/356
[15]: https://github.com/rust-num/num-bigint/pull/15


# Prior releases

No prior release notes were kept.  Thanks all the same to the many
contributors that have made this crate what it is!

