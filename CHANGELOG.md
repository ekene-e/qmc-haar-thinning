# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project uses
[Semantic Versioning](https://semver.org/).

Because the library promises bit-for-bit replicable output for a given
`Options` and seed, any change that alters the points produced from the same
inputs is a **behavioural change** and is called out explicitly.

## [0.1.0] - 2026-09-11

### Added
- Header-only C++23 library `haar::Thinner` implementing online
  (1 + ε)-thinning with sign feedback (Haar-thinning, DFG+19 / Theorem 1.2 of
  the paper) and linear feedback (Theorem 1.3), an optional uniform random
  shift, hyperbolic-cross and box scale sets, fixed or adaptive resolution,
  and dense or sparse counter storage.
- Exact integer decision rule: identical output on every platform for the
  same options, seed and sample stream.
- `haar::thin(n, options)` convenience, `Thinner::offer` online interface,
  `Thinner::run` with any `std::uniform_random_bit_generator`.
- Model-based automatic choice of the linear-feedback bound B.
- Quality measures: exact star discrepancy (d ≤ 2), Warnock's L2-star
  discrepancy (any d), sample means.
- `haar-thin` command-line tool (generation, stream thinning of external
  samples, CSV/TSV/binary output, statistics, discrepancy report).
- `haar-bench` reproducing the discrepancy, integration-error, throughput and
  saturation tables of the README.
- Unit tests including a brute-force reference implementation of the
  algorithm, and golden digests for replicability.
- CMake package (`find_package(haar)`), install rules, examples.
