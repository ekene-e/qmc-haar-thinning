# haar-thinning

Online **Haar-thinning**: turn a stream of (1 + ε) n i.i.d. uniform samples in
[0,1)^d into n points with polylogarithmic discrepancy and quasi-Monte Carlo
integration error, one sample at a time.

A header-only C++23 library, a command-line tool, tests and benchmarks
implementing

> Ekene Ezeunala, Agastya Vibhuti Jha, Haotian Jiang.
> **Quasi-Monte Carlo Beyond Hardy–Krause II: (1 + ε) n Samples Suffice.**

which builds on the Haar-thinning strategy of Dwivedi, Feldheim,
Gurel-Gurevich and Ramdas (*The power of online thinning in reducing
discrepancy*, PTRF 2019).

## What it does

Monte Carlo integration with n random points has error σ(f)/√n. Quasi-Monte
Carlo point sets reach Õ(1/n) but need a deterministic construction and a
function of bounded Hardy–Krause variation. Online thinning sits in between:
samples arrive one by one, the algorithm keeps or discards each (never two in
a row), and after consuming about (1 + ε/2) n samples it has retained n points
such that

* **integration error** is Õ_d(σ_SO(f)/n), where the smoothed-out variation
  σ_SO(f) can be far smaller than the Hardy–Krause variation
  (Theorem 1.2 — *uniformly-shifted Haar-thinning*), and
* **star discrepancy** is O_d(log^{d+1} n) (Theorem 1.3 — *linear-feedback
  Haar-thinning*).

Each decision costs Õ_d(1): the algorithm keeps a counter φ_t(H) for every
Haar function H up to a resolution, and rejects a candidate x with a
probability that depends on the counters of the Haar functions nonzero at x.

## Quick start

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build              # 5 test executables, ~10 s

build/haar-thin -n 100000 -d 2 > points.csv          # 100000 points in [0,1)^2
build/haar-thin -n 4096 -d 1 --report -q -o /dev/null # star discrepancy of the output
cat my_samples.txt | build/haar-thin --stream -d 3 -n 5000 > kept.txt
build/haar-bench --quick                              # reproduce the tables below
```

Requirements: a C++23 compiler (tested with GCC 13; `std::format` and
`std::expected` are used) and CMake ≥ 3.23. No dependencies.

### Library

```cpp
#include <haar/haar.hpp>

haar::Options opt;
opt.dim = 2;
std::vector<double> pts = haar::thin(100000, opt);   // 100000 x 2 doubles, row-major

double estimate = haar::mean({pts, opt.dim}, [](std::span<const double> x) {
    return std::sin(2 * std::numbers::pi * 32 * x[0]) * x[1];
});
```

Online, with your own sample source (the thinner never needs to know n in
advance; without `expected_n` it raises its resolution as points accumulate,
so every prefix of the output is a good point set):

```cpp
haar::Thinner thinner(opt);
std::vector<double> x(opt.dim);
while (thinner.size() < n) {
    my_source.fill(x);              // i.i.d. uniform in [0,1)^d
    if (thinner.offer(x)) use(x);   // true = retained; after a rejection the next x is always retained
}
```

`Thinner::run(n, gen)` accepts any `std::uniform_random_bit_generator`;
`Thinner::run(n)` uses the internal seeded generator. See `examples/`.

### Integrating into your build

* **FetchContent / add_subdirectory**: `add_subdirectory(haar-thinning)` then
  `target_link_libraries(app PRIVATE haar::haar)`. The tests, CLI and benchmarks
  are only built when this is the top-level project.
* **Installed package**: `cmake --install build --prefix /opt/haar`, then
  `find_package(haar 0.1 REQUIRED)` and link `haar::haar`.
* **Copy the headers**: everything is under `include/haar/`.

## Options

| field | default | meaning |
|---|---|---|
| `dim` | 1 | dimension d |
| `epsilon` | 0.5 | each sample is rejected with probability ≤ ε; about (1 + ε/2) n samples are consumed for n points |
| `feedback` | `linear` | `sign` = Haar-thinning of DFG+19 (eq. 3.1, Theorem 1.2); `linear` = linear-feedback Haar-thinning (eq. 4.1, Theorem 1.3) |
| `shift` | `random` | uniform random shift of the Haar system (Section 3.1); `none` or `fixed` (`shift_vector`) |
| `scale_set` | `hyperbolic` | Haar scale vectors used: `hyperbolic` (Σ j_i ≤ level, as in DFG+19) or `box` (max j_i ≤ level, the paper's Π_{≤ℓ}) |
| `level` | automatic | resolution L (or ℓ); automatic = ⌈log₂ n⌉ + `level_offset` from `expected_n`, or adaptive |
| `expected_n` | unset | intended n; fixes the level up front (`haar::thin` sets it) |
| `saturation_bound` | automatic | the bound B of the linear rule; smaller = stronger feedback |
| `on_saturation` | `clamp` | what to do when the linear field exceeds B: clamp the density (still a valid strategy) or `fail` (throw) |
| `max_functions_per_sample` | 2048 | work cap for the automatic level rule |
| `dense_budget_bytes` | 512 MiB | memory cap for the counter table (automatic level rule and dense/sparse choice) |
| `storage` | `automatic` | dense array if it fits the budget, otherwise a hash table |
| `seed` | `0x5EED` | seed for coins, shift and generated samples |

`Thinner::describe()` prints the effective configuration;
`Thinner::stats()` reports samples offered/retained/rejected, saturations,
rebuilds, and the largest field seen.

### Which settings to use

* **Integration of smooth functions** (the beyond-Hardy–Krause regime):
  the defaults (`linear` + `random` shift).
* **Smallest star discrepancy**: `shift = none`. The shift costs up to a
  factor 2^d in the discrepancy bound; it is what buys the σ_SO(f) integration
  guarantee.
* **Faithful Theorem 1.2 algorithm**: `feedback = sign`, `shift = random`.
  In practice its Haar discrepancies scale like N/ε and its star discrepancy
  is only at Monte Carlo level for n ≲ 10⁵; `linear` is 2–3× better.
* **Best empirical discrepancy**: `saturation_bound = 1` — the greedy limit
  (reject with probability ε whenever the field is positive). This is outside
  the paper's analysis but consistently wins in the benchmarks.
* **d ≥ 4**: work per sample is Θ(log^d n) in the uncapped algorithm; the
  automatic rule lowers the resolution until at most
  `max_functions_per_sample` Haar functions are visited per sample. Raise the
  cap (and the memory budget) if you can afford it.
* **Sequences**: leave `level` and `expected_n` unset and feed samples; the
  resolution grows with the retained count (footnote 7 of the paper).

## Results

All numbers from `haar-bench --quick` on one laptop core (GCC 13, `-O3`); the
full tables are in `docs/`. Discrepancies are unnormalised, as in the paper
(divide by n for the usual normalisation). *sign* and *linear* are the two
feedback rules; *shift* is the uniform random shift; *greedy* is
`saturation_bound = 1`.

**Star discrepancy** (mean of 2 runs). Monte Carlo grows like √n; thinning
grows polylogarithmically, so the gap widens with n.

| d | n | Monte Carlo | sign | sign + shift | linear | linear + shift (default) | greedy + shift |
|---|---:|---:|---:|---:|---:|---:|---:|
| 1 | 1 024 | 30.2 | 18.4 | 21.3 | 18.7 | 16.2 | 12.1 |
| 1 | 4 096 | 57.8 | 41.1 | 57.9 | 36.4 | 43.0 | 20.5 |
| 1 | 16 384 | 114.8 | 66.2 | 81.9 | 57.5 | 62.6 | 33.6 |
| 1 | 65 536 | ~222 | ~100 | — | ~63 | — | ~37 |
| 2 | 1 024 | 40.0 | 34.2 | 32.5 | 35.0 | 31.7 | 31.3 |
| 2 | 4 096 | 78.5 | 81.9 | 82.9 | 58.0 | 71.1 | 51.5 |
| 2 | 16 384 | 147.4 | 113.5 | 127.1 | 105.1 | 108.1 | 71.2 |

**L2-star discrepancy, d = 3** (mean of 2 runs):

| n | Monte Carlo | sign | sign + shift | linear | linear + shift (default) | greedy + shift |
|---:|---:|---:|---:|---:|---:|---:|
| 1 024 | 8.23 | 9.57 | 8.79 | 8.47 | 8.51 | 6.57 |
| 4 096 | 21.4 | 17.3 | 15.8 | 13.3 | 13.3 | 12.8 |
| 16 384 | 48.9 | 37.2 | 39.4 | 27.9 | 30.1 | 19.5 |

**RMS integration error, d = 2, n = 65 536** (8 runs each; full tables for
d = 1, 2, 3 and n = 2¹⁰…2¹⁶ in `docs/bench_integ_quick.md`):

| f | Monte Carlo | sign | sign + shift | linear | linear + shift (default) | greedy + shift |
|---|---:|---:|---:|---:|---:|---:|
| Σ sin(2π·16 x_i) | 4.4e-3 | 2.1e-3 | 2.5e-3 | 1.7e-3 | 2.4e-3 | 1.5e-3 |
| Π exp(−16 (x_i − ½)²) | 4.4e-4 | 9.7e-4 | 8.1e-4 | 4.4e-4 | 2.7e-4 | 1.3e-4 |
| Π √x_i | 1.1e-3 | 1.0e-3 | 1.1e-3 | 2.3e-4 | 4.0e-4 | 1.3e-4 |
| 1[Σ x_i < 1] | 2.2e-3 | 2.2e-3 | 2.0e-3 | 5.5e-4 | 7.8e-4 | 5.3e-4 |

Two honest remarks. First, the asymptotic advantage (polylog versus √n) is
real but the constants are not small: at n ≈ 10³–10⁴ the thinned sets are only
modestly better than random ones, and the sign rule can be worse. Second, for
high-frequency integrands (the sine with k = 16 in d = 3) thinning does not
help until n is large enough for the resolution to reach the frequency.

**Throughput**, linear feedback with shift, `expected_n = 200 000`:

| d | scale set | level | Haar functions / sample | table | points / s |
|---|---|---:|---:|---:|---:|
| 1 | hyperbolic | 19 | 19 | 2 MiB | ~8 000 000 |
| 2 | hyperbolic | 19 | 209 | 21 MiB | ~700 000 |
| 3 | hyperbolic | 19 | 1 539 | 126 MiB | ~40 000 |
| 3 | box | 7 | 511 | 8 MiB | ~200 000–340 000 |
| 4 | hyperbolic | 12 (capped from 19) | 1 819 | 1.5 MiB | ~50 000–90 000 |

About 4–6 ns per Haar function per candidate (timings on this laptop vary by
±30% between runs). The cost is inherent to the algorithm (Θ(N) counter reads
per sample); everything else is already integer arithmetic on cache-friendly
layouts.

## Replicability

The per-sample path uses no floating point beyond one exact subtraction and
one exact scaling: coordinates become dyadic fixed-point integers,
discrepancies are `int32` counters, ε is a 32-bit fixed-point number, and the
rejection coin is a 32-bit integer compared against an exact integer
threshold. The generator is xoshiro256++ with a specified seeding. Hence the
same `Options`, seed and sample stream give identical output on every
platform and compiler; `tests/test_replicability.cpp` pins golden digests, and
any change to them is recorded in `CHANGELOG.md`.

## Project layout

```
include/haar/         the library (header-only)
  thinner.hpp           Thinner: the online algorithm
  options.hpp           Options and validation
  discrepancy.hpp       star / L2-star discrepancy, sample means
  rng.hpp               xoshiro256++
  detail/geometry.hpp   scale sets and counter layout
  detail/table.hpp      dense and sparse counter tables
src/haar_thin.cpp     command-line tool
tests/                dependency-free unit tests (ctest)
bench/bench.cpp       benchmarks (haar-bench)
examples/             basic.cpp, stream.cpp
docs/ALGORITHM.md     paper-to-code mapping, constants, complexity
```

## Testing and benchmarking

```sh
ctest --test-dir build --output-on-failure
build/haar-bench --quick                 # a few minutes
build/haar-bench --section bound         # saturation calibration of the automatic B
```

`tests/test_thinner.cpp` contains a deliberately naive reference
implementation (every Haar function evaluated on every retained point at every
step, no tables) and checks that the optimised thinner makes the identical
retain/reject decision at every step for all combinations of feedback rule,
scale set, storage, shift and level policy.

## Citing

Please cite the paper (see `CITATION.cff`):

```bibtex
@article{EJJ26,
  title  = {Quasi-Monte Carlo Beyond Hardy--Krause II: $(1+\varepsilon)n$ Samples Suffice},
  author = {Ezeunala, Ekene and Jha, Agastya Vibhuti and Jiang, Haotian},
  year   = {2026}
}
```

## License

MIT. See `LICENSE`.
