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
