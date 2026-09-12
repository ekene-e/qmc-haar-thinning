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

MIT. See `LICENSE`. This repository was prepared with the aid of Claude Fable 5.1.
