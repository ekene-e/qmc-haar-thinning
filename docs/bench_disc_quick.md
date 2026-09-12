# haar-thinning benchmarks (quick)

### Star discrepancy, d = 1 (unnormalised, mean of 2 runs)

| n | Monte Carlo | sign (DFG+19) | sign + shift (Thm 1.2) | linear (Thm 1.3) | linear + shift (default) | greedy B=1 + shift |
|---|---:|---:|---:|---:|---:|---:|
| 1024 | 30.17 | 18.42 | 21.31 | 18.66 | 16.16 | 12.05 |
| 4096 | 57.76 | 41.14 | 57.88 | 36.39 | 43.01 | 20.45 |
| 16384 | 114.81 | 66.22 | 81.85 | 57.52 | 62.55 | 33.61 |

### Star discrepancy, d = 2 (unnormalised, mean of 2 runs)

| n | Monte Carlo | sign (DFG+19) | sign + shift (Thm 1.2) | linear (Thm 1.3) | linear + shift (default) | greedy B=1 + shift |
|---|---:|---:|---:|---:|---:|---:|
| 1024 | 39.96 | 34.20 | 32.53 | 34.97 | 31.69 | 31.29 |
| 4096 | 78.51 | 81.90 | 82.87 | 57.96 | 71.14 | 51.52 |
| 16384 | 147.37 | 113.48 | 127.12 | 105.08 | 108.06 | 71.15 |

### L2-star discrepancy, d = 3 (unnormalised, mean of 2 runs)

| n | Monte Carlo | sign (DFG+19) | sign + shift (Thm 1.2) | linear (Thm 1.3) | linear + shift (default) | greedy B=1 + shift |
|---|---:|---:|---:|---:|---:|---:|
| 1024 | 8.234 | 9.569 | 8.792 | 8.465 | 8.514 | 6.565 |
| 4096 | 21.395 | 17.279 | 15.762 | 13.308 | 13.342 | 12.817 |
| 16384 | 48.924 | 37.245 | 39.432 | 27.920 | 30.111 | 19.471 |
