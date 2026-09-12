# haar-thinning benchmarks (quick)

### Throughput, n = 200000 points (linear feedback, shift on)

| d | scale set | level | Haar functions / sample | storage | table MiB | time (s) | points / s |
|---|---|---:|---:|---|---:|---:|---:|
| 1 | hyperbolic | 19 | 19 | dense | 2.0 | 0.03 | 7819831 |
| 1 | box | 19 | 19 | dense | 2.0 | 0.02 | 8192994 |
| 2 | hyperbolic | 19 | 209 | dense | 21.0 | 0.30 | 669247 |
| 2 | box | 10 | 120 | dense | 4.0 | 0.13 | 1593034 |
| 3 | hyperbolic | 19 | 1539 | dense | 125.5 | 4.42 | 45249 |
| 3 | box | 7 | 511 | dense | 8.0 | 0.59 | 340762 |
| 4 | hyperbolic | 12 | 1819 | dense | 1.5 | 2.31 | 86658 |
| 4 | box | 5 | 1295 | dense | 4.0 | 1.64 | 122006 |
