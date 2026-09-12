# haar-thinning benchmarks (quick)

### Linear feedback: automatic bound B versus the field |Phi_t(x)|, n = 65536

Saturation events are decisions where |Phi_t(x)| > B and the density was clamped.

| d | eps | N | level | B (auto) | max abs Phi | saturations | max abs phi(H) | rejection rate |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 0.1 | 18 | 17 | 955 | 652 | 0 | 200 | 0.0472 |
| 1 | 0.5 | 18 | 17 | 344 | 264 | 0 | 81 | 0.1993 |
| 1 | 1.0 | 18 | 17 | 207 | 184 | 0 | 41 | 0.3319 |
| 2 | 0.1 | 171 | 17 | 2826 | 1832 | 0 | 282 | 0.0459 |
| 2 | 0.5 | 171 | 17 | 1508 | 1299 | 0 | 174 | 0.1990 |
| 2 | 1.0 | 171 | 17 | 1059 | 868 | 0 | 110 | 0.3336 |
| 3 | 0.1 | 1140 | 17 | 5779 | 4425 | 0 | 458 | 0.0478 |
| 3 | 0.5 | 1140 | 17 | 4189 | 3279 | 0 | 233 | 0.2014 |
| 3 | 1.0 | 1140 | 17 | 3266 | 2786 | 0 | 177 | 0.3319 |
