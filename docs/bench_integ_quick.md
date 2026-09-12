# haar-thinning benchmarks (quick)

### RMS integration error, d = 1, f = sum_i sin(2 pi 16 x_i) (8 runs)

| n | Monte Carlo | sign (DFG+19) | sign + shift (Thm 1.2) | linear (Thm 1.3) | linear + shift (default) | greedy B=1 + shift |
|---|---:|---:|---:|---:|---:|---:|
| 1024 | 1.22e-02 | 9.46e-03 | 1.12e-02 | 2.20e-02 | 1.34e-02 | 9.54e-03 |
| 8192 | 7.59e-03 | 5.72e-03 | 4.89e-03 | 7.05e-03 | 6.70e-03 | 4.64e-03 |
| 65536 | 3.27e-03 | 1.43e-03 | 1.21e-03 | 1.18e-03 | 9.85e-04 | 7.68e-04 |

### RMS integration error, d = 1, f = prod_i exp(-16 (x_i - 1/2)^2) (8 runs)

| n | Monte Carlo | sign (DFG+19) | sign + shift (Thm 1.2) | linear (Thm 1.3) | linear + shift (default) | greedy B=1 + shift |
|---|---:|---:|---:|---:|---:|---:|
| 1024 | 1.09e-02 | 1.02e-02 | 9.68e-03 | 7.91e-03 | 8.03e-03 | 2.38e-03 |
| 8192 | 3.48e-03 | 2.20e-03 | 1.88e-03 | 1.07e-03 | 1.26e-03 | 3.18e-04 |
| 65536 | 1.13e-03 | 4.48e-04 | 2.30e-04 | 2.05e-04 | 1.58e-04 | 6.92e-05 |

### RMS integration error, d = 1, f = prod_i sqrt(x_i) (8 runs)

| n | Monte Carlo | sign (DFG+19) | sign + shift (Thm 1.2) | linear (Thm 1.3) | linear + shift (default) | greedy B=1 + shift |
|---|---:|---:|---:|---:|---:|---:|
| 1024 | 7.75e-03 | 6.46e-03 | 6.95e-03 | 4.79e-03 | 5.66e-03 | 3.79e-03 |
| 8192 | 2.28e-03 | 1.31e-03 | 1.31e-03 | 6.48e-04 | 1.40e-03 | 5.12e-04 |
| 65536 | 6.92e-04 | 2.36e-04 | 2.44e-04 | 9.89e-05 | 1.54e-04 | 1.26e-04 |

### RMS integration error, d = 1, f = 1[ sum_i x_i < d/2 ] (8 runs)

| n | Monte Carlo | sign (DFG+19) | sign + shift (Thm 1.2) | linear (Thm 1.3) | linear + shift (default) | greedy B=1 + shift |
|---|---:|---:|---:|---:|---:|---:|
| 1024 | 1.71e-02 | 1.32e-02 | 1.38e-02 | 7.68e-03 | 1.33e-02 | 6.84e-03 |
| 8192 | 5.45e-03 | 1.81e-03 | 1.73e-03 | 1.03e-03 | 1.20e-03 | 1.27e-03 |
| 65536 | 2.34e-03 | 3.38e-04 | 4.54e-04 | 1.61e-04 | 2.99e-04 | 2.02e-04 |

### RMS integration error, d = 2, f = sum_i sin(2 pi 16 x_i) (8 runs)

| n | Monte Carlo | sign (DFG+19) | sign + shift (Thm 1.2) | linear (Thm 1.3) | linear + shift (default) | greedy B=1 + shift |
|---|---:|---:|---:|---:|---:|---:|
| 1024 | 2.34e-02 | 2.52e-02 | 2.03e-02 | 2.73e-02 | 2.72e-02 | 2.58e-02 |
| 8192 | 1.08e-02 | 9.08e-03 | 9.12e-03 | 1.04e-02 | 9.91e-03 | 9.07e-03 |
| 65536 | 4.37e-03 | 2.06e-03 | 2.50e-03 | 1.71e-03 | 2.39e-03 | 1.50e-03 |

### RMS integration error, d = 2, f = prod_i exp(-16 (x_i - 1/2)^2) (8 runs)

| n | Monte Carlo | sign (DFG+19) | sign + shift (Thm 1.2) | linear (Thm 1.3) | linear + shift (default) | greedy B=1 + shift |
|---|---:|---:|---:|---:|---:|---:|
| 1024 | 1.25e-02 | 4.25e-03 | 4.08e-03 | 4.74e-03 | 4.39e-03 | 1.89e-03 |
| 8192 | 3.24e-03 | 2.53e-03 | 2.40e-03 | 2.00e-03 | 1.55e-03 | 6.57e-04 |
| 65536 | 4.39e-04 | 9.67e-04 | 8.13e-04 | 4.39e-04 | 2.69e-04 | 1.26e-04 |

### RMS integration error, d = 2, f = prod_i sqrt(x_i) (8 runs)

| n | Monte Carlo | sign (DFG+19) | sign + shift (Thm 1.2) | linear (Thm 1.3) | linear + shift (default) | greedy B=1 + shift |
|---|---:|---:|---:|---:|---:|---:|
| 1024 | 8.07e-03 | 6.43e-03 | 7.21e-03 | 5.37e-03 | 6.65e-03 | 3.91e-03 |
| 8192 | 2.99e-03 | 2.71e-03 | 3.43e-03 | 1.40e-03 | 1.82e-03 | 1.01e-03 |
| 65536 | 1.06e-03 | 1.01e-03 | 1.07e-03 | 2.28e-04 | 4.01e-04 | 1.34e-04 |

### RMS integration error, d = 2, f = 1[ sum_i x_i < d/2 ] (8 runs)

| n | Monte Carlo | sign (DFG+19) | sign + shift (Thm 1.2) | linear (Thm 1.3) | linear + shift (default) | greedy B=1 + shift |
|---|---:|---:|---:|---:|---:|---:|
| 1024 | 1.74e-02 | 7.60e-03 | 6.96e-03 | 7.08e-03 | 1.04e-02 | 1.19e-02 |
| 8192 | 6.60e-03 | 4.30e-03 | 4.80e-03 | 2.37e-03 | 2.25e-03 | 2.34e-03 |
| 65536 | 2.18e-03 | 2.20e-03 | 1.96e-03 | 5.54e-04 | 7.82e-04 | 5.28e-04 |

### RMS integration error, d = 3, f = sum_i sin(2 pi 16 x_i) (8 runs)

| n | Monte Carlo | sign (DFG+19) | sign + shift (Thm 1.2) | linear (Thm 1.3) | linear + shift (default) | greedy B=1 + shift |
|---|---:|---:|---:|---:|---:|---:|
| 1024 | 2.55e-02 | 4.05e-02 | 3.83e-02 | 4.44e-02 | 3.45e-02 | 3.46e-02 |
| 8192 | 9.96e-03 | 1.43e-02 | 1.32e-02 | 1.20e-02 | 1.13e-02 | 1.10e-02 |
| 65536 | 4.41e-03 | 3.04e-03 | 2.73e-03 | 1.85e-03 | 4.11e-03 | 2.94e-03 |

### RMS integration error, d = 3, f = prod_i exp(-16 (x_i - 1/2)^2) (8 runs)

| n | Monte Carlo | sign (DFG+19) | sign + shift (Thm 1.2) | linear (Thm 1.3) | linear + shift (default) | greedy B=1 + shift |
|---|---:|---:|---:|---:|---:|---:|
| 1024 | 5.00e-03 | 4.43e-03 | 4.52e-03 | 4.20e-03 | 4.62e-03 | 2.53e-03 |
| 8192 | 1.61e-03 | 1.27e-03 | 1.29e-03 | 1.61e-03 | 1.03e-03 | 7.05e-04 |
| 65536 | 5.89e-04 | 4.67e-04 | 4.65e-04 | 3.29e-04 | 1.93e-04 | 1.55e-04 |

### RMS integration error, d = 3, f = prod_i sqrt(x_i) (8 runs)

| n | Monte Carlo | sign (DFG+19) | sign + shift (Thm 1.2) | linear (Thm 1.3) | linear + shift (default) | greedy B=1 + shift |
|---|---:|---:|---:|---:|---:|---:|
| 1024 | 6.73e-03 | 5.21e-03 | 5.21e-03 | 5.00e-03 | 4.78e-03 | 3.73e-03 |
| 8192 | 1.89e-03 | 2.84e-03 | 2.91e-03 | 1.96e-03 | 2.14e-03 | 1.34e-03 |
| 65536 | 9.51e-04 | 6.88e-04 | 7.19e-04 | 3.20e-04 | 5.55e-04 | 3.93e-04 |

### RMS integration error, d = 3, f = 1[ sum_i x_i < d/2 ] (8 runs)

| n | Monte Carlo | sign (DFG+19) | sign + shift (Thm 1.2) | linear (Thm 1.3) | linear + shift (default) | greedy B=1 + shift |
|---|---:|---:|---:|---:|---:|---:|
| 1024 | 1.54e-02 | 1.93e-02 | 2.01e-02 | 1.55e-02 | 1.64e-02 | 5.36e-03 |
| 8192 | 5.75e-03 | 8.25e-03 | 8.15e-03 | 6.58e-03 | 7.24e-03 | 4.91e-03 |
| 65536 | 2.69e-03 | 2.23e-03 | 2.04e-03 | 9.63e-04 | 1.26e-03 | 5.83e-04 |
