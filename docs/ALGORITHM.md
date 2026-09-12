# From the paper to the code

This note maps the algorithms of

> E. Ezeunala, A. V. Jha, H. Jiang. *Quasi-Monte Carlo Beyond Hardy–Krause II:
> (1 + ε) n Samples Suffice.*

onto the implementation, states precisely what the code computes, and records
the places where a practical implementation must choose constants the paper
leaves as "sufficiently large".

## 1. The (1 + ε)-thinning framework (paper §2.1)

At step *t* the strategy holds the retained set *A_t* and a target density
μ_t on [0,1)^d with μ_t(x) ∈ [1 − ε/2, 1 + ε/2]. It draws a sample x_t and a
coin ξ_t ~ Unif[0,1); if ξ_t ≤ μ_t(x_t) − ε/2 it retains x_t, otherwise it
draws y_t and retains that. So each sample is rejected with probability
r_t(x) = 1 + ε/2 − μ_t(x) ∈ [0, ε], and two consecutive samples are never
rejected.

`Thinner::offer(x)` is exactly one presentation of a sample: it returns
`true` if x was retained. After a rejection an internal flag forces the next
call to retain. `Thinner::run(n)` (or `haar::thin`) is the loop "draw x, offer
x" until *n* points are retained. Feeding samples from your own source works
the same way (`examples/stream.cpp`).

**Sample usage.** Because the mean of μ_t is 1, the expected rejection
probability of a *decided* sample is exactly ε/2, and the forced retentions
never consume a coin. A run of *n* points therefore uses about (1 + ε/2) n
samples in expectation and at most (1 + ε) n + o(n) with high probability
(Proposition 2.1). `Stats::offered / Stats::retained` reports the realised ratio.

## 2. Haar discrepancies (paper §2)

For a scale vector **j** ∈ ℕ^d and a point x, exactly one Haar function
H_{**j**,**k**} is nonzero at x. Writing u_i = ⌊x_i 2^R⌋ for the fixed-point
coordinate at the finest level *R*,

* the position in coordinate i is k_i = ⌊x_i 2^{j_i − 1}⌋ = u_i ≫ (R + 1 − j_i)
  for j_i ≥ 1 (and 0 for j_i = 0);
* the sign of the factor h_{j_i,k_i}(x_i) is +1 if bit R − j_i of u_i is 0
  and −1 otherwise.

`Thinner::prepare_point` computes these per coordinate and level,
`Thinner::enumerate_keys` combines them over all scale vectors of the scale
set. Every Haar function owns one `int32` counter holding
φ_t(H) = Σ_{z∈A_t} H(z) exactly. Counters are laid out scale-major
(`detail::Geometry`): all functions with the same **j** form one contiguous
block of 2^{Σ max(j_i − 1, 0)} cells, coarse blocks first. The table is a
plain array when it fits `Options::dense_budget_bytes`, otherwise an
open-addressing hash map (`detail/table.hpp`).

## 3. Scale sets: Π_{≤ℓ} and the hyperbolic cross

The paper's analysis uses Π_{≤ℓ}, all non-constant H_{**j**,**k**} with
max_i j_i ≤ ℓ and ℓ = ⌈10 log₂ n⌉. Only N_ℓ = (ℓ + 1)^d of them are nonzero at
any point, but the set has 2^{ℓ d} members, so the paper stores φ_t in a
dictionary. Both the number of functions visited per sample and the number of
counters ever touched then grow like n log^d n, which is too much beyond
tiny *n* in *d* ≥ 2.

The implementation offers two scale sets (`Options::scale_set`):

| `ScaleSet`   | scale vectors                        | functions per sample | counters                    |
|--------------|--------------------------------------|----------------------|-----------------------------|
| `box`        | max_i j_i ≤ ℓ (the paper's Π_{≤ℓ})    | (ℓ + 1)^d − 1        | 2^{ℓ d}                     |
| `hyperbolic` | Σ_i j_i ≤ L (dyadic boxes of volume ≳ 2^{−L}) | C(L + d, d) − 1 | ≈ 2^L L^{d−1} / (d − 1)!   |

The hyperbolic set is the resolution used by Dwivedi, Feldheim,
Gurel-Gurevich and Ramdas, and is the default: with L = ⌈log₂ n⌉ + 1 the
finest cells contain about half a point, memory is O(n log^{d−1} n), and the
work per sample is O(log^d n / d!).

*Why the analysis is unaffected.* The uncorrelation lemma (Lemma 3.1) never
uses the shape of the scale set: the involution ρ it constructs maps each
Haar function to one with the same scale vector, so it maps any scale set to
itself. The second-moment bound (Lemma 3.2) only uses the number of functions
nonzero at a point, and the high-order bound (Lemma 3.4) only needs every
excluded function to have ‖H‖₂² ≤ 2^{d−L}, which holds for both sets. The
constants change, the theorems do not. (This is our reading of the proofs; the
paper states the results for Π_{≤ℓ}.)

**Automatic level.** With `expected_n = n` set (which `haar::thin` does for
you) the level is fixed from *n*; with neither `level` nor `expected_n` the
level grows with the retained count, L(t) = ⌈log₂(t + 1)⌉ + `level_offset`,
and the table is rebuilt from the stored points at each increase (the point
*sequence* variant of footnote 7; the rebuilds cost O(n log^d n) in total).
Because C(L + d, d) explodes with *d*, the automatic rule additionally lowers
the level until the number of scale vectors is at most
`max_functions_per_sample` (default 2048) and the dense table fits the memory
budget. `Thinner::describe()` shows both the requested and the effective level.
An explicit `Options::level` is honoured as given.

## 4. Sign feedback (paper eq. 3.1, Theorem 1.2)

μ_t(x) = 1 + (ε / 2N) Σ_H sgn(−φ_t(H)) H(x), with N the number of scale
vectors (N_ℓ in the paper). Writing V = Σ_H sgn(φ_t(H)) H(x), the rejection
probability is r = (ε/2)(1 + V/N) ∈ [0, ε].

`Feedback::sign` computes V in one pass over the counters and rejects when
2 N X < E (N + V), where X is a uniform 32-bit integer and E = round(ε 2^32).
This is the same law as the paper's test, evaluated without floating point.

## 5. Uniform shift (paper §3.1)

With `Shift::random` (the default) the constructor draws s ~ Unif[0,1)^d from
the seeded generator, and every sample is processed as x − s mod 1 while the
original x is stored and returned. This is the uniformly-shifted
Haar-thinning algorithm of Theorem 1.2: combined with `Feedback::sign` it gives
the beyond-Hardy–Krause integration error Õ(σ_SO(f)/n); the paper notes that
combined with `Feedback::linear` the hidden polylog improves to log^d n.

The shift slightly *increases* star discrepancy (anchored boxes no longer
align with the shifted dyadic grid; the bound degrades by at most 2^d), so
turn it off when the goal is the smallest possible star discrepancy rather
than integration of smooth functions.

## 6. Linear feedback (paper eq. 4.1, Theorem 1.3)

μ_t(x) = 1 − (ε / 2B) Φ_t(x), Φ_t(x) = Σ_H φ_t(H) H(x). The rejection
probability is r = (ε/2)(1 + Φ/B), valid as long as |Φ_t(x)| ≤ B. The paper
proves that with B = O_{d,ε}(log^{d+1} n) for a sufficiently large constant,
|Φ_t| never exceeds B with high probability (Lemma 4.7), and stops the process
if it does.

The implementation computes Φ exactly and rejects when 2 B X < E (B + Φ).
When |Φ| > B it *clamps* Φ to ±B, counts the event in `Stats::saturated`, and
continues; `OnSaturation::fail` makes it throw instead, reproducing the
paper's stopping rule. Clamping keeps μ_t within [1 − ε/2, 1 + ε/2], so the run
remains a valid (1 + ε)-thinning strategy in every case; only the discrepancy
proof of §4 assumes the clamp is never active.

### Choosing B

B is the one constant the paper leaves open, and it matters: the equilibrium
size of each φ_t(H) is about √(B/ε), so smaller B means smaller discrepancy,
until the field starts saturating. The automatic rule
(`Thinner::automatic_bound`) models the variance of φ_t(H) as
min(B/ε, t‖H‖₂²) — coarse functions have reached the equilibrium of the
restoring drift, fine ones are still random walks — sums it over the
functions nonzero at a point,

V(B) = Σ_H min(B/ε, t ‖H‖₂²),

and picks the fixed point B = 4.5 √V(B). The result is Θ(log^d n / ε) up to
logarithms; the paper's extra log factor comes from a union bound that the
constant 4.5 replaces in practice. Measured on 65 536 points in d = 1, 2, 3
and ε ∈ {0.1, 0.5, 1}, the largest |Φ_t| observed was between 65% and 90% of
B and no saturation occurred (`haar-bench --section bound`).

Setting `saturation_bound` by hand trades discrepancy against saturations.
The extreme `saturation_bound = 1` rejects with probability ε whenever
Φ_t(x) > 0 and never when Φ_t(x) < 0: a greedy rule in the spirit of
Conjecture 2 of DFG+19, which is outside the paper's analysis but gives the
smallest discrepancies in our benchmarks.

## 7. Exactness and replicability

Everything that influences a decision is an integer: fixed-point
coordinates, counters, N, B, E and the 32-bit coin. The only floating-point
operations are x − s (one IEEE subtraction) and y · 2^R (exact scaling). No
result depends on FMA contraction, `long double`, or a library `sin`.
Consequently `Options` + seed + sample stream determine the output bit for
bit on every platform; `tests/test_replicability.cpp` pins three digests.

ε itself is quantised to E/2^32; the deviation from the nominal rejection
probability is below 2.4 × 10⁻¹⁰.

## 8. Complexity

Per sample: one pass computing (ℓ + 1) d fixed-point quantities, then N − 1
table reads (and N − 1 writes if retained), where N is the number of scale
vectors. With the defaults N = C(L + d, d) with L ≈ log₂ n, capped at 2048.
Measured on one core of a laptop (GCC 13, `-O3`): about 4–6 ns per Haar
function when the table is cache-resident, i.e. roughly

| d | N (n = 2^20) | points / s |
|---|--------------|------------|
| 1 | 22           | ~8 000 000 |
| 2 | 253          | ~700 000   |
| 3 | 2024 (capped)| ~40 000    |
| 4 | 1820 (capped at L = 12) | ~50 000 |

Memory: 4 bytes per counter, ≈ 2^L L^{d−1}/(d−1)! counters for the
hyperbolic set (d = 3, n = 2^20: about 265 MB; the default budget is 512 MB),
plus 8 d bytes per retained point.
