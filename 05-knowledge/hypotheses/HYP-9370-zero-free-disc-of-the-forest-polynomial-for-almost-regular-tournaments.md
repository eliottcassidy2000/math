---
id: HYP-9370
title: "Zero-free disc of the forest polynomial: for every regular or almost-regular tournament T, Phi_T(z) = sum_L N_L(T) (z/n)^L has no zero in |z| <= 11/10, and its zeros of least modulus are a purely imaginary pair +-i t_0(T); for almost-regular T with n >= 8 even, no zero in |z| <= 6/5. min_T t_0 grows with n from n = 8 on (almost-regular 1.222, 1.247, 1.300, 1.349, 1.396 at n = 8..16; regular 1.179, 1.219 at n = 9, 11; all 13 doubly regular 23-tournaments 1.494; minimum 1.118 at C_5). The zeros are deformed Szego zeros of the truncated one-loop exponential, so the z-series of the cluster sum log Phi_T - log rho_T converges at z = 1 with ratio at most 1/1.118."
status: >
  OPEN (opus-2026-10-22-S54). FINITE-EXACT for all 85 classes at n = 8 and all 13,333 classes at n = 10 (complete
  censuses: labelled masses 3,230,080 and 48,251,508,480 = the numbers of labelled regular 9- and 11-tournaments,
  count_regular.py); sampled (MCMC by score-preserving 3-cycle reversals, exact forest numbers) at n = 12 (3000),
  14 (3000), 16 (1500). Every smallest zero lies on the negative w = z^2 axis. The zero-free disc FAILS for unrestricted
  tournaments: uniformly random 12-tournaments have a zero in |z| <= 1 with frequency 0.003 (min 0.945).
session: opus-2026-10-22-S54
source: 01-canon/theorems/THM-4757-linear-forest-loop-expansion-of-hamiltonian-paths.md
scripts:
  - 04-computation/experiments/loop_nonregular_20261022/phi_zeros.py (+ phi_zeros.out): zeros, argument, spectral radius
  - 04-computation/experiments/loop_nonregular_20261022/phi_zeros_model.py (+ phi_zeros_model.out): truncated one-loop and Gaussian models
  - 04-computation/experiments/loop_nonregular_20261022/ar_tool.c: classes, samples, exact forest numbers c_L = (n - L)! N_L
related:
  - 05-knowledge/hypotheses/HYP-9369-almost-regular-loop-predictor-at-scale-complete-n10-census-and-the-joint-remainder.md
  - 05-knowledge/hypotheses/HYP-9366-cosine-twin-and-the-alternating-end-of-the-forest-series.md
  - 05-knowledge/results/session30_nonregular_loops_20261022/SYNTHESIS.md
---

# HYP-9370 — the zero-free disc of the forest polynomial

## Objects

`Φ_T(z) = Σ_L N_L (z/n)^L` is a polynomial of degree at most `n − 1`. Its odd coefficients vanish, so it is a polynomial `P_T(w)` in `w = z²` of degree `⌊(n−1)/2⌋`. By THM-4757,
```
log Φ_T(z) = X0(z) + (sum of all cluster classes),    X0(z) = log ρ_T(z),
```
and `X0` is analytic for `|z| < 1/max|x_j|`, where `±i x_j` are the eigenvalues of `S/n`. So the cluster sum `log Φ_T − X0` has a convergent power series in `z` exactly up to the nearest of:
* the zeros of `Φ_T`;
* the spectral singularities at `|z| = 1/max|x_j|`.

Session30 proved the coefficientwise bound `[z^L](log Φ − X0 − A1 − A2) = O_L(n^(−3))` (nonregular_forest_extension.md §8) and listed evaluation at `z = 1` as the open step. A zero of `Φ_T` inside the closed unit disc would make the `z`-series diverge at `z = 1`. Its absence is the necessary condition checked here.

## Statement

For every regular tournament (odd `n`) and every almost-regular tournament (even `n`, scores `n/2` and `n/2 − 1`):
* **(a)** `Φ_T(z) ≠ 0` for `|z| ≤ 11/10`. For almost-regular `T` with `n ≥ 8`, `Φ_T(z) ≠ 0` for `|z| ≤ 6/5`.
* **(b)** The zeros of least modulus are a purely imaginary pair `±i t_0(T)`, i.e. the smallest root of `P_T` is negative real.
* **(c)** Within each family, `min_T t_0(T)` increases with `n` from `n = 8` on, and tends to infinity.

## Evidence (FINITE-EXACT at n = 8, 10; sampled at n = 12, 14, 16)

| n | tournaments | min `t_0` | median `t_0` | max `t_0` | smallest root negative real | min `1/max|x_j|` |
|---|---|---|---|---|---|---|
| 8 | all 85 classes | 1.2222 | 1.2603 | 1.2950 | 85/85 | 1.591 |
| 10 | all 13,333 classes | 1.2467 | 1.2915 | 1.3431 | 13,333/13,333 | 1.584 |
| 12 | 3000 samples | 1.3000 | 1.3301 | 1.3624 | all | 1.710 |
| 14 | 3000 samples | 1.3491 | 1.3735 | 1.4003 | all | 1.819 |
| 16 | 1500 samples | 1.3955 | 1.4180 | 1.4442 | all | 2.042 |

The binding singularity is always the zero of `Φ_T`. It is not the spectrum: `1/max|x_j|` is larger in every row. The spread of `t_0` within one order is small. The minima increase steadily with `n`.

**Regular tournaments and small orders** (`phi_zeros_small_regular.out`). The regular `(n+1)`-tournaments are obtained as the unique completions of the complete almost-regular lists, once per vertex orbit (`complete_regular.py`).

| family | n | classes | min `t_0` | median | all smallest zeros imaginary |
|---|---|---|---|---|---|
| regular | 3 | 1 (`C_3`) | 1.2247 = √(3/2) | — | yes |
| almost-regular | 4 | 1 | 1.4142 = √2 | — | yes |
| regular | 5 | 1 (`C_5`) | 1.1180 = √5/2 | — | yes |
| almost-regular | 6 | 5 | 1.2335 | 1.2735 | yes |
| regular | 7 | 3 | 1.1311 | 1.1445 | yes |
| regular | 9 | 15 | 1.1790 | 1.1865 | yes |
| regular | 11 | 1223 | 1.2194 | 1.2391 | yes |
| doubly regular | 23 | 13 (incl. `P_23`) | 1.49384 | 1.49390 | yes |

The `n = 23` row uses mac-mini's exact forest numbers (`../sierpinski_ocf_20261022/drt23_forest_exact.json`; `phi_zeros_drt23.out`). All 13 classes have their nearest zero at `w = −2.2316…−2.2318`, so `t_0` agrees across classes to `10^(−4)`.

Two features stand out:
* Regular orders sit below the neighbouring almost-regular ones. Their degree in `w` is the same as for `n + 1`.
* Below `n = 8` the almost-regular minima are not monotone: 1.414, 1.234, 1.222 at `n = 4, 6, 8`.

The global minimum over everything tested is `C_5`, with `t_0 = √5/2`. There `N_2 = 20` and `N_4 = 0` (forest numbers `c = (120, 0, 120, 0, 0)`). So `Φ = 1 + (4/5)w`, with its single root at `w = −5/4`.

The almost-regular 4-tournament has `Φ = 1 + w/2`, root `w = −2`.

**Unrestricted tournaments are different.** For uniform random 12-tournaments (1500 samples), 0.33% have a zero in `|z| ≤ 1` (minimum modulus 0.945). The smallest zeros sit anywhere: 43% on the positive `w`-axis, 38% on the negative `w`-axis, 19% complex. At `n = 14` (800 samples) the minimum is 1.050. Score irregularity of order `√n` destroys the disc. This matches Session30's power counting: the tree level is `O(1)` at that scale.

## Mechanism: deformed Szegő zeros

For almost-regular `T` the low coefficients of `P_T` are positive and close to those of the one-loop exponential. At `n = 8` (first class):
* `P_T(w) = 1 + 0.75w + 0.0703w² + 0.00049w³`, while `ρ_T`'s series is `1 + 0.875w + 0.383w² + 0.127w³ + …`;
* `N_2/n² = 1 − 2/n` for every almost-regular `T`.

The higher coefficients are strongly suppressed: a forest with `L` arcs and `p` paths needs `L + p` of the `n` vertices. Along the real `z`-axis, `P_T(w)` with `w > 0` is a positive sum. Along the imaginary axis, `P_T(−t²)` alternates, and its first sign change is the nearest zero. This is the Szegő picture for the partial sums `Σ_(k≤m) w^k/k!` of the exponential: their zeros lie near the scaled Szegő curve, closest to the origin on the negative axis.

Quantitatively (`phi_zeros_model.out`), the smallest zero of the exact `P_T` sits at 0.95 (`n = 8`) to 0.79 (`n = 16`) of the smallest-root modulus of the same-degree truncation of `ρ_T(z)`. That truncated root is in turn close to the root of the truncated Gaussian `exp((n−1)z²/n)`: `|z|` = 1.351, 1.470, 1.542, 1.643, 1.716 for `n = 8, 10, 12, 14, 16`. The exact zeros recede more slowly than the free model, but they recede.

**Why (c) should hold.** If `Φ_T(z) → exp(z²)` locally uniformly for flat almost-regular families (each weight-`w` class is `O(n^(−w))` on compact sets), Hurwitz's theorem pushes every zero to infinity. Note that the free model's nearest zero grows like `√(0.139 n)`. This argument is circular as stated, since local uniform convergence of the cluster sum is what needs zero-freeness. So (c) is conjectural, with the monotone minima as evidence.

## What this does and does not give

* (a) makes the `z`-series of `log Φ_T − X0` converge geometrically at `z = 1`, with ratio at most `1/1.22`.
  * It does NOT bound the remainder of the weight-graded truncation at `z = 1`.
  * HYP-9369 shows that this remainder, taken separately for the cluster part `log Φ_T(1) − log ρ_T`, is numerically large at small `n`.
  * The useful object is the joint remainder of `log κ`.
* A uniform zero-free disc of radius `r > 1`, plus uniform bounds `|[z^L] (log Φ − X0 − A1 − A2)| ≤ C n^(−3) a^L` with `a < r`, would give `log Φ_T(1) = X0 + A1 + A2 + O(n^(−3))`. The required analytic input is a Lee–Yang / Kotecký–Preiss statement for the resummed core gas.
  * The plain absolute polymer bound fails at `|z| = 1`: it ignores the sign cancellation of the `±1` arc weights.
* **Refutation targets.** An almost-regular tournament with a zero of `Φ_T` in `|z| ≤ 6/5`. Or a non-imaginary nearest zero at some `n ≥ 8`. Or `min_T t_0(T)` decreasing between consecutive even orders.

## Relation to the earlier forest-series results

HYP-9366 concerns the alternating end `F(R_n, −1)` of the forest series of rotational tournaments. mac-mini's four-loop commit 98cbf5a233 (HYP-9371 lane) reports a zero of `Φ` at `z² = −2.232` for doubly regular 23-tournaments. The table above confirms it for all 13 classes: `t_0 = 1.494`, on the same negative-`w` axis and on the trend of the regular rows (1.118, 1.131, 1.179, 1.219 at `n = 5, 7, 9, 11`). Circulant families beyond `n = 11` were not tested.
