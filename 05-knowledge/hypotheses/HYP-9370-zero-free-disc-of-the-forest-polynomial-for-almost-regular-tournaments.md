---
id: HYP-9370
title: "Zero-free disc of the forest polynomial: for every regular or almost-regular tournament T, Phi_T(z) = sum_L N_L(T) (z/n)^L has no zero in |z| <= 11/10, and its zeros of least modulus are a purely imaginary pair +-i t_0(T); for almost-regular T with n >= 8 even, no zero in |z| <= 6/5. min_T t_0 grows with n from n = 8 on: almost-regular 1.222, 1.247 at n = 8, 10 (complete censuses) and at most 1.279, 1.310, 1.336, 1.358 at n = 12, 14, 16, 18, attained by the deleted rotational tournament R_(n+1) - v, which is the exact minimiser at n = 6, 8, 10; regular 1.179, 1.219 at n = 9, 11; all 13 doubly regular 23-tournaments 1.494; global minimum C_5 with sqrt(5)/2. Their modulus is comparable (0.79-0.95 in w) to that of the truncated one-loop exponential. Given (a), the z-series of the cluster sum log Phi_T - log rho_T converges at z = 1 with ratio at most 10/11 (5/6 for almost-regular n >= 8)."
status: >
  OPEN (opus-2026-10-22-S54; independently audited 2026-10-22,
  04-computation/experiments/loop_nonregular_20261022/audit_S54/REPORT.md). The audit confirmed every data point with its
  own forest-number code. It found that the n = 12-16 minima first reported were MCMC sample minima, above the deleted
  rotational values; it struck an unsound mechanism paragraph; and it added circulant families n = 13-17 and R_19, with no
  counterexample. FINITE-EXACT for all classes at n = 4, 6, 8, 10 (complete censuses: labelled masses equal the numbers
  of labelled regular 5-, 7-, 9-, 11-tournaments, count_regular.py), for all regular classes at n = 3, 5, 7, 9, 11, and
  for all 13 doubly regular 23-tournaments. Sampled (MCMC by score-preserving 3-cycle reversals, exact forest numbers) at
  n = 12 (3000), 14 (3000), 16 (1500). Every smallest zero lies on the negative w = z^2 axis. The disc FAILS for
  unrestricted tournaments: uniformly random 12-tournaments have a zero in |z| <= 1 with frequency about 0.1-0.3% (5/1500
  here, 21/20,000 in the audit, min 0.854; 1/3000 at n = 14), and all such zeros are real.
session: opus-2026-10-22-S54
source: 01-canon/theorems/THM-4757-linear-forest-loop-expansion-of-hamiltonian-paths.md
scripts:
  - 04-computation/experiments/loop_nonregular_20261022/phi_zeros.py (+ phi_zeros.out, phi_zeros_small_regular.out): zeros, argument, spectral radius
  - 04-computation/experiments/loop_nonregular_20261022/phi_zeros_model.py (+ phi_zeros_model.out): truncated one-loop and Gaussian models
  - 04-computation/experiments/loop_nonregular_20261022/phi_zeros_drt23.py (+ .out): doubly regular 23-tournaments
  - 04-computation/experiments/loop_nonregular_20261022/ar_tool.c: classes, samples, exact forest numbers c_L = (n - L)! N_L
  - 04-computation/experiments/loop_nonregular_20261022/audit_S54/ (zeros_small, zb_*, zeros_structured, rotational_check, szego_check): independent audit
related:
  - 05-knowledge/hypotheses/HYP-9369-almost-regular-loop-predictor-at-scale-complete-n10-census-and-the-joint-remainder.md
  - 05-knowledge/hypotheses/HYP-9366-cosine-twin-and-the-alternating-end-of-the-forest-series.md
  - 05-knowledge/results/session30_nonregular_loops_20261022/SYNTHESIS.md
---

# HYP-9370 — the zero-free disc of the forest polynomial

## Objects

`Φ_T(z) = Σ_L N_L (z/n)^L` is a polynomial of degree at most `n − 1`. Its odd coefficients vanish, so it is a polynomial `P_T(w)` in `w = z²` of degree `⌊(n−1)/2⌋`. By THM-4757,
```
log Φ_T(z) = X0(z) + (sum of all cluster classes),    X0(z) = log ρ_T(z).
```
`X0` is analytic for `|z| < 1/max|x_j|`, where `±i x_j` are the eigenvalues of `S/n`.
* Session30's eigenvalue bound `x_j² ≤ (n−1)/(2n)` gives `1/max|x_j| ≥ √(2n/(n−1)) > √2` (PROVED). So the spectrum can never bind below `√2`.
* The cluster sum `log Φ_T − X0` therefore has a convergent power series in `z` exactly up to the nearest zero of `Φ_T`, or up to `1/max|x_j|` if that is nearer.

Session30 proved the coefficientwise bound `[z^L](log Φ − X0 − A1 − A2) = O_L(n^(−3))` (nonregular_forest_extension.md §8). They listed evaluation at `z = 1` as the open step. A zero of `Φ_T` inside the open unit disc would make the `z`-series diverge at `z = 1`. Its absence is the necessary condition checked here.

## Statement

For every regular tournament (odd `n`) and every almost-regular tournament (even `n`, scores `n/2` and `n/2 − 1`):
* **(a)** `Φ_T(z) ≠ 0` for `|z| ≤ 11/10`. For almost-regular `T` with `n ≥ 8`, `Φ_T(z) ≠ 0` for `|z| ≤ 6/5`.
* **(b)** The zeros of least modulus are a purely imaginary pair `±i t_0(T)`; that is, the smallest root of `P_T` is negative real.
* **(c)** Within each family, `min_T t_0(T)` increases with `n` from `n = 8` on, and tends to infinity.

## Evidence

**Almost-regular tournaments.** FINITE-EXACT at `n = 8, 10`; sample statistics at `n = 12, 14, 16`. The column "`R_(n+1) − v`" is the deleted rotational tournament (audit, `rotational_check.out`). It is the exact minimiser of the complete censuses at `n = 6, 8, 10`. At `n = 12, 14, 16` it lies below every sample, and the sampling chain starts there.

| n | tournaments | min `t_0` | median | max | `R_(n+1) − v` | smallest root negative real | min `1/max|x_j|` |
|---|---|---|---|---|---|---|---|
| 8 | all 85 classes | 1.2222 | 1.2603 | 1.2950 | 1.2222 | 85/85 | 1.591 |
| 10 | all 13,333 classes | 1.2467 | 1.2915 | 1.3431 | 1.2467 | 13,333/13,333 | 1.584 |
| 12 | 3000 samples | 1.3000 (sample) | 1.3301 | 1.3624 | 1.2789 | all | 1.710 |
| 14 | 3000 samples | 1.3491 (sample) | 1.3735 | 1.4003 | 1.3095 | all | 1.819 |
| 16 | 1500 samples | 1.3955 (sample) | 1.4180 | 1.4442 | 1.3357 | all | 2.042 |
| 18 | `R_19 − v` only | — | — | — | 1.3576 | yes | — |

The binding singularity is always the zero of `Φ_T`; `1/max|x_j|` is larger in every row. The best available minima increase: 1.2222 < 1.2467 < 1.2789 < 1.3095 < 1.3357 < 1.3576 (`n = 8, …, 18`).

**Regular tournaments, small orders and structured families.** The regular `(n+1)`-tournaments are the unique completions of the complete almost-regular lists, once per vertex orbit (`complete_regular.py`). Medians are taken over completions; the per-class medians are 1.1414 at `n = 7` and 1.1876 at `n = 9` (audit).

| family | n | classes | min `t_0` | median (completions) | all smallest zeros imaginary |
|---|---|---|---|---|---|
| regular | 3 | 1 (`C_3`) | 1.2247 = √(3/2) | — | yes |
| almost-regular | 4 | 1 | 1.4142 = √2 | — | yes |
| regular | 5 | 1 (`C_5`) | 1.1180 = √5/2 | — | yes |
| almost-regular | 6 | 5 | 1.2335 | 1.2735 | yes |
| regular | 7 | 3 | 1.1311 (`P_7`) | 1.1445 | yes |
| regular | 9 | 15 | 1.1790 | 1.1865 | yes |
| regular | 11 | 1223 | 1.2194 | 1.2391 | yes |
| circulant regular (audit) | 13, 15, 17 | all | 1.2763, 1.3100, 1.3370 | — | yes |
| their vertex deletions (audit) | 12, 14, 16 | all | 1.2789, 1.3095, 1.3357 | — | yes |
| rotational `R_19` (audit) | 19 | 1 | 1.3589 | — | yes |
| doubly regular | 23 | 13 (incl. `P_23`) | 1.49384 | 1.49390 | yes |

The `n = 23` row uses mac-mini's exact forest numbers (`../sierpinski_ocf_20261022/drt23_forest_exact.json`; `phi_zeros_drt23.out`; recomputed independently by the audit). All 13 classes have their nearest zero at `w = −2.2316…−2.2318`. The doubly regular tournaments are not family minimisers.

The audit also tested further structured families with no counterexample: `R[K_2]` for all regular 5- and 7-tournaments, `C_3[C_5]`, `C_5[C_3]` and their deletions, `C_3[C_3[K_2]]` and `C_9[K_2]`. All have `t_0 ≥ 1.31`, with negative-real nearest zeros (`zeros_structured.out`).

Further observations:
* Regular orders sit below the neighbouring almost-regular ones only for `n ≤ 11`. From `q = 15` on, `R_q − v` lies below `R_q`.
* Below `n = 8` the almost-regular minima are not monotone: 1.414, 1.234, 1.222 at `n = 4, 6, 8`.
* The global minimum over everything tested is `C_5`, with `t_0 = √5/2`. There `N_2 = 20` and `N_4 = 0` (forest numbers `c = (120, 0, 120, 0, 0)`), so `Φ = 1 + (4/5)w`, with single root `w = −5/4`.
* The almost-regular 4-tournament has `Φ = 1 + w/2`, root `w = −2`.

**Unrestricted tournaments are different.** For uniform random 12-tournaments, about 0.1–0.3% have a zero in `|z| ≤ 1`:
* 5/1500 in this lane's sample, minimum 0.945;
* 21/20,000 in the audit's sample, minimum 0.854;
* 1/3000 at `n = 14`, minimum 0.903.

Every such zero is real (`w > 0`). It comes from a very negative `N_2`, i.e. `W > n(n−1)`, which gives a positive root of `P_T`. The nearest zeros overall split as 42% positive-real, 39% negative-real and 19% complex. Score irregularity of order `√n` destroys the disc. This matches Session30's power counting: the tree level is `O(1)` at that scale (THM-4761 §2).

## What is and is not understood about (b)

For almost-regular `T` the lowest coefficients of `P_T` are close to those of the one-loop exponential:
* `N_2/n² = 1 − 2/n` for every almost-regular `T`;
* at `n = 8` (first class), `P_T(w) = 1 + 0.75w + 0.0703w² + 0.00049w³`, while `ρ_T`'s series is `1 + 0.875w + 0.383w² + 0.127w³ + …`.

The higher coefficients are strongly suppressed by vertex exclusion: a forest with `L` arcs and `p` paths needs `L + p` of the `n` vertices.

The moduli are comparable to those of truncated exponentials (Szegő 1924, partial sums of `e^w`). In `w`, the smallest zero of the exact `P_T` is 0.95 (`n = 8`) to 0.79 (`n = 16`) of the smallest-root modulus of the same-degree truncation of `ρ_T(z)` (`phi_zeros_model.out`). That truncated root is close to the root of the truncated Gaussian `exp((n−1)z²/n)`: `|z|` = 1.351, 1.470, 1.542, 1.643, 1.716 for `n = 8, 10, 12, 14, 16`.

Statement (b) is NOT explained by this picture, for two reasons (audit, `szego_check.out`):
* At even `w`-degree (`n = 10, 14`) the truncated models have complex nearest zeros (`arg w ≈ 2.66–2.80`), while every exact `P_T` has a negative real one.
* `P_T` is not a positive polynomial. Negative coefficients and positive real roots occur in 22 of 85 almost-regular 8-classes and 5379 of 13,333 almost-regular 10-classes, though never as the nearest zero.

Why the nearest zero is always on the negative `w`-axis is open.

**Why (c) should hold.** Suppose `Φ_T(z) → exp(z²)` locally uniformly for flat almost-regular families, with each weight-`w` class `O(n^(−w))` on compact sets. Then Hurwitz's theorem pushes every zero to infinity; the free model's nearest zero grows like `√(0.139 n)`. This argument is circular as stated, since local uniform convergence of the cluster sum is what needs zero-freeness. So (c) is conjectural, with the increasing minima as evidence.

## What this does and does not give

* (a) makes the `z`-series of `log Φ_T − X0` converge geometrically at `z = 1`.
  * The guaranteed ratio is `10/11` in general and `5/6` for almost-regular `n ≥ 8`. The observed nearest zeros give smaller ratios (`1/1.22` and below).
  * It does NOT bound the remainder of the weight-graded truncation at `z = 1`.
  * HYP-9369 shows that this remainder, taken separately for the cluster part `log Φ_T(1) − log ρ_T`, is about `250/n³` at `n ≤ 16`. The useful object is the joint remainder of `log κ`.
* Suppose there were a uniform zero-free disc of radius `r > 1`, plus uniform bounds `|[z^L] (log Φ − X0 − A1 − A2)| ≤ C n^(−3) a^(−L)` with `1 < a < r`. Then `log Φ_T(1) = X0 + A1 + A2 + O(n^(−3))`.
  * The required analytic input is a Lee–Yang / Kotecký–Preiss statement for the resummed core gas.
  * The plain absolute polymer bound fails at `|z| = 1`, because it ignores the sign cancellation of the `±1` arc weights.
* **Refutation targets.**
  * A regular or almost-regular tournament with a zero of `Φ_T` in `|z| ≤ 11/10`, or an almost-regular one with `n ≥ 8` and a zero in `|z| ≤ 6/5`.
  * A nearest zero that is not purely imaginary.
  * `min_T t_0(T)` decreasing between consecutive orders of one family from `n = 8` on. The deleted rotational tournaments are the natural first test.

## Relation to the earlier forest-series results

HYP-9366 concerns the alternating end `F(R_n, −1)` of the forest series of rotational tournaments. mac-mini's four-loop commit 98cbf5a233 (HYP-9371 lane) reports a zero of `Φ` at `z² = −2.232` for doubly regular 23-tournaments. The table above confirms it for all 13 classes: `t_0 = 1.494`, on the same negative-`w` axis.

No prior work on zeros of linear-forest or path-cover polynomials of tournaments was found (audit). The nearest objects are the Rédei–Berge function (Grinberg–Stanley, arXiv:2307.05569, 2402.07606) and the cycle-path indicator polynomial of D'Antona–Munarini (Adv. Appl. Math. 2000).
