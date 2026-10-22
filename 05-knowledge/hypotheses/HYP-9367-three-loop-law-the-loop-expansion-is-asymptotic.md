---
id: HYP-9367
title: "Three-loop law: the linear-forest loop expansion of Hamiltonian paths is an asymptotic expansion in 1/n, one power per loop. For regular T, kappa(T)/rho(T) = 1 + Y_1 + Y_2 + O(n^-3), with Y_1 = X_1 + E_0[a]/n (the two-loop law, HYP-9365) and Y_2 = X_2 + X_1^2/2 + X_1 E_0[a]/n + delta/n + (E_0[a^2]/2 + E_0[b])/n^2, where X_2 is the sum over the 15 three-loop cores (loop_engine.py), a = C(L,2), b = L(L-1)(2L-1)/12, E_0 the one-loop measure and delta the two-loop shift of E[a]. The second order gives the universal doubly regular constant exactly: kappa = e(1 - 3/(2n) + 25/(8n^2) + O(n^-3)), from q^2 X_2 -> 29."
status: >
  OPEN, pre-audit (mac-mini-2026-10-22-sierp). The three-loop cluster sum is computed by a general engine (any excess) that
  reproduces the two-loop closed forms to 12 digits and the brute-force excess-2 cluster sums (all set partitions of
  run positions, R_15 and P_19 at L = 4, 6) exactly. The asymptotic statement is the open part. Evidence (FINITE-EXACT):
  after the second order the residual scales like n^-3 on Paley (q = 19..31) and rotational (n = 11..35) tournaments,
  and its mean and spread over all circulant classes at n = 21, 25 drop 4-6 fold.
session: mac-mini-2026-10-22-sierp
source: 01-canon/theorems/THM-4757-linear-forest-loop-expansion-of-hamiltonian-paths.md
scripts:
  - 04-computation/experiments/sierpinski_ocf_20261022/loop_engine.py: cores of any excess, decompositions, propagators
  - 04-computation/experiments/sierpinski_ocf_20261022/loop_engine_fast.py: tensorized 5-matrix basis, K_4 and doubled-square contractions
  - 04-computation/experiments/sierpinski_ocf_20261022/excess2_bruteforce_check.py (+ .out): engine against brute force
  - 04-computation/experiments/sierpinski_ocf_20261022/second_order.py (+ .out), second_order_circulants.py (+ .out)
related:
  - 05-knowledge/hypotheses/HYP-9365-two-loop-law-for-hamiltonian-paths-in-regular-tournaments.md
  - 05-knowledge/hypotheses/HYP-9366-cosine-twin-and-the-alternating-end-of-the-forest-series.md
---

# HYP-9367 — the three-loop law

## The engine

Cores of excess `e` are connected multigraphs, loops allowed, of minimum degree ≥ 3 with `E − V = e`. There are 3 cores at excess 1 and 15 at excess 2. The 15 include:
* the bouquet of three loops (`|G| = 48`);
* the 4-fold theta (`|G| = 48`);
* the doubled square (`|G| = 16`);
* the tetrahedron `K_4` (`|G| = 24`);
* eleven cores with loops or bridges.

Run decompositions are summed in three layers:
* half-edge pairings at core vertices, with Möbius weights;
* path types: one segment, or two or more;
* orientations of the runs through core vertices. Free runs are pre-summed, and closed runs are excluded.

Every path propagator is a combination of five matrices: `R_± = ±Ŝ(I ∓ Ŝ)^(−1)` and `Q_(s_1 s_2) = −R_(s_1) R_(s_2) (I − Ŝ²)(I + Ŝ²)^(−1)`. A core's sum is one tensor contraction of a `5^E` weight tensor with these matrices, divided by `|G|`.

**Validation.**

| check | engine | independent |
|---|---|---|
| figure-eight, theta, dumbbell on a non-vertex-transitive regular 11-tournament | 12 digits | closed forms (THM-4757) |
| excess 2, `R_15`, coefficient of `z⁴`, `z⁶` | +1.2444e-2, +1.0050e-1 | brute force +1.244e-2, +1.005e-1 |
| excess 2, `P_19`, coefficient of `z⁴`, `z⁶` | +7.8729e-3, +6.0330e-2 | brute force +7.873e-3, +6.033e-2 |

## The second order

`κ = Σ_L φ_L n^L/(n)_L`, with `log(n^L/(n)_L) = a/n + b/n² + O(L⁴/n³)`, `a = C(L,2)`, `b = L(L−1)(2L−1)/12`, and `log Φ = X_0 + X_1 + X_2 + …`. Expanding to second order:

    κ/ρ = 1 + [X_1 + E_0 a/n] + [X_2 + X_1²/2 + X_1 E_0 a/n + δ/n + (E_0 a²/2 + E_0 b)/n²] + O(n^(−3)).

* `E_0` uses the cumulants `D_k = (z∂_z)^k X_0` at `z = 1`.
* `δ = (D_2 X_1 + 2 D_1 X_0 · D_1 X_1 − D_1 X_1)/2`.

## Evidence

| tournament | residual after 1st order | after 2nd order | `n³ ×` 2nd-order residual |
|---|---|---|---|
| `P_19` | 4.30e-3 | 9.73e-4 | 6.67 |
| `P_23` | 2.89e-3 | 5.22e-4 | 6.35 |
| `P_31` | 1.55e-3 | 1.91e-4 | 5.69 |
| `R_19` | 1.88e-3 | 4.99e-4 | 3.42 |
| `R_27` | 1.24e-3 | 1.83e-4 | 3.61 |
| `R_35` | 8.54e-4 | 8.61e-5 | 3.69 |

On `R_11 … R_35`, `n³ ×` residual = 3.19, 3.30, 3.42, 3.53, 3.61, 3.66, 3.69. The first-order residual times `n²` drifts from 0.30 to 1.05 over the same range, which was the slow drift seen in HYP-9365.

| circulant classes | after 1st order: mean, sd | after 2nd order: mean, sd |
|---|---|---|
| `n = 21` (88) | 1.95e-3, 6.8e-4 | 3.8e-4, 1.8e-4 |
| `n = 25` (208) | 1.50e-3, 3.6e-4 | 2.3e-4, 7.1e-5 |

**Towers and all regular tournaments** (`second_order_regular.out`).
* Ternary towers: the residual goes from `+2.38·10^(−4)` to `−3.2·10^(−6)` at `S_3` (27 vertices, 75 times smaller), and from `+9.32·10^(−5)` to `+1.9·10^(−6)` at `S_4` (81 vertices, 49 times smaller).
* All 1223 regular 11-tournaments, where the bridged cores contribute: the residual mean drops from `5.5·10^(−3)` to `2.8·10^(−3)`, and the sd from `2.1·10^(−3)` to `1.1·10^(−3)`.
* The parameter-free ranking of `κ` there improves from Spearman 0.8956 to 0.9738, and still finds `P_11`.

**The doubly regular constant.** The predicted `z(q) = q²(κ_pred/e − 1 + 3/(2q))` for Paley `P_q` is 2.783, 2.845, 2.919, 2.978, 2.991, 3.018 at `q = 19, 23, 31, 43, 47, 59`.
* Extrapolating in `1/q` gives `c_2 ≈ 3.13`. The exact-data value is 3.07–3.10 (`drt_universal_curve.out`).
* The gap `z_exact − z_pred` is 0.334, 0.265, 0.178 at `q = 19, 23, 31`, i.e. `≈ 6/q`, the expected third order.
* The three-loop clusters alone are large, `q² X_2 → 28.9`. The finite-size and cross terms of the second order cancel most of this, leaving `c_2 − 5/3 ≈ 1.46`.

## The doubly regular constant is 25/8

For doubly regular tournaments every second-order term except `X_2` has an exact limit. The one-loop measure is `L = 2R` with `R ~ Poisson(1)`, and `n X_1(z) → −(7/2) z⁴`. Hence:

| term (×`n²`) | limit |
|---|---|
| `X_1²/2` | 49/8 |
| `X_1 E_0 a / n` | −21/2 |
| `δ/n` | −49 |
| `E_0 a²/2` | 21 |
| `E_0 b` | 29/6 |

The checks at `q = 43, 59, 83` approach these limits like `1/q` (`drt_c2_terms.out`).

* **The three-loop sum.** `q² X_2` for Paley `P_q` at ten primes `q = 19..83` (`drt_x2_limit.out`) fits

      q² X_2 = 29 − 21/q + 189/q² − 322.3/q³ + …

  The limit is 29.000000, and the `1/q` coefficient is −21.0000 in fits of degree 4–5 (residuals below `10^(−10)`).
* **The constant.** So `n² Y_2 → 29 − 27.5417 = 35/24`, and the first order contributes `5/3` (from `ρ/e` and the exact `y_pred = (3n² − 8n + 1)/(6(n−1)²)`). Therefore

      κ_DRT(n) = H 2^(n−1)/n! = e (1 − 3/(2n) + 25/(8n²) + O(n^(−3)))

  for every doubly regular tournament (DERIVED; the limit 29 is FINITE-EXACT to 7 digits).
* **Data.** The exact values `n²(κ/e − 1 + 3/(2n)) = 3.19, 3.12/3.13, 3.10–3.11, 3.09, 3.10` at `n = 15..31` scatter around `25/8 = 3.125`, with `O(1/n)` third-order corrections.
* That the constants 25, 29, 21, 189 echo the prompt's numbers is NUMEROLOGY.

## Status of the expansion

The coefficients grow with the order: excess-2 clusters carry Möbius weights up to `(s−1)!`, and `q² X_2 ≈ 29` against `c_2 ≈ 3`. So the loop expansion is presumably asymptotic, not convergent. The non-perturbative alternating end (HYP-9366) is exponentially small and invisible to it.
