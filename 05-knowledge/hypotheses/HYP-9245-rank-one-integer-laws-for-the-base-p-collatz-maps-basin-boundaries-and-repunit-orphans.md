---
id: HYP-9245
title: "Rank-one integer laws for the base-p Collatz maps C_p: (a) when C_p has several positive-integer cycles, the density of basin boundaries (n, n+1 in different basins) up to N is of order q_p(T_N), the Haar no-merge probability at the typical orbit length T_N = ln N/|Lambda_p|, so it tends to 0 like (ln N)^(-1/2) on a logarithmic clock; (b) on the repunit line R_K = p^K - 1 the orphan fraction (partner-class bottoms) decays like K^(-1/2), the rank-one first-return exponent, in every base (p = 2 is HYP-9242)"
status: >
  OPEN.
  PROVED around it (THM-4610): Haar coalescence for C_p (every odd prime p); basin boundaries have density 0; P(no merge by T) >= c T^(-1/2).
  NUMERICAL, part (a): p = 3, 11, 17, 23 at N = 10^7.
    - Second-basin shares 0.967, 0.018, 0.292, 0.234, stable from 10^6.
    - Boundary densities 0.0366, 0.0315, 0.3925, 0.3439.
    - In units of 2s(1-s) (the value for independent landing): 0.577, 0.883, 0.948, 0.959.
    - Haar tails at T_N (92, 116, 143, 169 steps): 0.505, 0.842, 0.902, 0.923. Ratios 1.14, 1.05, 1.05, 1.04, consistently above 1, as expected when smaller n have shorter orbits.
    - For p = 3, boundary x sqrt(ln N) settles near 0.147-0.150 from 10^6 to 10^7. For p >= 11 the boundary is still flat at these sizes, as the lazy skeleton (flip rate 2/p) predicts.
    - Also for p = 3: 0.459 of consecutive pairs n <= 10^6 merge at equal time above the cycles, against the Haar prediction 0.498 at T = 79 (S3).
  NUMERICAL, part (b): repunit lines to K = 8000.
    - p = 7: orphan fraction x sqrt(K) = 0.589, 0.565, 0.547, 0.581 on [800,1600) .. [6400,8000] (supports 1/2).
    - p = 11: 0.799, 0.773, 0.673, 0.581.
    - p = 13: 0.883, 0.684, 0.652, 0.634.
    - For p = 11 and 13 the local exponents run 0.55-0.87, so the exponent is not resolved there.
  HEURISTIC: the mechanism. A pair still unmerged when its orbits reach the cycles lands in independent-looking basins. The Haar tail is q_p(T) ~ C_p T^(-1/2), with a constant that grows with p because flips happen at rate 2/p.
source: opus-2026-10-08-S22, 05-knowledge/results/zp_rank_one_coalescence_20261008.md (sections 3-5)
related:
  - 01-canon/theorems/THM-4610-rank-one-haar-coalescence-for-the-base-p-collatz-maps-on-z-p.md
  - 01-canon/theorems/THM-4590-collatz-classes-are-equidistributed-residues-windows-slowly-varying-density.md (p = 2: o(x) cuts; the negative-integer basins of -1, -5, -17)
  - 05-knowledge/hypotheses/HYP-9242-orphan-law-deletion-orphans-decay-like-inverse-square-root-of-log-n.md (p = 2 orphan law; fitted 0.59 [0.42, 0.75])
  - 05-knowledge/hypotheses/HYP-9243-diffusive-coalescence-of-the-mersenne-line-and-the-limit-of-deletion-descent.md
  - 05-knowledge/hypotheses/HYP-9244-polya-trichotomy-for-haar-coalescence-the-debt-lattice-rank.md
scripts:
  - 04-computation/experiments/zp_rank_one_20261008/zp_basins_survey.py (+ .out) (S1)-(S3)
  - 04-computation/experiments/zp_rank_one_20261008/zp_tails.py (+ .out)
  - 04-computation/experiments/zp_rank_one_20261008/zp_repunit_lines.py (+ .out) (R4)
---

# HYP-9245 — rank-one integer laws for the base-p Collatz maps

## Statement

Let `C_p(x) = x/p` if `p | x` and `⌈(p+1)x/p⌉` otherwise, on the positive integers. Let `Λ_p = (1/p) ln(1/p) + ((p−1)/p) ln((p+1)/p)` and `T_N = ln N/|Λ_p|`. Let `q_p(T)` be the Haar probability that `y` and `y + 1` have not merged by time `T` (THM-4610).

* **(a) Basin boundaries.** Suppose `C_p` has at least two cycles on the positive integers, with limiting basin shares `s` and `1 − s` (or several shares `s_i`). Then

      #{n < N : n and n+1 in different basins}/N  ≍  q_p(T_N),

  with the constant close to `Σ_(i≠j) s_i s_j` (`2s(1 − s)` for two basins). Since `q_p(T) ≍ T^(−1/2)` is the rank-one law (lower bound PROVED in THM-4610 5(c)), the boundary density is `≍ (ln N)^(−1/2)`.
* **(b) Repunit orphans.** On `R_K = p^K − 1`, the fraction of `K` in `[X, 2X)` that are the least member of their partner class is `≍ X^(−1/2)`. For `p = 2` this is HYP-9242.

## Why

* Integers below `N` take about `T_N` steps to reach the cycles.
* Two neighbours that have not merged by then behave like independent draws among the basins.
* By THM-4610 the probability of not having merged is the Haar tail at that time. That tail is `≍ T^(−1/2)` because the debt skeleton is a fair simple random walk.
* The same first-return law governs which repunits find a smaller partner before the horizon.

## What would refute it

* **(a)** Take `p = 3` and `N` up to `10^10` (segmented sieve), or `p = 11`, `17`, `23` with `N` up to `10^9`. A boundary density whose ratio to `2s(1 − s) q_p(T_N)` drifts outside `[0.7, 1.5]`, or whose decay is faster than every power of `ln N`, refutes (a).
* **(b)** Take repunit lines to `K = 3·10^4` for `p = 7, 11, 13`. A sustained local exponent outside `[0.35, 0.65]` over the last two dyadic windows refutes (b) for that base.
