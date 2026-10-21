---
id: HYP-9245
title: "Rank-one integer laws for the base-p Collatz maps C_p: (a) if C_p has several positive-integer cycles whose basins have limiting densities, the density of basin boundaries (n, n+1 in different basins) below N is of the order of q_p(T_N), the Haar probability that y and y+1 have not merged by the typical orbit length T_N = ln N/|Lambda_p| (THM-4590's cut heuristic for 3x + 1, transplanted); with the rank-one tail q_p(T) of order T^(-1/2) this is (ln N)^(-1/2); (b) on the repunit line R_K = p^K - 1 the fraction of exponents in [X, 2X) that are least members of their partner class (exponents whose orbits first reach a cycle at the same element and the same height n - K) is of order X^(-1/2), the base-p form of HYP-9243's Mersenne-line prediction"
status: >
  OPEN. Restated after the session's independent audit (MISTAKE-649): attributions to THM-4590 and HYP-9243, the
  density assumption made explicit, 'partner class' defined, refutation test (b) replaced by a rule simulated under its
  null (zp_orphan_fit.py (O4)), and test (a) restated as a pre-registered ratio band (not null-checked).
  PROVED around it (THM-4610): Haar coalescence for C_p (every odd prime p); basin boundaries have density 0;
  q_p(T) >= c T^(-1/2).
  NUMERICAL, part (a) (zp_boundary_ratio.py; q_p(T_N) from 10^5 Haar paths, standard error <= 0.0016):
    - R(N) = boundary/(2s(1-s) q_p(T_N)) at N = 10^3..10^7: p = 3: 1.65, 1.61, 1.40, 1.26, 1.19; p = 11: 1.14, 1.12,
      1.06, 1.06, 1.05; p = 17: 1.12, 1.13, 1.08, 1.07, 1.06; p = 23: 1.34, 1.12, 1.06, 1.04, 1.04.
    - The ratio stays between 1.04 and 1.65 with a net downward drift at every p. Whether it tends to 1 is not resolved.
    - Second-basin shares 0.967, 0.018, 0.292, 0.234, stable from 10^6 (zp_basins_survey.py (S2)).
    - p = 3: 0.459 of consecutive pairs n <= 10^6 merge at equal time above the cycles, against the Haar prediction
      0.498 +- 0.003 at T = 79 (S3).
  NUMERICAL, part (b) (zp_orphan_fit.py on zp_repunit_lines.py's data, K <= 8000): Poisson maximum-likelihood exponents
    on K in [800, 8000], orphan events treated as independent: p = 7: 0.57 [0.21, 0.93]; p = 11: 0.69 [0.37, 1.01];
    p = 13: 0.71 [0.39, 1.03] (95% profile-likelihood intervals). All three contain 1/2; none resolves the exponent.
  HEURISTIC: the mechanism (section Why).
source: opus-2026-10-08-S22, 05-knowledge/results/zp_rank_one_coalescence_20261008.md (sections 3-5)
related:
  - 01-canon/theorems/THM-4610-rank-one-haar-coalescence-for-the-base-p-collatz-maps-on-z-p.md
  - 01-canon/theorems/THM-4590-collatz-classes-are-equidistributed-residues-windows-slowly-varying-density.md (p = 2: o(x) cuts; the cut heuristic (1 - sum d_i^2) q(4.8 log2 x), whose measured level is about 1.13 times the prediction; the negative-integer basins of -1, -5, -17)
  - 05-knowledge/hypotheses/HYP-9242-orphan-law-deletion-orphans-decay-like-inverse-square-root-of-log-n.md (p = 2 orphan law; fitted 0.59 [0.42, 0.75]; leaves faster decay open)
  - 05-knowledge/hypotheses/HYP-9243-diffusive-coalescence-of-the-mersenne-line-and-the-limit-of-deletion-descent.md (predicts the Mersenne orphan fraction K^(-1/2); counts orphans among odd K)
  - 05-knowledge/hypotheses/HYP-9244-polya-trichotomy-for-haar-coalescence-the-debt-lattice-rank.md
  - Goncalves-Greenfeld-Madrid arXiv:2111.06170, Remark 1.7 (numerically, the basins of the three cycles of 3x - 1 have densities of about 1/3 each)
scripts:
  - 04-computation/experiments/zp_rank_one_20261008/zp_basins_survey.py (+ .out) (S1)-(S4)
  - 04-computation/experiments/zp_rank_one_20261008/zp_boundary_ratio.py (+ .out)
  - 04-computation/experiments/zp_rank_one_20261008/zp_tails.py (+ .out)
  - 04-computation/experiments/zp_rank_one_20261008/zp_repunit_lines.py (+ .out) (R4)
  - 04-computation/experiments/zp_rank_one_20261008/zp_orphan_fit.py (+ .out)
---

# HYP-9245 — rank-one integer laws for the base-p Collatz maps

## Statement

Let `C_p(x) = x/p` if `p | x` and `⌈(p+1)x/p⌉` otherwise, on the positive integers. Let `Λ_p = (1/p) ln(1/p) + ((p−1)/p) ln((p+1)/p)` and `T_N = ln N/|Λ_p|`. Let `q_p(T)` be the Haar probability that `y` and `y + 1` have not merged at equal time by time `T`: the pair chain of THM-4610 from `(0, 1)`.

* **(a) Basin boundaries.**
  * *Assumptions.* `C_p` has cycles `c_1, …, c_r` (`r ≥ 2`) on the positive integers. Every positive integer reaches one of them. The basins `B_i` have natural densities `s_i`. The last is an assumption; for `3x − 1` compare Gonçalves–Greenfeld–Madrid, Remark 1.7.
  * *Claim.* Let `b(N) = #{n < N : n and n + 1 in different basins}/N`. Then

        b(N) ≍ q_p(T_N),

    that is, `R(N) = b(N)/(Σ_(i≠j) s_i s_j · q_p(T_N))` stays between two positive constants.
  * *Consequence.* `q_p(T) ≍ T^(−1/2)` is the rank-one tail law. Its lower bound is PROVED in THM-4610 5(c); its upper bound is at the sketch level of THM-4581 (4). With it, `b(N) ≍ (ln N)^(−1/2)`.
* **(b) Repunit orphans.**
  * *Partners.* For `R_K = p^K − 1`, let `n(K)` be the number of steps to the first value on a cycle, and `e(K)` that value. Call `K` and `K′` *partners* if `n(K) − K = n(K′) − K′` and `e(K) = e(K′)`. Equivalently, `R_K` and `R_K′` merge with time shift `K − K′` no later than their entry into the cycle: the S21 exit `R_K ⇝ R_K′`, completed by cycle entry. For `K ≤ 8000` this pair also determines the multiplication count and the cycle (`zp_orphan_fit.py` (O1)).
  * *Orphans.* `K` is an *orphan* if no smaller exponent is its partner.
  * *Claim.* The fraction of orphans in `[X, 2X)` is `≍ X^(−1/2)`.
  * *The case `p = 2`.* This is HYP-9243's prediction for the Mersenne line; HYP-9242 leaves faster decay open. Both count orphans among odd `K` only, since Mersenne class bottoms with `K ≥ 4` are odd. Per exponent their fraction is half as large.

## Why

* Integers below `N` take about `T_N` steps to reach the cycles. This is the base-`p` analogue of the standard heuristic `(2/ln(4/3)) ln n` for the Terras map; for the repunit lines it is confirmed to 0.2% in `zp_repunit_lines.py` (R3).
* Two neighbours that have not merged by then land in the basins roughly independently. This is an assumption. THM-4590 makes the same assumption for `3x + 1`, and its measured cut level is about 1.13 times the prediction; here `R(N)` is 1.04–1.65.
* By THM-4610 the probability of not having merged is the Haar tail at that time. That tail is `≍ T^(−1/2)` because the debt's flip skeleton is a fair simple random walk.
* The same first-return law governs which repunits find a smaller partner before the horizon. That is HYP-9243's mechanism on the Mersenne line.

## What would refute it

* **(a)** Assume the densities exist. Compute `b(N)` and `R(N)`, with `q_p(T_N)` from at least `10^5` Haar paths:
  * for `p = 3` up to `N = 10^10`;
  * for `p = 11, 17, 23` up to `N = 10^9`.

  A segmented basin sieve is needed. (a) is refuted for that base if `R(N)` leaves `[1/2, 2]` at any such `N`.
  * The `(ln N)^(−1/2)` form also needs the upper tail bound.
  * A slope test of `b(N)` against `ln N` is not used. At accessible `N`, `q_p(T_N)` is still pre-asymptotic for `p ≥ 11`: it falls only from 0.849 to 0.839 between `10^6` and `10^7` for `p = 11`.
* **(b)** Compute the orbits of `R_K` for every `K ≤ 10^5`, for `p = 7, 11, 13`.
  * This costs about 2000 times the `K ≤ 8000` run, which took 18 minutes on 12 cores. It needs a faster integer library or many more cores.
  * Fit `c K^(−α)` by Poisson maximum likelihood on `K ∈ [10^4, 10^5]`. (b) is refuted for that base if the 99% profile-likelihood interval for `α` excludes 1/2.
  * *Validity.* Under independent Bernoulli(`c K^(−1/2)`) orphans, this rule rejects a true law in 0.8% of simulated trials. It rejects `α = 1` in 99.4% of trials and `α = 1/4` in 53% (`zp_orphan_fit.py` (O4)).
  * *Correlations.* Orphan events may be correlated along the line, so a block bootstrap over `K` should confirm the interval's width before a rejection is accepted.
  * *Why the original rule was dropped.* It said: refute if the local exponents of the last two dyadic windows both fall outside `[0.35, 0.65]`. That rule rejects a true `K^(−1/2)` law in about 40% of trials (O4).
