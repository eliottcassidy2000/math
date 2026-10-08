---
id: HYP-9243
title: "Diffusive coalescence of the Mersenne line: the deletion clusters grow with exponent 1/2 in Terras time, both for the size-biased cluster N_E(T) of deletions absorbed into a generic odd source by time T (from T = 2^12 at E ~ 1e11) and for the per-cluster mean size at the orbit horizon T = sigma_T(x_K) ~ 7.64 K, where the clusters are the exact stopping-time level sets (THM-4605 statement 8) and the statement predicts that the Mersenne orphan fraction (= 2/(mean cluster size) per odd K) decays like K^(-1/2), one of the readings HYP-9242 leaves open"
status: >
  OPEN.
  NUMERICAL:
  (i) Six fresh odd E in [1e10, 1e11], deletions D <= 2048: upper median N(T) = 28, 138, 332 at T = 2^12, 2^17, 2^20;
  median N(T)/sqrt(T) between 0.23 and 0.52; 1-4 unabsorbed groups at 2^20. Six sources only; no trend is resolved.
  (ii) Exact horizon clusters for K <= 12800 from mac-mini's table (clusters = level sets of the stopping-time invariant,
  PROVED in THM-4605 statement 8). Per-cluster mean: 2/(mean size) x sqrt(K) = 1.49, 1.02, 1.05, 1.13, 0.90 on the dyadic windows
  [400, 800) .. [6400, 12800], which is HYP-9242's orphan fraction x sqrt(K) (1.55, 1.09, 1.01, 1.14, 0.89) up to window edges;
  13-30 clusters per window. Size-biased: median (cluster size below K)/sqrt(sigma_T(x_K)) = 0.25, 0.36, 0.40, 0.34, 0.44,
  a local exponent of about 0.6-0.7 in this range.
  The contiguous extent S(T) is not sqrt(T)-like at the horizon (median S/sqrt(T) rises from 0.08 to 0.27), so the hypothesis is
  stated for cluster sizes.
  HEURISTIC: the mechanism, a symmetric level walk between clusters. It is PROVED only in the Haar model (THM-4581 (1b));
  over exponents the coins are Haar-distributed in the limit (note section 4.4); for one integer source it is an assumption.
  Scope against THM-4593: that theorem pins the exponent 1/2 of the probability that a source has NO exit (finite lag sets,
  including Mersenne lags D <= 61; its 0.69 is a coalescence transient). This hypothesis concerns the GROWTH of the cluster,
  with D_max >= C sqrt(T), and the horizon level sets.
  The checkpoint version's "consequence" (deletion descent cannot ground giant Mersennes) is now PROVED as THM-4605
  statement 8 and is no longer part of this hypothesis. Its band agreement with HYP-9242 was a bookkeeping artifact (MISTAKE-588).
source: opus-2026-10-07-S21, 05-knowledge/results/mersenne_line_barriers_20261007.md (sections 4.2, 4.3, 4.5)
related:
  - 01-canon/theorems/THM-4605-mersenne-exits-by-pair-chain-absorption-the-t23-fan-escapes-at-depth-1.9e7.md (statement 8)
  - 01-canon/theorems/THM-4593-partner-coalescence-pins-the-diffusive-exponent-and-all-lags-merge-exponentially.md (the no-exit tail; one coin per step)
  - 05-knowledge/hypotheses/HYP-9242-orphan-law-deletion-orphans-decay-like-inverse-square-root-of-log-n.md
  - 01-canon/theorems/THM-4581-haar-coalescence-affinely-related-collatz-orbits-merge-almost-surely.md ((1b), (4), (7))
scripts:
  - 04-computation/experiments/mersenne_line_coalescence_20261007.py (+ .out) (S), (H2)
  - 04-computation/experiments/mersenne_line_barriers_20261007.py (+ .out) (F), (M)
---

# HYP-9243 — diffusive coalescence of the Mersenne line

## Statement

* **The size-biased cluster.** For a source exponent E and a Terras time T counted from `x_E`, let `N_E(T)` be the number of deletions `D ≤ D_max` whose source-reference chains (THM-4605) are absorbed by time T. This is the size of the source's cluster below it.
* **The horizon clusters.** At the horizon `T = σ_T(x_K)` the clusters are exactly the level sets of `(σ_T(M_K) − K, o(K))` (PROVED, THM-4605 statement 8).
* **Conjecture.** For generic odd E, `D_max ≥ C√T` and `1 ≪ T ≤ σ_T(x_E)`:
  * `N_E(T) = T^(1/2 + o(1))` in distribution;
  * the per-cluster mean size of the horizon clusters near K is `K^(1/2 + o(1))`.
* The second statement says that the Mersenne orphan fraction decays like `K^(−1/2+o(1))`. That fraction is 2/(mean cluster size) per odd K, since every cluster bottom with K ≥ 4 is odd (reset rule).
* HYP-9242, as corrected after its audit C, fits the exponent as α̂ = 0.59 with 95% CI [0.42, 0.75]. The value 1/2 is inside, so the two hypotheses are consistent. HYP-9242 leaves faster decay open; this hypothesis predicts 1/2.
* The constants at finite T and at the horizon are not claimed to be equal.
* In this range the exact horizon data support the per-cluster statement better than the size-biased one, whose local exponent is about 0.6–0.7.

## Mechanism (HEURISTIC)

* Two clusters are related by `y = 3^k x + e`.
* All relations to a common reference move in the same direction when they move (THM-4593 (1), "one coin per step"). So the relative level of two clusters changes only when exactly one of them flips.
* In the Haar model that relative level is itself a THM-4581 chain, a symmetric recurrent walk at flips (THM-4581 (1b)). For one integer source this is an assumption.
* Exponents D apart start at level −D and need about `D²` steps to meet (THM-4581 (7), heuristic).
* A merge also needs the offset to match, about ten visits to level 0 on average. So walks pass through one another and clusters interleave. At the fan bottom, the clusters with D ranges 1911–2902 and 2895–4000 overlap, and the group of D = 1 is not an interval.
* Within a lag set the partners coalesce among themselves first (THM-4593's Arratia-type reading). This is why the no-exit tail has exponent 1/2 rather than a larger one.
* **Scaling analogy.** One-dimensional coalescing random walks (Arratia 1979, 1981; Bramson–Griffeath 1980) have density about `(πt)^(−1/2)`. Order is not preserved here, so the analogy covers the scaling only.
* **Normalisation.** S19's debt L has variance 4 per odd step (`L ≈ 2k`). The chain level k has variance about 1/2 per Terras step.
* **Prior art.** LaDue (arXiv:1709.02979, 2017) studies clusters of integers with equal total stopping times; the clusters here are the Mersenne-line case, with exact stopping-time invariants.

## Evidence

* See the status. The two regimes are measured independently:
  * chains at E ≈ 10^11 for T ≤ `2^20`;
  * exact orbits at the horizon for K ≤ 12800.
* The t = 23 fan bottom is not evidence. It is one selected source with N(T) = 0 until 19,000,765, followed by one heavy-tailed jump to 1916.

## What would refute it

* **The size-biased cluster.** Take at least 20 fresh odd E near 10^11 with `D_max = 2^13`. Refuted if the median `N_E(T)` has a local exponent outside [0.3, 0.7], sustained over `T = 2^14` to `2^24`.
* **The per-cluster mean.** Compute exact horizon clusters for K up to 10^5. Refuted if (mean cluster size)/`√K` drifts monotonically by more than a factor 2 per decade of K. Equivalently: the orphan fraction times `√K`.
