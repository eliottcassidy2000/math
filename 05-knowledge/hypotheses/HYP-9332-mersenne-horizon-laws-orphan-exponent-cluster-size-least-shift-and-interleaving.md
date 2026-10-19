---
id: HYP-9332
title: "Mersenne horizon laws: (a) the orphan fraction among odd exponents is K^(-alpha + o(1)) with alpha in [0.55, 0.65] (the constant fraction x sqrt(K) drifts from 1.0 to 0.8 between K = 5000 and 50000), (b) the mean number of cluster members below K equals (0.5-0.6) sqrt(sigma_T(x_K)), the SRW-skeleton envelope 1/sqrt(pi) = 0.564 of THM-4712 (5), and the mean cluster size is (2.0-2.4) sqrt(K), (c) the least shift to a smaller cluster-mate has P(D = d) ~ c d^(-2) for odd d >= 3 (95.5% of odd exponents have d = 1), (d) the interleaving depth of clusters stays bounded (<= 3 to K = 50000) while the hole fraction stays near 0.15"
status: >
  OPEN. NUMERICAL evidence: the exact horizon of THM-4712 (6) (every 2^K - 1, K <= 50000). (a) refines HYP-9242 (whose
  fit was 0.59 with CI [0.42, 0.75] at K <= 12800): the least-squares exponent over dyadic windows 400-50000 is 0.63, and
  fraction x sqrt(K) still decreases at the top of the range. (b) is THM-4712 (5)'s comparison. (c) is a fit on 20
  odd d with counts 492 ... 2 (exponent -2.01). (d) is the census of THM-4712 (6e). HEURISTIC for (b): linearity plus
  the SRW skeleton of THM-4581; no mechanism is proposed for (c) or (d). What would refute: (a) the exponent settling
  outside [0.5, 0.7] with 4x more data; (b) the ratio leaving [0.4, 0.7]; (c) the d^-2 law failing on the d >= 43 tail
  with more exponents; (d) an exponent covered by 4 cluster spans.
session: mac-mini-2026-10-19-mersenne
source: 05-knowledge/results/mersenne_horizon_20261019_clusters_orphans_reanchoring.md
related:
  - THM-4712 (the data and the Haar-model mean-cluster law), THM-4581 (pair chains, rate 1/2), THM-4556 (vi)
  - HYP-9242 (orphan law, K <= 12800), HYP-9243 (diffusive coalescence, O(sqrt E) cluster span)
  - 04-computation/experiments/mersenne_horizon_20261019/ (horizon.c, analyze_horizon.py, analyze_horizon_extra.py, data .gz)
---

# HYP-9332 — Mersenne horizon laws

## Statements

* **(a) Orphan exponent.** The fraction of odd `K` in a dyadic window around `K` with no deletion child merging before 1 is `K^(−α + o(1))`, with `α` between `0.55` and `0.65`. Data: `fraction · sqrt(K_mid) = 0.92, 0.99, 0.84, 0.87, 0.76` on `6401–12800, 12801–19200, 19201–25600, 25601–38400, 38401–50000`; HYP-9242's `1/sqrt K` reading needs the constant to stop drifting.
* **(b) Members below `K` and cluster size.** `E[#cluster members below K] = κ · sqrt(σ_T(x_K))` with `κ ∈ [0.5, 0.6]` for large `K`, against the SRW envelope `1/sqrt π = 0.564` of THM-4712 (5); the mean cluster size (clusters with minimum in a window around `K`) is `(2.0–2.4) sqrt K`, the orphan fraction being `2/(mean size)` per odd `K` (barriers note §4.2).
* **(c) Least-shift law.** For odd `K` with a smaller cluster-mate, `P(least shift = d) ≈ c d^(−2)` for odd `d ≥ 3`; `P(d = 1) = 0.955`. Counts to 50000: `492, 204, 129, 75, 50, 35, 16, 13, 9, 11, 10` for `d = 3, …, 23`.
* **(d) Shallow interleaving.** Clusters are not intervals (23 of 231 are), but at most three cluster spans cover any exponent, and the hole fraction of the total span is `0.147`.

## Why these are plausible, and what is missing

* (a)–(b): in the Haar model the `D`-chains are symmetric level walks started at `D` (THM-4581); the mean number absorbed by time `T` is `Θ(sqrt T)` up to a root-log (THM-4712 (3)–(4)), and `T = σ_T(x_K) ≈ 7.6 K`. The orphan event `N_T = 0` is not controlled by the mean; its exponent exceeding `1/2` would mean that the first absorption (over all `D`) has a lighter tail than the single `D = 1` chain, which is what the least-shift law suggests.
* (c): no derivation. Under independent chains the least shift would decay super-exponentially; the power law shows that the chains fail together: when `D = 1` fails, the pair `(K−1, K−2)` belongs to a foreign cluster and the next candidates `K−3, K−4` are reached through the same mechanism one block further down.
* (d): the bounded depth suggests a one-dimensional coalescent with sparse defects rather than an ultrametric hierarchy; a proof would have to control how two clusters can cross.
