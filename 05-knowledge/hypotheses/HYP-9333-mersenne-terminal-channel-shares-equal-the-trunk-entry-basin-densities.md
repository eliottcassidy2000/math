---
id: HYP-9333
title: "The Mersenne line enters 1 only through 5, 85 and 341, and the shares of these three terminal channels (per cluster, 215/8/8 of 231 to K = 50000; per exponent 0.924/0.038/0.038) converge to the trunk-entry basin densities of THM-4517, e_2 = 0.93796, e_4 = 0.0236, e_5 = 0.0378, with e_i = 0 for 3 | i and e_7, e_8, ... negligible"
status: >
  OPEN. FINITE-EXACT: the channel census of THM-4712 (6f) (every K <= 50000; no cluster mixes channels; 5461 and
  21845 never occur). NUMERICAL: the per-cluster shares 0.931 / 0.035 / 0.035 against 0.938 / 0.024 / 0.038, within the
  Poisson noise of 8 +- 3 clusters. Not claimed: that the Mersenne line (a null set of sources) must follow the Haar
  basin densities; the channel of a cluster is decided by a single orbit, so the effective sample is the number of
  clusters, not of exponents. What would refute: a fourth channel (5461 or 21845) on the Mersenne line, or the 85-share
  per cluster settling away from 0.024 with ten times more clusters.
session: mac-mini-2026-10-19-mersenne
source: 05-knowledge/results/mersenne_horizon_20261019_clusters_orphans_reanchoring.md
related:
  - THM-4517 (no root-uniform positive density of Collatz basins; the trunk-entry spectrum e_i), THM-4712 (6f), THM-4556 (ii)
---

# HYP-9333 — terminal channels of the Mersenne line

* Every orbit reaches 1 through an odd predecessor `a_m = (4^m − 1)/3`; `a_m` is reachable from above only when `3 ∤ a_m`, i.e. `3 ∤ m`. THM-4517 measured the basin densities of these entry points among all integers below `2^32`: `e_2 = 0.93796` (through 5), `e_4 = 0.023647` (85), `e_5 = 0.037789` (341), `e_i = 0` to seven decimals for `3 | i`, the rest negligible.
* On the Mersenne line to `K = 50000` only the first three channels occur. The hypothesis is that the cluster-level shares converge to `(e_2, e_4, e_5)` as `K → ∞` and that `5461 = a_7` and `21845 = a_8` eventually occur with their tiny densities.
* The orbit of a cluster's bottom decides the channel of the whole cluster (a cluster shares one future), so the 215/8/8 split is a sample of 231 trials: the 85-channel count 8 is `1.4` standard deviations above `231 · 0.0236 = 5.5`.
