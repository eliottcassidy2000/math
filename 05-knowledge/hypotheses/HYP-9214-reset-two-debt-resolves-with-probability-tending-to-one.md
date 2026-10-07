---
id: HYP-9214
title: "Reset-2 debt resolves almost surely for large sources: for a uniformly random first-reset-2 odd n of B bits, P(n and (n-1)/2 merge at equal odd-step time before 1) -> 1 as B -> infinity (empirically 0.19, 0.33, 0.41, 0.50, 0.61, 0.70, 0.79, 0.82 for B = 16, 32, ..., 2048; 0.88 at 4096, 0.91 at 8192)"
status: >
  OPEN. NUMERICAL (seeded Monte Carlo, 1000 samples per size, 600 at B = 2048; script sevens_20261006_debt_resolution_trend.py);
  the merges are collision switches in the sense of THM-4555 (independently audited 2026-10-06: every sampled merge up to
  B = 2048 is a collision; corrections in MISTAKE-572).
source: mac-mini-2026-10-06-mod1819, 05-knowledge/results/seven_twentyone_mersenne_openai_math_20261006.md, section 1
related:
  - 01-canon/theorems/THM-4555-uniform-switches-are-collisions-at-minus-one-trailing-ones-deletion.md
  - 05-knowledge/results/checked_switch_phase19_20261004.md (debt states (4)-(5); 239 residual seeds)
---

# HYP-9214 — the reset-2 debt resolves almost surely

**Conjecture.** Let `n` be a uniformly random odd `B`-bit number whose first reset exponent is 2. Then the probability that `n` and `(n − 1)/2` merge at equal odd-step time before reaching 1 tends to 1 as `B -> ∞`.

**Evidence (NUMERICAL, seed 2026).**

| `B` | 16 | 32 | 64 | 128 | 256 | 512 | 1024 | 2048 |
|---|---|---|---|---|---|---|---|---|
| merge rate | 0.189 | 0.330 | 0.409 | 0.500 | 0.605 | 0.695 | 0.785 | 0.817 |
| median lag beyond the run (noisy) | 8 | 15 | 23 | 30 | 57 | 73 | 128 | 140 |

The median-lag row is noisy. At `B = 128` the audit's 40 replicate runs gave medians from 33 to 50 (mean 39.6), so the entry 30 is not reproduced. The merge rates reproduce within about 1.5σ.

A second run (seed 99, 400 samples each) extends the table: 0.838 at `B = 2048`, 0.880 at 4096, and 0.910 at 8192. The non-merge share falls by a factor of about 0.74 per doubling of `B`, roughly `B^(−0.4)`.

* The fraction of sources that merge within 10 steps beyond the run stays near 0.12 (0.11–0.14) at every size. As a share of merges it falls from about 0.66 to 0.15. So the growth comes from late coalescence, not from short collisions.
* Even debt heights resolve more often than odd ones. The committed script gives 0.652 vs 0.580 at `B = 256` and 0.432 vs 0.398 at `B = 64`; an earlier seeded run gave 0.684 vs 0.574 and 0.453 vs 0.368.
  * The gap narrows at larger sizes (audit: 0.813 vs 0.771 at `B = 1024`, 0.861 vs 0.825 at `B = 2048`).
  * This is consistent with `checked_switch_phase19`, where 220 of the 239 residual seeds have odd normalized debt height.

**Model.**
* Write `x_j = 2^(L_j) y_j + Δ_j` for the two orbits. Then `L' = L + (b − a)` and `Δ' = (3Δ + 1 − 2^L)/2^a`.
* `L_j` performs (heuristically) a symmetric walk with increments `b − a`, the difference of the two exponents. The audit measured a mean increment of −0.004.
* The step before a first merge is a sibling pair `x = 4^k y + (4^k − 1)/3` with `k ∈ Z \ {0}`. In about 2/3 of merges `k < 0`. In walk variables this is `L = 2k`, `Δ = (4^k − 1)/3` while `|Δ| << y`. All 301 merges the audit checked pass through such a configuration.
* Recurrence of `L` makes repeated attempts likely. The open question is whether the per-visit success probability stays bounded below while `|Δ|` fluctuates.

**Consequence if true.** Combined with THM-4555 (iv), the trailing-ones switch `n => (n−1)/2` would certify almost all large reset-2 sources by collisions. The obstruction to a total rewrite family would then be a density-zero set, not a positive fraction.
