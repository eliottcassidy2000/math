---
id: HYP-9214
title: "Reset-2 debt resolves almost surely for large sources: for a uniformly random first-reset-2 odd n of B bits, P(n and (n-1)/2 merge at equal odd-step time before 1) -> 1 as B -> infinity (empirically 0.19, 0.33, 0.41, 0.50, 0.61, 0.70, 0.79, 0.82 for B = 16, 32, ..., 2048)"
status: >
  OPEN. NUMERICAL (seeded Monte Carlo, 1000 samples per size, 600 at B = 2048); the merges are collision switches
  in the sense of THM-4555 (all equal-length merges with (n+1)/2^D - 1 found below 2*10^4 are collisions).
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
| median lag beyond the run | 8 | 15 | 23 | 30 | 57 | 73 | 128 | 140 |

* The share of merges within 10 steps stays near 0.12 at every size. So the growth comes from late coalescence, not from short collisions.
* Even debt heights resolve more often than odd ones: 0.684 vs 0.574 at `B = 256`, and 0.453 vs 0.368 at `B = 64`. This is consistent with `checked_switch_phase19`, where 220 of the 239 residual seeds have odd normalized debt height.

**Model.**
* Write `x_j = 2^(L_j) y_j + Δ_j` for the two orbits.
* `L_j` performs a symmetric walk with increments `b − a`, the difference of the two exponents.
* A merge requires a visit to a sibling configuration `L = 2k`, `Δ = (4^k − 1)/3`.
* Recurrence of `L` makes repeated attempts likely. The open question is whether the per-visit success probability stays bounded below while `|Δ|` fluctuates.

**Consequence if true.** Combined with THM-4555 (iv), the trailing-ones switch `n => (n−1)/2` would certify almost all large reset-2 sources by collisions. The obstruction to a total rewrite family would then be a density-zero set, not a positive fraction.
