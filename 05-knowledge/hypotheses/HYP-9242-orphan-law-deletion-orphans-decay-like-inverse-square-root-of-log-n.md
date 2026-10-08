---
id: HYP-9242
title: "Orphan law: the proportion of residual first-reset-2 sources n with no deletion partner (no child (n+1)/2^D - 1 merging at equal time before 1) decays like c/sqrt(log n); on the Mersenne line the orphan fraction among odd K in [X, 2X] is about 1/sqrt(X) (fraction*sqrt(K) = 1.09, 1.01, 1.14, 0.89 on the windows [800,1600), ..., [6400,12800]); orphans are enriched 10x among exponents with long normalized orbits, are not explained by single long post-run ones-runs or by the leading digits of 3^K, and are rescued by shallow 3-adic predecessors at the base rate (about 40%)"
status: >
  OPEN. UPPER BOUND at the sketch level of THM-4581 (4), for Haar-random t: P(orphan) <= C L^(-1/2) (log L)^2 (see 'Partial result').
  NUMERICAL: Mersenne orbits computed exactly for every K <= 12800 (odd-step count and Terras length to 1; partner = smaller
  exponent with equal odd count and Terras difference K - K', which THM-4556 (ii) shows equivalent to an equal-time merge before 1 for K <= 6000);
  general residual sources n = 2^K t - 1 (K = 9, 17; t of 25..800 random bits; 300 per row). HEURISTIC: the pair chain's debt is a
  martingale with non-absorption tail T^(-1/2) (THM-4581) and the time available before the child reaches 1 is proportional to log n.
source: mac-mini-2026-10-07-twoanchor (continuation), 05-knowledge/results/runcompress_orphans_cayley_20261007.md (section 3)
related:
  - 01-canon/theorems/THM-4556-the-mersenne-line-is-a-chain-of-debt-states-odd-shift-distance-2-adic-periodicity.md ((vi): 50 of 2500 odd a in [10^3, 6*10^3])
  - 01-canon/theorems/THM-4581-haar-coalescence-affinely-related-collatz-orbits-merge-almost-surely.md (the T^(-1/2) rate)
  - 01-canon/theorems/THM-4601-two-anchor-reduction-of-the-residual-first-reset-two-branch.md, THM-4603 (drift)
scripts:
  - 04-computation/experiments/runcompress_20261007/ (mersenne_sigma*.py, orphan_law_mersenne.py (+ .out), orphan_law_general.py (+ .out),
    mersenne_partners_by_shell.py, orphan_structure.py, orphan_time.py, orphan_rescue.py; data mersenne_sigma_12800.txt)
---

# HYP-9242 — the orphan law

## Statement

* A residual source n (first reset letter 2) is a **deletion orphan** if no child `h_D = (n+1)/2^D − 1`, `1 ≤ D ≤ K−1`, merges with it at equal time before 1. Equivalently, no D has `o(h_D) = o(n)` and `σ_T(h_D) = σ_T(n) − D`.
* Conjecture: among residual sources of size about 2^L, the orphan proportion is `≍ L^(−1/2)`. The constant depends on how many children are available.
* On the Mersenne line (n = 2^K − 1, K odd) the orphan fraction is about `1/√K`.

## Evidence (NUMERICAL)

**Mersenne windows:**

| K window | odd K | orphans | fraction | fraction × √K_mid |
|---|---|---|---|---|
| [200, 400) | 100 | 6 | 0.060 | 1.01 |
| [400, 800) | 200 | 13 | 0.065 | 1.55 |
| [800, 1600) | 400 | 13 | 0.0325 | 1.09 |
| [1600, 3200) | 800 | 17 | 0.0213 | 1.01 |
| [3200, 6400) | 1600 | 27 | 0.0169 | 1.14 |
| [6400, 12800] | 3200 | 30 | 0.0094 | 0.89 |

* Over the last three doublings the local decay exponent is about 0.6. That is 1/2 within the counting noise (30 ± 5.5 in the last window).
* **General residual sources:** at K = 9 the orphan fraction falls from 0.557 (B = 25) to 0.087 (B = 800); at K = 17 from 0.427 to 0.070. Scaled by √(K + B) this is about 3 and about 2.2 respectively.
* **Long orbits.** Split odd K in [1000, 12800] into thirds by normalised post-run length `(σ_T − K)/K`. The orphan rates are 0.31%, 0.92% and 3.0%.
  * A long orbit means excess odd density in the source's tail: a positive debt drift (THM-4603 (1)) in aggregate.
  * It is not a single long ones-run. The orphans' longest post-run ones-runs match the non-orphans' (Mann–Whitney z = −0.6, K in [1000, 6400]).
* **No digit effect.** The orphan rate is 1.9–2.3% in each quarter of `frac(K log2 3)`.
* **Nearest partner.** It is K−1 (D = 1) for 92% of non-orphans with K ≥ 1000. In each shell v2(K−1), the reset pair succeeds in 82–91% of cases and h_3/h_4 in 82–88% (K in [600, 6400]).
* **3-adic rescue.** A smaller node in the Terras backward tree within depth 36 certifies n. This holds for 37.7% of the 53 orphans in [1000, 6400], against 40% for all odd K. So the rescue looks independent of orphan status.
* **Doubly uncovered.** 33 of 2700 odd K in [1000, 6400] (1.2%) have neither a deletion partner nor a shallow predecessor.

## Partial result: the upper bound for random t

**Claim** (at the sketch level of THM-4581 (4)). Fix K, and let t be uniform among the odd L-bit integers in the residual class. Then

    P(n = 2^K t − 1 is a deletion orphan) ≤ C(K) L^(−1/2) (log L)^2.

**Reasons.**
* Every Terras step at most halves, so the child h_1 = (n−1)/2 cannot reach 1 before Terras time `log2 h_1 ≥ L − 1`.
* Hence an absorption of the D = 1 chain by time `L − K − 2` is an equal-time merge at a value ≥ 2. That is a deletion certificate.
* The chain's driving bits up to that time are the child's first parity bits. These are exactly uniform over t in the class.
* The two-run length J is geometric, so the post-run state (J, 1) (THM-4601 (v)) has a debt with finite mean.
* THM-4581 (4) gives non-absorption probability ≤ C T^(−1/2) (log T)^2 from a state of bounded debt. Averaging over J gives the claim.

**Open.** The matching lower bound, and the Mersenne line, whose post-run bits are equidistributed in K only to depth about log2 K (THM-4601 (vi)).

## Why

* The deletion certificate is the absorption of the pair chain's debt walk.
* That walk has zero drift on average and a `T^(−1/2)` non-absorption tail (THM-4581).
* The time available before the child reaches 1 is about 4.8 times the child's bit length.
* So P(orphan) ≈ c/√(log n), with a constant that decreases as more children are available.
* Excess odd steps in the source's tail push the debt up, which explains the enrichment in long orbits.

## Open

* A proof of the exponent: a quantitative THM-4581 together with control of the available time.
* The constant, and its dependence on K for general t.
* A certificate family for the doubly uncovered exponents other than descent.
