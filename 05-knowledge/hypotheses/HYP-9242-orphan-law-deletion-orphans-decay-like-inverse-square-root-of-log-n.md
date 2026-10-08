---
id: HYP-9242
title: "Orphan law: the proportion of residual first-reset-2 sources n with no deletion partner (no child (n+1)/2^D - 1 merging at equal time before 1) tends to 0 like a negative power of log n (upper bound (log n)^(-1/2) up to logs at sketch level, from the child D = 1 alone); the effective exponent at log2 n <= 3200 grows with the number of children (about 0.42 for K = 2, 3; 0.54-0.63 for K = 9-33; on the Mersenne line 0.59, 95% CI [0.42, 0.75], with fraction*sqrt(K_mid) = 1.09, 1.01, 1.14, 0.89 on [800,1600), ..., [6400,12800]); orphans' post-run orbits are about one standard deviation longer than random orbits of the same size, and they are rescued by shallow 3-adic predecessors at the population rate (about 40%)"
status: >
  OPEN. UPPER BOUND at the sketch level of THM-4581 (4), for t uniform among the odd L-bit integers in the residual class:
  P(orphan) <= C L^(-1/2) (log L)^2 (see 'Partial result'). NUMERICAL: Mersenne orbits computed exactly for every K <= 12800 (odd-step
  count and Terras length to 1; reproduced byte-for-byte by audit C). The partner criterion (smaller exponent with equal odd count and
  Terras difference K - K') is, for orbits that reach 1, exactly an equal-time zero-debt merge at a value >= 8 for every K: aligned orbits
  that first reach 1 together agree from 8 on (checked directly by audit C on 173 partner and 155 non-partner pairs, K <= 12800, and
  exhaustively for the orphans 1039, 1081, 1113). General residual sources n = 2^K t - 1 (session: K = 9, 17, 300 per row; audit C:
  K = 2..33, 2000 per row, B <= 3200). HEURISTIC: the pair chain's debt is a zero-drift walk with non-absorption tail T^(-1/2)
  (THM-4581), and the time available before the child reaches 1 is proportional to log n. Exponent and long-orbit readings corrected
  after audit C (MISTAKE-587).
source: mac-mini-2026-10-07-twoanchor (continuation), 05-knowledge/results/runcompress_orphans_cayley_20261007.md (section 3)
related:
  - 01-canon/theorems/THM-4556-the-mersenne-line-is-a-chain-of-debt-states-odd-shift-distance-2-adic-periodicity.md ((vi): 50 of 2500 odd a in [10^3, 6*10^3])
  - 01-canon/theorems/THM-4581-haar-coalescence-affinely-related-collatz-orbits-merge-almost-surely.md (the T^(-1/2) rate)
  - 01-canon/theorems/THM-4601-two-anchor-reduction-of-the-residual-first-reset-two-branch.md, THM-4603 (drift)
scripts:
  - 04-computation/experiments/runcompress_20261007/ (mersenne_sigma*.py, orphan_law_mersenne.py (+ .out), orphan_law_general.py (+ .out),
    mersenne_partners_by_shell.py, orphan_structure.py, orphan_time.py, orphan_rescue.py; data mersenne_sigma_12800.txt)
  - 04-computation/experiments/runcompress_20261007/audit_C/ (c4_*.py: orbits, partner merges, null model, exponent fits, rescue; c5_partial_result.py)
  - 04-computation/experiments/depthspec_20261007/ (depth spectrum of the D = 1, 3 certificates)
---

# HYP-9242 — the orphan law

## Statement

* A residual source n (first reset letter 2) is a **deletion orphan** if no child `h_D = (n+1)/2^D − 1`, `1 ≤ D ≤ K−1`, merges with it at equal time before 1.
* For orbits that reach 1 this is equivalent to: no D has `o(h_D) = o(n)` and `σ_T(h_D) = σ_T(n) − D`. The reason is that aligned orbits that first reach 1 together already agree from 8 on.
* **Conjecture.** Among residual sources of size about 2^L, the orphan proportion is ≍ L^(−1/2) up to logarithms when small K dominate. For fixed larger K and on the Mersenne line the decay may be faster:
  * the effective exponents for K = 9–33 at L ≤ 3200 are 0.54–0.63, which excludes 1/2;
  * on the Mersenne line the orphan fraction is about K^(−α) with α̂ = 0.59 (95% CI [0.42, 0.75]). 1/√K fits as well as K^(−0.7) or log K / K^0.7.

## Partial result: the upper bound for random t

**Claim** (at the sketch level of THM-4581 (4)). Fix K, and let t be uniform among the odd L-bit integers in the residual class. Then

    P(n = 2^K t − 1 is a deletion orphan) ≤ C L^(−1/2) (log L)^2.

**Reasons.**
* Every Terras step at most halves, so the child h_1 cannot reach 1 before Terras time log2 h_1 ≥ L − 1.
* An absorption of the D = 1 chain by chain time L − K − 2 is therefore an equal-time merge at a value ≥ 2^(K+2): a deletion certificate.
* The chain's driving bits are a forced prefix, (0, 1) or (0, 0, 1), followed by exactly uniform bits, even conditional on the two-run data (J, r) (audit C, exhaustive for L = 16).
* The post-run state is (J, 1) with P(J ≥ m) = 4^(1−m).
* THM-4581 (4) gives non-absorption ≤ C T^(−1/2) (log T)^2, with a constant at most linear in the starting debt (its sketch). Averaging over J gives the claim. The D = 1 chain's law does not depend on K, so C can be taken uniform for K ≤ L/2.
* Numerically, P(not absorbed by L−K−2)·√L = 7.1 → 14.3 for L = 100 → 1600: pre-asymptotic, inside the (log L)^2 envelope.

**Open.**
* The matching lower bound.
* Whether the exponent with many children exceeds 1/2.
* The Mersenne line, whose post-run bits are equidistributed in K only to depth about log2 K (THM-4601 (vi)).

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

* There are 135 orphans among odd K ≤ 12800.
* Over the last three doublings the local exponent is 0.60, with a 95% interval of about [0.28, 0.91]. A fit over odd K in [200, 12800] gives 0.59 [0.42, 0.75]. The data do not distinguish 1/2 from 0.6–0.7 or from log-corrected laws (audit C).

**General residual sources** (audit C, 2000 per row, B ≤ 3200). The fitted exponents are:
* about 0.42 for K = 2, 3;
* 0.54 for K = 9, 0.605 for K = 17, 0.62 for K = 33.

For K ≥ 9, 1/2 is rejected. For K = 2 the orphan event is exactly the failure of the single D = 1 chain (exponent 1/2 by THM-4581), yet it measures 0.42. So effective exponents at these sizes are biased downward. The √(K+B)-scaled fractions are not constant: 3.0 → 2.5 for K = 9 and 2.7 → 1.7 for K = 17.

**Long orbits.**
* Partners share σ_T − K exactly. So non-orphans inherit their post-run length from a smaller exponent, while an orphan's is fresh.
* The tercile rates 0.31%, 0.92%, 3.0% (odd K in [1000, 12800]) are therefore not evidence. A null model with the observed class structure and random-orbit lengths gives 0.25–0.56%, 0.66–1.22%, 2.75–3.20%.
* The genuine signal: measured against random orbits of the same size, orphans' post-run lengths are +1.0 standard deviations long on average (median +0.94, 86% positive, n = 83), and non-orphans' −0.01. This is consistent with a positive aggregate debt drift.
* It is not a single long ones-run: orphans' longest post-run ones-runs match the non-orphans' (Mann–Whitney z = −0.38).

**Other checks.**
* *Leading digits.* For odd K in [600, 6400] the orphan rate is 1.9–2.3% in each quarter of frac(K log2 3), and 1.6–2.2% on [1000, 6400]; the differences are within Poisson noise.
* *Nearest partner.* It is K−1 (D = 1) for 92% of non-orphans with K ≥ 1000. In each shell v2(K−1), the reset pair succeeds in 82–91% of cases and h_3/h_4 in 82–88% (K in [600, 6400]).
* *3-adic rescue.* A smaller node in the Terras backward tree within depth 36 certifies n. This holds for 20 of the 53 orphans in [1000, 6400] (0.377), against 0.403 for all 2700 odd K (hypergeometric p = 0.41). Depth-3 rescue happens exactly when K ≡ 5 mod 6.
* *Doubly uncovered.* 33 of the 2700 odd K (1.2%) have neither a deletion partner nor a shallow predecessor.

## Why

* The deletion certificate is the absorption of the pair chain's debt walk, which has zero drift and a T^(−1/2) non-absorption tail (THM-4581).
* The available time is about 4.8 times the child's bit length. This is the single-chain prediction.
* With many children the observed decay is faster (effective exponent ≈ 0.6), plausibly because several weakly correlated debt walks must all fail.

## Sheet diagnostic (2026-10-08, mac-mini-2026-10-08-reframes; NUMERICAL)

* Atlas P1 requires 2-adic statements to be sheet-blind. The conjugate of the Mersenne line under 3x − 1 is `P_K = 2^K + 1`.
  * Its deletion children are `P_(K−D)`, and its run ends are `2·3^(K−1) + 1`.
  * Its orbits end in the three positive 3x − 1 cycles, with shares 1037/1096/866 for `K ≤ 3000`.
  * Its residual (first-letter-2) class is the even `K`; the odd `K` always merge with `K − 1` three steps after the run.
* In the residual class, `K ≤ 3000`, the orphan counts are 75 (3x − 1) against 77 (3x + 1). The window fractions are:

  | window | 3x − 1 | 3x + 1 |
  |---|---|---|
  | [100, 400) | 0.107 | 0.093 |
  | [400, 1600) | 0.040 | 0.043 |
  | [1600, 3000] | 0.020 | 0.023 |

* So the orphan law is sheet-blind within noise. Script: `04-computation/experiments/reframes_20261007/sheet_diagnostic.py` (+ `.out`).
