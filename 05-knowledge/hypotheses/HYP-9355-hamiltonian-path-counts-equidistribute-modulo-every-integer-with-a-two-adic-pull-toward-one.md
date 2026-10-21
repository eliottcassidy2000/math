---
id: HYP-9355
title: "The Hamiltonian-path count of a uniformly random labelled n-tournament equidistributes modulo every integer as n -> infinity (over all residues for odd moduli, over the odd residues for powers of 2; H is always odd by Redei), and at finite n the leading 2-adic deviation is a radial pull toward 1: the frequency of H = r (mod 2^k) decreases with the 2-adic distance |r - 1|_2"
status: >
  OPEN. FINITE-EXACT for n <= 9, labelled counts from the class census (h_patterns.py (Q4), (Q4b)). For p = 3, 5, 7, 11, 13 the fraction with p | H is within 1.6% of 1/p at n = 9.
  At n = 7, 8, 9 the relative mass of H = 1 (mod 8) is 1.179, 1.064 and 1.015; the mod-8 masses at n = 9 order
  1 > 5 > 7 > 3. Mod 16 at n = 7 the masses fall from residues 1 and 9 (1.205, 1.154) to 7 and 15 (0.824, 0.790);
  the residues at 2-adic distance 1/2 from 1 lie in between but are not monotone (h_patterns.py (Q4b)).
  The sampling predictions below are committed BEFORE the sampler is run.
source: opus-2026-10-21-S50, 05-knowledge/results/h_spectrum_multiplicities_20261021.md (to be written)
related:
  - 01-canon/theorems/THM-466-the-2-adic-digits-of-H-are-odd-cycle-collection-counts.md (H = sum alpha_k 2^k mod 2^m)
  - 01-canon/theorems/THM-478-digit-one-degree-law.md, THM-1980-the-two-adic-edge-of-the-redei-count.md (2-adic theory of H)
  - 05-knowledge/results/actual_hamilton_multiplicities_20261021/session24-multiplicities/PROOF.md (concurrent multiplicities)
scripts:
  - 04-computation/experiments/h_spectrum_multiplicities_20261021/h_patterns.py (+ .out)
  - 04-computation/experiments/h_spectrum_multiplicities_20261021/h_residue_sampling.py (to be written and run after this commit)
---

# HYP-9355 — residues of the Hamiltonian-path count

## Statement

* Let `T` be a uniformly random labelled tournament on `n` vertices, with each arc an independent fair coin, and let `H = H(T)`.
* **Equidistribution.** For every odd `M` and every residue `r`, `P(H ≡ r mod M) → 1/M`. For `M = 2^k`, `P(H ≡ r mod 2^k) → 2^(1−k)` for every odd `r`.
* **Two-adic pull toward 1.** At finite `n`, `P(H ≡ r mod 2^k)` is larger for `r` 2-adically closer to 1, and the excess at `r = 1` decays as `n` grows.
  * By THM-466, `H ≡ 1 + 2α_1 + 4α_2 + … (mod 2^k)`, where `α_j` counts sets of `j` disjoint odd cycles. So the pull says that `α_1 + 2α_2 + …` is 2-adically biased toward 0 at small `n`.

## Evidence (FINITE-EXACT, labelled, n ≤ 9)

| `p` | `n = 7` | `n = 8` | `n = 9` | `1/p` |
|---|---|---|---|---|
| 3 | 0.3147 | 0.3297 | 0.3388 | 0.3333 |
| 5 | 0.1681 | 0.1806 | 0.2005 | 0.2 |
| 7 | 0.1086 | 0.1377 | 0.1420 | 0.1429 |
| 11 | 0.0884 | 0.0882 | 0.0919 | 0.0909 |
| 13 | 0.0985 | 0.0756 | 0.0764 | 0.0769 |

* For `n ≤ 6` no tournament has `7 | H`. The fraction for 7 jumps from 0 to 0.1086 at `n = 7`; for 7, the approach to `1/p` starts late.
* Relative masses mod 8 (each divided by 1/4), for residues 1, 3, 5, 7:

  | `n` | 1 | 3 | 5 | 7 |
  |---|---|---|---|---|
  | 7 | 1.179 | 0.985 | 1.029 | 0.807 |
  | 8 | 1.064 | 0.979 | 1.019 | 0.938 |
  | 9 | 1.015 | 0.992 | 1.000 | 0.993 |

## Pre-registered sampling predictions

The sampler draws `4·10^6` labelled tournaments at `n = 10`, `2·10^6` at `n = 11` and `5·10^5` at `n = 12`, with fixed seeds. It counts paths exactly by subset DP. Error bars are binomial.
* **P1.** For `p ∈ {3, 5, 7, 11, 13}` and `n ∈ {10, 11, 12}`: `|p·freq(p | H) − 1| < 0.02`.
* **P2.** For `n ∈ {10, 11, 12}`: residue 1 is the most frequent odd residue mod 8.
* **P3.** At `n = 10`: `1 < 4·freq(H ≡ 1 mod 8) < 1.0153`, the latter being the `n = 9` value. The bias persists but is smaller.
* **P4 (informative).** The excess `4·freq(H ≡ 1 mod 8) − 1` decreases from `n = 10` to `n = 12`, within errors.

**Refutation.** P1, P2 or P3 fails by more than 3 standard errors. If P4 fails, the decay claim is weakened.
