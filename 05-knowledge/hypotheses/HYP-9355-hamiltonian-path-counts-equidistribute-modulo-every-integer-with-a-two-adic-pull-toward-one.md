---
id: HYP-9355
title: "The Hamiltonian-path count of a uniformly random labelled n-tournament equidistributes modulo every integer as n -> infinity (over all residues for odd moduli, over the odd residues for powers of 2; H is odd by Redei). The small-n excess of H = 1 (mod 8), clear for n <= 10, is not a persistent law: a pooled 10^7-sample replication at n = 11 shows no excess"
status: >
  OPEN (equidistribution). REFUTED as a persistent law: the 'two-adic pull toward 1' of the original title, after the
  S50 audit's replication (see the results section). FINITE-EXACT for n <= 9 (labelled counts from the class census,
  h_patterns.py (Q4), (Q4b)): for p = 3, 5, 7, 11, 13 the fraction with p | H is within 1.64% of 1/p at n = 9; the
  relative mass of H = 1 (mod 8) is 1.179, 1.064, 1.015 at n = 7, 8, 9. Mod 16 the masses are monotone in the 2-adic
  distance |r - 1|_2 only at n = 9 (at n = 7, residue 11 at distance 1/2 has 1.102, above residue 5 at distance 1/4 with
  1.000; at n = 8, residue 13 at 1/4 has 1.034, above residue 9 at 1/8 with 1.027).
  Pre-registered sampling (fad6a40581, before the sampler existed; seeds not pre-registered):
  - P1 PASSED (largest deviation 0.79%, p = 7 at n = 12);
  - P2 NOT CONFIRMED at n = 11: residue 5 leads in the session's two runs and in the audit's two independent runs;
    pooled over 10^7 samples, residue 1's excess is -0.0004 +- 0.0005;
  - P3 PASSED (n = 10: 1.0030 +- 0.0009, about 3.4 s.e. above 1);
  - P4 not monotone.
  Related, not overlapping: Session25 proves that every parity-admissible residue modulo every M is ATTAINED by a
  strong tournament (session25_spectrum_wiring_20261021); this hypothesis concerns the DISTRIBUTION.
source: opus-2026-10-21-S50, 05-knowledge/results/h_spectrum_multiplicities_20261021.md (section 3)
related:
  - 01-canon/theorems/THM-466-the-2-adic-digits-of-H-are-odd-cycle-collection-counts.md (H = sum alpha_k 2^k mod 2^m)
  - 01-canon/theorems/THM-478-digit-one-degree-law.md, THM-1980-the-two-adic-edge-of-the-redei-count.md (2-adic theory of H)
  - 05-knowledge/results/session25_spectrum_wiring_20261021/SYNTHESIS.md (residue attainment)
scripts:
  - 04-computation/experiments/h_spectrum_multiplicities_20261021/h_patterns.py (+ .out)
  - 04-computation/experiments/h_spectrum_multiplicities_20261021/h_residue_sampling.py (+ .out; the first attempt crashed at n = 12, its partial output kept as h_residue_sampling_run1_crashed_at_n12.out)
  - 04-computation/experiments/h_spectrum_multiplicities_20261021/audit_S50/a7_residue_replication.py, a7b, a7c (+ .out; independent replication)
---

# HYP-9355 — residues of the Hamiltonian-path count

## Statement

* Let `T` be a uniformly random labelled tournament on `n` vertices, with each arc an independent fair coin, and let `H = H(T)`.
* **Conjecture (equidistribution).** For every odd `M` and every residue `r`, `P(H ≡ r mod M) → 1/M`. For `M = 2^k`, `P(H ≡ r mod 2^k) → 2^(1−k)` for every odd `r`.
* **Withdrawn (the 2-adic pull).** The original version also conjectured that, at every finite `n`, residues 2-adically closer to 1 are more frequent. The exact data contradict monotonicity at `n = 7` and 8, and the excess at residue 1 is absent at `n = 11` (below). What survives is a finite statement: at `n ≤ 10`, `H ≡ 1 (mod 8)` is over-represented.

## Evidence (FINITE-EXACT, labelled, n ≤ 9)

| `p` | `n = 7` | `n = 8` | `n = 9` | `1/p` |
|---|---|---|---|---|
| 3 | 0.3147 | 0.3297 | 0.3388 | 0.3333 |
| 5 | 0.1681 | 0.1806 | 0.2005 | 0.2 |
| 7 | 0.1086 | 0.1377 | 0.1420 | 0.1429 |
| 11 | 0.0884 | 0.0882 | 0.0919 | 0.0909 |
| 13 | 0.0985 | 0.0756 | 0.0764 | 0.0769 |

* For `n ≤ 6` no tournament has `7 | H`.
* Relative masses mod 8, residues 1, 3, 5, 7:

  | `n` | 1 | 3 | 5 | 7 |
  |---|---|---|---|---|
  | 7 | 1.179 | 0.985 | 1.029 | 0.807 |
  | 8 | 1.064 | 0.979 | 1.019 | 0.938 |
  | 9 | 1.015 | 0.992 | 1.000 | 0.993 |

## Results of the pre-registered test, and the audit's replication

| | `n = 10` | `n = 11` | `n = 12` |
|---|---|---|---|
| `p·freq(p | H)`, `p = 3` | 1.0034 ± 0.0007 | 1.0010 ± 0.0010 | 0.9957 ± 0.0020 |
| `p = 7` | 0.9992 ± 0.0012 | 1.0002 ± 0.0017 | 0.9921 ± 0.0035 |
| mod 8 (session), residues 1/3/5/7 | 1.0030 / 0.9961 / 1.0020 / 0.9989 | 0.9982 / 1.0003 / 1.0010 / 1.0005 | 1.0032 / 0.9990 / 1.0002 / 0.9977 |
| mod 8, pooled with the audit | — | 0.9996 / 0.9991 / 1.0019 / 0.9994 (`N = 10^7`) | 1.0027 / 0.9991 / 0.9999 / 0.9982 (`N = 1.5·10^6`) |

* **P1 passed.** The largest deviation is 0.79%.
* **P2 not confirmed.**
  * At `n = 11`, residue 5 leads in all four runs: the session's two, and the audit's own two with an independent generator and DP.
  * One audit run fails P2 by the session's own criterion.
  * Pooled, residue 5 minus residue 1 is +0.0023 (2.6 s.e., with the multinomial covariance), and residue 1's excess is −0.0004 ± 0.0005.
* **P3 passed.**
* **P4 was not monotone.**
* **The excess of `3 | H`** decays: 1.016, 1.0034, 1.0010, 0.996 at `n = 9, 10, 11, 12`.

**Verdict.** Equidistribution is supported through `n = 12`. The 2-adic pull is a small-`n` effect: resolved through `n = 10`, absent at `n = 11`.
