---
id: HYP-9336
title: "Gilbreath depth tracks the maximal prime gap: the first row k(N) of the absolute-difference triangle of the first N primes that is 1 followed only by 0s and 2s (the absorbing sea, after which the leading 1 persists) satisfies k(N) = (0.8-1.05) x G(p_N), where G(p_N) is the largest gap between consecutive primes below p_N; data k = 35, 65, 97, 147, 175, 239 against G = 34, 72, 114, 154, 222, 248 for N = 10^3, ..., 10^7, 1.1 x 10^7 (p_N < 2 x 10^8), and Odlyzko's 635 rows for p < 10^13 against G = 674; equivalently k(N) ~ c (log p_N)^2 with c about 0.5"
status: >
  NARROW BAND REFUTED; broader prime-specific correlation OPEN. The
  exact all-prefix scan below 2 x 10^8 gives k/G=69/112 at N=31569,
  outside the title's 0.8-1.05 band. The wider [0.5,1.5] falsification
  test survives for N>=1000 in this finite range only. No all-N bound
  or asymptotic is proved. The proposed isolated-gap disappearance
  mechanism is false for general difference automata; and later
  safe-sea entry requires a separate check of all earlier left edges.
  An explicit nonprime increasing prime-parity family has fixed G=6,
  every positive left edge 1, and k(N)=N-3. Thus any prime-specific
  depth bound must use more than maximal gap and parity.
session: mac-mini-2026-10-19-duality
source: 05-knowledge/results/duality_20261019_multiplication_addition_golden_circle_gilbreath.md
related:
  - THM-4511 (wall theorem for Gilbreath's automaton), HYP-9163 (size-6 extinction law), 05-knowledge/results/collatz_oscillation_gilbreath_20260926.md (the 0/2 sea is XOR; fronts move left at speed 1; row 65 for p < 2 x 10^5), gilbreath_fermat_platonic_20260926.md
  - 04-computation/experiments/duality_20261019/gilbreath_depth.py (+ .out)
---

# HYP-9336 — Gilbreath depth tracks the maximal prime gap

> **Current-truth correction, 2026-10-19.** The title's narrow `0.8–1.05` ratio band is false even within the original numerical range: at `N=31,569`, `p_N=370,609`, the exact values are `k=69`, `G=112`, so `k/G≈0.616071`. An all-prefix scan of the 11,078,937 primes below 200 million gives this as the minimum for N≥1000, while the maximum is 35/34. The wider `[0.5,1.5]` finite-range test survives; no asymptotic conclusion follows from that.
>
> **Mechanism correction.** Cancellation of the central descendant of one large gap does not show that all its descendants die. A single 6 in an infinite background of 2s produces a persistent `{0,4}` XOR tail. More strongly, the increasing nonprime sequence `2,3,5,9,15,17,19,…` has maximal gap 6, every positive leading difference 1, and exact sea-entry depth `k(N)=N−3` for all N≥5. Primality and consecutiveness are missing hypotheses in any attempted transfer from general gap words.
>
> **Certification correction.** `k(N)<N` alone only certifies the edge after sea entry. Earlier leading entries must also be checked. The nonconsecutive prime word `(2,3,7,13)` has edge `(2,1,3,1)` and sea entry at depth 3, showing the distinction. The new scan explicitly checks all earlier leading entries. It does not furnish a counterexample to the consecutive-prime conjecture.
>
> See [the corrected depth analysis and exact controls](../results/gilbreath_depth_corrections_20261019.md). The original title, sample table and heuristic are retained below as historical context, with their status governed by this notice.

## Setting

Row 0 is the sequence of primes; row `r+1` is the absolute difference of consecutive entries of row `r`. Gilbreath's conjecture: every row `r ≥ 1` starts with 1. Once a row is `1` followed only by `0`s and `2`s on an initial segment of length `L`, every later row is too (`|1 − 0| = |1 − 2| = 1`, and the tail evolves by XOR, Pascal mod 2), so the leading 1 persists for `L` more rows. The **depth** `k(N)` is the least `r` such that row `r` has this form on its first `N − r` entries; Gilbreath's conjecture for the first `N` primes follows from `k(N) < N`.

## Data (`gilbreath_depth.out`)

| `N` primes | `p_N` | depth `k(N)` | max gap `G(p_N)` | `k/G` | `(log p_N)²` |
|---|---|---|---|---|---|
| 10³ | 7,919 | 35 | 34 | 1.03 | 81 |
| 10⁴ | 104,729 | 65 | 72 | 0.90 | 134 |
| 10⁵ | 1,299,709 | 97 | 114 | 0.85 | 198 |
| 10⁶ | 15,485,863 | 147 | 154 | 0.95 | 274 |
| 10⁷ | 179,424,673 | 175 | 222 | 0.79 | 361 |
| 11,078,937 | 199,999,991 | 239 | 248 | 0.96 | 365 |

The rightmost defect of each row sits at the row's right end (the frontier moves left by exactly one position per row), so the depth is decided by the longest-lived interior defects.

## Why plausible, what is missing

* The entries of row 1 are the gaps; a large gap `g` between typical neighbours `a, b` produces `(g − a, g − b)` and then `|a − b|`: three rows. The persistence must come from clusters of unusual gaps; the automaton of THM-4511 (size-4 walls that never cross a sea 2) suggests that a defect's lifetime is governed by the distance to the next wall, i.e. by the local gap statistics, whose extreme is the maximal gap.
* The `(log p)²` reading is Cramér's model for `G`; the data cannot separate `k ≈ G` from `k ≈ 0.5 (log p)²`.
* Odlyzko's 635 rows for `p < 10^13` (cited from memory) against `G = 674` fits the `0.8–1.05` band.
