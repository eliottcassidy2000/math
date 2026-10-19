---
id: HYP-9336
title: "Gilbreath depth tracks the maximal prime gap: the first row k(N) of the absolute-difference triangle of the first N primes that is 1 followed only by 0s and 2s (the absorbing sea, after which the leading 1 persists) satisfies k(N) = (0.8-1.05) x G(p_N), where G(p_N) is the largest gap between consecutive primes below p_N; data k = 35, 65, 97, 147, 175, 239 against G = 34, 72, 114, 154, 222, 248 for N = 10^3, ..., 10^7, 1.1 x 10^7 (p_N < 2 x 10^8), and Odlyzko's 635 rows for p < 10^13 against G = 674; equivalently k(N) ~ c (log p_N)^2 with c about 0.5"
status: >
  OPEN. NUMERICAL: the exact depths for p_N < 2 x 10^8 (gilbreath_depth.py, numpy, 11,078,937 primes, 239 rows); the ratio
  k/G = 1.03, 0.90, 0.85, 0.95, 0.79, 0.96. CITED (from memory, to be checked against Odlyzko, Math. Comp. 61 (1993)): all rows
  from about row 635 on are 0/2 after the leading 1 for the primes below 10^13; the largest gap below 10^13 is 674 (first after
  7,177,162,611,713; T5K maximal-gap table), ratio 0.94. The repo's S9 note recorded k = 65 for p < 200,000 (reproduced here at
  N = 10^4, p_N = 104,729). No mechanism is proved: an isolated large gap dies in three rows, so the long-lived structures are
  not single gaps; the S9/THM-4511 automaton (fronts moving left at speed 1, size-4 walls) is the place to look. What would
  refute: k(N)/G(p_N) leaving [0.5, 1.5] for larger N; what would prove: a bound k <= G + O(1) from the automaton.
session: mac-mini-2026-10-19-duality
source: 05-knowledge/results/duality_20261019_multiplication_addition_golden_circle_gilbreath.md
related:
  - THM-4511 (wall theorem for Gilbreath's automaton), HYP-9163 (size-6 extinction law), 05-knowledge/results/collatz_oscillation_gilbreath_20260926.md (the 0/2 sea is XOR; fronts move left at speed 1; row 65 for p < 2 x 10^5), gilbreath_fermat_platonic_20260926.md
  - 04-computation/experiments/duality_20261019/gilbreath_depth.py (+ .out)
---

# HYP-9336 — Gilbreath depth tracks the maximal prime gap

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
