# Pre-registered predictions (klein-2026-11-16-two3)

Written to `prediction_wheels.txt` at 2026-11-15 22:57 MST, before `brier_sheets.py` was run and before the base-2 table
beyond `N = 123` and the base-3 table beyond `N = 93` were read (the two `wheels2.py` runs were still in progress).

| # | prediction | verdict |
|---|---|---|
| P-W1 | for every `N ≤ 150` where the exact search decides, the exact answer equals the nested-dyadic criterion | **confirmed**: identical at all decided `N`; the five periods the first search left undecided (90, 100, 126, 132, 150) are settled in WHEELS.md §2 — three by the `a = 1` identity, 100 by the density bound, 132 by the two-arrangement capacity argument (it does not cover) |
| P-W2 | no period `N ≤ 150` admits a self-dual (two-sheet) cover by its even-order primes | **confirmed on every decided period** (24, 40, 48, 56, 60, 64, 72, 80, 84, 112, 128 fail); 96 and 140 fail; 120 undecided after `2·10^8` nodes |
| P-W3 | in base 3 the coverable periods `N ≤ 120` are the multiples of 24 that are `≥ 48`, plus 90 | **refuted by one period**: the coverable base-3 periods are 48, 72, 90, 96, 108, 120 — 108 = 4·27 also covers |

Not predicted (computed afterwards): the base-2 table itself (the naive guess "`12 | N`" was abandoned when 112 = 16·7
covered, before any prediction was written), Corollary C and the least periods 80 and 112, and the classical-set
identification.
