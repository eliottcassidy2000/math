# Pre-registration (kind-pasteur-2026-11-13-S27)

Written before the order-11 run below. An earlier launch of the same run was stopped unread: its output file was
empty when it was killed. Reserved identifiers: THM-5152, HYP-9679 (next free on origin/main at 3194e2073f).

## Object

Session 312 (CLUSTER_CHART.md) proves that a tournament is a 3-adic determinant orientation chi(det(v_i, v_j)),
v_i in Q_3^2, iff some vertex switch of it lies in the grammar G generated from a vertex by ordinal joins X => Y and
cyclic joins C3[X,Y,Z]. Call these the *chart class*. The question is which odd numbers are Hamiltonian-path counts
of chart-class tournaments.

## Already computed (reported, not predicted)

- G up to isomorphism, by canonical modular-decomposition trees (linear nodes flattened, C3 nodes up to rotation):
  1, 1, 2, 4, 9, 24, 65, 182, 528, 1597 classes at orders 1..10.
- Every vertex switch of every G member through order 10. The odd values below 400 not attained are 7, 21, 31, 35,
  53, 113, 125, 127, 149, 173, 181, 203, 215, 217, 235, 251, 263, 267, 277, 287, 305, 341, 343.
- Busch's bound gives min H over strong tournaments of order n = 3, 5, 9, 15, 25, 45, 75, 125, 225, 375 for
  n = 3..12 (THM-4744).
- The chart class is hereditary (Session 312), and H is multiplicative over ordinal sums.
- Hence a value below f(11) = 225 that is missed through order 10 and is not a product of attained values is missed
  at every order: 31, 35, 53, 113, 127, 149, 173, 181, 203, 217 (besides 7, 21).
- The bare grammar G, without switching, misses 13, 19, 29, ... (406 odd values below 1200 through order 13).

## Decisive run P1 (order 11, every switch, values below 3000)

Prediction:
- the permanent misses below f(12) = 375 (missed through order 11, and not a product of attained values) are the
  ten listed above together with at least one value in (225, 375) — candidates 251, 263, 277, 287, 341, 343;
- 31 and 127 stay missed (Mersenne primes);
- confidence about 0.65 for "at least one more".

## Correction appended after registration (before the P1 result was read)

- The order-10 count of G classes is **1,558**, not 1,597. The file g10.txt has 2,374 lines over orders 1-10, and
  orders 1-9 account for 816 of them. The read-only sweep regenerated 1,558 independently and matched Session 268.
- **Credit.** The spectrum of the unswitched grammar G is KNOWN. It is Session 268's composable-tournament census
  (session268_composable_tournaments_20261111/modular, n <= 16, 74,149 distinct H at order 16), which also names
  13 as outside the family. Only the switching closure, the chart class proper, is new here.

## Outcome

- **P1: CONFIRMED**, 2026-11-13. Order 11 has 4,674 classes and took 48 min 45 s.
- The permanent omissions below 375, besides 7, are 21, 31, 35, 53, 113, 127, 149, 173, 181, 203, 217, 251, 263,
  277, 287, 341 and 343.
- All six candidates in (225, 375) are omitted at every order, and 31 and 127 stay omitted.
