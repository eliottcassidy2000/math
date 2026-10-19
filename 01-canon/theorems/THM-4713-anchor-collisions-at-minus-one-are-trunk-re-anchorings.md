---
id: THM-4713
title: "Anchor collisions at -1 are trunk re-anchorings: a canonical head h of total A shares the numerator 2^a - 3 of the single-letter head (a+1) iff f_h(-1) = (4^k - 1)/3 - 2^(1-A) with a = A - 1 + 2k, i.e. h lands one unit of its own denominator below the k-th predecessor (4^k-1)/3 of 1; the k = 1 heads are exactly those whose numerator is a Mersenne number 2^(A-1) - 1; THM-4640's rigidity list {7, 13, 17, 21, 28, 31} is explained (five k = 1 heads, one k = 2 head below 5) and is complete to a <= 44; each gives the uniform switch (a+1, c) ~ (h, c+2k) of distance |h| - 1, and 125 = 2^7 - 3 is the only perfect power among 2^a - 3, 3 <= a <= 400"
status: >
  PROVED: statements 1-3 (a direct corollary of THM-4640 (1) and (3), the ladder relation specialised to a single-letter
  head; the proof is two lines of arithmetic). FINITE-EXACT: statement 4 (complete backward search for every a <= 44 and
  every k; the perfect-power check to a = 400). Not claimed: that these are new to the literature (THM-4640 cites
  Bohm-Sontacchi, Lagarias, Blecksmith-McCallum-Selfridge for the numerators), or anything about the Mersenne exponent
  classes beyond the displayed switch. Pre-audit.
session: mac-mini-2026-10-19-mersenne
source: 05-knowledge/results/mersenne_horizon_20261019_clusters_orphans_reanchoring.md
scripts:
  - 04-computation/experiments/mersenne_horizon_20261019/anchor_reanchoring.py (+ anchor_reanchoring_64.out if finished; the a <= 44 run is reproduced in the note)
related:
  - THM-4640 (the repair algebra at -1: numerators K(s), positivity, anchor rigidity list, ladders, Gersonides carries)
  - THM-4555 (uniform switches are collisions at -1; (vi): (8, c) ~ (4, 1, 1, c+2))
  - THM-4600 (anchors c_a = 1/(2^a - 3)), THM-4601 (ladder compiler at +1; (iii) the clearing identity with (4,1,1))
  - THM-4517 (trunk numbers (4^k-1)/3 as entry points of 1), THM-4556, THM-4712
---

# THM-4713 — anchor collisions at −1 are trunk re-anchorings

## Setting (THM-4640)

* `f_c(x) = (3x+1)/2^c`; a reduced head `h = (a_1, …, a_m)` has `a_1 ≥ 2`, total `A = Σ a_i`, and anchor image `f_h(−1) = N/2^(A−1)` with `N` odd (THM-4555 (ii)). Its **numerator** is `K(h) = 3N + 2^(A−1)`, the numerator of `3 f_h(−1) + 1` over `2^(A−1)`; two words `(h, c)`, `(h', c')` collide at `−1` iff `K(h) = K(h')` and `A + c = A' + c'` (THM-4640 (1)).
* **Canonical** heads have `a_1 ≥ 3`. The single-letter head `(a+1)` has `f = −2^(−a)`, `A = a + 1`, `N = −1`, `K = 2^a − 3`: the reciprocal anchor of THM-4600.
* `a_k = (4^k − 1)/3 = 1, 5, 21, 85, …` are the odd predecessors of 1 (`U(a_k) = 1`).

## Statements

1. **Trunk re-anchoring (PROVED).** Let `h ≠ (a+1)` be a reduced head with total `A` and odd numerator `N`. Then

       K(h) = 2^a − 3     iff     a − A + 1 = 2k with k ≥ 1   and   f_h(−1) = a_k − 2^(1−A).

   In words: `h` collides with the single letter exactly when it lands one unit of its own denominator `2^(A−1)` below a trunk number. For `k = 1` this reads `N = 2^(A−1) − 1`: the heads whose numerator is a Mersenne number are exactly the collisions with the anchor at `+1`.

2. **The switch (PROVED).** For such `h` and every `c ≥ 1`, `(a+1, c) ~ (h, c + 2k)` is a collision at `−1` (equal totals `a + 1 + c`), and by THM-4555 (iv) it is the uniform switch of distance `D = |h| − 1`: every odd `n` whose word starts with `1^r (a+1, c)`, `r ≥ D`, merges at equal odd-step time with `(n+1)/2^D − 1`, whose word starts with `1^(r−D) (h, c + 2k)`.

3. **Ladder form (PROVED).** Statement 1 is THM-4640 (3) with the short head `(a+1)`: `f_h(−1) = 4^k f_(a+1)(−1) + (4^k − 1)/3`.

4. **The complete list to `a ≤ 44` (FINITE-EXACT).** The canonical heads with numerator `2^a − 3` other than `(a+1)` are

   | `a` | `k` | head `h` | total `A` | `N` | distance `D` |
   |---|---|---|---|---|---|
   | 7 | 1 | (4, 1, 1) | 6 | 31 = 2^5 − 1 | 2 |
   | 13 | 1 | (5, 1, 2, 3, 1) | 12 | 2047 = 2^11 − 1 | 4 |
   | 17 | 1 | (4, 3, 2, 2, 1, 1, 3) | 16 | 2^15 − 1 | 6 |
   | 21 | 1 | (6, 1, 3, 1, 1, 1, 2, 2, 3) | 20 | 2^19 − 1 | 8 |
   | 28 | 2 | (5, 2, 3, 1, 1, 1, 4, 1, 1, 2, 1, 1, 2) | 25 | 5·2^24 − 1 | 12 |
   | 31 | 1 | (4, 1, 3, 6, 2, 3, 3, 2, 2, 3, 1) | 30 | 2^29 − 1 | 10 |

   and there are no others for `a ≤ 44`. This reproduces THM-4640 (1)'s rigidity list `{7, 13, 17, 21, 28, 31}` (`a ≤ 31`, one extra representation each) and explains it: five re-anchorings to `+1` and one, at `a = 28`, to `5`. The first row is THM-4555 (vi)'s `(8, c) ~ (4, 1, 1, c+2)`, the Catalan carry `9 = 8 + 1` of THM-4640 (5); `(4, 1, 1)` is also the block of THM-4601 (iii)'s clearing identity. Among `2^a − 3`, `3 ≤ a ≤ 400`, the only perfect power is `125 = 5^3` at `a = 7`.

## Proofs

**1.** `K(h) = 3N + 2^(A−1) = 2^a − 3` gives `3(N + 1) = 2^a − 2^(A−1) = 2^(A−1)(2^(a−A+1) − 1)`. The right side is divisible by 3 iff `a − A + 1` is even, say `2k`; `k ≥ 1` because `N + 1 > 0` (`N ≥ 1` for canonical `h` by THM-4640's positivity, and `N = −1` only for the single-letter head). Then `N + 1 = 2^(A−1) a_k`, i.e. `f_h(−1) = N/2^(A−1) = a_k − 2^(1−A)`. Conversely, `f_h(−1) = a_k − 2^(1−A)` gives `N = 2^(A−1) a_k − 1` and `K(h) = 3·2^(A−1) a_k − 3 + 2^(A−1) = 2^(A−1)(3a_k + 1) − 3 = 2^(A−1+2k) − 3`. ∎

**2.** Equal numerators and equal totals `a + 1 + c = A + (c + 2k)` make the two words collide (THM-4640 (1)); THM-4555 (iv) with `|u'| − |u| = (|h| + 1) − 2` gives the switch. ∎

**3.** `4^k·(−2^(−a)) + a_k = −2^(2k−a) + a_k = a_k − 2^(1−A)` since `a − 2k = A − 1`. ∎

**4.** For each `a ≤ 44` and each `k` with `A = a + 1 − 2k ≥ 2`, the heads with anchor image `(2^(A−1) a_k − 1)/2^(A−1)` are enumerated backward: the last letter `c` of a head with numerator `N` and total `A` leads to the numerator `(N − 2^(A−c−1))/3`, which must be an odd integer, or to the single-letter head when `N = 2^(A−2) − 3`; the recursion is memoised. Every listed head is checked against its exact anchor image, and `3N + 2^(A−1) = 2^a − 3` is verified. The perfect-power check uses `sympy.perfect_power`. ∎

## Remark

The six exceptional anchors `2^a − 3 ∈ {125, 8189, 131069, 2097149, 268435453, 2147483645}` are the only reciprocal anchors `1/c_a` of THM-4600 (`a ≤ 44`) reached by a second canonical path from `−1`. On the Mersenne line the single letter `a + 1` is the first post-run letter of `M_E` iff `E` is even with `v_2(E) = a − 2`, so these are the even-exponent classes `E ≡ 2^(a−2) (mod 2^(a−1))` that admit, beyond the reset switch, a direct switch of distance `|h| − 1` to `M_(E−|h|+1)` whenever the child's post-run word begins with `h`.
