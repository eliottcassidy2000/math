---
id: THM-4735
title: "Debt cones covering 63: valuation-dependent last-step cones, optimal terminal splits, and the debt-two obstruction to universal cone coverage. In the affine ledger of an odd Collatz word w = (a_1..a_m) (z = n-1, z_end = lambda z + P - D, charge_palettes_boundary_20261020g/charge.md): (A) the payment theorem holds for every cone constant c <= c_max(a) := (2^a + 4)/(2^a - 3) when the final actual valuation is a >= 3 (and c < 8 for a = 2): a prefix u with D_u <= c (lambda_u - 1) followed by a descending step of valuation a has D_total < 2 <= n - h, so positive charge pays; the sharpened cone covers source 63 (prefix ratio 1.0573 <= 20/13 at its final valuation 4) and 48,741 of the 50,000 first exits below 100,001 (487,405 of 500,000 below 1,000,001), against 47,424 (473,943) for the c = 1 cone. (B) Optimal terminal splits: with a terminal packet v of multiplier lambda_v < 1 the certificate value is V(u,v) = (D_u/(lambda_u - 1))(1 - lambda_v) + D_v, and the best split covers 49,997 / 50,000 (499,958 / 500,000). (C) Obstruction: V(u,v) >= D_w for every split (because lambda_u lambda_v < 1), so a first-exit word with D_w >= 2 admits no cone/terminal-packet certificate at the uniform threshold 2, whatever c and split: there are exactly 3 such sources below 100,001 (4575, 42703, 65567; debts 2.16, 2.53, 2.13; drops 110, 482, 8722) and 42 below 1,000,001; universal coverage by threshold-two cones is FALSE. (D) Every one of the 500,000 first exits pays directly (D_w < n - h); the first-exit payment conjecture is recorded in HYP-9349"
status: >
  PROVED (pre-audit): (A) the sharpened payment theorem (the packet's proof with the constant tracked), (C) the
  inequality V(u,v) >= D_w and the resulting obstruction. FINITE-EXACT (cone_splits.out, cone_splits_1e6.out, exact
  rational arithmetic): all counts above; the packet's 47,424 / 2,576 split reproduced exactly (its "last-step cone"
  counts single-step words as trivially in the cone); 63's word, ratios and best split displayed; the 42 debt-two
  sources listed with words. CITED: charge_palettes_boundary_20261020g/charge.md (ledger, cone closure under
  concatenation, payment theorem, first-reset clock, the (1,1,1,1,3) family, source 63 as the open hostile),
  first_exit_complement_games_20261020 and adaptive_budgets_smooth_interfaces_20261020 (the decomposition).
  NOT PROVED: that every first exit pays (D_w < n - h) for all sources; any statement about sources without a first
  exit (Collatz); any ROOT-current consequence.
source: mac-mini-2026-10-20-debtcones (owner's prompt: debt cones covering 63)
depends_on:
  - 05-knowledge/results/charge_palettes_boundary_20261020g/charge.md
related:
  - 05-knowledge/results/first_exit_complement_games_20261020/README.md
  - 05-knowledge/hypotheses/HYP-9349-first-exit-payment-annulus-increments-and-the-idoneal-dictionary-of-105.md
scripts:
  - 04-computation/experiments/debt_cones_20261020/cone_splits.py
  - 04-computation/experiments/debt_cones_20261020/cone_obstruction_check.py
outputs:
  - 04-computation/experiments/debt_cones_20261020/cone_obstruction_check_1e6.out
  - 04-computation/experiments/debt_cones_20261020/cone_splits.out
  - 04-computation/experiments/debt_cones_20261020/cone_splits_1e6.out
  - 04-computation/experiments/debt_cones_20261020/cone_splits.json
---

# THM-4735 -- which cones cover 63, and which first exits no cone can cover

**PROVED (pre-audit) + FINITE-EXACT** as itemized in the status.

## 0. The ledger (charge.md, retained verbatim)

For an odd source `n > 1` let `w = (a_1, ..., a_m)` be the valuations `a_i = v_2(3x_i + 1)`
of the odd steps until the first odd iterate `h < n` (the first exit). With `z = n - 1`,

```text
z_end = h - 1 = lambda z + P - D,   lambda = 3^m / 2^K,  K = sum a_i,
P = sum_{a_i = 1} lambda_(i+1:m),    D = sum_{a_i >= 3} (1 - 4/2^{a_i}) lambda_(i+1:m),
```

(`lambda_(i+1:m)` the product of the later multipliers `3/2^{a_j}`). Concatenation:
`lambda_uv = lambda_u lambda_v`, `D_uv = lambda_v D_u + D_v`. The **cone** `0 <= D <= c (lambda - 1)`,
`lambda >= 1`, is closed under concatenation. **Positive charge pays** when
`h - 1 + D_w < n - 1`, i.e. `D_w < n - h`; since `n - h` is a positive even integer, `D_w < 2`
suffices. The packet's payment theorem: prefix `u` in the cone with `c <= 1`, then one actual
descending step of valuation `a`, gives `D_total < 2`. It left source 63 open: 63's first exit
(34 odd steps, 5 resets, final valuation 4, `h = 61`) pays, but its prefix has `D_u > lambda_u - 1`.

## 1. Theorem A: the admissible cone constant depends on the final valuation

Let `u` lie in the cone with constant `c` and let the final step have valuation `a >= 2`,
`mu = 3/2^a`, descending (`h < n`). The positive intercept gives `mu lambda_u < 1`, hence

```text
D_total = mu D_u + D_a <= c mu (lambda_u - 1) + D_a < c (1 - mu) + D_a,    D_a = [a >= 3] (1 - 4/2^a).
```

So `D_total < 2` whenever `c (1 - mu) + D_a <= 2`, i.e.

```text
c <= c_max(a) := (2 - D_a)/(1 - mu) = (2^a + 4)/(2^a - 3)   (a >= 3),      c < 8   (a = 2):
c_max(3) = 12/5 = 2.4,  c_max(4) = 20/13 = 1.538...,  c_max(5) = 36/29,  c_max(a) -> 1.
```

For source 63 the prefix of 33 steps has `D_u/(lambda_u - 1) = 1.0573...` and the final
valuation is `4`: `1.0573 <= 20/13`, so the valuation-4 cone covers 63. The same family covers
every first exit whose prefix ratio is at most `c_max(a_m)`:

```text
sources 3..100001 (50,000):   c = 1 cone 47,424 (= packet)   valuation cones 48,741   unresolved 1,259
sources 3..1000001 (500,000): c = 1 cone 473,943             valuation cones 487,405  unresolved 12,595
```

(25,000 / 250,000 of these are single-step words `(a)`, in every cone with the empty prefix.)
The first sources the valuation cones miss are `47, 79, 191, 223, 251, 283, 319, 463, 487, ...`.

## 2. Theorem B: optimal terminal splits

The packet's terminal-packet form (prefix `u` in the cone, packet `v` with `lambda_v < 1`)
pays whenever `c (1 - lambda_v) + D_v <= 2`. Taking the smallest admissible `c`,
`c = D_u/(lambda_u - 1)`, the certificate value of a split `w = u v` is

```text
V(u, v) = (D_u/(lambda_u - 1)) (1 - lambda_v) + D_v        (V(empty, w) = D_w),
```

and `w` is covered iff some split has `V <= 2` (strictly: `D` is non-integral when a reset
exists). Minimising over splits covers `49,997 / 50,000` and `499,958 / 500,000`; the
histogram of the minimal value is concentrated at `0` (no reset), `0.5` (one reset of
valuation 3 at the end), `0.8-1.0`, with a thin tail to `3.125` (source 437,247).

## 3. Theorem C: the debt-two obstruction

For every split, `V(u, v) >= D_w`, with equality only when the prefix carries no debt.
*Proof.* `D_w = lambda_v D_u + D_v`. For a first exit `lambda_w = lambda_u lambda_v < 1`
(if `lambda_w >= 1` the exit `lambda z + P - D < z` would need `D > P + (lambda - 1) z`;
exact check: no first-exit word below `1,000,001` has `lambda_w >= 1`, cone_obstruction_check.out
and cone_obstruction_check_1e6.out). Hence `lambda_v < 1/lambda_u`,
`1 - lambda_v > (lambda_u - 1)/lambda_u`, and `(D_u/(lambda_u - 1))(1 - lambda_v) >= D_u/lambda_u
>= lambda_v D_u`, strictly when `D_u > 0`. The check also confirms `V >= D_w` on all 123,470
admissible splits below `100,001` and all 1,241,405 below `1,000,001` (minimum gap `0`, at
debt-free prefixes) and the identity `pays <=> P_w < (1 - lambda_w) z` on every source. `#`

Hence a first exit with `D_w >= 2` has **no** cone / terminal-packet certificate at the
uniform threshold `2`, for any cone constant and any split; it pays only through the actual
drop `n - h`. Exactly three sources below `100,001` are of this kind,

```text
4575  : D = 2.1621, n - h = 110,  m = 29, 5 resets
42703 : D = 2.5269, n - h = 482,  m = 41, 7 resets
65567 : D = 2.1327, n - h = 8722, m = 27, 4 resets
```

and 42 below `1,000,001` (largest debt `3.125`, source `437,247`, `92` odd steps). So the
packet's "universal coverage by such controllers" is false as stated; the correct universal
target is the direct payment `D_w < n - h`, which holds for every first exit tested
(all 500,000) and is recorded as a conjecture in HYP-9349 together with what a proof needs
(`P_w < (1 - lambda_w) z`: the descent is never caused by the debt alone).

## 4. The cone covering 63, concretely

```text
63 -> 95 -> 143 -> 215 -> 323 -> 485 -> 91 -> 137 -> 103 -> 155 -> 233 -> 175 -> 263 -> 395 -> 593
   -> 445 -> 167 -> 251 -> 377 -> 283 -> 425 -> 319 -> 479 -> 719 -> 1079 -> 1619 -> 2429 -> 911
   -> 1367 -> 2051 -> 3077 -> 577 -> 433 -> 325 -> 61
word (1,1,1,1,1,4,1,2,1,1,2,1,1,1,2,3,1,1,2,1,2,1,1,1,1,1,3,1,1,1,4,2,2,4),  lambda = 0.9258, D = 1.5305,
prefix ratio 1.0573 (fails c = 1, passes c_max(4) = 20/13);  best split value 1.5305 (the whole word as packet).
```

The minimal split value for 63 is attained by the empty prefix: 63 is covered most cheaply by
the "terminal packet = whole word" certificate `D_w = 1.53 < 2`, and the valuation-4 cone is
the smallest *prefix* cone covering it.

## 5. Files

`04-computation/experiments/debt_cones_20261020/cone_splits.py` (exact `Fraction` arithmetic;
`cone_splits.out` to `100,001`, `cone_splits_1e6.out` to `1,000,001`, `cone_splits.json`).
