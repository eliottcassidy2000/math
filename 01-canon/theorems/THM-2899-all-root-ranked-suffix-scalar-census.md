---
id: THM-2899
title: "All-root ranked-suffix scalar census and uniform five-slot parity entry"
status: >
  PROVED + FINITE-EXACT + VERIFIED.  The globally sealed scalar test closes
  26,609 of 41,415 marked suffix branches and five whole roots.  Every one
  of the 14,806 scalar-hard branches satisfies the THM-2895 five-slot
  parity-entry inequality, with an exact finite-core cutoff at most 2,782.
source: codex-lrc-rank-impossible-overlap-2026-07-29
depends_on:
  - THM-735-bonferroni-simultaneous-multi-peel-defeats-the-clustered-non-isolated-wall
  - THM-2892-eight-body-five-slot-heavy-triangle-closure
  - THM-2893-complement-cap-finite-core-flag-lemma
  - THM-2895-singleton-complement-parity-descent-and-four-root-j6-closure
  - THM-2896-seven-body-adaptive-six-cover-hitting-gate-atlas
related:
  - THM-2897-partition-cap-tropical-convolution-and-alternating-pair-ladder
  - THM-2898-unique-max-gate-five-parity-matching-closure
  - THM-2900-flag-conditioned-rank-selective-partition-closure
verification:
  - 04-computation/lrc14_j6_all_root_ranked_suffix_scalar_census_codex_20260729.py
  - 05-knowledge/results/lrc14_j6_all_root_ranked_suffix_scalar_census_codex_20260729.out
  - 05-knowledge/results/lrc14_j6_all_root_ranked_suffix_scalar_hard_ledger_codex_20260729.out
---

# THM-2899 -- all-root ranked-suffix scalar census and uniform five-slot parity entry

**PROVED + FINITE-EXACT + VERIFIED.**

## 1. Marked suffix universe

Let

```text
E in C({1,...,14},7)
```

and let `G_E` be its lonely-time carrier.  THM-2896 assigns `E` an ordered
external six-cover hitting gate

```text
A_E=(a_1,...,a_K),                 K<=25.
```

If a putative six-cover has earliest gate member `a_j`, THM-2893's marked
suffix refinement replaces the root by the literal carrier

```text
C_(E,j)=G_E minus D_(a_j),         h_(E,j)=|C_(E,j)|,
```

and forbids the complete prefix

```text
P_(E,j)={a_1,...,a_j}
```

from the five remaining labels.  This gives exactly

```text
3432 roots,                 41415 marked apex branches.     (1)
```

The branch rank, literal carrier, and excluded prefix are theorem-bearing
data.  They are not reconstructed from the unmarked residual.

## 2. Globally sealed singleton order statistics

For an allowed external label `w`, put

```text
c_(E,j)(w)=|C_(E,j) intersect D_w|.
```

Let

```text
q_1>=q_2>=...>=q_5
```

be the five largest allowed singleton coverages.  They are genuine global
order statistics, not finite-window proposals.  If the carrier has `r`
interval components, THM-735(ii) gives

```text
c_(E,j)(w)<h_(E,j)/7+(99/70)r/(7w).                       (2)
```

The verifier scans every allowed `15<=w<tail_first`, where

```text
tail_first=max(
  1601,
  ceil((99/70)r/(7(q_5-h_(E,j)/7)))
).                                                        (3)
```

It checks `q_5>h/7`, so `(3)` is defined, and strictness in `(2)` seals
every omitted label below or at the retained fifth rank.

For every five allowed labels `Q`,

```text
|C_(E,j) intersect union_(w in Q)D_w|
 <=sum_(w in Q)c_(E,j)(w)
 <=q_1+...+q_5.                                          (4)
```

Thus

```text
q_1+...+q_5<h_(E,j)                                      (5)
```

excludes every five-cover on that marked branch.

## 3. Exact scalar census

The exact result is

```text
branches                                         41415
strict scalar closures                           26609
scalar-hard branches                             14806
whole roots closed scalarly                          5
nonterminal roots under the scalar test           3427
roots with a closed rank before a later open rank   912.  (6)
```

By THM-2896 stratum:

```text
stratum       roots    branches    closed    hard    terminal
low             792       9447       6394    3053       0
one            1848      21986      14133    7853       4
both            792       9982       6082    3900       1. (7)
```

The scalar-hard rank distribution is

```text
rank:   1    2    3    4    5   6   7   8  9  10  11  12  13
hard:3417 3308 2853 2194 1446 844 444 180 78  24  14   3   1. (8)
```

In particular, hardness is early but not an initial segment: `912` roots
are genuinely nonmonotone in the coverage order.

## 4. Five new whole-root closures

The five scalar-terminal bodies are

```text
(1,2,3,4,5,6,13),
(1,2,3,4,6,7,14),
(1,2,3,4,6,11,13),
(1,2,3,4,6,12,13),
(7,8,9,11,12,13,14).                                    (9)
```

For each body in `(9)` and every six-element set
`Q subset {15,16,...}`, the hitting gate gives a unique earliest apex, and
`(5)` closes its marked suffix.  Hence

```text
|G_(E union Q)|>0.                                       (10)
```

If an added label lies in `{1,...,14} minus E`, the resulting eight-body
subfamily is covered by THM-2892, which closes the other five labels.
Therefore `(10)` holds for every legal six-label extension, not only the
external sector.

The five bodies in `(9)` are disjoint from the four THM-2895 bodies and
from THM-2898's unique max-gate body
`(1,8,10,11,12,13,14)`.  Thus these are five genuinely new root closures.
The three proved theorems close ten pairwise distinct roots and leave

```text
3432-4-1-5=3422                                           (10a)
```

seven-body roots in the official residual.

## 5. Uniform five-slot parity entry

On every one of the `14806` scalar-hard carriers, the exact global
singleton cap satisfies

```text
q_1<3h/7.                                                (11)
```

There are zero failures in every stratum and in every branch rank
`1,...,13`.  The global minimum strict margin is

```text
3h/7-q_1 = 369209/35315280                               (12)
```

at

```text
E=(2,4,9,10,12,13,14),   rank=1,   apex=22,
h=62558/315315,          q_1=c(16)=75245/1009008.
```

Apply THM-2895 with `p=5`.  For each hard branch, the finite high core is

```text
H_4={w allowed : c(w)>=(h-q_1)/4}.                        (13)
```

Every hypothetical five-cover contains at least two labels of `H_4`.
Consequently every hard branch reduces to literal three-cover residuals
behind the unordered pairs of its finite `H_4`.

This is the uniform theorem-level change in the frontier:

```text
scalar layer closes, or the same finite 5 -> 3 -> 1
parity grammar applies, on every one of the 41415 branches. (14)
```

For `mu=3h/7-q_1`, the exact discrepancy cutoff for `(13)` is

```text
N_4=ceil(4(99/70)r/(7mu))-1.                              (14a)
```

Its nearest-rank quantiles on all `14806` hard branches are

```text
percent       0   25   50   75   90   95    99   100
N_4         215  427  508  612  729  816  1071  2782.    (14b)
```

The maximum `2782` occurs on the same minimum-margin carrier in `(12)`.
Thus the uniform branch descent has a moderate exact cutoff even though
the competing root-level cutoff in Section 7 has a multimillion outlier.

Statement `(14)` does not assert that the child residuals close.  THM-2900
proves the strongest cheap child cap on four scoped roots only; its
all-root census remains open.

## 6. Adaptive ranked complement flags

The same retained top five give four nested cheap complement caps.  For
`1<=d<=4`, put

```text
B_d^sc=q_1+...+q_d,          s=5-d.                        (15)
```

Subadditivity makes `B_d^sc` a global `d`-set cap.  THM-2893 applies when

```text
B_d^sc<(d+2)h/7.                                           (16)
```

It then forces at least `d+1` labels of every putative five-cover into the
finite high core

```text
H_s={w:c(w)>=(h-B_d^sc)/s}.                                (17)
```

The exact eligibility counts on the `14806` scalar-hard branches are

```text
d   finite-core grammar                     eligible    failures
1   H4 pair -> three-label residual           14806          0
2   heavy H3 triple -> two-label residual      14555        251
3   heavy-pair K4 -> one-label residual        11699       3107
4   five H1 labels -> finite full-cover check   6180       8626. (18)
```

Choosing the largest eligible `d` branch by branch gives

```text
strongest d       1       2       3       4
branches        251    2856    5519    6180.               (19)
```

Eligibility is strictly nested here because every retained `q_i>h/7`.
The `251` rows at `d=1` are exactly the scalar-cap boundary where the H4
pair grammar is the only ranked trigger.  At the other extreme, all
survivors at ranks `11`, `12`, and `13` satisfy the `d=4` trigger.

These are sufficient scalar caps, not exact union caps.  In particular,
failure of `q_1+q_2<4h/7` does not imply failure of the stronger exact
pair-union condition `B_2<4h/7`; an exact all-root pair-cap census can
repair some or all of those `251` rows.

## 7. Competing root-level parity route

The already computed THM-2896 root profiles also give a zero-cost test of
THM-2895 at `p=6`.  Exactly

```text
3200/3432
```

roots satisfy

```text
q_1(E)<2|G_E|/7,                                         (20)
```

while `232` fail.  The eligible counts by stratum are

```text
low 772/792,       one 1728/1848,       both 700/792.     (21)
```

For a positive margin `eta=2h/7-q_1`, the exact discrepancy cutoff for the
five-high-label core is

```text
N=ceil(5(99/70)r/(7 eta))-1.                             (22)
```

Its nearest-rank quantiles on the `3200` eligible roots are

```text
percent       0    25    50    75    90    95     99      100
N           302  1042  1531  2507  4890  9371  38437  7775459. (23)
```

Thus the root `6->4->2` grammar is available on most roots, but its extreme
finite core can be millions wide.  The marked suffix `5->3->1` route is
uniform and should be compared by actual proof cost rather than by seed or
core cardinality alone.

## 8. Boundary rows

The smallest positive scalar margin is

```text
5557/3988104120
```

at `E=(3,4,5,7,8,13,14)`, rank `4`, apex `17`.  The closest failure is

```text
-17/119819700
```

at `E=(3,5,6,7,10,12,13)`, rank `3`, apex `22`.

The canonical hostile is the even arithmetic progression

```text
E=(2,4,6,8,10,12,14),    rank=1,    apex=22,
scalar margin=-27077/270270.                              (24)
```

At the opposite depth extreme, the unique rank-13 survivor is

```text
E=(1,5,8,9,11,12,13),    K=17,      rank=13,    apex=42,
prefix=(21,17,23,19,34,35,50,28,31,49,40,46,42),
h=8683/42042,             r=34,
margin=-2378297/11503321830.                              (25)
```

Its exact top five are

```text
37:823/18648,
20:3713/90090,
29:4195/102312,
51:116317/2858856,
30:4769/120120.                                          (26)
```

The full hard ledger preserves the corresponding data on all `14806`
survivors; `(24)` and `(25)` are canonical hostile controls, not proposed
counterexamples.

## 9. Scope

The scalar computation proves the five root closures in `(9)` and the
uniform parity-entry statement `(11)`--`(14)`.  It does not close the other
`3427` roots under its own scalar test.  After composing with THM-2895 and
THM-2898, `3422` roots remain.  The theorem does not enumerate their H4
pair complexes, prove the all-root THM-2900 child cap, close the seven-body
rung, or prove LRC(14).

## 10. Verification

The verifier contains no Python `assert`; all theorem-bearing guards use
explicit exceptions and therefore remain active under optimization.  It
hash-pins the THM-2896 gate engine, its transcript, and the literal-residual
engine; reconstructs all `3432` roots and all `41415` apex carriers; checks
`165660` scalar/vector controls and `41415` independent literal/direct
carrier controls; and locks the full root-and-branch digest

```text
9770a2ef0e7ae063a482e475f47047de08344cd1fc24880f742651d4c71d167d.
```

Ordinary and optimized full replays are byte-identical, with SHA-256

```text
e6261bb6c4dd13624ddccba32cac7a90f7973929c1c0dc4881743f4aec5d4206.
```

The canonical artifacts are

```text
04-computation/lrc14_j6_all_root_ranked_suffix_scalar_census_codex_20260729.py
SHA-256 e0ff69252870f194549bba61289c1c5b15bef451e37a72836d71f9e71b1016e9

05-knowledge/results/lrc14_j6_all_root_ranked_suffix_scalar_census_codex_20260729.out
SHA-256 dbd3dc5a8c44a55957a6e1ce660ca0e89fcd70e6c0d06d5ba47dc3a22f40c680

05-knowledge/results/lrc14_j6_all_root_ranked_suffix_scalar_hard_ledger_codex_20260729.out
SHA-256 6be9a6c9218f3b42b2eea733c9050f5d35160664af0f19390337b3c5be57cb37
```
