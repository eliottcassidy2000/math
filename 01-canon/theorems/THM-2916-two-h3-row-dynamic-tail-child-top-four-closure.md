---
id: THM-2916
title: "Two-H3-row dynamic-tail child top-four closure"
status: >
  PROVED + FINITE-EXACT + VERIFIED.  After THM-2913 the exact residual
  has no one-H3-row roots; its next stratum consists of 690 roots with
  two ordinary H3 rows.  A per-child globally sealed dynamic tail closes
  5,251 of 5,618 open ordered-centre children and every child on 394
  roots, all additive.  The proved union becomes 745 and the residual
  2,687.  The 367 failed children on 296 roots remain open to stronger
  pair/Hunter/toothpick certificates.  This is not LRC(14).
source: codex-post2913-residual-scout-2026-07-29
depends_on:
  - THM-735-bonferroni-simultaneous-multi-peel-defeats-the-clustered-non-isolated-wall
  - THM-2904-all-root-ranked-h1-hunter-closure
  - THM-2907-paircap-exception-h4-heavy-link-child-census
  - THM-2911-all-root-finite-ranked-h1-hunter-closure
  - THM-2912-one-h3-row-exact-child-top-four-closure
  - THM-2913-one-h3-row-pair-hunter-toothpick-closure
related:
  - THM-2915-all-open-centre-exact-child-top-four-closure
verification:
  - 04-computation/lrc14_j6_two_h3_row_child_top4_census_codex_20260729.py
  - 05-knowledge/results/lrc14_j6_two_h3_row_child_top4_census_codex_20260729.out
  - 05-knowledge/results/lrc14_j6_two_h3_row_child_top4_census_codex_20260729.ledger.out
---

# THM-2916 -- two-H3-row dynamic-tail child top-four closure

**PROVED + FINITE-EXACT + VERIFIED.**

## 1. Exact next residual stratum

THM-2911 proves a union of `181` seven-body roots.  THM-2912 takes the
union to `314`, and THM-2913 closes the whole one-H3-row stratum, taking
it to

```text
proved union                                             351
finite residual                            3,432-351=3,081.           (1)
```

Reconstruct all `11,842` THM-2904 parent rows and discard its closed
branches together with THM-2907's `52` pair-cap exceptions.  The remaining
ordinary H3 map has `11,511` rows on `3,401` roots.  Exact set difference
against the `351`-root union in `(1)` gives the residual histogram

```text
ordinary H3 rows/root       2    3    4    5    6   7   8  9  10
residual roots            690  929  751  449  181  60  13  7   1.    (2)
```

In particular, there is no remaining one-row root.  This theorem treats
exactly the first bin of `(2)`: `690` roots carrying `1,380` ordinary H3
rows.

## 2. The marked ordered child

For one parent row retain its full theorem-bearing state

```text
E       seven-body root,
a       marked apex,
P       inherited excluded gate prefix,
C       literal parent carrier,
x_i     an open ordered hostile centre,
A_i     the hostile centres earlier than x_i.            (3)
```

THM-2904 orders hostile centres by decreasing parent coverage and then
increasing label.  A hypothetical five-cover assigned to `x_i` leaves
four distinct labels covering

```text
R_i=C minus D_(x_i),                                    (4)
```

and those labels avoid

```text
F_i=P union {x_i} union A_i.                            (5)
```

The earlier-centre exclusion is forced by the earliest-maximum allocation;
later centres remain allowed.  Across the `1,380` rows in `(2)`, there are
exactly

```text
open ordered-centre children                           5,618.          (6)
```

The verifier reconstructs `(4)` both by literal subtraction and directly
from `E union {a,x_i}`.  It checks the interval list, component count, and
mass on both paths.  Thus neither an unmarked residual nor a quotient
carrier replaces the child.

## 3. Per-child global tail seal

For every label `w>=15` outside `(5)`, put

```text
c_i(w)=mu(R_i intersect D_w),       h_i=mu(R_i).
```

If `R_i` has `r_i` interval components, THM-735 gives the strict bound

```text
c_i(w)<h_i/7 + (99/70) r_i/(7w).                        (7)
```

Rather than impose one large common horizon, the verifier starts at
`N=127` and successively replaces `N` by `2N+1`.  At each stage it
retains the exact four largest scanned coverages

```text
q_(i,1)>=q_(i,2)>=q_(i,3)>=q_(i,4)
```

and stops only when

```text
q_(i,4)>=h_i/7 + (99/70)r_i/[7(N+1)].                   (8)
```

Every omitted label has `w>=N+1`, so the strict inequality `(7)` makes
its coverage strictly smaller than the right side of `(8)`.  Hence the
retained four values are the genuine global top four, including the
possible equality boundary in `(8)`.

All `5,618` children seal.  Their exact stopping-horizon histogram is

```text
N              255    511   1023  2047
children        19   3,086  2,501    12.                (9)
```

The computation scans `4,060,613` allowed labels in total.  Child
component counts range from `18` to `62`, and the minimum exact
fourth-rank tail gap is

```text
3833/4116752640>0.                                      (10)
```

Thus no equality in the tail seal is actually attained.

## 4. Strict closure and root recomposition

For any four allowed labels `Y`, subadditivity and the global rank seal
give

```text
mu(R_i intersect union_(y in Y)D_y)
 <=sum_(y in Y)c_i(y)
 <=q_(i,1)+q_(i,2)+q_(i,3)+q_(i,4).                    (11)
```

Therefore the strict deficit

```text
Delta_i=h_i-sum_(j=1)^4 q_(i,j)>0                       (12)
```

closes the entire ordered-centre child.  The exact census is

```text
globally tail-sealed children                           5,618
strict top-four closures                                5,251
zero margins                                                0
negative margins                                          367
two-row roots with every open child closed                394
roots retaining at least one failed child                 296.        (13)
```

The smallest positive margin is

```text
22664/350525175
E=(1,2,3,5,7,9,11), rank=1, apex=24, centre=26.         (14)
```

The closest failure is

```text
-1/5045040
E=(1,2,4,6,7,9,11), rank=1, apex=20, centre=16.         (15)
```

Every selected root lay outside the union `(1)`, so all `394` closed
roots are additive.  Both ordinary H3 rows of each such root have every
open pivot closed by `(12)`; their other branches were already discharged
by the inherited THM-2904/2907 partition.  Hence

```text
proved union                                      351+394=745
finite residual                             3,432-745=2,687.          (16)
```

The output contains the explicit sorted `394`-root and `296`-root lists.
The SHA-256 digest of the sorted closed tuple, which is also the additive
tuple, is

```text
c68d09676683f6204df3b04353a3b3107ebbb4285d13a3b6001446372e351e1b.
```

## 5. Boundary and scope

All promoted inequalities in `(12)` are strict.  The `367` negative
children are certificate failures, not four-cover witnesses and not
counterexamples.  Pair overlap, Hunter trees, and recursive second-centre
geometry are deliberately absent from this theorem; they form a logically
separate repair lane.

The dynamic tail is a computational stopping law, not an assumed uniform
period bound.  The order in `(3)` allocates proof branches and is not a
tournament on runners.  This theorem treats only the two-row bin in `(2)`,
does not close its `296` surviving roots, does not treat higher-row bins,
and does not prove LRC(14).

## 6. Verification

The verifier hash-pins the THM-2912 and THM-2913 sources and outputs,
reconstructs the complete THM-2904 parent universe, derives the `351`-root
union from explicit root lists, and checks every histogram and set
difference above.  Repository-text hashes normalize CRLF to LF and reject
lone carriage returns.

Ordinary and optimized eight-worker replays produce byte-identical summary
and `5,618`-line ledgers.  The semantic ledger digest is

```text
5ae96776a6ffe025c3246dbb3e5381c040d7ad2f4f99315115651b70622252c6.
```

LF-normalized SHA-256 values are

```text
source  7d3c4e82abb0ce3af13c43c1e03f09d4be1ee3dbfd496ca05273672cd2a462a6
output  6b31abaadd4f089a9f98a9eea49845c0ed0123a810bb4bf4f3c0309c519e7df9
ledger  a202807a1da472ed6c19109425478922b2218c48f2afdb39c67a0b9f2f5fb7c4.
```

Reproduction:

```bash
python3 04-computation/lrc14_j6_two_h3_row_child_top4_census_codex_20260729.py --workers 8
python3 -O 04-computation/lrc14_j6_two_h3_row_child_top4_census_codex_20260729.py --workers 8
```
