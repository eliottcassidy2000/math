---
id: THM-2904
title: "All-root ranked H1 Hunter closure"
status: >
  PROVED + FINITE-EXACT + VERIFIED.  Every G5-surviving hard row has a
  discrepancy-finite hostile-centre core of size at most 13.  An ordered
  Hunter-pivot bound closes 4,071 centre flags and 279 hard branches.
  Exact root recomposition adds six multi-hard roots, taking the proved
  union to 88 and the official residual to 3,344.
source: root-2026-07-29
depends_on:
  - THM-735-bonferroni-simultaneous-multi-peel-defeats-the-clustered-non-isolated-wall
  - THM-2892-eight-body-five-slot-heavy-triangle-closure
  - THM-2893-complement-cap-finite-core-flag-lemma
  - THM-2897-partition-cap-tropical-convolution-and-alternating-pair-ladder
  - THM-2899-all-root-ranked-suffix-scalar-census
  - THM-2901-all-hard-exact-global-pair-cap-and-route-partition
  - THM-2903-one-hard-actual-h3-link-and-bad-triangle-closure
  - THM-2905-all-hard-hunter-star-envelope-census
related:
  - THM-2902-omission-six-ranked-h1-depth-one-closure
  - THM-2907-paircap-exception-h4-heavy-link-child-census
verification:
  - 04-computation/lrc14_j6_all_hard_ranked_h1_hunter_pivot_census_codex_20260729.py
  - 05-knowledge/results/lrc14_j6_all_hard_ranked_h1_hunter_pivot_census_codex_20260729.out
---

# THM-2904 -- all-root ranked H1 Hunter closure

**PROVED + FINITE-EXACT + VERIFIED.**

## 1. Statement

On every one of the `11,842` scalar-hard marked suffixes not closed by
THM-2905's strict `G_5<h` test, a hypothetical five-cover has a canonical
maximum singleton inside an exact finite core of size at most `13`.
Ordering that core and retaining the exact pair cap closes `4,071` of the
`55,293` possible centre pivots and all pivots on `279` hard branches.

After composing with THM-2899, THM-2905, and the already proved root union,
six new multi-hard seven-body roots close.  Thus the proved union has size

```text
88
```

and the official residual is

```text
3432-88=3344.                                             (1)
```

## 2. The hostile-centre threshold

Fix a THM-2905-surviving marked carrier `C`, with mass `h`, excluded prefix
`P`, exact singleton ranks `q_1>=...>=q_5`, and exact pair-union cap
`beta_2`.  Define

```text
g(a)=a+sum_(r=2)^5 min(a,q_r,beta_2-a),                  (2)
```

on `0<=a<=min(q_1,beta_2)`.  Since the row was not closed by the strict
`G_5<h` test, the maximum of `g` is at least `h`.  Each summand in `(2)` is
the minimum of affine functions, so `g` is continuous and concave.  Its
hostile superlevel is therefore a closed interval.  Put

```text
lambda=min{a:g(a)>=h}.                                   (3)
```

Suppose five allowed labels cover `C`, and order their singleton coverages
as `a_1>=...>=a_5`.  The Hunter-star proof in THM-2897 gives

```text
h<=U_C(Q)<=g(a_1).                                       (4)
```

Consequently `a_1>=lambda`.  Thus every possible maximum singleton belongs
to the actual hostile-centre core

```text
H_1^star={w allowed:c_C(w)>=lambda}.                     (5)
```

This is stronger than merely knowing that some high singleton exists:
`(3)` remembers the entire failed star envelope and selects only labels
capable of serving as its centre.

## 3. Exact finite cores

Here “allowed” retains THM-2899's exact marked-suffix meaning:
`w>=15` and `w` is outside the excluded prefix.  The piecewise-linear
breakpoints of `(2)` are

```text
0, min(q_1,beta_2), beta_2/2, q_r, beta_2-q_r
                                      for r=2,3,4,5,     (6)
```

clipped to the legal domain.  The verifier solves the first affine crossing
of height `h` exactly and checks affine midpoint identities on every cell.
On all `11,842` rows it finds the strict finite-tail gap

```text
lambda-h/7>0.                                            (7)
```

If `C` has `r` components, THM-735(ii) gives

```text
c_C(w)<h/7+gamma/w,             gamma=(99/70)r/7.        (8)
```

Hence every member of `(5)` satisfies

```text
w<=ceil(gamma/(lambda-h/7))-1.                           (9)
```

The complete exact census is

```text
G_5-surviving hard rows                              11,842
pair-cap exceptions included                            52
allowed labels scanned                           4,797,677
cutoff range                                       149..1013
actual hostile-centre labels                         55,293
actual core-size range                                1..13
empty actual cores                                        0
median actual core size                                   5. (10)
```

The analytic bin is already modest, but actual membership removes another
factor of about `87`.  The resulting object is a list of at most thirteen
candidate maxima per hard row.

## 4. Ordered pivot bound

Order the actual core as

```text
x_1 prec x_2 prec ... prec x_m
```

by decreasing parent coverage, breaking equal-coverage ties by increasing
label.  Write

```text
a_i=c_C(x_i),             C_i=C minus D_(x_i),
h_i=h-a_i,                rho_i=beta_2-a_i.             (11)
```

Assign a hypothetical cover to the earliest one of its maximum-coverage
labels.  If its centre is `x_i`, none of `x_1,...,x_(i-1)` may occur among
the four leaves.

For a later core label `x_j`, retain its exact parent coverage `a_j`.  Every
label outside the core has parent coverage strictly below `lambda`.  Also,
for every allowed leaf `y`,

```text
c_(C_i)(y)
 =U_C({x_i,y})-a_i
 <=beta_2-a_i=rho_i.                                    (12)
```

Form the multiset

```text
M_i={min(a_j,rho_i):j>i}
       union {min(lambda,rho_i)} repeated four times,    (13)
```

and let `S_i` be the sum of its four largest members.  Four formal copies
of the noncore bound suffice because the child has four slots.  Equations
`(11)`--`(13)` show that the sum of the four child singleton coverages is
at most `S_i`.  Therefore the strict inequality

```text
S_i<h_i                                                   (14)
```

excludes the entire pivot branch.

There is one strict refinement at a zero displayed margin.  If
`rho_i>=lambda` and no four later-core bounds attain `S_i`, then every
four-leaf choice either has bound sum below `S_i` or uses a noncore label.
In the latter case its parent coverage, and hence its child coverage, is
strictly below `lambda`.  Therefore

```text
S_i=h_i, rho_i>=lambda, and no four later-core bounds attain S_i
  implies the pivot branch is excluded.                  (14')
```

This is the ordered form of the Hunter star.  The parent order restores
compatibility that `G_5` discarded: exact later-core ranks replace four
independently maximized global ranks, while `(12)` retains the centre--leaf
overlap credit.

The exact pivot census gives

```text
possible centre pivots                                  55,293
certified pivot closures                                 4,071
pivots surviving (14)/(14')                             51,222
raw equality pivots                                          4
equality pivots closed by strict noncore bounds               4
open equality pivots                                         0
hard branches with every pivot closed                      279
hard branches surviving this test                       11,563. (15)
```

None of the `52` pair-cap exception branches closes by `(14)`.  Their H4
heavy-link recursion remains a separate live lane.

## 5. Whole-root recomposition

Exactly ten bodies have all their `G_5`-surviving hard rows closed by
`(14)`:

```text
(1,2,3,4,5,6,11),
(1,2,3,4,5,8,9),
(1,2,3,5,7,10,14),
(1,2,4,8,11,12,13),
(1,2,5,6,10,12,13),
(2,4,6,8,11,12,13),
(3,4,6,7,9,11,12),
(5,6,7,8,10,11,12),
(5,7,9,10,11,12,14),
(6,7,8,9,11,12,13).                                    (16)
```

Four of `(16)` are one-hard roots already in THM-2903:

```text
(1,2,3,5,7,10,14),
(1,2,4,8,11,12,13),
(1,2,5,6,10,12,13),
(2,4,6,8,11,12,13).                                    (17)
```

The other six are new:

```text
(1,2,3,4,5,6,11),
(1,2,3,4,5,8,9),
(3,4,6,7,9,11,12),
(5,6,7,8,10,11,12),
(5,7,9,10,11,12,14),
(6,7,8,9,11,12,13).                                    (18)
```

For every body in `(16)`, hard rows absent from this census already satisfy
THM-2905's `G_5<h` test, and scalar-soft rows satisfy THM-2899.  THM-2892
supplies the eight-body chamber.  Therefore every body in `(16)` is a
whole-root closure.

The exact proved union through THM-2905 has size `82`; its intersection
with `(16)` is precisely `(17)`.  Hence `(18)` adds six roots and proves
`(1)`.

## 6. Strict equality and failure boundary

The smallest positive pivot margin is

```text
1/336336
```

at

```text
E=(2,3,4,5,7,12,13), rank 2, apex 18, centre 16.        (19)
```

The four raw equality pivots are

```text
E=(1,2,3,4,7,9,14),  rank 1, apex 24, centre 20;
E=(1,2,3,4,9,12,14), rank 4, apex 44, centre 15;
E=(1,2,6,8,9,11,12), rank 1, apex 20, centre 52;
E=(1,3,4,6,8,9,13),  rank 2, apex 21, centre 20.        (20)
```

In every row in `(20)`, the displayed centre is the last label of its
ordered core and `rho_i>=lambda`.  Thus all four leaves are noncore and
each has child coverage strictly below `lambda`.  Equation `(14')` closes
all four pivots.  There are no open equality pivots.

The closest genuinely failed pivot has raw margin

```text
-1/2828280
```

at

```text
E=(1,6,7,8,10,11,13), rank 1, apex 19, centre 37.       (21)
```

The order on `H_1^star` is an acyclic allocation of proof branches, not a
tournament on runners.  Parent coverages may tie, and the label tie-break
has no intrinsic orientation content.  The faithful object is a ranked
core with a symmetric pair-cap sidecar.  It forgets actual child
intersections beyond the inherited centre--leaf cap, so failure of `(14)`
is not evidence for a cover.

## 7. Verification and scope

The verifier hash-pins the THM-2899 and THM-2901 ledgers, both exact
interval engines, and the seven canonical theorem outputs composing the
proved `82`-root union through THM-2905.  It artifact-derives every prior
root group, THM-2903's overlap and proof digest, and THM-2905's additive
set difference.  It joins rows by body, rank, apex, and excluded prefix;
reconstructs every literal carrier; checks scalar/vector agreement; solves
every threshold in `(3)`; proves every cutoff in `(9)`; retains every exact
core coverage; and locks every pivot bound in `(13)`--`(14')`.  All guards
are explicit and survive optimized Python.  Repository-text hashes
normalize CRLF to LF and reject lone carriage returns; the imported
residual/vector engines are scoped to the same canonical reader.

The complete core-and-pivot semantic digest is

```text
ec878244b922ba5f48633614a86a1f9706c1fbdd0ebd6c61f020291cfd737bab.
```

Canonical artifacts:

```text
04-computation/lrc14_j6_all_hard_ranked_h1_hunter_pivot_census_codex_20260729.py
SHA-256 644104b0de90654466e75c6531109736b0445aadb357eee2413e8787ac3a53fa

05-knowledge/results/lrc14_j6_all_hard_ranked_h1_hunter_pivot_census_codex_20260729.out
SHA-256 0933c67a108b6d588e36737fb2b17b325ca36146976cfb035bebe036a6234036
```

Locked ordinary and optimized replays are byte-identical:

```bash
python3 04-computation/lrc14_j6_all_hard_ranked_h1_hunter_pivot_census_codex_20260729.py
python3 -O 04-computation/lrc14_j6_all_hard_ranked_h1_hunter_pivot_census_codex_20260729.py
```

This theorem closes the six roots in `(18)`.  It does not close the
`11,563` hard rows surviving its ordered-pivot test, any pair-cap exception
branch, the remaining `3,344` roots, or LRC(14). ∎
