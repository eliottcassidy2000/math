---
id: THM-2912
title: "One-H3-row exact child top-four closure"
status: >
  PROVED + FINITE-EXACT + VERIFIED.  On the 210 post-THM-2904/2907
  residual seven-body roots having exactly one ordinary H3 row, exact
  ordered-centre deletion and a globally sealed child top-four census
  close 765 of 807 open pivots and all pivots on 172 roots.  Against
  THM-2911's current 181-root union, 39 overlap and 133 are additive,
  taking the proved union to 314 and the finite residual to 3,118.  This
  is not LRC(14).
source: codex/hunter-star-resume-2026-07-29
depends_on:
  - THM-735-bonferroni-simultaneous-multi-peel-defeats-the-clustered-non-isolated-wall
  - THM-2899-all-root-ranked-suffix-scalar-census
  - THM-2901-all-hard-exact-global-pair-cap-and-route-partition
  - THM-2904-all-root-ranked-h1-hunter-closure
  - THM-2907-paircap-exception-h4-heavy-link-child-census
  - THM-2911-all-root-finite-ranked-h1-hunter-closure
related:
  - THM-2893-complement-cap-finite-core-flag-lemma
  - THM-2897-partition-cap-tropical-convolution-and-alternating-pair-ladder
verification:
  - 04-computation/lrc14_j6_one_h3_row_child_top4_census_codex_20260729.py
  - 05-knowledge/results/lrc14_j6_one_h3_row_child_top4_census_codex_20260729.out
---

# THM-2912 -- one-H3-row exact child top-four closure

**PROVED + FINITE-EXACT + VERIFIED.**

## 1. Inherited finite universe

Retain THM-2899's marked suffix data: a seven-body root `E`, branch rank,
marked apex, forbidden prefix `P`, literal good carrier `C`, and mass
`h=|C|`.  THM-2904 applies the Hunter-star test to every one of the
`11,842` rows that survived the scalar `G_5<h` test.  Its exact ordered
hostile-centre core is

```text
H_1^star(C)=(x_1,...,x_s),       1<=s<=13,                (1)
```

ordered by decreasing singleton coverage and then increasing label.
THM-2904 already closes all pivots on `279` rows.  THM-2907 separately
discharges the `52` pair-cap-exception H4 rows.  The post-THM-2904/2907
ordinary H3 universe therefore has

```text
ordinary H3 rows                                      11,511
rows on the previously proved 88 roots                    76
rows on that checkpoint's residual roots                11,435
checkpoint residual roots                                3,344. (2)
```

The number of residual ordinary H3 rows on a root has exact histogram

```text
rows/root       1    2    3    4    5    6   7   8  9  10
roots          210  726  943  754  449  181  60  13  7   1. (3)
```

This theorem treats exactly the first bin of `(3)`.  Thus every selected
root has one remaining ordinary H3 row; all its other marked branches have
already been discharged by the inherited partition.

## 2. Exact ordered-centre children

Suppose a hypothetical five-cover of `C` is assigned, as in THM-2904, to
the earliest one of its maximum-coverage labels.  If the assigned centre is
`x_i`, put

```text
R_i=C minus D_(x_i),              h_i=|R_i|.              (4)
```

The four leaves are distinct allowed labels outside

```text
P union {x_i} union {x_1,...,x_(i-1)}.                   (5)
```

The first exclusion is the marked-suffix rule, the second prevents reuse
of the centre, and the third is forced by earliest-maximum allocation:
an earlier hostile-core label has parent coverage at least that of `x_i`
and would receive the cover first.  Later core labels remain allowed.

THM-2904 has already excluded every pivot satisfying its parent
Hunter-star bound.  On the 210 rows in `(3)`, exactly

```text
open ordered-centre children                              807             (6)
```

remain.  The verifier reconstructs `C` from the body and apex, deletes the
selected danger comb literally, and checks

```text
|R_i|=h-c_C(x_i)>0.                                      (7)
```

No quotient carrier or unmarked residual is substituted for `(4)`.

## 3. Globally sealed child top four

For every label `w` allowed by `(5)`, write

```text
c_i(w)=|R_i intersect D_w|.
```

The verifier evaluates `c_i(w)` exactly for every allowed
`15<=w<=2000`.  If `R_i` has `r_i` interval components, the discrepancy
bound used in THM-2899/2904 gives, for every omitted `w>=2001`,

```text
c_i(w)<h_i/7 + (99/70)r_i/(7w)
      <=h_i/7 + (99/70)r_i/(7*2001).                     (8)
```

Here the exact child component range is

```text
20<=r_i<=56.
```

Let

```text
q_(i,1)>=q_(i,2)>=q_(i,3)>=q_(i,4)
```

be the four largest scanned singleton coverages.  On all `807` children,
the exact fourth-rank tail gap is positive; its minimum is

```text
q_(i,4)-[h_i/7+(99/70)r_i/(7*2001)]
 >=12590819/2860285428 >0.                               (9)
```

Thus the four scanned ranks are the genuine global top four among all
labels allowed by `(5)`.  For any four leaves `Y`,

```text
|R_i intersect union_(y in Y)D_y|
 <=sum_(y in Y)c_i(y)
 <=q_(i,1)+q_(i,2)+q_(i,3)+q_(i,4).                     (10)
```

Consequently the strict exact deficit

```text
Delta_i=h_i-sum_(k=1)^4 q_(i,k)>0                        (11)
```

excludes the entire centre child.  This uses only the union bound on the
literal child; it does not assume independence or discard overlaps.

## 4. Exact census and whole-root recomposition

The finite result is

```text
ordered-centre children                                  807
globally tail-sealed children                            807
strict top-four closures                                 765
zero margins                                               0
negative margins                                          42
one-row roots with every pivot closed                    172
one-row roots surviving this test                         38. (12)
```

The smallest positive margin in `(11)` is

```text
Delta=43/2751840
E=(1,2,4,8,10,11,12), rank=1, apex=28, centre=52.        (13)
```

The closest failure is already negative:

```text
Delta=-3/140140
E=(1,2,3,5,6,11,12), rank=1, apex=27, centre=21.         (14)
```

The exact sorted list of all `172` closed roots is the explicit
`closed_roots=` record in the pinned output artifact.  Its independent
list digest is

```text
fcb9e88dd0b5e8aa855e0164636a17f42b4f1620407b9a0a780e50fd48336217.
```

The verifier reconstructs the previous `88`-root union from the
hash-pinned THM-2904 dependencies.  It then imports THM-2911's explicit
`93` additive roots, checks that their union with the baseline has size
`181`, and takes exact set difference.  Of the `172` roots above, `39`
already lie in the THM-2911 union and `133` are additive.  Since each
selected root had exactly one remaining ordinary H3 row at checkpoint
`(2)`, closing every pivot on that row closes the whole root.  Therefore
the current composition is

```text
proved root baseline from THM-2911                  181
THM-2912 overlap                                     39
THM-2912 additive roots                             133
proved root union                           181+133=314
finite root residual                     3432-314=3118. (15)
```

## 5. Boundary and scope

There are no equality closures hidden in `(12)`: all promoted children
have the strict margin `(11)`, and all `42` unclosed children have negative
margin.  Inequality `(9)` is also strict, so no omitted label can tie the
retained fourth rank.  The 38 listed failed roots are not declared
counterexamples; they only survive this sufficient top-four certificate.

This theorem is confined to the one-row bin in `(3)`.  It does not treat
roots outside that checkpoint stratum, prove that any surviving H3 branch
has a five-cover, classify equality families, or prove LRC(14).

## 6. Verification

The verifier hash-pins the final THM-2904 source/output, the THM-2907
census and endpoint output artifacts, and THM-2911's locked composition
output.  It reconstructs all `11,842` parent rows, checks the complete
filter/histogram controls `(2)`--`(3)`, and uses exact rational interval
arithmetic throughout.  Repository-text dependency hashes explicitly
normalize CRLF to LF and reject lone carriage returns.  The child ledger
semantic digest is

```text
39f3a9fb8ec6447baf96512bee3ee174e2390639ab3bbf0a6a36dcb5cdf0274e.
```

Single-worker, eight-worker, and optimized eight-worker replays are
byte-identical.  Final SHA-256 values are

```text
source  4859162ec9d03f7f711de99e97122c916e708adfdb8b86581e9e9534c735fa9f
output  454d87c8beeb81405b031cce4b40bdda0d385cfcd9c48e6fcf4eb810cfc00c5a.
```

Reproduction:

```bash
python3 04-computation/lrc14_j6_one_h3_row_child_top4_census_codex_20260729.py --workers 1
python3 -O 04-computation/lrc14_j6_one_h3_row_child_top4_census_codex_20260729.py --workers 8
```
