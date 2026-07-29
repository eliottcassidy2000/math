---
id: THM-2915
title: "All-open-centre exact child top-four closure"
status: >
  PROOF CANDIDATE / CLEAN REPLAY PENDING.  Exact dynamic-tail top-four
  closure on all 51,222 open THM-2904 centre children closes 46,356
  children and 8,112 parent rows.  Branchwise composition with the proved
  E/H/P/Q routes gives 964 roots additive over the current 351-root
  baseline, for candidate union 1,315 and residual 2,117.
source: codex-thm2915-all-open-centre-child-top4-2026-07-29
depends_on:
  - THM-735-bonferroni-simultaneous-multi-peel-defeats-the-clustered-non-isolated-wall
  - THM-2904-all-root-ranked-h1-hunter-closure
  - THM-2907-paircap-exception-h4-heavy-link-child-census
  - THM-2911-all-root-finite-ranked-h1-hunter-closure
  - THM-2912-one-h3-row-exact-child-top-four-closure
  - THM-2913-one-h3-row-pair-hunter-toothpick-closure
related:
  - THM-2893-complement-cap-finite-core-flag-lemma
  - THM-2897-partition-cap-tropical-convolution-and-alternating-pair-ladder
verification:
  - 04-computation/lrc14_j6_all_open_centre_child_top4_closure_thm2915.py
  - 05-knowledge/results/lrc14_j6_all_open_centre_child_top4_closure_thm2915.out
  - 05-knowledge/results/lrc14_j6_all_open_centre_child_top4_closure_thm2915.ledger.out
---

# THM-2915 -- all-open-centre exact child top-four closure

**PROOF CANDIDATE / CLEAN REPLAY PENDING.**

Ordinary and optimized theorem-grade computations have completed with
every exact guard passing and byte-identical artifacts.  Promotion awaits
the final fresh-checkout replay and final artifact hashes.

## 1. Candidate statement

For every one of the `51,222` ordered hostile-centre pivots left open by
THM-2904, retain its literal marked child and compute the exact global four
largest allowed singleton coverages.  The strict top-four deficit closes

```text
open ordered centre children                         51,222
strict child top-four closures                       46,356
child failures                                        4,866
equality children                                         0
globally tail-sealed children                        51,222.          (1)
```

At parent-row level this closes `8,112` of the `11,563` rows having at
least one open THM-2904 centre.  Composing these rows with the previously
proved pair-exception, finite-H1, fully pivot-closed, and one-H3 toothpick
routes closes `1,285` route roots.  Exact set difference against the
current `351`-root union through THM-2913 gives

```text
route roots intersecting the current union              321
new roots                                                964
proved union candidate                                 1,315
official residual candidate                3,432-1,315=2,117.        (2)
```

This is a finite exact closure candidate, not a proof of LRC(14).

## 2. The inherited marked child

Fix a THM-2904-surviving marked carrier.  Its theorem-bearing data are

```text
C          literal post-apex carrier,
h=mu(C)    carrier mass,
P          inherited excluded gate prefix,
x_1,...,x_m
           hostile centres, ordered by decreasing coverage and then label.
```

THM-2904 allocates a hypothetical five-cover to the earliest one of its
maximum-coverage labels.  On the branch with centre `x=x_i`, put

```text
R=C minus D_x,
h_R=mu(R),
F=P union {x} union {x_1,...,x_(i-1)}.                  (3)
```

The four remaining labels are distinct integers at least `15`, avoid
`F`, and must cover `R`.  The earlier-centre sidecar in `(3)` is
load-bearing: allowing an earlier centre would merge two branches of the
allocation and invalidate the root projection.

Every one of the `51,222` children is reconstructed in two independent
forms:

```text
R = C minus D_x
  = G_(E union {a,x}),                                  (4)
```

where `E` is the seven-body root and `a` is the marked apex.  The ordered
interval lists, component counts, and exact masses agree in every case.

## 3. Dynamic finite tail with no arbitrary cap

Scan every allowed label through

```text
M_0=2,000
```

and let `q_4(M_0)` be the fourth retained singleton coverage on `R`.
If `r_R` is the number of components, THM-735 gives for every label `w`

```text
c_R(w)<h_R/7+gamma/w,        gamma=(99/70)r_R/7.         (5)
```

The computation finds

```text
q_4(M_0)-h_R/7>0
```

on every child.  Therefore the exact integer

```text
N_0=max(2,000,
        ceil(gamma/(q_4(M_0)-h_R/7))-1)                  (6)
```

is a sufficient head horizon.  The verifier scans literally through
`N_0`, recomputes `q_4`, and checks

```text
q_4 >= h_R/7 + gamma/(N_0+1).                           (7)
```

For any omitted `w>N_0`, `(5)` is strict and its right side is at most the
right side of `(7)`.  Thus the omitted value is strictly below `q_4`,
even if equality holds in `(7)`.  The four retained ranks are consequently
the exact global allowed top four.  Equality is deliberately accepted in
the tail seal; closure itself remains strict.

Exactly one child extends past `2,000`, to `2,134`:

```text
E=(1,6,8,9,10,12,13), rank 4, apex 44, centre 156.
```

Its tail gap, also the global minimum, is

```text
4163/4039185150.                                        (8)
```

There is no discovery cap such as `100,000` in the theorem computation:
`(6)` supplies the complete finite universe separately for each child.

## 4. Exact child closure and failure boundary

Let the exact allowed singleton ranks be

```text
q_1>=q_2>=q_3>=q_4.
```

Subadditivity gives for any four remaining labels

```text
mu(R intersect D_Q)<=q_1+q_2+q_3+q_4.                  (9)
```

Hence the strict inequality

```text
q_1+q_2+q_3+q_4<h_R                                    (10)
```

closes the centre branch.  Equation `(10)` holds on `46,356` children and
fails on `4,866`; no child is an equality case.

The smallest positive margin is

```text
412297/346686319980
```

at

```text
E=(1,2,5,7,8,11,13), rank 3, apex 36, centre 50.        (11)
```

The closest failure is

```text
-1/5045040
```

at

```text
E=(1,2,4,6,7,9,11), rank 1, apex 20, centre 16.         (12)
```

Failure of `(10)` means only that this sufficient four-singleton invoice
does not close the child.  It is not evidence for a four-cover.

## 5. Faithful branchwise composition

Use the following sets of THM-2904 survivor-row keys
`(E,rank,apex,P)`:

```text
E   52 pair-cap exception rows discharged by THM-2907,
H   3,090 finite-H1 rows discharged by THM-2911,
P   279 rows whose every hostile-centre pivot closes in THM-2904,
T   8,112 rows whose every remaining centre child closes by (10),
Q   38 one-H3 rows discharged by THM-2913.              (13)
```

Before adjoining `Q`, the exact atoms on the `11,842` survivor rows are

```text
neither  3,375       E       52       H       24
HP         215       HT   2,851       P       64
T        5,261.                                            (14)
```

The only nonzero intersections are

```text
|H intersect P|=215,          |H intersect T|=2,851.      (15)
```

In particular, every exception row has at least one failed top-four child:
`E intersect T` is empty.  This explains why the older routes can add
genuine strength.  The predicate defining `T` requires every open centre
child to pass `(10)`, whereas `E` and `H` close the entire parent row by
different information and can absorb a failed child.

Adjoining `Q` replaces `37` neither-atoms by `Q` and one `H` atom by
`HQ`; there is no `Q` intersection with `E`, `P`, or `T`.  Thus
THM-2913 contributes its expected `37` roots and no additional hidden
cross-route synergy.

## 6. Root quantifiers and the twelve synergy roots

The `3,432` seven-body roots split first as

```text
no G5-surviving row                                      21
all THM-2904 pivots closed                               10
at least one open-centre child                        3,401.          (16)
```

The live `3,401`-root universe is the correct denominator for the pure
child route.  Exactly

```text
1,219 close,          2,182 fail.                       (17)
```

Adding the ten vacuous `P` roots gives `1,229` pure `P union T` roots.
Against the canonical union of `314` through THM-2912, this pure route has

```text
overlap 277, additive 952, union 1,266, residual 2,166. (18)
```

Full `E union H union P union T` composition closes `1,248` roots.  It
gains nineteen over the pure route; seven were already known, while
twelve are genuinely additive.  The six exception synergies are

```text
(1,2,7,8,9,10,11)
(1,6,7,8,9,10,11)
(2,3,4,6,10,13,14)
(2,3,4,7,8,10,14)
(2,4,6,7,10,13,14)
(2,7,8,9,10,11,12).                                    (19)
```

Each root in `(19)` has exactly one pure-route-missing exception row, and
all six are among THM-2907's literal H4 closures rather than its
endpoint-repair case.  The six finite-H1 synergies are

```text
(1,2,3,4,6,7,9)
(1,2,3,6,9,11,12)
(1,2,4,11,12,13,14)
(1,2,5,9,11,12,13)
(1,3,4,8,9,10,11)
(1,3,4,9,10,11,12).                                    (20)
```

Each root in `(20)` has exactly one pure-route-missing finite-H1 row.
Those six rows are genuine literal-search closures in the locked THM-2911
shards, not a relabeling of the top-four test.  No root needs a mixed
`E`-and-`H` synergy.

Finally, adjoining `Q` gives `1,285` branch-composed route roots.  Their
intersection with the live `351`-root union is `321`, and their additive
part is the `964` roots in `(2)`.

## 7. Inherited-slice and hostile controls

On the pair-exception slice the exact census is

```text
exception rows                                           52
children                                                400
top-four closed / failed                            272 / 128
fully T-completed exception rows                           0.          (21)
```

The histogram of failed children per exception row is

```text
1:5, 2:26, 3:14, 4:6, 5:1.                             (22)
```

Thus the six root gains in `(19)` really use THM-2907; they are not
silently certified by `(10)`.

On the one-H3 slice, the general verifier reproduces THM-2912 exactly:

```text
ordinary parent rows                                    210
ordered children                                        807
closed / failed children                           765 / 42
slice-closed roots                                       172.         (23)
```

Only `170` of those `172` are pure whole-route closures.  The two
difference bodies

```text
(2,3,6,10,12,13,14),       (2,6,7,8,9,10,11)            (24)
```

have an additional pair-exception row with failed top-four children.
The `E` route closes those extra rows, so all `172` belong to the full
composition.  The `42` failed children and `38` parent keys are exactly
the `Q` slice subsequently closed by THM-2913.

## 8. Exact lists, digests, and verification

The complete `1,219` live-route root list, `964` additive list, and
`1,315` current union are printed explicitly in the locked main output.
Selected body-list SHA-256 digests are

```text
live route 1,219
  b25f51db595f01bd6c24db6f6d25ed9a230f6c7abf10841a1ef4aacdb793f371
additive 964
  4b2dcf6945aa80e8896f22115df5096f028dce5d9936b3148b0083c259657254
current union 1,315
  47c0b23646ae5744f5354d6475aa23283c275341e6665ba5532640cddea0c41f.
```

The exact open-pivot key digest is

```text
7da02bd0669380e5e9e2554b4f38118fefd99db6f5058c59b449d360816185c9.
```

The full ledger records the parent mass and component count, centre rank
and coverage, full earlier-centre prefix, parent margin, direct child mass
and components, `q_4(2000)`, initial and final required horizons, scanned
label count, exact tail bound and gap, global top four, strict margin, and
closure bit for every child.  Its semantic digest is

```text
274114167a4e173d242fd0d62b980593df5fffb6adfc9cf116ec171fe628b1ff.
```

The first full replay produced ledger file SHA-256

```text
798cd660ab60e2021b28074a1390af3f6b1367c99f2d0ab63a581513f7871071.
```

Final source/main hashes and the ordinary/optimized replay attestation
will be inserted after the locked clean replay.  All guards are explicit
and remain active under optimized Python.  The row order is deterministic
and independent of worker scheduling.

The proof object is a ranked allocation tree with symmetric coverage
sidecars.  Its tie-break is a gauge, not an intrinsic tournament on
runners.  The theorem claims only the `964` additive roots in `(2)`; it
does not close any of the `4,866` failed children, the remaining `2,117`
roots, or LRC(14). ∎
