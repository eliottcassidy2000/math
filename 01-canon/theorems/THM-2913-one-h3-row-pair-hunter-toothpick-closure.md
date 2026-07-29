---
id: THM-2913
title: "One-H3-row pair-Hunter toothpick closure"
status: >
  PROVED + FINITE-EXACT + VERIFIED.  The 38 roots and 42 ordered children
  left open by the THM-2912 child top-four census all close: 22 by exact
  pair partition, 10 more by the four-slot Hunter envelope, and the final
  10 through 29/29 strict terminal pivots on finite hostile
  second-centre cores.  Exact recomposition adds 37 roots over the live
  THM-2911/2912 union, taking the proved union to 351 and the residual to
  3,081.
source: h4-exception-children-2026-07-29
depends_on:
  - THM-735-bonferroni-simultaneous-multi-peel-defeats-the-clustered-non-isolated-wall
  - THM-2893-complement-cap-finite-core-flag-lemma
  - THM-2897-partition-cap-tropical-convolution-and-alternating-pair-ladder
  - THM-2904-all-root-ranked-h1-hunter-closure
  - THM-2911-all-root-finite-ranked-h1-hunter-closure
  - THM-2912-one-h3-row-exact-child-top-four-closure
related:
  - THM-2901-all-hard-exact-global-pair-cap-and-route-partition
  - THM-2907-paircap-exception-h4-heavy-link-child-census
verification:
  - 04-computation/lrc14_j6_one_h3_row_pair_hunter_toothpick_closure_codex_20260729.py
  - 05-knowledge/results/lrc14_j6_one_h3_row_pair_hunter_toothpick_closure_codex_20260729.out
  - 05-knowledge/results/lrc14_j6_one_h3_row_pair_hunter_toothpick_closure_codex_20260729.ledger.out
  - 05-knowledge/results/lrc14_j6_one_h3_row_pair_hunter_toothpick_closure_codex_20260729.replay.out
---

# THM-2913 -- one-H3-row pair-Hunter toothpick closure

**PROVED + FINITE-EXACT + VERIFIED.**

## 1. Statement

THM-2912 leaves exactly `42` ordered four-cover children, lying on `38`
of its `210` one-H3-row bodies.  Every one of those `42` children closes:

```text
exact child pair partition, 2 beta_2<h                    22
additional four-slot Hunter-star closures                  10
finite hostile second-centre children                      10
second-centre pivots, all top-three closed               29/29.       (1)
```

At whole-root level the first route closes `19` of the `38` bodies, its
union with the Hunter route closes `28`, and the recursive route closes
the final `10`.  Thus all `38` THM-2912-open bodies close.  Exactly one,

```text
(1,2,3,5,6,8,13),                                         (2)
```

already belongs to the THM-2911 proved union.  The other `37` are
additive over the live THM-2911 plus THM-2912 union.  Consequently the
exact recomposition is

```text
proved union through THM-2911                             181
after the 172 THM-2912 roots                              314
after the 37 additive roots here                          351
official residual                            3,432-351=3,081.          (3)
```

Equivalently, the complete `210`-root one-H3-row stratum meets the
THM-2911 union in `40` roots and adds `170`.

## 2. The inherited marked child

Fix one THM-2912-open branch.  Its theorem-bearing object is not merely an
interval carrier.  It retains

```text
C       the THM-2904 parent carrier,
x       its allocated first hostile centre,
P       its inherited excluded gate prefix,
A_x     all earlier hostile centres,
R=C minus D_x,
h=mu(R),
F=P union {x} union A_x.                                  (4)
```

Every remaining label is at least `15` and lies outside `F`.  THM-2912
already seals its exact singleton ranks

```text
q_1>=q_2>=q_3>=q_4.                                       (5)
```

The ordered sidecar `F` is load-bearing: reintroducing an earlier centre
would merge different proof branches.

## 3. Exact global child pair caps

Let `r` be the number of interval components of `R`.  THM-735 gives, for
every allowed label `w`,

```text
c_R(w)<h/7+gamma/w,             gamma=(99/70)r/7.          (6)
```

The verifier evaluates every allowed label through `2,500`.  On all `42`
rows the fourth finite singleton rank is at least the omitted-label
majorant

```text
tau=h/7+gamma/2501,                                       (7)
```

so `(5)`, and in particular `q_1`, is globally sealed.  It then computes
the largest literal union `H_2` of two distinct finite-head danger sets,
using singleton-sum branch-and-bound followed by exact residual
subtraction.  Any pair with an endpoint beyond `2,500` has union strictly
below

```text
q_1+tau.                                                  (8)
```

Every finite winner satisfies `H_2>q_1+tau`; hence

```text
beta_2(R)=H_2                                             (9)
```

is the attained exact global allowed pair cap.  The computation scans
`104,320` singleton labels and pays only `342` exact pair unions.  Its
smallest self-sealing gap `H_2-(q_1+tau)` is

```text
871202831/264970545840.                                   (10)
```

Partitioning four putative cover labels into two pairs now gives

```text
mu(R intersect D_Q)<=2 beta_2(R).                         (11)
```

The strict inequality `2 beta_2<h` closes `22` children.  Its smallest
positive margin is

```text
17047/29099070
```

at body `(1,3,6,8,9,11,13)`, apex `21`, first centre `20`.
The closest pair-partition failure is

```text
-103/144144
```

at body `(2,3,4,5,6,10,12)`, apex `22`, first centre `16`.

## 4. The four-slot Hunter envelope

The pair partition forgets that two separately maximizing pairs may be
incompatible.  Let a putative four-cover be ordered by decreasing
singleton coverage, and write `a` for its maximum.  For any other member
of rank `j`, its contribution outside the maximum-centre danger set is at
most each of

```text
a,             q_j,             beta_2-a.                (12)
```

The first bound uses maximality of `a`, the second uses the exact global
rank, and the third uses the exact pair cap.  Therefore Hunter's star
inequality gives

```text
G_4=max_(0<=a<=min(q_1,beta_2))
       [a+sum_(j=2)^4 min(a,q_j,beta_2-a)],               (13)
```

and `G_4<h` excludes a four-cover.  The objective is continuous,
piecewise affine, and concave.  Its complete breakpoint set is

```text
0, min(q_1,beta_2), beta_2/2,
q_j, beta_2-q_j                         (j=2,3,4),        (14)
```

clipped to the legal interval, so `(13)` is maximized exactly.

Among the `20` pair-partition failures, `(13)` closes `10`.  The smallest
positive Hunter margin is

```text
17/21021
```

at body `(1,2,3,5,6,9,10)`, apex `24`, centre `22`.  The
closest remaining failure is

```text
-49/163020
```

at body `(1,2,3,5,7,9,10)`, apex `24`, centre `39`.
Failure is only failure of a sufficient envelope, not evidence for a
cover.

## 5. Hostile second centres: the toothpick recursion

On each of the `10` remaining children, let

```text
lambda=min{a:G_4(a)>=h},                                  (15)
```

where `G_4(a)` denotes the bracketed expression in `(13)`.  The verifier
solves the first affine crossing exactly and finds

```text
lambda-h/7>0                                              (16)
```

on every row.  Hence the maximum singleton in any hypothetical four-cover
belongs to the finite core

```text
H_1(R)={w not in F:c_R(w)>=lambda}.                        (17)
```

By `(6)`, every such label satisfies

```text
w<=ceil(gamma/(lambda-h/7))-1.                            (18)
```

The exact cutoffs range from `241` to `409`; the smallest gap in `(16)`
is

```text
1335417/90250160.                                         (19)
```

The ten actual core sizes are

```text
2,2,2,2,3,3,5,4,4,2,                                    (20)
```

for `29` second-centre pivots in total.

Order each core by decreasing coverage, breaking ties by increasing
label.  Allocate a hypothetical cover to the earliest one of its
maximum-coverage labels `y`.  After replacing `R` by

```text
S=R minus D_y,                                            (21)
```

the other three labels avoid `F`, `y`, and every earlier second-centre
label.  This is the same ordered-centre operation as the preceding layer,
now applied one tooth lower.  The unmarked carrier `S` alone would not
remember these exclusions.

## 6. Exact three-slot terminal

For every one of the `29` marked grandchildren `(21)`, the verifier scans
all allowed labels through `2,500`.  The THM-735 tail lies at or below the
third retained singleton rank, with minimum gap

```text
19740543491/2998793637840.                                (22)
```

Thus the exact global ranks `p_1>=p_2>=p_3` are known.  Subadditivity
gives, for any remaining three labels,

```text
mu(S intersect D_Q)<=p_1+p_2+p_3.                         (23)
```

Every pivot satisfies the strict terminal inequality

```text
p_1+p_2+p_3<mu(S).                                        (24)
```

The computation scans `71,974` singleton labels at this stage.  The
smallest margin in `(24)` is

```text
127333/19339320
```

at body `(1,3,9,10,11,12,13)`, apex `21`, first centre
`23`, second centre `24`.  Thus all `29` pivots, all `10` recursive
children, and finally all `42` THM-2912-open children close.

## 7. Explicit closed-root list

The `38` closed bodies are

```text
(1,2,3,4,5,7,9)       (1,2,3,4,5,9,10)
(1,2,3,4,5,9,14)      (1,2,3,4,7,9,11)
(1,2,3,4,9,10,11)     (1,2,3,5,6,8,13)
(1,2,3,5,6,9,10)      (1,2,3,5,6,9,11)
(1,2,3,5,6,11,12)     (1,2,3,5,7,9,10)
(1,2,3,6,7,8,14)      (1,2,3,6,7,12,13)
(1,2,3,7,9,11,14)     (1,2,4,5,6,7,10)
(1,2,4,8,9,11,12)     (1,2,5,6,7,8,10)
(1,2,5,6,8,9,10)      (1,2,5,6,8,10,12)
(1,2,5,6,8,10,13)     (1,2,5,6,9,10,12)
(1,2,5,7,8,10,11)     (1,3,4,5,6,9,11)
(1,3,4,5,8,9,11)      (1,3,4,5,9,12,13)
(1,3,4,6,9,12,13)     (1,3,5,6,7,10,11)
(1,3,5,6,9,10,12)     (1,3,5,6,9,12,13)
(1,3,5,9,10,11,14)    (1,3,6,8,9,11,13)
(1,3,9,10,11,12,13)   (1,4,5,6,9,10,13)
(2,3,4,5,6,10,12)     (2,3,5,6,8,13,14)
(3,4,5,6,9,10,11)     (3,4,5,6,11,12,13)
(3,4,5,7,9,10,14)     (3,5,6,9,10,11,12).              (25)
```

The SHA-256 digest of the sorted tuple in `(25)` is

```text
62f177cca38cfa9dc8a1aa33c9842efa61ab7c9a4270cc024f248cb3ee5fd37b.
```

Removing the overlap `(2)` leaves the explicit `37`-root additive list
printed in the locked output, with digest

```text
35dd2a9391ec0ba4e4cc95821382dd1ca8f08dc2571ca1c9ae6d25855beb8f54.
```

## 8. Equality, validity, and scope

There are no equality cases at any theorem-bearing terminal:

```text
2 beta_2=h                 0
G_4=h                      0
p_1+p_2+p_3=mu(S)          0.                             (26)
```

Every child and grandchild is reconstructed both by literal subtraction
and directly from its full family; pair winners receive a repeated literal
residual check.  Every forbidden prefix is carried forward.  Pair endpoints
and cover labels are distinct.  The global pair cap and both singleton rank
banks include equality at their finite heads and use the strict discrepancy
estimate only to seal omitted labels.  Strict inequalities, rather than
decimal comparisons, certify `(1)`.

The relation in `(12)` is an undirected pair-overlap sidecar.  The order
used in `(17)` is a branch-allocation gauge and does not define a
tournament on runners.  “Toothpick” names the self-similar marked
operation

```text
(carrier, forbidden prefix, k slots)
  -> (literal residual, enlarged forbidden prefix, k-1 slots),        (27)
```

not a new combinatorial object or an unproved all-depth principle.

This theorem concerns only the `42` children and `38` roots in
`(25)`.  Even after `(3)`, `3,081` seven-body roots remain outside the
proved union.  It does not prove LRC(14).

## 9. Verification

The verifier hash-pins the THM-2912 child source/output and the THM-2911
verifier/output, explicitly normalizing CRLF to LF and rejecting lone
carriage returns so checkout policy cannot alter a text dependency.  The
final dependency pins are

```text
THM-2912 source
  d2810560a7d002d7eeadecc6a50a7733c90585527295aa5e85e72775739b839b
THM-2912 output
  454d87c8beeb81405b031cce4b40bdda0d385cfcd9c48e6fcf4eb810cfc00c5a
THM-2911 verifier
  e0ac67539f7ff09376645a62beef0a9ac7d0766a2e749666f94d1fd4d6487b15
THM-2911 output
  e5c58cc2eb325928c00839c2593450ea7cce8945b3835898ec83c6c5f42fac9b.
```

It reconstructs the full
`11,842`-row THM-2904 universe, the `210` one-H3-row bodies, all `807`
ordered children, and then the exact `42`-child residual.  The semantic
ledger digest is

```text
236df90a5497eb85b6536b9cbfebc1e8a7b30990fe9bd24d3a890d748c3bc428.
```

The final run was performed in a forced-CRLF worktree.  Ordinary
eight-worker, optimized eight-worker, and ordinary serial executions
produced byte-identical main and ledger outputs.  LF-normalized SHA-256
hashes are

```text
source
  14e56e124197cd1bdae841efa195a1e7c282e7ea368a610e5f4d56509431858b
main output
  3604644a9691b13e7fa245249b68c9586ec2775996834f7761f32eb0b89f3e64
full ledger
  5419b87511bbf51c43a2bed9647e82cb5178ad99ac8f667d0e66318caa049632.
```
