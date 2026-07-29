---
id: THM-2920
title: "Two-H3-row pair-Hunter recursive-toothpick closure"
status: >
  PROVED + FINITE-EXACT + VERIFIED.  All 367 ordered children left by
  THM-2916 close: 149 by exact pair partition, 111 more by the four-slot
  Hunter envelope, 102 by finite second-centre top-three recursion, and
  the final five by a lawful three-slot pair/Hunter endpoint.  All 296
  surviving two-H3-row roots are additive.  The proved union becomes
  1,041, the residual becomes 2,391, and the complete 690-root two-H3
  stratum is closed.  This is not LRC(14).
source: codex-post2913-residual-scout-2026-07-29
depends_on:
  - THM-735-bonferroni-simultaneous-multi-peel-defeats-the-clustered-non-isolated-wall
  - THM-2893-complement-cap-finite-core-flag-lemma
  - THM-2897-partition-cap-tropical-convolution-and-alternating-pair-ladder
  - THM-2913-one-h3-row-pair-hunter-toothpick-closure
  - THM-2916-two-h3-row-dynamic-tail-child-top-four-closure
related:
  - THM-2915-all-open-centre-exact-child-top-four-closure
verification:
  - 04-computation/lrc14_j6_two_h3_row_pair_hunter_recursive_toothpick_closure_codex_20260729.py
  - 05-knowledge/results/lrc14_j6_two_h3_row_pair_hunter_recursive_toothpick_closure_codex_20260729.out
  - 05-knowledge/results/lrc14_j6_two_h3_row_pair_hunter_recursive_toothpick_closure_codex_20260729.ledger.out
---

# THM-2920 -- two-H3-row pair-Hunter recursive-toothpick closure

**PROVED + FINITE-EXACT + VERIFIED.**

## 1. Statement

THM-2916 closes every ordered first-centre child on `394` of the `690`
roots having exactly two ordinary H3 rows.  It leaves `367` children on
the other `296` roots.  Every one of those children closes:

```text
exact pair partition                                      149
additional four-slot Hunter envelope                      111
finite second-centre top-three recursion                  102
three-slot recursive pair/Hunter endpoint                   5.        (1)
```

At whole-root level the cumulative route counts are

```text
pair only                  102
pair or Hunter             195
through first recursion    292
through three-slot G3      296.                           (2)
```

All `296` roots lie outside the `745`-root union through THM-2916.
Consequently

```text
proved union                                      745+296=1,041
finite residual                           3,432-1,041=2,391.          (3)
```

Together with the `394` THM-2916 roots, this closes the complete
`690`-root two-H3-row stratum.

## 2. Inherited child and exact pair cap

For each THM-2916 failure retain the full ordered state

```text
E       seven-body root,
a       marked apex,
P       inherited excluded prefix,
x       allocated first hostile centre,
A_x     all earlier first centres,
R       C minus D_x,
F       P union {x} union A_x.                            (4)
```

The verifier reconstructs `R` both by literal subtraction from the parent
carrier `C` and directly from `E union {a,x}`.  It checks mass, components,
prefix, and ordered exclusions before using the child.

Let `h=mu(R)` and let `beta_2(R,F)` be the exact largest union of two
distinct allowed danger sets.  As in THM-2913, the computation scans
allowed labels through `2,500`, seals the singleton ranks by THM-735, and
evaluates every finite pair whose singleton-sum upper bound can beat the
current winner.  A finite winner strictly beats

```text
q_1 + [h/7+(99/70)r/(7*2501)],                            (5)
```

so any pair containing an omitted label is smaller.  Thus the reported
pair cap is global, not a search-horizon maximum.

Partitioning four possible remaining labels into two pairs shows that

```text
2 beta_2(R,F)<h                                            (6)
```

closes `149` children.  The computation scans `911,490` allowed singleton
labels, pays `3,011` exact pair unions, and has minimum strict pair-tail
gap

```text
31927313083/13373652271980.                               (7)
```

The smallest positive margin in `(6)` is `1/40040`.

## 3. Four-slot Hunter repair and first toothpick

Write the globally sealed singleton ranks as
`q_1>=q_2>=q_3>=q_4`.  The four-slot Hunter envelope is

```text
G4=max_(0<=u<=min(q_1,beta_2))
   [u+sum_(j=2)^4 min(u,q_j,beta_2-u)].                   (8)
```

Its exact finite breakpoint set is the one in THM-2913.  The strict
inequality `G4<h` closes another `111` children; its smallest positive
margin is

```text
307/19684665.                                             (9)
```

For every remaining child, the first crossing of the objective in `(8)`
gives a hostile singleton threshold `lambda>h/7`.  THM-735 therefore
reduces all possible second centres to a finite exact core.  The verifier
scans `29,548` labels to build those cores.  They have sizes `1` through
`5`, with histogram

```text
size             1   2   3   4   5
children        17  33  32  21   4.                      (10)
```

The earliest-maximum allocation of a second centre `y` enlarges the
forbidden sidecar to

```text
F_y=F union {y} union A_y,                                (11)
```

where `A_y` contains all earlier centres in its finite core.  Direct and
sequential reconstructions of `R minus D_y` agree.  Its global top three
singleton coverages are sealed at horizon `2,500`.  Among `283` ordered
second pivots, `278` close strictly by their top-three sum.  This closes
`102` whole children.  The smallest positive terminal margin is

```text
643567/367447080,                                         (12)
```

and the minimum second-level tail gap is

```text
2110444943/375374939940.                                  (13)
```

## 4. Lawful three-slot endpoint

Exactly five children retain one open second pivot.  At this point two
centres have been fixed, so **three** labels remain.  A bare top-two sum is
therefore not a terminal certificate; computing a pair statistic does not
consume a label.  This is the slot-count correction recorded in
`MISTAKE-324`.

For each of the five twice-subtracted carriers, the verifier recomputes
the exact global pair cap `beta_2`, the global top three
`q_1>=q_2>=q_3`, and the three-slot Hunter envelope

```text
G3=max_(0<=u<=min(q_1,beta_2))
   [u+sum_(j=2)^3 min(u,q_j,beta_2-u)].                   (14)
```

The simpler `q_1+beta_2<h` closes three tips.  It fails on two, with
margins `-543079/116035920` and `-191/504504`.  The lawful envelope
`G3<h` closes all five, with exact strict margins

```text
4211/280280, 2459243/348107760, 52873/2802800,
143831/12892880, 7009/420420.                             (15)
```

Thus three tips use the pair-plus-singleton specialization and two require
genuine Hunter compatibility.  The minimum deep pair-tail gap is

```text
27321799/12267154900.                                     (16)
```

No terminal equality occurs at any of the pair, `G4`, first-toothpick,
pair-plus-singleton, or `G3` stages.

## 5. Identity, recomposition, and scope

The verifier rebuilds all `11,842` parent rows, selects exactly `1,380`
rows on the `690` two-H3 roots, and reconstructs all `5,618` open
ordered-centre children.  Its `5,251/367` THM-2916 partition agrees
exactly with the pinned source artifacts.  The sorted additive
`296`-root tuple has digest

```text
e3045198e08804c78025bd532111377309882911e08bc50604aa7119ac266c71,
```

and the complete sorted `690`-root two-H3 tuple has digest

```text
772f67e5711cb009012a9c4abeb1b9a288195126f382244231f1f93362b63efc.
```

The prefix and earlier-centre data in `(4)` and `(11)` are proof-bearing
sidecars, not quotient decoration.  Every inequality is sufficient for
noncoverage; a failed intermediate route would not be a cover witness.
This theorem closes one exact residual stratum.  It does not treat roots
with three or more ordinary H3 rows and does not prove LRC(14).

## 6. Verification

The script hash-pins THM-2913 and THM-2916 sources and outputs, checks the
complete identity joins, performs direct carrier reconstruction at every
depth, seals every infinite tail, records all exact witnesses and margins,
and recomposes roots from full parent keys.

Ordinary and optimized eight-worker runs produce byte-identical summaries
and ledgers.  The semantic ledger digest is

```text
e1434329c46118991b1fe357be5f87a01a22b81a814d760166ba3582de79e83a.
```

LF-normalized SHA-256 values are

```text
source  049be1e331fb0c4fc46e703ffb37d61be2b3ec3b781835f480e55ba89bdf894e
output  1a38fd441dfd77a4f5d30d45d3160febc33d2d4eeb6247b223f10a1e31a8aefb
ledger  9c96a24c90c07c69d96b86f29355f9a86599c9d9174e0260dd952aa357f7d7f1.
```

Reproduction:

```bash
python3 04-computation/lrc14_j6_two_h3_row_pair_hunter_recursive_toothpick_closure_codex_20260729.py --workers 8
python3 -O 04-computation/lrc14_j6_two_h3_row_pair_hunter_recursive_toothpick_closure_codex_20260729.py --workers 8
```
