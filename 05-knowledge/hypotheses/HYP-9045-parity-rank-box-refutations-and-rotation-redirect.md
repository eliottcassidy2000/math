---
id: HYP-9045
title: "Parity and rank-box refutations, with rotation and archimedean redirects"
status: >
  RESOLVED NEGATIVES / OPEN REDIRECTS. Two proposed mechanisms are ruled out
  at their stated scope by elementary algebra and exact finite referees. The
  replacement rotation and Euler-characteristic programs remain hypotheses.
source: mac-mini-2026-07-27-S143
related:
  - HYP-9040-cancellation-completeness-synthesis
  - THM-2052
  - THM-2337
  - THM-2356
  - LEM-020-redei-involution-parity-layers
---

# HYP-9045 -- parity/rank-box negative verdicts and redirects

## 1. Mod-two forcing is the wrong coefficient prime

The mixed `91` line lives on odd-order `7 x 13` character data. Multiplication
by two is invertible on this coefficient group. Hence the usual averaging
contraction kills positive-degree `C_2` cohomology: a mirror or palindrome
involution cannot force this line by a mod-two Bockstein or a two-torsion
fixed-layer class. This is a structural obstruction to the mechanism proposed
in HYP-9040, not a proof that every possible parity-flavoured count vanishes.

**Open redirect.** Use the native `C_7 x C_13` rotations instead. Translations
by `1/7` and `1/13` act freely off degenerate words, so a correctly typed
equivariant count should reduce to its fixed/degenerate contribution modulo
seven or thirteen. The decisive test is the rotation-equivariance of
THM-2356's target--jet pairing with all word and ancestry labels retained.

## 2. The rank-twelve box is archimedean, not a second CRT channel

For an all-small height vector `v`, the admissible pair relations
`v_j e_i-v_i e_j` span `v^perp`; exact AP, geometric-progression, deep-well,
and random controls reproduce rank twelve. Modulo `p=7,13`, zero coordinates
do not create the hoped-for sparse-rank drop: a support-three recovery restores
the missing basis direction, and the recorded hostile families retain rank
twelve. Reproduce with:

```bash
python3 04-computation/lrc14_rank12_box_anatomy_macmini_S143.py
```

The matching stored output is
`05-knowledge/results/lrc14_rank12_box_anatomy_macmini_S143.out`. Thus the
rank-twelve box is the all-small height regime and is invisible to the proposed
CRT rank test. Its shallow shells belong to the existing plateau/strata
analysis; at depth it rejoins the spectral `91`-line problem.

## 3. Euler-ledger correction

The first fast fibre stratification proposed under HYP-9040 is **REFUTED** at
degenerate strata: it undercounts the `(0,1,0)` control relative to the exact
Groebner computation. The previously checked discriminant-cube and individual
Euler-characteristic calculations survive, but the displayed global balance
does not until an exact multiplicity-aware count map is supplied.

## Scope

These are negative mechanism verdicts and redirects, not LRC(14) or planar-JC
closures. No live LRC row is removed, and no Euler argument proves or refutes
the planar Jacobian conjecture.
