---
id: HYP-3845
title: The polygon-partition theorem (owner's statement) = Davenport-Mirsky-Newman-Rado on Z/N -- PROVED and FORMALIZED sorry-free in Lean 4 / mathlib style; the discrete shadow of the tower floor's Fraenkel rigidity
status: CONFIRMED + FORMALIZED (lake build clean: 0 errors, 0 warnings, 0 sorries)
source: klein-2026-07-01-S89
script: 04-computation/lean/TournamentH7/TournamentH7/PolygonPartitionDMNR.lean
related:
  - mac-mini THM-594(C) (continuous Mirsky-Newman -- the R/Z twin)
  - opus HYP-3902 F1 (tiling rigidity), HYP-3901 (D_7(k/7)=0 tiling times)
  - mac-mini HYP-3850(a) (Mirsky-Newman floor)
  - 03-artifacts/drafts/LRC14-formalization-status-2026-06-22.md (the Lean corpus this joins)
---

# HYP-3845: the polygon-partition theorem, machine-checked

## Statement (owner's form)

If the vertices of a regular N-gon are colored so that every color class itself forms the
vertices of a regular polygon, then two color classes are congruent.

## Decoding = DMNR

A regular m-gon inscribed with vertices among the N-gon's vertices is exactly a residue
class {k < N : k = a mod q}, q = N/m | N. "All classes pairwise non-congruent" = a
partition of Z/N into residue classes with PAIRWISE DISTINCT moduli -- forbidden by the
Davenport-Mirsky-Newman-Rado theorem. Proof: evaluate sum_{k<N} w^k = 0 at w = a primitive
root of unity of order the LARGEST modulus q_max (>= 2 once there are >= 2 classes and
moduli are distinct); each smaller-modulus class sums to a full geometric cycle = 0
(w^{q_i} != 1 since q_max does not divide q_i < q_max, and (w^{q_i})^{N/q_i} = w^N = 1);
the largest-modulus class contributes w^{a} * (N/q_max) != 0. Contradiction.

## Formalization (the mathlib seed)

`PolygonPartition.exists_eq_modulus` in PolygonPartitionDMNR.lean: for q i | N, 0 < q i,
2 <= card, and the exact-cover hypothesis (∀ k < N, ∃! i, k % q i = a i % q i), there
exist i != j with q i = q j. Supporting lemmas: `residueClass_eq_image` (the class is the
image of range (N/q) under j -> a%q + j*q) and `residueClass_expSum` (geometric
evaluation). **Builds sorry-free against mathlib v4.30.0** (~170 lines). Notes for an
actual mathlib PR: trim `import Mathlib` to the needed modules; consider restating over
`ZMod N`; the exact-cover hypothesis could be packaged as `Finpartition`.

## Why it matters here

This is the discrete Z/N shadow of the tower floor's key input: the cluster density
vanishes exactly at arc-TILINGS of the circle (opus D_7(k/7) = 0), and Mirsky-Newman
rigidity says distinct-speed systems cannot tile except on the degenerate locus --
continuous version proved by mac-mini (THM-594(C)), Z/N version now MACHINE-CHECKED here.
The pair gives the tower floor's rigidity input at both ends of the discretization, and
the formalized version is the first piece of the LRC14 corpus in submission-ready shape.
