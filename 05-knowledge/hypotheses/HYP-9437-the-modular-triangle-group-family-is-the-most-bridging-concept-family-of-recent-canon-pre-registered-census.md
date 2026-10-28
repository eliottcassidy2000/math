---
id: HYP-9437
title: "Pre-registered keyword census: in recent canon (THM-4400..4839) the modular/triangle-group family (PSL_2(Z), Farey, Berggren, Pythagorean, Pell, modular forms, triangle groups) is the most bridging of ten concept families against a degree-preserving null model"
status: TESTED 2026-10-28 - P1 CONFIRMED, P2 REFUTED, P3 REFUTED (the modular family is pervasive but is not the bridging hub)
session: opus-2026-10-28-S68
date: 2026-10-28
related: [THM-4839, THM-4837, THM-4838, THM-4808, THM-4825]
note: 05-knowledge/results/attractors_modular_triangle_world_abc_genealogy_20261028.md
script: 04-computation/experiments/attractor_census_canon_20261028.py
---

# HYP-9437 - which concept family bridges the repository?

## Why

The owner asked which "larger unifying theories lurk as attractors" for the repository's novel results. Concurrent
lanes propose methodological unifiers:
- session84: guarded observers;
- codex-incidence: guarded attachment algebras;
- session72: operation-relative representation.

They explicitly decline any dynamical "attractor". This session's mathematical thesis is that the main **object-level**
attractor is the world of the modular group PSL_2(Z) = Delta(2,3,oo) with:
- its congruence and triangle-group quotients (PSL_2(F_p) = Delta(2,3,p) quotients; the Euclidean groups Delta(2,4,4) and
  Delta(2,3,6) of the Gaussian and Eisenstein grids);
- its trees (Farey, Stern-Brocot, Berggren via Gamma(2));
- its Euler characteristics, which are Bernoulli values (THM-4838).

A thesis like this is easy to tell as a story. The census below is a falsifiable check of one consequence: if the
family is a hub, theorems that mention it should co-occur with many other families more often than a length-controlled
null model allows.

## Protocol (fixed before running)

- **Corpus:** every `01-canon/theorems/THM-<id>-*.md` with 4400 <= id <= 4839, excluding RESERVED stubs. The window
  covers roughly September-October 2026. Each file is read for its first 6000 characters.
- **Families:** ten regex families: MODULAR, QUADRATIC, FINITE_GEOMETRY, BINARY_ALGEBRA, DIVISOR, COLLATZ, TOURNAMENT,
  ZETA, ABC, OBSERVER. The exact patterns are in the script, committed in this same commit.
- **Null model:** curveball randomisation of the file-family incidence, which preserves both degree sequences. 2000
  samples with seed 20261028 and 5N swaps between samples.
- **Statistic:** z_ij for each of the 45 pairs. The lift-degree of a family is the number of partners with z >= 2.

## Predictions

- **(P1)** MODULAR touches at least 20% of the files in the window.
- **(P2)** MODULAR has the largest lift-degree of the ten families; ties count as success.
- **(P3)** At least one of the three highest-z pairs lies inside {MODULAR, QUADRATIC, FINITE_GEOMETRY}, the
  arithmetic-geometry core.

## What the outcome means

- P2 failing means that, on this crude lexical measure, some other family (for instance BINARY_ALGEBRA or OBSERVER)
  bridges more. The note would then report that family as the census attractor and weaken its own thesis accordingly.
- Passing is weak evidence only. Keyword co-occurrence is not mathematical dependence; the decisive content of the
  thesis lies in the proved bridges of THM-4839 and the cited dictionary in the note.

The full canon (all ids) is reported as an exploratory, non-registered comparison.

## Results (confirmatory run, 2026-10-28, after the pre-registration commit 593105aedd)

`python3 attractor_census_canon_20261028.py 4400 4839 2000` gives 320 files, 2.82 families per file on average, and 81.6%
of files touching at least two families.

- **(P1) CONFIRMED.** MODULAR touches 53.8% of the files. Only BINARY_ALGEBRA, at 57.5%, is more frequent.
- **(P2) REFUTED.** Lift-degrees (partners with z >= 2): TOURNAMENT 2; FINITE_GEOMETRY, BINARY_ALGEBRA, MODULAR and
  COLLATZ 1 each; the rest 0.
- **(P3) REFUTED.** The top three pairs are FINITE_GEOMETRY-TOURNAMENT (z = 7.29), COLLATZ-MODULAR (3.31) and
  BINARY_ALGEBRA-TOURNAMENT (2.06). Inside the arithmetic-geometry core the best pair is FINITE_GEOMETRY-MODULAR at
  only z = 1.20.
- **Strongest negative pairs:** COLLATZ-OBSERVER (-5.78), COLLATZ-TOURNAMENT (-3.25) and MODULAR-OBSERVER (-3.10).
  OBSERVER has the least positive-z strength (1.19).

**Exploratory full-canon comparison** (not registered; 3957 files, 300 null samples;
`attractor_census_canon_20261028_fullcanon.out`):
- FINITE_GEOMETRY-TOURNAMENT is again the strongest pair (z = 13.36).
- Lift-degree leaders are BINARY_ALGEBRA 5, DIVISOR 4, MODULAR 3 and COLLATZ 3.
- OBSERVER is again the least bridging.

**Reading.** The session's first thesis, the modular family as the hub, is refuted on this measure. What survives is
two object-level hubs:
- finite geometry with tournaments;
- Collatz with modular arithmetic.

They share a common binary language and sit beside an observer layer that is lexically isolated. THM-4839 supplies
explicit morphisms across the weak FINITE_GEOMETRY-MODULAR gap. The census itself is lexical evidence only.
