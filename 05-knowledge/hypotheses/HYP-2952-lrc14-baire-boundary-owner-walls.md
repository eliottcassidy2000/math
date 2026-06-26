---
id: HYP-2952
title: Baire boundary witnesses are C27 owner walls in path-planning form
status: COMPUTATIONAL SCOUT / boundary-witness transfer; not a proof
source: codex-2026-06-24
tags: [lrc14, baire, boundary, c27, owner-wall, path-planning, farey, unital]
related:
  - HYP-2951
  - HYP-2950
  - HYP-2949
  - HYP-2947
  - HYP-2943
  - HYP-2937
  - HYP-2248
results:
  - 04-computation/lrc14_borel_baire_haar_path_codex.py
  - 05-knowledge/results/lrc14_borel_baire_haar_path_codex.out
external:
  - https://en.wikipedia.org/wiki/Baire_measure
  - https://en.wikipedia.org/wiki/Haar_measure
---

# HYP-2952: Baire boundary witnesses as owner walls

HYP-2950 separates robust Baire interior from boundary-only witnesses.
HYP-2951 turns that split into a sixth any-angle search carrier.  HYP-2952 is
the direct LRC14 transfer:

```text
Baire boundary cells are the path-planning analogue of C27 owner walls.
```

## Claim

In a grid planner, a line-of-sight interval changes only when it crosses a
named obstacle tangency.  In the LRC14 branch atlas, a witness interval changes
only when it crosses a named arithmetic wall:

```text
Farey parent/child wall,
divisibility threshold,
C27 shell owner wall,
AP/Goddyn-Wong transfer wall,
q=3 unital pair-reuse wall,
K33 near-miss wall,
PH bad-child rank wall.
```

So the proof should not treat boundary witnesses as leftovers after measure
estimates.  It should promote them to first-class labelled walls.

## Finite Toy

The scout compares two size-20 subsets of `C_80`:

```text
contiguous arc: robust=18, boundary=2
dust:           robust=0,  boundary=20
```

They have identical Haar mass.  The difference is entirely categorical and
boundary-theoretic.

For LRC14 this is the warning:

```text
two residual packets can have the same invariant mass
while only one has a robust perturbation class.
```

The packet with no robust core is not automatically irrelevant.  It may be the
exact equality wall that decides tightness.

## Relation To C27 And q=3 Unitals

At `n=14`, the recurrent shell constant is:

```text
C = 2n - 1 = 27 = 3^3.
```

That makes C27 owner labels and q=3 unital blocks natural boundary charts, but
only after AP/Goddyn-Wong labels are attached.  HYP-2949 already warns that one
unital chart cannot freely recombine incompatible branch packets because
`lambda=1` pair uniqueness creates a reuse obstruction.

HYP-2952 restates this as a boundary rule:

```text
Do not cross a C27 owner wall without paying a chart-change cost.
Do not merge boundary-only packets unless their Borel branch labels agree.
Do not discard a measure-zero wall until its exact equality case is checked.
```

## Candidate Lemma

```text
After AP/Goddyn-Wong and C27 labels are attached, every boundary-only LRC14
survivor is either
  (a) an exact equality witness in the tight branch,
  (b) a Farey-child near-miss such as the known det=-1 migration pattern,
  (c) blocked by q=3 unital pair uniqueness,
  (d) or forced to drop PH-style bad-child rank.
```

This is not yet an LRC14 proof.  It is a proposed local theorem schema for
turning "measure-zero exception" into a finite list of named walls.
