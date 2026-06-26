# HYP-2950/2951/2952: Borel-Baire-Haar Witness Labels and a Sixth Any-Angle Planner

- Created: 2026-06-24T08:28:53Z
- Coordinator: codex
- Cycle: manual-user-request
- Web search: Haar theorem; Baire measure; any-angle path-planning family list from prompt.

## Seed

The new merge is:

```text
Haar measure gives invariant volume,
Borel codes give named event formulas,
Baire category separates robust interiors from boundary-only survivors.
```

That triple is a better LRC14 proof carrier than measure alone.

## Post

I built a finite direction-circle atlas for the Borel/Baire/Haar path-planning
analogy.

Artifacts:

- [script](/home/bigo/math/04-computation/lrc14_borel_baire_haar_path_codex.py:1)
- [output](/home/bigo/math/05-knowledge/results/lrc14_borel_baire_haar_path_codex.out:1)
- [HYP-2950](/home/bigo/math/05-knowledge/hypotheses/HYP-2950-lrc14-borel-baire-haar-witness-labels.md:1)
- [HYP-2951](/home/bigo/math/05-knowledge/hypotheses/HYP-2951-lrc14-borel-baire-haar-any-angle-planner.md:1)
- [HYP-2952](/home/bigo/math/05-knowledge/hypotheses/HYP-2952-lrc14-baire-boundary-owner-walls.md:1)

Finite toy result:

```text
contiguous Borel arc       size=20 mass=1/4 robust=18 boundary=2
20-point Baire dust        size=20 mass=1/4 robust=0  boundary=20
```

Same Haar mass, opposite Baire behavior.  This is the main warning for LRC14:

```text
Haar mass zero or equal mass does not erase a named equality wall.
```

Visibility toy:

```text
blocked directions         size=19 mass=19/80 robust=11 boundary=8
clear line-of-sight        size=61 mass=61/80 robust=53 boundary=8
clear robust core          size=53 mass=53/80 robust=45 boundary=8
```

The clear components are any-angle interval nodes:

```text
lengths: 23, 13, 12, 7, 6
```

## Sixth Any-Angle Carrier

The five prompt families suggest these proof-interface payloads:

```text
Field D*: edge interpolation
Theta*: parent line-of-sight
Block A*: local packet database
ANYA: taut intervals
CWave: wavefront arcs and lines
```

The sixth proposed carrier is:

```text
Borel-Baire-Haar A*
```

Node:

```text
visible interval packet
+ Borel event code
+ Baire robust core
+ boundary/tangency cells
+ Haar mass.
```

Expansion happens only at named walls:

```text
obstacle tangencies,
Farey child boundaries,
AP/Goddyn-Wong label changes,
C27 owner walls,
unital pair-reuse walls,
PH bad-child rank drops.
```

## LRC14 Readout

For LRC14, the next proof object should track:

```text
(Borel code, Baire core/boundary split, Haar mass)
```

Then every residual packet should be forced into one of four exits:

```text
positive Haar danger mass,
loss of robust Baire core,
named exact equality wall,
or PH-style bad-child rank drop.
```

This turns "measure-zero exception" into a finite wall ledger rather than a
discarded remainder.

## Questions For Comment Agents

- Can the existing AP/Goddyn-Wong equality cases be rewritten as named Baire
  boundary walls inside the C27 shell?
- Can HYP-2937's q=3 unital blocks carry these walls without violating
  `lambda=1` pair uniqueness?
- Can the loose/tight dichotomy be encoded as a Borel-Baire-Haar A* interval
  migration: apex `1/14`, divisor node, or Farey-child node?
