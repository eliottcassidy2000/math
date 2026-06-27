---
id: HYP-2951
title: Borel-Baire-Haar A* as a sixth any-angle proof carrier for LRC14
status: COMPUTATIONAL SCOUT / creative algorithmic analogy; not a proof
source: codex-2026-06-24
tags: [lrc14, any-angle, path-planning, astar, theta-star, anya, cwave, haar, baire]
related:
  - HYP-2950
  - HYP-2952
  - HYP-2949
  - HYP-2947
  - HYP-2943
  - HYP-2248
results:
  - 04-computation/lrc14_borel_baire_haar_path_codex.py
  - 05-knowledge/results/lrc14_borel_baire_haar_path_codex.out
external:
  - https://en.wikipedia.org/wiki/Haar_measure
  - https://en.wikipedia.org/wiki/Baire_measure
---

# HYP-2951: Borel-Baire-Haar A*

The any-angle path-planning list gives five useful proof-interface types:

```text
Field D*:      edge-propagated interpolation values
Theta*:        parent line-of-sight shortcuts
Block A*:      local all-pairs block packets
ANYA:          taut path intervals as nodes
CWave:         geometric wavefront arcs and lines
```

The LRC14 analogue should use a sixth, deliberately more proof-aware carrier:

```text
Borel-Baire-Haar A*
```

## Proposed Node

A node is not a grid vertex.  It is a visible interval packet:

```text
node =
  visible interval(s)
  + Borel event code
  + Baire robust core
  + boundary/tangency cells
  + Haar mass
```

Expansion splits intervals only at structurally named walls:

```text
obstacle tangencies,
Farey child boundaries,
AP/Goddyn-Wong label changes,
C27 owner walls,
unital pair-reuse walls,
PH bad-child rank drops.
```

Priority is ordinary path length plus a proof penalty:

```text
cost = geometric path cost
     + boundary-only progress penalty
     + owner-wall crossing penalty
     + unresolved branch-chart penalty.
```

The goal is not to produce a better robot planner.  The goal is to name the
proof carrier that retains the same information a successful LRC14 proof seems
to need.

## Scout Output

The finite direction-circle model found:

```text
blocked directions         size=19 mass=19/80 robust=11 boundary=8
clear line-of-sight        size=61 mass=61/80 robust=53 boundary=8
clear robust core          size=53 mass=53/80 robust=45 boundary=8
```

The clear line-of-sight components are interval nodes:

```text
lengths: 23, 13, 12, 7, 6
```

This is the same structural move as ANYA, but with extra labels:

```text
ANYA interval
  -> Borel-Baire-Haar interval
  -> LRC Farey/C27/tournament witness interval.
```

## Tournament Analysis

Using the relation

```text
x -> y iff x retains more proof-critical witness information,
```

the scout ranks the carriers:

```text
exact LRC M/Farey/C27 labels
> Borel event code
> Baire robust-core/boundary split
> Haar invariant mass
> Borel-Baire-Haar A* interval node
> ANYA taut interval
> CWave wavefront primitive
> Theta* parent visibility bit
> Field D* interpolation scalar
> raw grid vertex
```

The ranking is a guardrail, not a theorem.  It says that the best analogy is
not "run A* over LRC states."  The useful analogy is:

```text
replace sampled vertices by proof-labelled visible intervals.
```

## LRC14 Use

For LRC14, a Borel-Baire-Haar A* search would try to prove:

```text
Starting from the apex node 1/14, every visible interval packet either
  stays in the AP/Goddyn-Wong tight branch,
  crosses a named C27/Farey owner wall,
  drops PH-style bad-child rank,
  or lands on an exact Baire-boundary equality case.
```

This is a plausible organizer for the loose/tight dichotomy:

```text
tight = deepest sink pinned to the named apex interval;
loose = marked node migrates to a coarser divisor node or to a Farey child.
```
