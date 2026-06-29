# HYP-3478: LRC14 Small-Touch / No-Hard Geometry

**Status:** RESERVED STUB / computation in progress; not a proof.

**Claimed by:** codex-2026-06-29  
**Tangent:** T1438  
**Technique handles:** LTI-438, LTT-338  
**Open question:** OPEN-Q-108

## Question

HYP-3476 leaves six non-AP boundary-current exceptions that have no HYP-3475
hard mirror orbit and no E/branch pair-current edge support:

```text
random_covering_001
random_covering_039
random_covering_062
random_covering_074
random_covering_086
random_covering_101
```

This hypothesis reserves a geometry audit for those six rows.  The goal is to
describe what their zero-edge singleton dead-cover projections actually are:
component intervals, mirror pairing, blocker-owner pairs, adjacent touching
gate types, unit-delta versus cover-delta sidecar status, and any local
owner-current/two-adic/signed-SPEC proof payload suggested by that geometry.

## Planned Test

Use the HYP-3450/HYP-3451/HYP-3453/HYP-3438 row bank and the HYP-3472/HYP-3476
current utilities to compute, for each of the six rows:

- dead component intervals and minimal branch-cover labels;
- mirror pairing of singleton dead components and owner-pair residues;
- baseline dead-cover projection size and edge support;
- touching low-rank E/branch gates, endpoint kinds, branch masks, adjacency,
  and cover-delta vectors;
- the intersection with the colored gate unit-delta split
  (`random_covering_039`, `random_covering_074` expected as cover-delta
  sidecar rows; the other four expected as unit-delta singleton-current rows).

## Assumption Challenge

Do not assume the tournament vertices are runners, arcs, or raw row names.
This audit will compare alternate vertex sets: singleton dead components,
mirror-paired dead components, blocker-owner pairs, fixed circle intervals,
section boundaries, touching gate events, residues, cover arcs, Fourier/color
sidecars, and proof obligations.  The quotient should preserve the LRC
predicate "which terminal discharge packet remains after the low-rank
E/branch producer" while explicitly recording what it destroys: interval
order, branch orientation, endpoint wall labels, and owner-current locality.

## Current Pull

If the six rows all reduce to balanced mirror-paired singleton components,
then the next proof target is a finite singleton-current lemma: mirror-balanced
dead components with only small adjacent E/branch touches must discharge by
bounded owner-current, endpoint-spine, two-adic descent, exact-period,
signed-SPEC/Rprime, or state-lift debt.  If the two cover-delta rows
(`random_covering_039`, `random_covering_074`) have a distinct geometry, split
them from the four clean unit-delta rows before formalizing the terminal
packet.

## Pointers

HYP-3478, HYP-3477, HYP-3476, HYP-3475, HYP-3472, HYP-3471, HYP-3455,
HYP-3453, HYP-3451, HYP-3450, HYP-3438, THM-523, T1438, LTI-438, LTT-338,
OPEN-Q-108.
