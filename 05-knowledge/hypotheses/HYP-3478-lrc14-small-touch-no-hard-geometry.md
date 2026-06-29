# HYP-3478: LRC14 Small-Touch / No-Hard Geometry

**Status:** SUPERSEDED RESERVATION / filled by the HYP-3478 geometry atlas;
not a proof.

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

This hypothesis reserved a geometry audit for those six rows.  The active
evidence file is now
`05-knowledge/hypotheses/HYP-3478-lrc14-small-touch-geometry-atlas.md`, with
script
`04-computation/lrc14_small_touch_geometry_atlas_codex_20260629.py` and stored
output
`05-knowledge/results/lrc14_small_touch_geometry_atlas_codex_20260629.out`.
The goal was to describe what their zero-edge singleton dead-cover projections
actually are:
component intervals, mirror pairing, blocker-owner pairs, adjacent touching
gate types, unit-delta versus cover-delta sidecar status, and any local
owner-current/two-adic/signed-SPEC proof payload suggested by that geometry.

## Filled Test

The HYP-3478 atlas used the HYP-3450/HYP-3451/HYP-3453/HYP-3438 row bank and
the HYP-3472/HYP-3476 current utilities to compute, for each of the six rows:

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

The atlas confirms that the six rows reduce to mirror-paired singleton
components.  Five rows have one mirror pair and `random_covering_001` has two.
All dead-cover projections have zero edges, all hard-orbit counts are zero, and
every pocket has at least two complete E/branch gate touches.  The next proof
target is a finite singleton-current lemma for the four branch-unit rows
(`001`, `062`, `086`, `101`), followed by the cover-delta sidecar clause for
`039` and `074`.

## Pointers

HYP-3478, HYP-3477, HYP-3476, HYP-3475, HYP-3472, HYP-3471, HYP-3455,
HYP-3453, HYP-3451, HYP-3450, HYP-3438, THM-523, T1438, LTI-438, LTT-338,
OPEN-Q-108.
