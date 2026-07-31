---
id: THM-2980
title: "Projected k2 band 1600--1679 exact closure"
status: >
  PROVED + FINITE-EXACT + VERIFIED + INDEPENDENTLY AUDITED.  In the
  projected k=2 scalar necessary sector inherited from THM-2941, every row
  with 1600<=z_1<=1679 is empty.  A fresh monolithic 3,003-body by 80-height
  atlas has 68 candidate rows, exactly the disjoint set closed by twelve
  typed component referees; there are no missing, extra, or duplicate rows.
  Consequently the projected k=2 first-drift cap is z_1<=1599.  This does
  not close the sector below 1600, any other k sector, the six-body/seven-tail
  rung, or LRC(14).
source: codex-lrc14-k2-band-1600-1679-closure-2026-07-30
depends_on:
  - THM-2941-critical-seven-slot-scalar-wall-and-balanced-boundary
  - THM-2970-located-low-cell-torsion-tail-closure
  - THM-2972-all-four-high-punctured-packet-torsion-closure
  - THM-2984-projected-k3-signed-ray-attainment-and-unit-phase-gate
related:
  - THM-2995-projected-k2-band-1580-1599-translated-capacity-closure
verification:
  - 04-computation/lrc14_j7_k2_1600_1679_integrated_promotion_audit_thm2980.py
  - 05-knowledge/results/lrc14_j7_k2_1600_1679_integrated_promotion_audit_thm2980.out
---

# THM-2980 -- projected k2 band 1600--1679 exact closure

**PROVED + FINITE-EXACT + VERIFIED + INDEPENDENTLY AUDITED.**

## Statement and universe

Use the lossless projected `k=2` scalar necessary sector and first-drift
coordinate `z_1` of THM-2941.  The body universe is the `3,003` six-subsets
of `{1,...,14}`.  Later nonaligned labels remain ordered and distinct; no
finite label horizon is introduced.

Every scalar row in

```text
1600 <= z_1 <= 1679
```

is empty.  Therefore the first drift in this projected necessary sector
satisfies

```text
z_1 <= 1599.                                                   (1)
```

This is not a physical-cover census.  In particular, (1) does not close
`z_1<=1599`, another alignment sector, the full six-body/seven-tail rung, or
LRC(14).

## Complete atlas

The promotion referee reruns the scalar atlas as one `80`-height block on
all bodies.  Its universe has exactly

```text
3,003 * 80 = 240,240
```

body-height rows.  Exactly `68` survive the scalar screens.  Their complete
height distribution is

```text
1604:1, 1608:2, 1612:19, 1616:1, 1618:1, 1620:8,
1625:1, 1630:1, 1631:1, 1642:1, 1645:1, 1650:12,
1656:5, 1660:1, 1665:1, 1668:9, 1670:2, 1672:1.        (2)
```

All other heights in the band are empty already at the atlas stage.  Of the
rows in (2), `28` carry the forced-high-tail type; those rows are retained,
not silently passed through an all-label routine.

## Exact closure and set join

Twelve component referees close the rows in descending handoff order.  Their
case counts are

```text
3, 10, 1, 5, 12, 2, 2, 1, 8, 2, 19, 3,
```

whose sum is `68`.  Their typed terminal partition is

```text
EMPTY                     40
SCALAR-EMPTY              22
LITERAL-CLOSURE-BELOW       1
LITERAL-TORSION-EMPTY       1
SECTION-CLOSURE-BELOW       4.                         (3)
```

The integrated audit pins the source and output hash of every component,
parses only certified terminal `CASE` records of the types in (3), and joins
their `(z_1,E)` keys.  It then compares that join with the independently
rerun monolithic atlas.  The exact comparison is

```text
monolithic atlas rows       68
component terminal rows     68
duplicates                   0
missing                      0
extra                        0
upper/lower overlap          0.                        (4)
```

Thus every possible row is either absent from (2) or occurs exactly once in
a proved terminal component.  This proves the band statement and (1).

## Hostile controls and audit boundary

The component transcripts retain positive and zero forced-high scalar
controls.  A row enters (4) only after its declared scalar, literal packet,
torsion, projected-cell, or complete-section certificate has succeeded.
This prevents a positive scalar near miss from being reclassified as empty
merely because a later handoff cap passed it.

The normal and optimized promotion transcripts are byte-identical.  The
monolithic atlas profile, atlas-survivor set, component key set, component
manifest, and semantic payload have separate frozen SHA-256 digests in the
verification output.  This audit proves completeness of the finite band; it
does not promote any assertion about the newly reserved band below `1600`.
