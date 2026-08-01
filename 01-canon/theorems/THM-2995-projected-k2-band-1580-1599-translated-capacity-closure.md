---
id: THM-2995
title: "Projected k2 band 1580--1599 translated-capacity closure"
status: >
  PROVED + FINITE-EXACT + VERIFIED + INDEPENDENTLY AUDITED.  In the
  projected k=2 scalar necessary sector inherited from THM-2941, every row
  with 1580<=z_1<=1599 is empty.  A fresh monolithic 3,003-body by 20-height
  atlas has 26 candidate rows, exactly the set closed by five typed component
  referees; there are no missing, extra, or duplicate rows.  Together with
  THM-2980, the projected k=2 first-drift cap is z_1<=1579.  This does not
  close the sector below 1580, any other k sector, the six-body/seven-tail
  rung, or LRC(14).
source: codex-lrc14-k2-band-1580-1599-closure-2026-07-31
depends_on:
  - THM-2941-critical-seven-slot-scalar-wall-and-balanced-boundary
  - THM-2980-projected-k2-1600-1679-zero-high-and-high-order-phase-closure
  - THM-2984-projected-k3-signed-ray-attainment-and-unit-phase-gate
related:
  - THM-2980-projected-k2-1600-1679-zero-high-and-high-order-phase-closure
verification:
  - 04-computation/lrc14_j7_k2_1580_1599_integrated_promotion_audit_thm2995.py
  - 05-knowledge/results/lrc14_j7_k2_1580_1599_integrated_promotion_audit_thm2995.out
---

# THM-2995 -- projected k2 band 1580--1599 translated-capacity closure

**PROVED + FINITE-EXACT + VERIFIED + INDEPENDENTLY AUDITED.**

## Statement and universe

Use the lossless projected `k=2` scalar necessary sector and first-drift
coordinate `z_1` of THM-2941.  The body universe is the `3,003` six-subsets
of `{1,...,14}`.  Later nonaligned labels remain ordered and distinct; no
finite label horizon is introduced.

Every scalar row in

```text
1580 <= z_1 <= 1599
```

is empty.  Combining this band with THM-2980 gives

```text
z_1 <= 1579.                                                   (1)
```

This is a closure of a projected necessary sector, not a physical-cover
census.  In particular, (1) does not close `z_1<=1579`, another alignment
sector, the full six-body/seven-tail rung, or LRC(14).

## Complete atlas

The promotion referee reruns the scalar atlas as one `20`-height block on
all bodies.  Its universe has exactly

```text
3,003 * 20 = 60,060
```

body-height rows.  Exactly `26` survive the scalar screens, with complete
height distribution

```text
1581:3, 1586:11, 1588:1, 1590:3,
1594:2, 1595:2, 1599:4.                              (2)
```

All other heights in the band are empty already at the atlas stage.  Of the
rows in (2), `17` carry the forced-high-tail type; those rows are retained
rather than passed through an inapplicable all-label routine.

## Exact closure and set join

Five component referees close the rows in descending handoff order.  Their
case counts are

```text
4, 4, 3, 12, 3,
```

whose sum is `26`.  Their typed terminal partition is

```text
EMPTY                       9
SCALAR-EMPTY               13
SECTION-CLOSURE-BELOW       4.                       (3)
```

The four section closures are exact translated complete-section witnesses,
not positive scalar rows relabelled as empty.  The integrated audit pins the
source and output hash of every component, parses only certified terminal
`CASE` records of the types in (3), and joins their `(z_1,E)` keys.  It then
compares that join with both the pinned atlas and the freshly rerun monolithic
atlas.  The exact comparison is

```text
monolithic atlas rows       26
pinned atlas rows           26
component terminal rows     26
duplicates                   0
missing                      0
extra                        0.                       (4)
```

The component handoff caps are

```text
1598, 1593, 1589, 1585, 1579.                       (5)
```

Thus every possible row is either absent from (2) or occurs exactly once in
a proved terminal component.  Equations (4)--(5), together with THM-2980,
prove (1).

## Hostile controls and audit boundary

The four nonpositive exact one-high cases are matched one-for-one with the
four zero-failure translated complete-section witnesses.  Every other case
record has a certified `EMPTY` or `SCALAR-EMPTY` conclusion.  This guards the
key failure mode in which a zero scalar supremum is mistaken for packet
emptiness without checking translated capacity.

The monolithic atlas profile, atlas-survivor set, component manifest,
component key set, and semantic payload have separate frozen SHA-256 digests
in the verification output.  The audit is independent at the
atlas-versus-component set-join level; it is not claimed to be an independent
implementation of every component argument.
