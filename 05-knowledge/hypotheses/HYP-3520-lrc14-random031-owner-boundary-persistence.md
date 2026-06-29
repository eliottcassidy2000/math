---
id: HYP-3520
title: LRC14 random031 owner-boundary persistence
status: RESERVED / in-progress owner-cobordism certificate; not an LRC14 proof
source: codex-2026-06-29 owner-boundary persistence pass after HYP-3494/HYP-3493/HYP-3510/HYP-3511
tangent: T1520
technique: LTI-520
tournament_technique: LTT-420
script: 04-computation/lrc14_random031_owner_boundary_persistence_codex_20260629.py
result: 05-knowledge/results/lrc14_random031_owner_boundary_persistence_codex_20260629.out
reflection: 07-reflections/lrc14-random031-owner-boundary-persistence-codex-20260629.md
related:
  - HYP-3494
  - HYP-3493
  - HYP-3511
  - HYP-3510
  - HYP-3490
  - HYP-3486
  - HYP-3485
  - HYP-3483
  - HYP-3482
  - HYP-3481
  - HYP-3402
  - HYP-3034
  - OPEN-Q-108
---

# HYP-3520: Random031 Owner-Boundary Persistence

## Claim Reservation

This packet reserves the owner-boundary persistence lane requested after
HYP-3494.  The current exact random031 data already says:

```text
hard_components = (43,54)
forbidden_seam_owners = (23,45,93,113,147,169,173)
pure_bypass_flow_owners = (23,93,113)
owner_boundary_debt = (45,147,169,173)
pure_bypass_cells = 12
```

The intended certificate is not a proof that the seam is harmless.  It is the
finite statement that the lower-delta bypass is a flow carrier on the same two
hard components, while the four missing labels persist as a boundary word.  In
this reading, the hard pair is a forbidden seam and the complement carries the
phase flow; owner labels are relative boundary charge.

## Known

HYP-3493 already gives a relative seam-sheaf table with `79` legal
horizontal+mirror components, no mixed/debt stalks, and one pure bypass stalk.
That stalk has size `12`, branch split `6/6`, hit components `(43,54)`, owners
`(23,93,113)`, and seam debt `(45,147,169,173)`.

HYP-3510 gives the coarse connected phase-carrier side: all `282` phase
witnesses live in the seam-complement transport object after the max-delta seam
is deleted, and the lower-delta bypass is hit exactly `12` times.

HYP-3511 removes the main free-hole ambiguity by bracketing all `40` free-hole
cells against ordinary endpoint-rank-2 packets.

## Missing

The missing finite object is a quotient-price matrix showing which sidecars are
necessary to reconstruct the owner-boundary word.  A scalar count, component
pair, phase block, or bypass owner word alone is expected to forget the
seam-only owners.  The proposed exact sidecar is the owner-current cobordism
word:

```text
forbidden_seam_owner_word - pure_bypass_flow_owner_word
  = (45,147,169,173).
```

## Why Reserve This Namespace

HYP-3494 made owner-boundary persistence the live proof target.  This packet
will either promote it to a deterministic finite certificate or name the first
missing sidecar/debt.  The tournament vertices should be proof sidecars and
quotient choices, not runners or arcs.
