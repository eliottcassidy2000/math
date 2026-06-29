---
id: HYP-3523
title: LRC14 random031 certificate spigot
status: RESERVED / spigot-inspired terminal-dispatch experiment; not an LRC14 proof
source: codex-2026-06-29 spigot-algorithm reframe after HYP-3521/HYP-3520
tangent: T1523
technique: LTI-523
tournament_technique: LTT-423
script: 04-computation/lrc14_random031_certificate_spigot_codex_20260629.py
result: 05-knowledge/results/lrc14_random031_certificate_spigot_codex_20260629.out
reflection: 07-reflections/lrc14-random031-certificate-spigot-codex-20260629.md
related:
  - HYP-3521
  - HYP-3520
  - HYP-3512
  - HYP-3511
  - HYP-3510
  - HYP-3494
  - HYP-3493
  - HYP-3490
  - HYP-3486
  - HYP-3483
  - HYP-3034
  - OPEN-Q-108
---

# HYP-3523: Random031 Certificate Spigot

## Claim Reservation

The spigot-algorithm prompt suggests a proof discipline for
`random_covering_031`:

```text
mixed-radix state     -> terminal certificate ledger
reduce-and-carry      -> discharge plus named debt propagation
predigit buffer       -> terminal packet waiting for possible boundary carry
safe digit emission   -> proof certificate whose future debt is bounded
tail error            -> unprocessed owner/free-hole/quotient debt
```

This packet reserves an executable test of that analogy.  The target is a
certificate spigot over HYP-3521's `77` terminal certificates: ordinary
rank-2 routes should emit immediately, free-hole packets should emit after
their bracket/doublet carry is resolved, and the pure bypass should emit only
with the HYP-3520 owner-boundary sidecar
`45:+1,147:+1,169:+1,173:+1`.

## Missing

The missing finite object is a streaming audit that records, at each emitted
certificate, the current carry buffer:

```text
free_hole_carry
owner_boundary_carry
private_firewall_guard
vertical_halfturn_guard
```

The experiment should test whether HYP-3521 already gives a no-backtracking
terminal spigot, or whether some certificate remains a predigit that cannot be
safely emitted without a new sidecar.

## Assumption Challenge

Tournament vertices should be proof-output states and carry buffers, not
runners, arcs, or scalar counts.  The preserved predicate is terminal
certificate emission without later invalidation.  Destroyed information should
be named as carry/debt rather than silently dropped.
