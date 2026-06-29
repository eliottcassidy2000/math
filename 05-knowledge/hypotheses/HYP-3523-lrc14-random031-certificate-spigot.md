---
id: HYP-3523
title: LRC14 random031 certificate spigot
status: EVIDENCE / streaming terminal-certificate carry audit; not an LRC14 proof
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

## Claim

The spigot-algorithm prompt suggests a proof discipline for
`random_covering_031`:

```text
mixed-radix state     -> terminal certificate ledger
reduce-and-carry      -> discharge plus named debt propagation
predigit buffer       -> terminal packet waiting for possible boundary carry
safe digit emission   -> proof certificate whose future debt is bounded
tail error            -> unprocessed owner/free-hole/quotient debt
```

The executable test supports the analogy.  HYP-3521's `79` legal component
events stream in branch/u order to `77` emitted terminal certificates covering
all `282` cells.  Ordinary rank-2 route components emit immediately.  The
`10` free-hole singles emit after ordinary bracket carry.  The two free-hole
doublets each hold one predigit event and then emit one collapsed certificate
when the mate arrives.  The unique pure bypass emits only with HYP-3520 and
HYP-3522 owner-filtration sidecars.

Key readout:

```text
component_events=79
emitted_certificates=77
component_cells_covered=282
emitted_certificate_cells=282
held_predigit_events=2
held_predigit_buffer_cells=4
untyped_pending_carry=0
typed_residual_carry=(45,173)
```

The bypass carry is:

```text
seam_owners=(23,45,93,113,147,169,173)
transport_word=(23,93,113)
branch_boundary_owners=(23,93,147,169)
bracket_lift=(147,169)
residual_after_branch_boundary=(45,173)
```

## Interpretation

This reframe changes the random031 terminal packet from a static ledger into a
no-backtracking proof procedure.  The hard pair is not treated as a wall.  It
is a forbidden seam whose complement carries phase flow, while the certificate
stream records the delayed proof carry required to emit across nearby
features.

There is no untyped terminal carry left in the stream.  What remains is typed:
formalize ordinary route emissions, HYP-3511 doublet buffering, HYP-3522
owner-filtration carry through transport `(23,93,113)` and bracket lift
`(147,169)`, and then discharge the residual pair `(45,173)` under the
HYP-3490/HYP-3513 private-firewall route sidecar.  Every emitted certificate
also carries the vertical-halfturn guard; the vertical address projection is
not a legal sheet gluing.

## Micro Release Test

A second, smaller quotient test is useful for future proof compression.  The
safe bypass packet is the only carrier in the tested family that both emits
the certified prefix and keeps all typed carry:

```text
full_spigot_packet              score=87
owner_filtration_without_route  score=68  illegal_loss=HYP-3490 route
transport_plus_bracket          score=26  illegal_loss=residual_pair,route
HYP3520_owner_current_only      score=21  illegal_loss=branch_lift,route
terminal_class_only             score=-1  illegal_loss=owner,residual,route
raw_counts                      score=-21 illegal_loss=owner layers,route,free-hole type
raw_seven_owner_shadow          score=-28 illegal_loss=transport,bracket,residual,route
```

This is not a second proof object; it is a compression guardrail.  It says a
candidate quotient may forget raw phase, runner, or count head data only after
the emitted proof packet is certified and the residual/route tail remains
named.

## Proof Pull

1. Formalize ordinary route certificates as immediate spigot digits.
2. Formalize HYP-3511 doublet buffering: first half is a predigit, second half
   emits the doublet certificate.
3. Formalize bypass emission as owner-filtration carry: transport
   `(23,93,113)`, bracket lift `(147,169)`, residual `(45,173)`.
4. Attach HYP-3513 route sidecar `R` and the vertical guard to every emitted
   certificate before any quotient compression.

## Assumption Challenge

Tournament vertices are proof-output states and carry buffers, not runners,
arcs, raw witnesses, or scalar counts.  The preserved predicate is terminal
certificate emission without later invalidation.  Destroyed information is
named as carry/debt rather than silently dropped: branch/u order, terminal
class, bracket carry, owner carry, route guard, and vertical guard must be
present before scalar compression.
