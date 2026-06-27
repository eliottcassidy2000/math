---
id: HYP-3119
title: LRC14 tournament edge witness recursion
status: RESERVED / edge-witness synthesis and executable scout pending; not a proof
source: codex-2026-06-27-S268
tangent: T1194
technique: LTI-255
tournament_technique: LTT-153
script: 04-computation/lrc14_tournament_edge_witness_recursion_codex_s268.py
result: 05-knowledge/results/lrc14_tournament_edge_witness_recursion_codex_s268.out
related:
  - HYP-3118
  - HYP-3117
  - HYP-3116
  - HYP-3115
  - HYP-3112
  - HYP-3108
  - HYP-3106
  - HYP-3105
  - HYP-3054
  - HYP-3050
  - HYP-3049
  - HYP-2963
  - OPEN-Q-108
---

# HYP-3119: LRC14 Tournament Edge Witness Recursion

## Reservation Claim

This lane reserves the S268 response to the edge-as-witness prompt.  The
working object is not a tournament vertex and not a scalar edge count.  It is
an oriented edge `tail -> tip` whose proof value is the pair of recursive
witness obligations seen from both ends:

```text
edge_witness(tail -> tip) =
  (tail_packet,
   tip_packet,
   outside_four_sector_deck,
   tail_deletion_child,
   tip_deletion_child,
   repaired_coordinate_or_named_debt)
```

The existing sources already point at this shape:

- HYP-3050 says an edge is naturally dualistic: it has a tail and a tip, and
  outside vertices split into the four sector words around that directed edge.
- HYP-3054 names `edge_tail_tip_sector_word` as an observer-extension cut
  coordinate.
- HYP-3106 turns directed-edge sectors into a controlled-forgetting functor.
- HYP-3112 and HYP-3115 expose proof-relevant edge boundaries through
  ear-payload edges and Lee-Yang/Ising domain-wall edges.
- HYP-3116/HYP-3118 say any such edge witness is legal only if the destroyed
  coordinate is either retained, resurrected by a sidecar, or routed to named
  residual debt.

Planned computation: enumerate small tournaments, compute each directed edge's
four-sector deck, recursively compare the tail-deletion and tip-deletion
children, and run Tournament Analysis on edge-witness reframes rather than
runners, roots, scalar edges, or proof-gate names.

Challenged assumption: an edge is not just a relation between two vertices.
For LRC-style proof search it is a bidirectional proof obligation: the tail
asks what survives after pushing forward, while the tip asks what survives
after pulling back.  A legal witness must make both recursions compatible.
