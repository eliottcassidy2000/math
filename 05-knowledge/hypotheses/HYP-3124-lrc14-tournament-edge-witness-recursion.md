---
id: HYP-3124
title: LRC14 tournament edge witness recursion
status: EVIDENCE / executable edge-witness recursion scout plus S271 class-deck stress supplement; not a proof
source: codex-2026-06-27-S268; extended by codex-2026-06-27-S271
tangent: T1198
technique: LTI-259
tournament_technique: LTT-157
script: 04-computation/lrc14_tournament_edge_witness_recursion_codex_s268.py
result: 05-knowledge/results/lrc14_tournament_edge_witness_recursion_codex_s268.out
related:
  - HYP-3123
  - HYP-3122
  - HYP-3121
  - HYP-3120
  - HYP-3119
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

# HYP-3124: LRC14 Tournament Edge Witness Recursion

## Claim

The S268 edge-witness lane was renumbered after a rebase showed that
HYP-3119/T1194/LTI-255/LTT-153 belongs to the S269 niche archive bridge ledger
and that HYP-3121 is occupied by an incoming covering-case synthesis.  A later
rebase also showed HYP-3122 is the phi4 cap-stabilizer lane and HYP-3123 is
the chiral/Cech proof-angle lane, so the completed edge-witness scout now
uses HYP-3124/T1198/LTI-259/LTT-157.  The working object here is not a
tournament vertex and not a scalar edge count.  It is an oriented edge
`tail -> tip` whose proof value is the pair of recursive witness obligations
seen from both ends:

```text
edge_witness(tail -> tip) =
  (endpoint_role_word,
   outside_four_sector_deck,
   tail_deletion_child_signature,
   tip_deletion_child_signature,
   recursive_tail_child_edge_deck,
   recursive_tip_child_edge_deck,
   observer_gluing_payload_orbit,
   missing_input_vector,
   coordinate_resurrection_sidecar_or_named_debt,
   terminal_exit)
```

In LRC14 language, a local edge quotient is legal only if the packet route is
constant on the quotient fiber, or the lost tail/tip coordinate is
reconstructed, dual-annihilated, recursed into a smaller child, or named as
residual debt.

## Evidence From The S268 Scout

The script
`04-computation/lrc14_tournament_edge_witness_recursion_codex_s268.py` exhaustively
enumerates labelled tournaments through `n=5` and records each directed edge's
four-sector word, the one-sided endpoint-deletion children, the paired
tail-deletion/tip-deletion child signature, and one recursive layer inside both
children.  Stored output:
`05-knowledge/results/lrc14_tournament_edge_witness_recursion_codex_s268.out`.

Exact census:

| n | labelled tournaments | directed-edge instances | sector words | sector + tail child | sector + tip child | sector + child pair | depth-1 signatures | sector groups split by child pair |
|---|----------------------|-------------------------|--------------|---------------------|--------------------|---------------------|--------------------|------------------------------------|
| 2 | 2 | 2 | 1 | 1 | 1 | 1 | 1 | 0 |
| 3 | 8 | 24 | 4 | 4 | 4 | 4 | 4 | 0 |
| 4 | 64 | 384 | 10 | 14 | 14 | 16 | 16 | 6 |
| 5 | 1024 | 10240 | 20 | 52 | 52 | 80 | 80 | 20 |

For `n=5`, every sector group is split by the paired endpoint-deletion child
object.  The largest refinements are sector `(1,1,1,0)` and `(0,1,1,1)`, each
with `7` child-pair signatures and `960` labelled edge instances.  The four
sector word is therefore a strong local observable but not a witness by
itself; the first natural witness carrier is the sector word plus both
endpoint-deletion children.

This also separates the user's three equality lenses:

- Equinumerosity: the four sector sizes count how many outside vertices lie in
  each tail/tip relation class.
- Equidistribution: sector-size profiles and edge-deck histograms compare how
  those local classes are distributed across a tournament family.
- Equidecomposability: the paired tail-deletion and tip-deletion child decks
  ask whether the local edge can be cut into two recursively compatible proof
  pieces.  The S268 census says this is the first level that sees hidden
  witness distinctions invisible to sector counts.

Tournament Analysis over proof reframes uses edge-witness schemas as vertices,
not runners, roots, spins, or scalar edge counts.  Pairwise observable:
majority retention across LRC predicate, tail/tip symmetry, recursive closure,
sector resolution, observer gluing, missing-input repair, domain-wall
alignment, formal readiness, and scalar guardrail axes.  Fingerprint:
`score_hist={0:1,1:1,2:1,3:1,4:1,5:1,7:3,9:1,10:1}`,
`directed_3cycles=1`, `scc_sizes=[3,1,1,1,1,1,1,1,1]`,
`hamiltonian_path_count=3`, and `7` edge flips against a locality-first gauge.
Selected Hamiltonian path:

```text
coordinate_resurrection_edge_sheaf
-> edge_witness_two_ended_packet
-> cross_sector_orientation_word
-> paired_tail_tip_deletion_recursion
-> proof_circuit_edge_gate
-> domain_wall_edge_classifier
-> ear_payload_edge_mass
-> outside_four_sector_deck
-> tail_deleted_one_sided_recursion
-> tip_deleted_one_sided_recursion
-> raw_edge_count_scalar
```

## S271 Class-Deck Stress Supplement

The S271 follow-up script
`04-computation/lrc14_tournament_edge_witness_recursion_codex_20260627.py`
keeps the same directed-edge carrier but shifts the finite test from labelled
edge instances through `n=5` to unlabelled class decks through `n=6`.  Stored output:
`05-knowledge/results/lrc14_tournament_edge_witness_recursion_codex_20260627.out`.

Class-deck readout:

| n | classes | sector counts | sector internal | roleless children | recursive children | full edge witness |
|---|---------|---------------|-----------------|-------------------|--------------------|-------------------|
| 3 | 2 | 2/2 | 2/2 | 2/2 | 2/2 | 2/2 |
| 4 | 4 | 4/4 | 4/4 | 4/4 | 4/4 | 4/4 |
| 5 | 12 | 12/12 | 12/12 | 12/12 | 12/12 | 12/12 |
| 6 | 56 | 55/56 | 55/56 | 56/56 | 56/56 | 56/56 |

The only `n=6` sector-count/internal collision is the converse pair
`344/345` with score sequence `(2,2,2,3,3,3)`, `c3=8`, and `H=43`.  Retaining
the paired endpoint children repairs that collision.  At the directed-edge
instance level, all `43` nontrivial `n=6` sector-internal fibers split by
tail/tip deletion children, and `16` recursive fibers split further by the
cross-sector orientation word.  This says the cross-sector word is not
cosmetic; it is the remaining covariance/phase payload inside recursive edge
fibers.

The covering-floor translation is now carried by HYP-3125/HYP-3127 rather than
by this HYP-3124 witness card.  Read this supplement as evidence that the
tail/tip deletion children and cross-sector word are genuine information
coordinates before an `R-safe packet -> Q-safe packet` is compressed into the
multi-far `Rprime` floor.

## Integration With Existing Threads

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

The immediate proof-route consequence is:

```text
edge_witness_certificate =
  four_sector_deck
  + paired_endpoint_deletion_recursion
  + repair_sidecar_or_named_debt
```

Next executable test: attach this packet to HYP-3115 one-swap/domain-wall
edges and ask which walls become legal observer-gluing discharges, which
recurse to smaller tail/tip children, and which remain named HYP-2963/HYP-3098
debt.  When the same witness is used for covering-floor rows, use HYP-3125's
`edge_floor_packet` fields.

Challenged assumption: an edge is not just a relation between two vertices.
For LRC-style proof search it is a bidirectional proof obligation: the tail
asks what survives after pushing forward, while the tip asks what survives
after pulling back.  A legal witness must make both recursions compatible.
