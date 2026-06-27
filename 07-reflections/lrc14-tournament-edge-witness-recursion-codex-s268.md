# LRC14 Tournament Edge-Witness Recursion Reflection

Anchors: HYP-3124, T1198, LTI-259, LTT-157, OPEN-Q-108.

The useful shift in this pass is that an edge is not just a relation between
two tournament vertices.  It is a two-ended proof obligation.  The tail asks
what survives after deleting or pushing past the tail endpoint; the tip asks
the dual question after deleting or pulling back from the tip endpoint.  A
legal witness should make those two recursive children compatible and name any
coordinate lost in the quotient.

The S268 census makes the distinction concrete.  Four-sector words are
equinumerosity/equidistribution data: they count and distribute outside
vertices around `tail -> tip`.  They are useful but too coarse.  Through
labelled tournaments `n<=5`, sector counts are `1,4,10,20`, while sector plus
paired endpoint-deletion child signatures are `1,4,16,80`; at `n=5`, every
sector group splits by child pair.  The paired child object is the
equidecomposability datum: it asks whether the local edge can be cut into two
recursively compatible proof pieces.

The resulting invariant is:

```text
edge_witness_certificate =
  four_sector_deck
  + paired_endpoint_deletion_recursion
  + repair_sidecar_or_named_debt
```

Tournament vertices in this pass were edge-witness reframes and repair
schemas, not runners, roots, spins, domain walls, or scalar edge counts.  The
selected path was:

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

The creative next test is not to rank more edge statistics.  It is to attach
the edge-witness packet to HYP-3115 one-swap/domain-wall edges, then classify
each wall as an observer-gluing discharge, a smaller tail/tip recursion, or
named HYP-2963/HYP-3098 residual debt.  That would merge the edge-as-witness
idea with the circuit missing-input and coordinate-resurrection guardrails
rather than adding another scalar analogy.
