# LRC14 Random031 Gate-Gluing Obstruction

HYP-3439 closed the AP/84m bridge in the expected direction: canonical rank
`6` is the `m=1` corridor-fence case, AP tails drop to the rank-`5` core, and
HYP-3452 gives the endpoint phase.  The broad-bank variant exposed the useful
negative result: rank `6` is not globally canonical.  `random_covering_031` is
the lone noncanonical rank-`6` exception in the critical bridge set.

HYP-3455 makes that exception finite.  The row has rescue core
`(23,45,93,113,147,169)`, a connected `15`-edge overlap graph, and exactly two
maximal survivor gates with total cover delta `7`.  Those gates are components
`43` and `54`; they are mirror partners, have word `B-S-B-S-B`, and have
delta `(3,4)` / `(4,3)`.  Their adjacent cover-owner union is
`(23,45,93,113,147,169,173)`, so all six rescue owners appear in the hard
gate-pair payload and `173` is the only extra owner.

This changes the proof obligation.  The noncanonical high-rank obstruction is
not a scalar family to average away.  It is a seven-owner local clause:

```text
connected rank-6 rescue graph
+ mirror max-delta survivor gate pair
+ branch-coloured component-cover escape router
```

The component-cover side is favorable: `random_covering_031` has `98`
components, only `4` dead-both components, `94` low-rank escapes, maximum
dead-pair rank `2`, and danger score `0.000868`.  A failed proof would have to
explain why the two hard gates can block every low-rank escape despite the
tiny dead-cover projection.

Proof routes that look strongest:

1. Menger/Green route.  Treat the two max-delta survivor gates as terminals in
   the branch-coloured blocker graph.  Prove any attempt to saturate them cuts
   off or exposes one of the `94` low-rank component escapes.
2. Schwarz-Christoffel route.  Treat the gate endpoints as labelled
   prevertices.  The mirror pair supplies boundary order; the cover deltas
   `(3,4)` and `(4,3)` are the accessory-parameter payload.
3. Owner-current route.  The extra owner `173` should become a signed current
   imbalance or a removable accessory owner.  The rescue owners themselves are
   already forced by the hard gate-pair union.
4. Signed-SPEC/BDH route.  Averaging can enter only after the seven-owner
   local clause is retained; raw rank, survivor mass, or dead fraction would
   forget the proof object.

Tournament Analysis: vertices are proof obligations and finite gate/gluing
carriers, not runners.  The pairwise observable retains rank-`6` graph
connectivity, max-delta gate payload, two-colour escape count, mirror
involution, and scalar-firewall risk.  The priority path is:

```text
max_delta_survivor_gate_pair
-> rank6_rescue_overlap_graph
-> two_color_component_escape_router
-> owner_delta_gluing_clause
-> mirror_involution_cut_word
-> endpoint_spine_wall_lift
-> bdh_signed_spec_sidecar
-> raw_rank6_scalar
```

Assumption challenge: considered runners, rescue owners, gaps, survivor gates,
fixed circle sections, section boundaries, wall events, residues, cover arcs,
endpoint labels, component-cover graph nodes, Fourier modes, matroid circuits,
and proof obligations.  Chosen vertices are the rank-`6` rescue graph, the two
max-delta survivor gates, and component-cover escape nodes.  This preserves
the branch-relocation predicate while destroying raw runner order and scalar
mass.  The challenged assumption is that noncanonical rank `6` means an
uncontrolled family; here it is a finite mirror-symmetric gate-gluing clause.
