# LRC14 Small-Touch / No-Hard Geometry

HYP-3478 makes the six HYP-3476 small-touch rows visible.  They are not a
vague random residual: they are mirror-balanced singleton dead-component
packets.

The exact audit finds `14` dead components across the six rows, with histogram
`{2:5, 4:1}`.  Every dead component has a singleton `B0`/`B1` owner pair, and
the interval mirror sends it to the component with the owners swapped.  Every
row has zero owner imbalance after summing `B0` minus `B1`, and every dead
projection has zero edges.

That explains why HYP-3476 pair-current failed: the projection is already a
set of isolated vertices, so no E/branch pair can remove projection edges.  It
also explains why the rows are not invisible to the gate machinery.  E/branch
gates touch dead labels in all six rows: `56`, `12`, `26`, `42`, `34`, and
`30` touching gates respectively.

The proof packet splits cleanly.  `random_covering_039` and
`random_covering_074` are the cover-delta sidecar rows, with minimum E/branch
gate deltas `(1,2)`.  The four clean unit-delta rows are
`random_covering_001`, `random_covering_062`, `random_covering_086`, and
`random_covering_101`.  The last has one nuance: its absolute shortest
unit-delta E/branch gate does not touch a dead label, but a slightly longer
unit-delta touching gate does, so it stays in the unit-delta singleton packet.

Proof-facing consequence: the remaining work is not a bigger pair-current
enumeration.  It is a local singleton-current theorem for mirror-paired
components, probably split into a two-row cover-delta clause and a four-row
unit-delta clause.  A raw global owner-current imbalance will not see it
because the imbalance cancels row-by-row; the current has to be local to
mirror pairs, touching gates, endpoint-spine data, two-adic descent,
exact-period/state-lift, or signed-SPEC/Rprime payload.

Tournament Analysis used proof carriers, not runners or raw row names:
singleton components, mirror-paired components, owner pairs, touching gate
events, unit-delta/cover-delta sidecars, and terminal proof obligations.  The
quotient preserves the zero-edge singleton-current discharge predicate.  It
destroys raw runner order, interval order if not retained, branch orientation,
endpoint wall labels, and owner-current locality; those fields must remain in
any formal packet statement.
