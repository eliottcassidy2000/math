# HYP-2975: LRC14 taut bridge graph curvature

**Status:** reserved proof-interface stub, codex-2026-06-24-S155.

This is a local refinement of the incoming HYP-2970 endpoint-credit winding
cycle dual, not a replacement for it.  HYP-2970 works on the open-cover
transition graph of danger arcs and uses endpoint credit
`K((s,m),(r,n)) = 14(rm-sn)+r+s`.  HYP-2975 asks what happens at the boundary
where that graph degenerates: the zero-length taut vertices, positive safe
bridges, and signed owner-current left by the exact endpoint arrangement.

For a 13-speed row at threshold `1/14`, sweep the open danger arcs and record:

- positive safe intervals as directed bridges from the endpoint owner that
  stops covering on the left to the endpoint owner that starts covering on the
  right;
- isolated safe points as zero-length taut vertices with positive cover depth
  on both sides and zero point-depth at the event;
- the signed owner-current and mod-14/gcd support of these transfers.

The proof target is not that this graph by itself proves LRC14.  The honest
target is narrower: a strict counterexample would have no positive bridge and
no taut vertex, while boundary-only equality atoms should have a forced
zero-curvature taut transfer graph.  If every non-AP/GW labelled packet either
creates a positive bridge, violates Toeplitz PSD (HYP-2974), or carries a named
K33/state-lift transfer, then the taut graph is the local combinatorial layer
between the Haar/Baire packet route and the endpoint/Fourier dual routes.

Evidence still missing in this stub:

- an exact endpoint-sweep audit over the named AP/GW/petal/K33/covering rows;
- a bank scan checking whether any non-AP/GW row has zero positive bridges and
  a nonempty taut graph;
- Tournament Analysis fingerprints using bridge/taut-state vertices rather
  than raw runners.

Assumption challenge: runners are deliberately not the chosen tournament
vertices.  The considered vertex sets are runners, raw endpoints, endpoint
owner labels, positive safe intervals, isolated taut points, boundary-current
states, missed-depth sectors, K33 state-lift obligations, and proof routes.
This stub chooses endpoint-owner transfer states because they preserve the LRC
predicate "open witness, boundary-only equality, or fully covered residual"
before scalarization.  It destroys row magnitude away from labelled endpoints,
so it cannot replace exact `M`, qdiv, Farey scale, endpoint credit, Toeplitz
PSD, or C27/K33 labels.

