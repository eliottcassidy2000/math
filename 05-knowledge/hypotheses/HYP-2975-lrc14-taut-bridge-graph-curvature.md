# HYP-2975: LRC14 taut bridge graph curvature

**Status:** proof-interface / exact bounded audit, codex-2026-06-24-S155.

This is a local refinement of the incoming HYP-2970 endpoint-credit winding
cycle dual, not a replacement for it.  HYP-2970 works on the open-cover
transition graph of danger arcs and uses endpoint credit
`K((s,m),(r,n)) = 14(rm-sn)+r+s`.  HYP-2975 asks what happens at the boundary
where that graph degenerates: the zero-length taut vertices, positive safe
bridges, and signed owner-current left by the exact endpoint arrangement.
After rebasing over HYP-2977, read this as the local endpoint-owner quotient
complementing the global spectral-shadow quotient: HYP-2975 keeps boundary
ownership and loses Fourier energy, while HYP-2977 keeps Fourier energy and
loses endpoint ownership.

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

Evidence added by S155:

- Named-row audit: AP and Goddyn-Wong `12->24` have safe mass `0`, no positive
  bridges, and exactly six isolated taut vertices.  Each taut vertex has one
  endpoint-owner pair with speed sum `0 mod 14`, and the total owner-current is
  zero.
- Named non-AP/GW rows expose positive bridges: `12->36` has safe mass
  `1/1260`, `10->20` has `1/980`, `13->26` has `1/182`, `P10+K33` has
  `4/2205`, `12->84` has `563/105105`, and `drop13 add182` has
  `4637/194040`.
- AP-neighborhood bank scan with one-swaps through added value `160` and
  two-swaps through added value `36`: `21645` primitive rows, `21644`
  positive-open rows, and one zero-open row, namely the GW single swap
  `12->24`.  The smallest positive safe mass is `1/1260` at `12->36`.
- Tournament Analysis over proof carriers has fingerprint
  `score_hist={0:1,2:2,3:1,4:2,6:1,7:1}`, `directed_3cycles=3`, and SCC sizes
  `[5,1,1,1]`.  The leading order is endpoint-credit winding,
  taut-vertex-current, Toeplitz PSD dual, C27/K33 labels, positive open bridge,
  exact M/Farey scale, multiplicity-moment dual, raw runner vertices.

Artifacts:

- `04-computation/lrc14_taut_bridge_graph_codex_s155.py`
- `05-knowledge/results/lrc14_taut_bridge_graph_codex_s155.out`
- `07-reflections/lrc14-taut-bridge-graph-curvature-codex-s155.md`

Evidence still missing:

- an unbounded theorem converting the finite taut audit into a global
  endpoint-current rigidity statement;
- a proof that every zero-positive-bridge residual either has AP/GW taut
  current, violates HYP-2974 Toeplitz PSD at bounded degree, or constructs the
  HYP-2908/THM-572 K33 state lift;
- comparison of the taut-current support with the actual HYP-2970
  endpoint-credit potential solver once that solver exists.
- comparison with the HYP-2977 spectral-shadow packets, especially whether
  high-frequency relation tails determine the same endpoint-owner current.

Assumption challenge: runners are deliberately not the chosen tournament
vertices.  The considered vertex sets are runners, raw endpoints, endpoint
owner labels, positive safe intervals, isolated taut points, boundary-current
states, missed-depth sectors, K33 state-lift obligations, and proof routes.
This stub chooses endpoint-owner transfer states because they preserve the LRC
predicate "open witness, boundary-only equality, or fully covered residual"
before scalarization.  It destroys row magnitude away from labelled endpoints,
so it cannot replace exact `M`, qdiv, Farey scale, endpoint credit, Toeplitz
PSD, or C27/K33 labels.

Readout: HYP-2975 gives a useful local falsifier shape.  A future strict
counterexample must have no positive bridge, no AP/GW zero-sum taut current,
and no harmonic or K33/state-lift exit.  The bounded audit found no such row.
