# LRC14 Route-Triple Center-Control Addendum

This is the S236 follow-up to the HYP-3069 medianized route-center gate.

HYP-3069 showed that full Boolean sidecars make all `220` named route triples
deterministic, while raw projections leave `122` ambiguous triples and create
`70` median-completion obligations.  The addendum isolates a smaller interface
test:

```text
raw route labels alone are a clique and should have no centers;
legal sidecar pages should turn the same route leaves into a median carrier.
```

The scout has `15` route leaves.  In the raw route-label clique, all
`455` leaf triples have empty median intersection.  In the legal sidecar tree,
all `455` have exactly one center.

Serious example centers:

```text
Q_WITNESS / COVERING_MOMENT / K33_STATE_LIFT
  -> positive_residual_router

AP_BOUNDARY / GW_BOUNDARY / Q_WITNESS
  -> boundary_router

FEJER_INTERVAL / RAMANUJAN_PROJECTOR / HAAR_ZETA
  -> harmonic_certificate_backend

MOSER_PARTIAL_CUBE / TOEPLITZ_SCALE_GATE / ROTH_MINKOWSKI_FENCE
  -> guardrail_sidecar_hub

C27_PETAL / K33_STATE_LIFT / F7_THM572_DEBT
  -> resonant_state_lift_router

Q_WITNESS / AP_TAIL_Q13 / COVERING_MOMENT
  -> primitive_period_router
```

The last row is the useful tiebreak: two primitive-clock legs center at the
primitive-period router before the owner-strip leg gets to decide the route.

New packet fields proposed:

```text
route_triple_center_control
raw_route_clique_center_status
legal_sidecar_tree_center_status
median_center_expected_page
center_page_depth
center_page_majority_reason
guardrail_sidecar_center
center_control_exit
```

Tournament Analysis vertices are proof-interface states and sidecar hubs, not
runners.  Candidate alternatives considered were runners, gaps, route labels,
certificates, sidecar columns, discharge exits, and proof obligations.  The
pairwise observable is predicate retention plus median uniqueness, sidecar
legality, first-missing-sidecar clarity, discharge namedness, and formal
checkability.  The fingerprint is transitive with one Hamiltonian path, which
is correct here: this tournament is an ordering discipline for legal proof
carriers, not a scalar theorem signal.

Next target:

```text
For every HYP-2963 coarse fiber, every route/status-changing triple should
either have a unique named center after legal sidecars are attached, expose
the first missing sidecar coordinate, or route to AP/GW boundary, primitive
clock, owner-strip descent, harmonic certificate, state lift, or THM-572/F7.
```
