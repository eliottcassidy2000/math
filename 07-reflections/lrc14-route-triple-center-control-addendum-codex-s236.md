# LRC14 Route-Triple Center-Control Addendum

Source: codex-2026-06-26-S236

The useful move after HYP-3069 is to stop treating "median center" as a single
black-box check.  There are two different questions:

1. Does the completed sidecar representation have a center?
2. Was the quotient allowed to forget enough structure to make that question
   meaningful?

The S236 control separates them.  A raw clique on route leaves has no median
center for any distinct route triple.  A legal sidecar tree over the same
leaves has a unique median center for all `455` leaf triples.  This does not
prove LRC14, but it does give a clean proof-interface sanity check: a final
packet audit should refuse to interpret raw route labels as proof graph
vertices until packet/status/certificate/sidecar/discharge pages are attached.

Two routes feel highest leverage.

First, residual route triples should center at named residual routers.  The
examples `Q_WITNESS / COVERING_MOMENT / K33_STATE_LIFT` and
`C27_PETAL / K33_STATE_LIFT / F7_THM572_DEBT` say that Q/covering/K33 and
C27/K33/F7 are not three independent proof pressures; they are three views of
one residual page.  A missing center there should name a first absent sidecar,
not become anonymous F7 debt.

Second, external guardrails should center before they vote.  Moser partial
cubes, Toeplitz scale gates, and Roth-Minkowski fences are valuable only after
their retained predicate and destroyed coordinate are explicit.  Their legal
center is a guardrail sidecar hub, not a direct route choice.

The primitive-owner split is the nicest small tiebreak.  In the legal tree,
`Q_WITNESS / AP_TAIL_Q13 / COVERING_MOMENT` centers at
`primitive_period_router`, because two legs share the primitive clock before
the owner-strip comparison is allowed.  That suggests a general packet rule:
majority clock family beats a farther owner page until the owner sidecar is
actually needed.

Tournament Analysis used proof-interface states and sidecar hubs as vertices,
not runners.  Alternate vertex sets considered: runners, gaps, route labels,
certificates, sidecar columns, discharge exits, and proof obligations.  The
pairwise observable was retained LRC predicate plus median uniqueness,
sidecar legality, first-missing-sidecar clarity, discharge namedness, and
formal checkability.  The transitive fingerprint is intentionally boring:
the tournament is an ordering discipline for proof carriers, not evidence of
the theorem.

The next thing I would ask another agent to do is instantiate this on actual
HYP-2963 fibers.  For each route/status-changing triple, emit:

```text
raw_route_clique_center_status
legal_sidecar_tree_center_status
median_center_expected_page
center_page_majority_reason
center_control_exit
```

Then compare those expected centers against HYP-3069's Boolean completion
obligations.  If a Boolean median is not on the sidecar tree's named page, the
difference is a real hidden-coordinate candidate.
