---
id: HYP-3070
status: SYNTHESIS / exact finite control scout; not a proof of LRC14
source: codex-2026-06-26-S236
tags: [lrc14, median-graphs, route-triples, sidecars, controlled-forgetting, tournament-analysis]
script: 04-computation/lrc14_route_triple_center_control_codex_s236.py
result: 05-knowledge/results/lrc14_route_triple_center_control_codex_s236.out
related:
  - HYP-3069
  - HYP-3068
  - HYP-3067
  - HYP-3066
  - HYP-3065
  - HYP-3064
  - HYP-3063
  - HYP-3062
  - HYP-3057
  - HYP-3056
  - HYP-3054
  - HYP-2963
  - THM-572
  - OPEN-Q-108
---

# HYP-3070: LRC14 Route-Triple Center-Control Addendum

The final LRC14 proof interface should not ask only whether a route triple has
a center in some completed Boolean cube.  It should first check whether the
route vocabulary itself is legal.

HYP-3069 gives the Boolean side of the story: full packet/route/certificate/
sidecar/discharge coordinates make all `220` serious route triples
deterministic, but raw projections leave `122/220` ambiguous triples and
produce `70` median-completion obligations.  HYP-3070 adds a sharper control:

```text
raw route labels alone are a clique and should be centerless;
legal sidecar pages should turn the same route leaves into a median carrier.
```

This is a quotient guardrail.  If a raw route-label projection already seems
to have a center, it is probably an accidental scalar or vocabulary collapse.
If legal sidecars do not provide a center, the missing coordinate must be named
before the triple becomes theorem debt.

## Exact Scout

The S236 script builds two finite models on the same `15` route leaves.

1. **Raw projection:** the route leaves form a clique.  For every distinct
   route triple, the pairwise intervals are just the endpoint edges, so the
   triple intersection is empty.

2. **Legal sidecar interface:** the route leaves attach to a packet/status/
   certificate/sidecar/discharge tree.  Trees are median graphs, so every leaf
   triple has one center.  The point is not that the final proof graph must be
   a tree; the point is that the legality sidecars must behave at least this
   cleanly on serious triples.

Main output:

```text
route_leaf_count=15
all_route_leaf_triples=455
raw_projection_unique_centers=0
raw_projection_empty_centers=455
legal_sidecar_unique_centers=455
legal_sidecar_empty_centers=0
```

Six named serious triples are checked as proof-interface examples:

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

The last row is the useful tiebreak.  Even though `COVERING_MOMENT` lives under
the owner-strip page, the two primitive-clock legs force the center at
`primitive_period_router` before owner-strip comparison.  This suggests a
local proof rule: when two legs share a clock family, center there before
using a more global owner or residual page.

## Proposed Packet Fields

Add these fields beside HYP-3069's Boolean median fields:

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

Suggested `center_control_exit` values:

```text
unique_named_center
first_missing_sidecar
raw_projection_false_center
sidecar_vocabulary_ambiguous
ap_gw_boundary_stop
primitive_clock_majority
owner_strip_descent
harmonic_certificate_backend
state_lift_debt
thm572_f7_debt
```

## Proof Angles Selected

### 1. Residual Route Triples

The center-control version of the residual problem is:

```text
Q / covering / K33 -> positive_residual_router
C27 / K33 / F7 -> resonant_state_lift_router
```

This isolates the open-residual proof surface before analytic or state-lift
certificates are invoked.  A future HYP-2963 audit should check whether every
route/status-changing residual triple lands in a named router, and if not,
which first sidecar is absent.

### 2. Guardrail Sidecar Triples

External proof motifs should not directly vote on a route:

```text
Moser partial cube / Toeplitz scale gate / Roth-Minkowski fence
  -> guardrail_sidecar_hub
```

This turns the wild analogy stack into a legal sidecar object.  The motif can
influence the proof only after it states the retained LRC predicate, destroyed
coordinate, sidecar vocabulary, and discharge rule.

## Tournament Analysis

Candidate vertex sets considered:

```text
runners, gaps, fixed circle sections, route labels, sidecar columns,
certificates, discharge exits, proof obligations, residual routers
```

Chosen vertices:

```text
packet/route/certificate/sidecar/discharge proof states
```

Pairwise observable:

```text
predicate retention, median uniqueness, sidecar legality,
first-missing-sidecar clarity, discharge namedness, formal checkability
```

Gauge:

```text
lexicographic retained-payload vector;
ties follow the printed Hamiltonian path
```

The scout tournament is transitive with `score_hist={0:1,...,11:1}`, no
directed 3-cycles, singleton SCCs, and one Hamiltonian path:

```text
labelled_packet_sheaf
> route_triple_center_control
> medianized_route_center_gate
> median_owner_root_spine
> desargues_median_lens
> boundary_status_gate
> positive_residual_router
> sidecar_observability_matrix
> harmonic_certificate_backend
> resonant_state_lift_router
> primitive_period_router
> raw_route_label_triangle
```

## Assumption Challenge

The challenged assumption is that route labels are the natural vertices of the
final proof graph.  They are not.  Raw route labels are a deliberately bad
quotient: as a clique they have no triple centers, and as scalar labels they
can also create false centers by forgetting sidecars.  The LRC predicate
preserved by the quotient is not the route name; it is the sidecar-completed
fact that a route/status-changing triple has a unique legal center, or names
the first missing sidecar.

Thus the next proof target becomes:

```text
In every HYP-2963 coarse fiber, every serious route triple either
has a unique named center after legal sidecars are attached,
exposes the first missing sidecar coordinate,
or routes to AP/GW boundary, primitive clock, owner-strip descent,
harmonic certificate, state-lift, or THM-572/F7 debt.
```
