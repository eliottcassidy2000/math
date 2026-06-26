---
id: HYP-3074
status: SYNTHESIS / executable proof-interface scout; not a proof of LRC14
source: codex-2026-06-26-S240
tags: [lrc14, median-graphs, partial-cubes, proof-interface, sidecars, route-triples, controlled-forgetting, tournament-analysis]
related:
  - HYP-3073
  - HYP-3072
  - HYP-3071
  - HYP-3070
  - HYP-3069
  - HYP-3068
  - HYP-3067
  - HYP-3066
  - HYP-3064
  - HYP-3063
  - HYP-3062
  - HYP-3061
  - HYP-3056
  - HYP-3054
  - HYP-3037
  - HYP-3027
  - HYP-2997
  - HYP-2995
  - HYP-2963
  - HYP-2887
  - THM-572
  - OPEN-Q-099
  - OPEN-Q-108
---

# HYP-3074: Route-state closure median interface for LRC14

HYP-3073/S239 reserves the renormalized polymer / Dirichlet bridge: signed
polymer blocks and Dirichlet boundary currents should become state coordinates
before an energy or Schur-complement quotient is trusted. HYP-3072/S238 frames
the adjacent cross-carrier pullback resonance problem: carrier rows, hidden
coordinates, endpoint-owner strips, primitive decks, analytic clocks,
automaton states, partial-cube Theta classes, and formal exits should be
tested as proof vertices before scalarizing them. HYP-3071/S237
supplies the first-tooth observability and cycle-class image side of the
current residual certificate interface. HYP-3069/S235 turns the
final route-center gate into a median completion over
packet/route/certificate/sidecar/discharge carriers, and HYP-3070/S236 checks
raw route labels against legal route-triple center control.  This addendum
asks the next legality question: after each witness is closed under forced
sidecars, is the coordinate-wise median still a legal proof state?

```text
state = packet / route / certificate / sidecar / discharge
```

The target interface is:

```text
serious route triples have a unique legal center after legal sidecars are attached.
```

Here "unique" is not the hard part.  Coordinate-wise majority gives a unique
median in the ambient Boolean cube.  The hard part is making sure the legal
LRC proof states form a gated enough subcomplex that the median center is
still legal.  If the median leaves the legal subcomplex, the failure must name
one of:

```text
missing gated partial-cube sidecar
missing cycle-class image
missing observer-cut repair
explicit F7/THM-572 residual debt
```

## Interface

Let each proof witness be a bit vector plus labels over packet coordinates.
The S240 scout uses these coordinate families:

```text
packet:
  lrc_predicate
  packet_scale_exact
  endpoint_owner
  topology_safe

route:
  q_witness
  ap_gw_boundary_atom
  route_word_automatic
  fibbinary_partial_cube
  moser_two_lane_split
  doubled_triangular_bridge
  relation_lattice_covolume
  low_height_wall_ledger

certificate:
  fejer_interval_certificate
  toeplitz_square_scale_gate
  haar_zeta_square
  certificate_cycle_generators
  cycle_generator_checked
  hodge_positive_shadow
  hodge_cycle_image

sidecar:
  native_transition_word
  magnitude_cocycle
  theta_class_word
  gated_subcube_status
  median_interval_status
  simplex_face_rank
  bridge_rank_line_id
  ordered_quad_collapse_mode
  positive_scale_noncollapse
  observer_cut_payload
  cross_sector_orientation
  deletion_fiber_profile
  rectangle_hourglass_residue
  cycle_class_image_status

discharge:
  residual_capacitor_cut
  residual_tooth
  state_lift_obligation
  f7_state_lift_target
  residual_debt_exit
  discharge_resolved
  discharge_named_debt
```

A legal-sidecar closure operator `cl` attaches forced coordinates.  Examples:

```text
route_word_automatic
  -> native_transition_word + magnitude_cocycle

fibbinary_partial_cube
  -> theta_class_word + gated_subcube_status + median_interval_status

doubled_triangular_bridge
  -> bridge_rank_line_id + simplex_face_rank

toeplitz_square_scale_gate
  -> ordered_quad_collapse_mode + positive_scale_noncollapse

cycle_generator_checked + certificate_cycle_generators + hodge_positive_shadow
  -> cycle_class_image_status + hodge_cycle_image

state_lift_obligation
  -> f7_state_lift_target + residual_debt_exit
```

For a route triple `(a,b,c)`, define:

```text
center(a,b,c) = majority(cl(a), cl(b), cl(c)).
```

The triple is serious only if `center(a,b,c)` is legal.  If it is not legal,
the center is not a proof object; it is a named missing-sidecar obligation.

## Two proof angles

### 1. Automatic / Moser / fibbinary / partial-cube bridge

Raw automatic shadows are too lossy.  The S240 scout tests:

```text
automatic_word_shadow
partial_cube_bridge_rank
q_witness_exact_scale
```

The raw median remembers only:

```text
fibbinary_partial_cube
```

and fails because the center lacks:

```text
theta_class_word
gated_subcube_status
median_interval_status
```

After legal sidecars are attached, the median passes.  The center keeps the
gated partial-cube coordinates plus `magnitude_cocycle`; the non-majority
witnesses still carry the native-transition and bridge-rank data required to
audit the closure.  This gives the right proof-carrier meaning to the Moser
and fibbinary material: they are not scalar sequence tests; they are
coordinates in a gated cube/bridge-rank proof state.

### 2. Hodge / Toeplitz / Fejer certificate center

The Hodge-cycle lens from HYP-3066 says positivity is not realization.  The
median interface turns that warning into a direct test.  The scout tests:

```text
hodge_positive_shadow
fejer_interval_certificate
toeplitz_square_scale
```

Even after closure, the median fails because it has a positive closed packet
cochain and certificate generators but no:

```text
hodge_cycle_image
residual_debt_exit
```

The repaired triple:

```text
hodge_generated_cycle
fejer_interval_cycle
toeplitz_square_cycle
```

passes only after checked cycle generators force:

```text
cycle_class_image_status
hodge_cycle_image
```

The named-debt triple:

```text
hodge_phantom_debt_state
k33_state_lift
f7_residual_debt
```

also passes, but only because the median explicitly contains
`residual_debt_exit`.  This is the intended use of F7/THM-572: a coordinate in
the discharge state, not an informal hard-case label.

## Observer-cut warning

The triple:

```text
observer_cut_sector
haar_zipper_square
residual_capacitor_cut
```

still fails after closure.  The median sees only `observer_cut_payload`, while
the three repair lanes disagree: cross-sector orientation, cycle-image, and
residual-debt exit each appear in only one witness.  This is useful negative
evidence.  An observer-cut proof route must either make one repair sidecar
majority-visible or split the triple before medianization.

## Exact S240 scout

Script:

```text
04-computation/lrc14_route_state_closure_median_s240.py
```

Stored result:

```text
05-knowledge/results/lrc14_route_state_closure_median_s240.out
```

Readout:

```text
45 state coordinates
17 named proof states
11 legal-sidecar closure rules
5 route triples

automatic_partial_cube_router:
  raw median FAIL
  closed median PASS

hodge_toeplitz_fejer_phantom:
  raw median FAIL
  closed median FAIL

hodge_toeplitz_fejer_repaired:
  raw median FAIL
  closed median PASS

named_residual_debt_exit:
  raw median PASS
  closed median PASS

observer_cut_incompatible_repairs:
  raw median FAIL
  closed median FAIL
```

## Tournament Analysis

Vertices are proof states, not runners.  The tested vertex alternatives were:

```text
packet states
route obligations
certificate generators
sidecar coordinates
discharge exits
```

Pairwise observable:

```text
weighted retained proof-state coordinates, with a legality bonus
```

Binary relation:

```text
A -> B when A has higher retained-coordinate score;
ties follow the declared Hamiltonian path.
```

Fingerprint from S240:

```text
raw_score_hist = {0:1, 1:1, ..., 16:1}
closed_score_hist = {0:1, 1:1, ..., 16:1}
raw_directed_3cycles = 0
closed_directed_3cycles = 0
raw_hamiltonian_path_count = 1
closed_hamiltonian_path_count = 1
edge_flips_after_legal_sidecars = 59
```

The edge flips are the point: once legal sidecars are attached, many raw
boundary/q-witness states are no longer the strongest carriers.  The proof
order changes when lost coordinates are restored.

## Assumption challenge

Do not assume the vertices are runners, arcs, or even rows.  For this
interface the vertices should be proof states and the quotient should preserve:

```text
the legality and uniqueness of the route-triple median center.
```

It destroys raw runner names, raw time geometry, and distinctions between
witnesses that agree on all retained proof coordinates.  The challenged
assumption is that discharge mode can be a label outside the state.  It cannot:
if discharge is not a coordinate, a median center can look legal while hiding
F7/THM-572 debt.

## Next theorem target

Build the HYP-2963 medianization ledger.  For each live packet route, emit:

```text
state_id
packet fields
route fields
certificate fields
sidecar closure fields
discharge fields
legal_closure_added
median_triples_entering_this_state
failed_median_reason
```

Then prove or refute:

```text
Every serious LRC14 route triple has a legal closed median, and every failed
closed median is one of the four named missing-sidecar/debt exits.
```
