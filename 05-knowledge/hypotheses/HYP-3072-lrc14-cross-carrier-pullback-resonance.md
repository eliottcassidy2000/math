---
id: HYP-3072
title: LRC14 cross-carrier pullback resonance
status: EVIDENCE / finite carrier-portfolio scout; not a proof of LRC14
source: codex-2026-06-26-S238
tangent: T1154
script: 04-computation/lrc14_cross_carrier_pullback_resonance_codex_s238.py
result: 05-knowledge/results/lrc14_cross_carrier_pullback_resonance_codex_s238.out
tags: [lrc14, carrier-pullback, controlled-forgetting, proof-obligations, squarefree-blindness, partial-cube, tournament-analysis]
related:
  - HYP-3071
  - HYP-3070
  - HYP-3069
  - HYP-3066
  - HYP-3065
  - HYP-3063
  - HYP-3058
  - HYP-3039
  - HYP-3032
  - HYP-3029
  - HYP-3026
  - HYP-2963
  - THM-572
  - OPEN-Q-108
---

# HYP-3072: LRC14 Cross-Carrier Pullback Resonance

This pass deliberately chooses a different remaining proof angle from the
recent median-center and cycle-class line.  Instead of asking for one more
scalar, route graph, or Hodge matrix, it treats the CPI pullback rows as a
finite carrier portfolio.

The useful question becomes:

```text
Which mixed portfolios cover the remaining proof obligations,
and which single carriers are dangerous because they preserve a signal while
destroying the coordinate that obligation needs?
```

This is a controlled-forgetting theorem shape.  A quotient is legal only after
it names the LRC predicate preserved, the coordinate destroyed, and the
sidecar, certificate, descent, AP/GW boundary stop, or THM-572/F7 residual
that pays for the forgetting.

## Scout

The S238 script:

```text
04-computation/lrc14_cross_carrier_pullback_resonance_codex_s238.py
05-knowledge/results/lrc14_cross_carrier_pullback_resonance_codex_s238.out
```

encodes `22` carrier pullbacks and `9` remaining proof obligations.  The
carriers include labelled packets, closed arc-Cech topology, safe-component
stalk/barcode/normal-fan geometry, endpoint-owner strips, Ramanujan primitive
decks, Mobius squarefree-blindness, Haar zeta, Fejer intervals, CRT/p-adic
trees, Moser/fibbinary partial cubes, observer cuts, rectangle/hourglass
flows, Hodge cycles, median route centers, Toeplitz scale gates,
Roth-Minkowski walls, K33 state lifts, carrier fusion, sidechannel repair,
Farey/Kpq scale, pair-good blocker grammar, and the hyperbolic
Fermat-Catalan power guard.

It records a duodecimal incident word over the core payload alphabet:

```text
status, route, exact_scale, topology, owner, period_deck,
analytic_certificate, automaton_partial_cube, crt_padic,
observer_cut, hodge_cycle, formal_exit
```

Example readouts:

```text
carrier_fusion_switchboard         word=111110001000
labelled_packet_sheaf              word=111000000000
median_route_center_control        word=010010000001
exact_farey_kpq_scale              word=011000000000
automatic_partial_cube_sidecar     word=000000010000
mobius_squarefree_blindness_report word=000000000000
```

The word alone is not enough; the lost-coordinate ledger is the second half of
the object.

## Main Finite Findings

The global cover result is the strongest negative readout:

```text
target_axes = 23 axes across the nine obligations
first full covers appear only at size 9
```

So no small universal scalar-like bundle appears in this finite atlas.  The
recurring skeleton is:

```text
packet/fusion
  + owner/topology/stalk
  + period/analytic
  + observer/formal/state-lift
```

Many individual obligations do have compact covers:

```text
automatic_partial_cube_route_purity:
  automatic_partial_cube_sidecar + carrier_fusion_switchboard

two_tail_ap_clock_search:
  endpoint_owner_strip_current + crt_padic_residual_tree

k33_f7_state_lift_exit:
  hodge_cycle_lift + carrier_fusion_switchboard
  or k33_state_lift_incidence + carrier_fusion_switchboard

multi_scale_binding_spectrum:
  observer_extension_cut_payload + exact_farey_kpq_scale
  or rectangle_hourglass_diagonal_flow + exact_farey_kpq_scale
```

The squarefree `q=23` obligation is a useful warning.  It requires exact
scale, Haar zeta, owner, period deck, squarefree-blindness, and route.  The
first size-3 covers are:

```text
mobius_squarefree_blindness_report
  + crt_padic_residual_tree
  + carrier_fusion_switchboard

ramanujan_primitive_period_deck
  + mobius_squarefree_blindness_report
  + carrier_fusion_switchboard
```

Thus `mu^2/phi` belongs in the proof, but only as a blindness report paired
with prime-power/period and packet-family carriers.

## Creative Proof Angle A: Blindness-Pair Proof

The first new route is to make every attractive quotient carry a paired
restoration obligation.

Examples:

```text
mu^2/phi capacity
  destroys prime powers / owner / route
  pair with Ramanujan primitive deck or CRT/p-adic tree

Moser/fibbinary partial cube
  destroys status / route / exact scale / owner
  pair with exact Farey scale plus topology/owner or fusion

observer cut
  destroys status / route / exact scale
  pair with formal exit plus exact scale

owner strip
  destroys primitive clock / exact scale
  pair with Ramanujan deck or CRT tree
```

The theorem target:

```text
No carrier may be used as a quotient unless every destroyed coordinate needed
by the active obligation has a named restoration, annihilation, descent,
AP/GW boundary stop, or THM-572/F7 exit.
```

This turns failures into proof data instead of bibliography.  A failed quotient
is useful precisely because it says which sidecar must be restored.

## Creative Proof Angle B: Resonance-Portfolio Proof

The second route is to stop trying to force every residual through one master
object.  Work with minimal portfolios per obligation, then glue the portfolios
by the formal exit ledger.

The top pair/triple resonances repeatedly involve:

```text
safe_stalk_barcode_normal_fan
  + rectangle_hourglass_diagonal_flow
  + exact_farey_kpq_scale

safe_stalk_barcode_normal_fan
  + observer_extension_cut_payload
  + exact_farey_kpq_scale

safe_stalk_barcode_normal_fan
  + automatic_partial_cube_sidecar
  + exact_farey_kpq_scale
```

These are not accidental.  They combine local owner/stalk geometry, an
observer/edge-sector payload, and exact binding-scale/Farey data.  That is the
same structural shape behind the earlier A000568 perspective defect, the S217
rectangle/hourglass residue, the Moser/fibbinary partial-cube warning, and the
multi-scale tournament-spectrum fix for magnitude blindness.

The theorem target:

```text
For each residual family, choose the minimal carrier portfolio whose union of
retained axes covers the family obligation.  A global proof glues these local
portfolios by legal exits, not by replacing them with a single scalar.
```

Incoming-mainline connection, 2026-06-26: the checked Lean moment-dual audit in
`07-reflections/lrc14-lean-moment-dual-audit-codex-20260626.md` is exactly a
`formal_exit` carrier refinement for this portfolio route.  Its verified
`p0_le_Ly` bridge says the cover-to-moment relaxation is no longer an informal
sidecar; the remaining proof payload is a finite rational dual certificate plus
the concrete event bridge from `shapeOf`/`witnessG2` into cover/safe/good sets.
In S238 language, this is not a scalar replacement for the portfolio.  It is a
legal-exit row that must name the retained `missCount` predicate, the destroyed
endpoint-owner/location coordinates, and the certificate column that restores
enough of that loss for Part A.

## Tournament Analysis

Vertices are proof-carrier pullbacks / CPI rows, not runners.  Candidate
vertex sets considered included runners, gaps, route labels, proof obligations,
CPI rows, hidden coordinates, endpoint-owner strips, primitive decks, analytic
clocks, automaton states, partial-cube Theta classes, CRT roots, observer
cuts, rectangle/hourglass residues, Hodge cochains, median-center pages,
state-lift sectors, and formal exits.

Pairwise observable:

```text
(full_obligation_count,
 weighted_axis_coverage,
 critical_axis_hits,
 payload_width,
 -destroyed_count,
 -cost)
```

Gauge: orient toward more retained proof payload; ties use the CPI pullback tie
path.  S238 reports:

```text
score_hist = {0:1, 1:1, ..., 21:1}
directed_3cycles = 0
scc_sizes = 22 singleton SCCs
hamiltonian_path_count = 1
```

Tie Hamiltonian path:

```text
carrier_fusion_switchboard
> labelled_packet_sheaf
> median_route_center_control
> exact_farey_kpq_scale
> sidechannel_repair_ladder
> ramanujan_primitive_period_deck
> endpoint_owner_strip_current
> fejer_interval_backend
> k33_state_lift_incidence
> crt_padic_residual_tree
> safe_stalk_barcode_normal_fan
> haar_zeta_square
> toeplitz_square_scale_gate
> hyperbolic_power_guard
> pair_good_blocker_grammar
> roth_minkowski_low_height_wall
> observer_extension_cut_payload
> hodge_cycle_lift
> closed_arc_cech_topology
> rectangle_hourglass_diagonal_flow
> automatic_partial_cube_sidecar
> mobius_squarefree_blindness_report
```

This ordering should not be read as proof rank.  It is a payload-retention
audit: broad fused and packet carriers lead, raw squarefree capacity comes last
unless restored by period/CRT/owner sidecars.

## Next Exact Job

Emit actual HYP-2963 packet rows with:

```text
carrier_pullback_row_id
core_incident_word
preserved_lrc_predicate
destroyed_coordinate
required_sidecar
blindness_pair_id
resonance_portfolio_id
status_mixing_result
route_mixing_result
legal_exit_status
```

Then test whether the listed portfolios make the residual coarse fibers
status-pure and route-pure before any new theorem debt is named.
