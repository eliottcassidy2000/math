# LRC14: Cross-Carrier Pullback Resonance

I tried a different remaining proof angle from the recent median-center and
cycle-class line: treat the CPI carrier-pullback rows themselves as a finite
proof object.

The scout is:

```text
04-computation/lrc14_cross_carrier_pullback_resonance_codex_s238.py
05-knowledge/results/lrc14_cross_carrier_pullback_resonance_codex_s238.out
HYP-3072
```

It encodes `22` carriers and `9` remaining proof obligations.  The core
payload alphabet is duodecimal:

```text
status, route, exact_scale, topology, owner, period_deck,
analytic_certificate, automaton_partial_cube, crt_padic,
observer_cut, hodge_cycle, formal_exit
```

The strongest readout is negative in the useful way:

```text
target_axes = 23
first full covers appear only at size 9
```

So no small universal scalar-like bundle is showing up.  The recurring full
cover skeleton is:

```text
packet/fusion
  + owner/topology/stalk
  + period/analytic
  + observer/formal/state-lift
```

Individual obligations are much smaller.  Examples:

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

The q=23 squarefree residual gives the clean warning.  `mu^2/phi` belongs in
the proof only as a blindness report.  The size-3 covers are:

```text
mobius_squarefree_blindness_report
  + crt_padic_residual_tree
  + carrier_fusion_switchboard

ramanujan_primitive_period_deck
  + mobius_squarefree_blindness_report
  + carrier_fusion_switchboard
```

So the proof role of `mu^2/phi` is not "capacity proves equality."  It is:

```text
this capacity meter is blind to repeated-prime / prime-power packet data;
now attach the period/CRT/owner sidecar that pays for that blindness.
```

Two proof routes seem worth developing.

First route: a blindness-pair theorem.  Every quotient must carry the carrier
that restores the coordinate it destroys.  Automata need exact Farey scale,
topology, owner, and route.  Observer cuts need formal exit plus exact scale.
Owner strips need primitive clocks.  Squarefree capacity needs Ramanujan/CRT
prime-power restoration.

Second route: a resonance-portfolio theorem.  Do not force all residuals
through one master quotient.  Pick the minimal portfolio per obligation, then
glue the local portfolios by the legal exit ledger.  The top recurring triples
combine:

```text
safe_stalk_barcode_normal_fan
observer or rectangle/hourglass cut payload
exact_farey_kpq_scale
```

This connects several older threads: the A000568 rootless perspective defect,
S217 rectangle/hourglass residues, Moser/fibbinary partial cubes, the
Fermat-Catalan hyperbolic guard, and multi-scale tournament spectra.  They are
all useful only as controlled-forgetting sidecars with exact packet payloads
attached.

Tournament Analysis used proof carriers as vertices, not runners:

```text
score_hist = {0:1, 1:1, ..., 21:1}
directed_3cycles = 0
scc_sizes = 22 singleton SCCs
hamiltonian_path_count = 1
```

Next exact job: emit actual HYP-2963 packet rows with
`carrier_pullback_row_id`, `core_incident_word`, `preserved_lrc_predicate`,
`destroyed_coordinate`, `required_sidecar`, `blindness_pair_id`,
`resonance_portfolio_id`, `status_mixing_result`, `route_mixing_result`, and
`legal_exit_status`, then test status/route purity before naming more theorem
debt.

