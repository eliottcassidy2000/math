---
id: HYP-3301
title: LRC14 sheaf exactness and Farey-cusp transfer proof angles
status: SYNTHESIS / executable proof-angle scout; not an LRC14 proof
source: codex-2026-06-28
tangent: T1356
technique: LTI-356
tournament_technique: LTT-256
script: 04-computation/lrc14_sheaf_cusp_transfer_angles_codex_20260628.py
result: 05-knowledge/results/lrc14_sheaf_cusp_transfer_angles_codex_20260628.out
reflection: 07-reflections/lrc14-sheaf-exactness-cusp-transfer-codex-20260628.md
related:
  - HYP-3300
  - HYP-3265
  - HYP-3257
  - HYP-3255
  - HYP-3253
  - HYP-3247
  - HYP-3246
  - HYP-3243
  - HYP-3242
  - HYP-3234
  - HYP-3231
  - HYP-3230
  - HYP-3102
  - HYP-2969
  - HYP-2963
  - HYP-2954
  - HYP-2704
  - THM-573
  - THM-523
  - OPEN-Q-108
---

# HYP-3301: LRC14 Sheaf Exactness And Farey-Cusp Transfer

## Claim

After HYP-3300, two remaining LRC14 proof angles look different enough to
deserve their own target theorem:

```text
Angle 1: first-obstruction sheaf exactness
Angle 2: Farey-cusp renormalization transfer
```

The first route is not another observability matrix.  It asks for the
quotient/observer overlap cochain complex.  Whenever a quotient forgets a
payload, the first forgotten payload is a cocycle.  The proof target is that
each such cocycle is:

```text
exact,
killed by zeta_7 contact holonomy,
lifted to an endpoint/finite chamber,
descended to a smaller family,
stopped at AP/Goddyn-Wong equality,
or routed to named K33/H7/state-lift debt.
```

The second route does not try to prove the covering branch by a single
denominator or scalar moment.  It treats `qdiv>14` rows as a scale-normal cusp
transfer:

```text
exact-period boundary
  -> boundary-moment image
  -> positive floor,
     impossible AP/GW kernel,
     K33/H7 named debt,
     or new zero-open kernel.
```

This reframes the remaining hard part as a theorem about the kernel of the
boundary-moment transfer, not a search for a new global extremal number.

## Executable Readout

The scout ranked ten proof carriers.  The top two were:

```text
B01 first_obstruction_sheaf_exactness       score 117
B02 farey_cusp_renormalization_transfer    score 112
```

The next tier was:

```text
B04 contact_graph_kernel_classifier         score 77
B05 qsqrt7_residue_to_magnitude_lift        score 77
B03 boundary_moment_chain_map               score 76
B06 schwarz_christoffel_contact_polygon     score 56
B09 proof_circuit_uniform_gate_test         score 52
B08 hyperbolic_reciprocal_defect            score 49
B07 loglog_mertens_scale_entropy            score 44
B10 raw_cross_discipline_analogy            score -33
```

Readout: HYP-3265's contact graph, HYP-3255/HYP-3257's residue/magnitude
split, and the old boundary-moment map are not rival approaches.  They are
supporting charts for the two larger target theorems.

## First-Obstruction / Sheaf Exactness Packet

The toy exactness packet uses these chart sidecars:

```text
unit_contact_matching
danger_nerve_hole
endpoint_arrangement_cell
zeta7_contact_holonomy
shell_lag_commutator
index_degree_packet
floor_margin_certificate
signed_address_word
first_obstruction_cocycle
repair_cover_rank
Qsqrt7_residue_orbit
magnitude_jump_anchor
unit_nullspace
exact_period_boundary
boundary_moment_image
Farey_parent_interval
cusp_principal_part
renormalization_depth
hensel_unit
K33_kernel_label
state_lift_H7
schwarz_christoffel_turning_word
accessory_parameter_debt
```

Rows are overlap failures between existing charts:

```text
unit_contact_to_danger_nerve
lag_residue_to_shell_magic
index_degree_to_floor_gap
residue_to_magnitude_split
chart_change_to_first_obstruction
exact_period_to_boundary_moment
covering_to_state_lift_kernel
farey_cusp_to_scale_normal
polygon_contact_to_endpoint_owner
renormalized_kernel_to_named_debt
```

The toy matrix has full row rank:

```text
rows = 10
cols = 23
GF(2) rank = 10
minimal chart-sidecar hitting-set size = 5
```

This does not prove the theorem.  It says the current chart vocabulary is
large enough to separate the named overlap failures without collapsing back
to raw scalar telemetry.

The proof-facing theorem should be:

```text
For every primitive residual packet and every legal quotient/observer edge,
the first obstruction cochain is exact, holonomy-repaired, endpoint-lifted,
family-descended, stopped at AP/GW equality, or routed to named debt.
```

This pulls HYP-3102's "first forgotten payload is a syndrome" rule into the
current HYP-3253/HYP-3265/HYP-3300 frontier.

## Farey-Cusp Renormalization Transfer

The second theorem target uses a small state potential:

```text
potential = (covering_depth, cusp_height, kernel_unknown, sidecar_debt)
```

States:

```text
F0_qdiv_witness:
  terminal loose witness with M >= 1/qdiv > 1/14

F1_q14_unit_contact_core:
  AP/GW equality or strict off-unit margin by contact classifier

F2_q14_killed_contact_covering:
  promote to Phi_14d or covering floor sidecar

F3_qgt14_exact_period_cover:
  exact-period boundary feeds boundary-moment map

F4_positive_boundary_moment_image:
  strict covering floor witness

F5_AP_GW_kernel_in_covering:
  forbidden by qdiv>14 unless scale labels were lost

F6_K33_H7_named_kernel:
  route to K33/state-lift theorem debt

F7_unknown_zero_open_kernel:
  the actual remaining target: prove empty or name new sidecar
```

Useful transitions either lower the packet to a witness, produce positive
boundary-moment slack, prove that the AP/GW kernel is illegal in `qdiv>14`, or
name K33/H7 debt.  The only bad transition is:

```text
F3_qgt14_exact_period_cover -> F7_unknown_zero_open_kernel
```

Thus the transfer theorem should be:

```text
Every qdiv>14 exact-period boundary has positive boundary-moment image,
or lies in a kernel already forbidden by qdiv>14,
or lies in named K33/H7 state-lift debt,
or supplies the first genuinely new zero-open kernel class.
```

This is the old boundary-moment adjunction made compatible with HYP-3230's
Farey/Stern-Brocot cap-kernel recursion, HYP-3231's scale-normal packet
ledger, HYP-3255's residue/magnitude warning, and HYP-3265's contact-graph
case split.

## Cross-Disciplinary Synthesis

The useful external ideas are not proof shortcuts.  They become exact only
when translated into packet fields:

```text
Cech/sheaf theory:
  chart overlaps and first obstruction cocycles.

Local systems:
  zeta_7 contact holonomy as a connection coordinate.

Modular/Farey geometry:
  qdiv and exact-period packets moving on a cusp tree.

Schwarz-Christoffel:
  contact points as polygon vertices; accessory parameter = off-unit
  endpoint-owner debt.

Mertens / loglog:
  obstruction-prime hierarchy as proof-complexity telemetry, not witness.

Fermat-Catalan / Hurwitz / Markov:
  sum reciprocal < 1 as hyperbolic triple-cusp defect after a concrete triple
  of packet widths is declared.
```

The creative point is that the same controlled-forgetting discipline governs
all of them: preserved predicate, destroyed coordinate, repair sidecar,
terminal exit.

## Tournament Analysis

Vertices are proof carriers and chart transitions, not runners or arcs.

Pairwise observable:

```text
which carrier preserves boundary/open status,
chart gluing,
covering exits,
scale normality,
owner/topology,
residue/magnitude split,
and formal packet schema
while destroying fewer coordinates.
```

Switch/gauge:

```text
A -> B iff A has larger weighted proof-carrier score;
ties prefer fewer destroyed coordinates and then stable key order.
```

Fingerprint:

```text
score_hist = {-33:1, 44:1, 49:1, 52:1, 56:1, 76:1, 77:2, 112:1, 117:1}
directed_3cycles = 0
scc_sizes = [1,1,1,1,1,1,1,1,1,1]
hamiltonian_path_count = 1
priority_path =
  first_obstruction_sheaf_exactness
  -> farey_cusp_renormalization_transfer
  -> contact_graph_kernel_classifier
  -> qsqrt7_residue_to_magnitude_lift
  -> boundary_moment_chain_map
  -> schwarz_christoffel_contact_polygon
  -> proof_circuit_uniform_gate_test
  -> hyperbolic_reciprocal_defect
  -> loglog_mertens_scale_entropy
  -> raw_cross_discipline_analogy
```

## Status

HYP-3301 is a proof-angle synthesis and executable scout, not a proof of
LRC14.  Its immediate value is to split two remaining targets:

```text
1. exactness of the chart-overlap obstruction complex;
2. kernel classification for the Farey-cusp boundary-moment transfer.
```

If either theorem fails, the failure should be useful: it names the first
unspanned cocycle or the first real zero-open cusp kernel rather than leaving
the proof stuck behind an anonymous residual bucket.
