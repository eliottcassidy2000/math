---
id: HYP-3300
title: LRC14 observability-matrix and discrete-Morse proof angles
status: SYNTHESIS / executable proof-angle scout; not an LRC14 proof
source: codex-2026-06-28
tangent: T1346
technique: LTI-346
tournament_technique: LTT-246
script: 04-computation/lrc14_observability_morse_angles_codex_20260628.py
result: 05-knowledge/results/lrc14_observability_morse_angles_codex_20260628.out
reflection: 07-reflections/lrc14-observability-morse-proof-angles-codex-20260628.md
related:
  - HYP-3249
  - HYP-3248
  - HYP-3247
  - HYP-3246
  - HYP-3245
  - HYP-3244
  - HYP-3243
  - HYP-3242
  - HYP-3241
  - HYP-3240
  - HYP-3239
  - HYP-3238
  - HYP-3237
  - HYP-3236
  - HYP-3235
  - HYP-3234
  - HYP-3233
  - HYP-3232
  - HYP-3231
  - HYP-3230
  - HYP-3229
  - HYP-3228
  - HYP-3227
  - HYP-3225
  - HYP-3224
  - HYP-3223
  - HYP-3219
  - HYP-3108
  - HYP-3069
  - HYP-3070
  - HYP-3048
  - OPEN-Q-108
---

# HYP-3300: LRC14 Observability-Matrix And Discrete-Morse Proof Angles

## Claim

The next useful LRC14 push should not be another single scalar extremality
search.  Two remaining proof angles look different enough from the recent
topology/geometry cover atlas to be worth isolating:

```text
Angle 1: sidecar observability matrix
Angle 2: finite chamber discrete-Morse descent
```

Both are controlled-forgetting statements.

The observability route says that every quotient used in the proof must retain
a small basis of hidden coordinates sufficient to separate each residual pair
whose LRC status, route status, or terminal exit can change.

The Morse route says that after those coordinates are retained, the finite
chamber complex should admit a vector-valued descent.  Every nonterminal
chamber should have a legal descending wall, and every critical chamber should
be one of:

```text
strict open safe cell
AP/Goddyn-Wong Phi14 equality
Phi14d dilation equality
finite Toeplitz/Green/root-motion discharge
state-lift H=7 contradiction
named residual debt
```

HYP-3243 supplies the chamber carriers.  HYP-3244 supplies the finite
lift/compress descent discipline.  HYP-3245 supplies lag/autocorrelation
transport as a scalar projection that the new matrix must police.  Incoming
HYP-3246/HYP-3247 supply the Chebyshev/unit-equioscillation index and the
three binding complement-pair word.  Incoming HYP-3248/HYP-3249 turn that into
a q-uniform Chebyshev frame and an index-prediction packet with analytic index
= topological degree, a Gauss-sum index word, and a named Borsuk-Ulam forcing
gap.  HYP-3300 treats these as boundary/core sidecars rather than as a
replacement for the finite descent proof.  Incoming HYP-3253 sharpens the
guardrail: the index packet describes the saddle, while floor/margin data,
contact holonomy, or a finite endpoint-chamber certificate still has to prove
the terminal step.  Incoming HYP-3255 adds the matching residue/magnitude
warning: `Q(sqrt-7)` and the unit-grid organize the residue layer, but exact
minimal-cover rigidity lives at the magnitude/equioscillation/census layer.
Incoming HYP-3257/HYP-3258 turn that warning into explicit columns: unit
rank-3 nullspace, blind residue/height ledger, strict safe-component atlas,
covering 14-multiple kill switch, and binding/covering census split.
The S289 coordination layer supplies two useful but nonterminal columns:
`Roth_Halasz_discrepancy_bound` for bulk chambers and
`Hensel_Krasner_valuation_unit` for scale-lift stability.

## Executable Readout

The scout ranked ten proof angles.  The top two were:

```text
A02 finite_chamber_discrete_morse_descent      score 96
A01 sidecar_observability_matrix              score 87
```

The supporting columns then ranked:

```text
A05 odd_negative_duality_resurrection          score 60
A03 green_toeplitz_thomson_energy              score 57
A08 lee_yang_ear_root_motion                   score 57
A06 tiling_half_tiling_descent_span            score 56
A04 lag_autocorrelation_transport              score 54
A09 median_route_center_closure                score 53
A07 roth_halasz_hensel_lift_packet             score 42
A10 raw_single_scalar_extremality              score -24
```

The important point is not the exact scoring.  The priority order says the
older routes should be promoted to coordinates inside two broader proof
programs:

```text
observability matrix decides what must be retained
Morse descent decides how the retained packet terminates
```

## Sidecar Observability Matrix

The proposed columns are:

```text
endpoint_owner
boundary_cocircuit
phi_witness_address
dilation_grid
toeplitz_normal_fan_slack
green_resistance_slack
odd_negative_payload
lee_yang_root_word
tiling_descent_packet
lag_transport_signature
unit_equioscillation_index
binding_complement_pair_word
analytic_topological_index_equalizer
gauss_sum_index_word
borsuk_ulam_forcing_gap
roth_halasz_discrepancy
hensel_krasner_unit
state_lift_H7
cech_betti_hole
ear_payload
```

The scout used thirteen proposed residual-pair rows:

```text
open_safe_vs_closed_boundary
AP_vs_GW_phi14_equality
chebyshev_units_vs_covering_floor
index_prediction_vs_naive_bu_gap
base_phi14_vs_dilation_phi14d
toeplitz_pair_plucker_trap
green_low_connectivity_trap
odd_negative_false_terminal
lee_yang_root_collision_vs_confinement
tiling_lift_vs_half_tiling_descent
bulk_discrepancy_vs_core_witness
state_lift_H7_residual
scale_lift_unit_vs_chart_drift
```

The toy matrix has `13` rows, `20` columns, and full row rank over `GF(2)`:

```text
gf2_rank = 13
```

This does not prove anything yet; it is a design test.  It says that the
current sidecar vocabulary can plausibly separate the named residual rows
without forcing all information into one scalar.

Minimal proposed hitting sets have size `5`; examples include:

```text
endpoint_owner, dilation_grid, toeplitz_normal_fan_slack,
unit_equioscillation_index, ear_payload

endpoint_owner, green_resistance_slack, lee_yang_root_word,
unit_equioscillation_index, hensel_krasner_unit

boundary_cocircuit, dilation_grid, green_resistance_slack,
cech_betti_hole, ear_payload
```

The next proof target is sharper than this hitting-set toy:

```text
For the actual residual-pair rows, find a sidecar column basis B such that
each status-changing pair is separated, reconstructed, dual-annihilated,
descended, stopped at an AP/GW/Phi14 boundary, or routed to named debt.
```

This is the HYP-3048 matrix-atlas theorem target made concrete for the
HYP-3243/HYP-3244/HYP-3245/HYP-3248/HYP-3249 frontier.

## Discrete Morse / Lyapunov Descent

The second proposed route uses a vector-valued energy rather than a scalar
extremal:

```text
tope_boundary_depth
toeplitz_deficit_order
green_resistance_excess
lag_transport_cost
unit_equioscillation_defect
index_equality_defect
gauss_sum_word_depth
bu_forcing_gap_status
root_collision_depth
odd_negative_debt
chart_change_depth
hensel_unit_failure
state_lift_distance
```

Descent cases:

```text
strict_open_tope
  terminal safe witness

chebyshev_unit_equioscillation
  unit-group boundary packet or covering-floor handoff

index_theorem_prediction_boundary
  descriptive analytic/topological index equality plus floor/margin, contact
  holonomy, endpoint-chamber certificate, or named Borsuk-Ulam forcing gap

AP_or_GW_phi14_boundary
  closed equality atom, not a trap

phi14d_dilation_witness
  constructed equality on finer witness grid

finite_toeplitz_green_trap
  descend by moment-cone or electrical slack

lag_transport_outward_mass
  push outward surplus into AP support or trap sidecar

odd_negative_false_terminal
  retain parity/sign debt before positive compression

bulk_discrepancy_chamber
  Roth-Halasz bulk descent unless core witness takes over

local_scale_lift_chamber
  Hensel-Krasner unit prevents scale-lift drift

state_lift_H7_exit
  forbidden tournament-state contradiction or named debt
```

The missing theorem is an acyclic matching theorem on finite proof packets:

```text
Every noncritical chamber has a matched legal wall of strictly smaller
lexicographic energy, after all selected observability columns are attached.
```

This reframes HYP-3225/HYP-3236 trap discharge and integrates the new
HYP-3246/HYP-3247/HYP-3248/HYP-3249 Chebyshev/index frame.  Green resistance
and Toeplitz slack are not just measurements; they are candidate Morse
coordinates.  The HYP-3245 outward-lag transport is not terminal; it is a way
to locate the wall where the descent should occur.  The unit-equioscillation
index, binding complement-pair word, analytic/topological equalizer,
Gauss-sum index word, and Borsuk-Ulam forcing-gap flag are boundary/core
coordinates: they distinguish the unit-group Chebyshev packet and index
prediction from the covering-floor handoff.  HYP-3253 turns this into a
non-overclaim rule: index equality is a boundary label until paired with
S-dependent margin, contact holonomy, unit-nullspace blind-residue data,
magnitude-level exact-cover data, or endpoint-chamber descent data.
HYP-3238/HYP-3239 odd-negative payloads are not decorations; they stop false
positive-even terminal cells from being declared critical.

## Cross-Disciplinary Synthesis

Control and coding theory supply the central model:

```text
residual pairs = error classes
sidecars = syndrome bits / observability columns
proof-safe quotient = full observability on status-changing fibers
```

Optimization supplies the legal certificate language.  Missing sidecars are
missing dual variables in a Farkas, KKT, Schur-complement, or proof-circuit
certificate.

Oriented matroids and topology supply the finite state space: topes,
cocircuits, Cech bars, cover holes, and normal-fan cells.

Electrical networks supply the energy: Toeplitz curvature, positive
conductance, Green resistance, Thomson currents, Rayleigh bottlenecks, and
negative-covariance leakage.

Arithmetic dynamics supplies scale discipline.  Roth-Halasz is a bulk
discrepancy descent candidate; Hensel-Krasner is a local unit sidecar that
prevents scale-lift drift at the `7 x 2` apex.  These are useful only when the
row-level finite chamber and odd/sign coordinates are still visible.

Circuit complexity supplies a guardrail: the proof circuit should expose
missing-input vectors.  The bounded core stays below the small-degree wall
only when gates retain the sidecar they would otherwise forget.

## Tournament Analysis

Tournament vertices are proof angles and sidecar carriers, not runners, arcs,
or sectors.

Pairwise observable:

```text
which angle preserves more LRC proof obligations while destroying fewer
scarce coordinates
```

Switch/gauge:

```text
A -> B iff A has larger weighted obligation score; ties prefer fewer
destroyed coordinates and then stable key order
```

Fingerprint:

```text
score_hist={-24:1, 42:1, 53:1, 54:1, 56:1, 57:2, 60:1, 87:1, 96:1}
directed_3cycles=0
scc_sizes=[1,1,1,1,1,1,1,1,1,1]
hamiltonian_path_count=2
priority_path=
  finite_chamber_discrete_morse_descent
  -> sidecar_observability_matrix
  -> odd_negative_duality_resurrection
  -> green_toeplitz_thomson_energy
  -> lee_yang_ear_root_motion
  -> tiling_half_tiling_descent_span
  -> lag_autocorrelation_transport
  -> median_route_center_closure
  -> roth_halasz_hensel_lift_packet
  -> raw_single_scalar_extremality
```

## Assumption Challenge

Alternate vertices considered:

```text
runners
speed gaps
fixed circle sections
endpoint boundaries
wall crossings
residues
cover arcs
Fourier modes
matroid circuits
proof obligations
sidecar columns
finite chambers
state-lift exits
```

Chosen vertices: proof angles plus sidecar carriers.

Preserved LRC predicate: strict-open safety, closed-boundary equality,
finite trap discharge, compression legality, bulk/core handoff, and state-lift
exit status.

Destroyed by the tempting runner/arc quotient:

```text
endpoint owners
odd/negative sign
finite chamber address
tiling descent packet
root-motion word
Hensel unit
H=7 state-lift obligation
```

Challenged assumption: the remaining proof target is a scalar extremal.  The
better next target is a column-basis theorem plus an acyclic chamber-descent
theorem.

## Next Tests

1. Replace the toy residual rows by actual rows sampled from HYP-2963,
   HYP-3202, HYP-3225, HYP-3236, HYP-3243, HYP-3244, HYP-3245, HYP-3248,
   HYP-3249, HYP-3253, HYP-3255, HYP-3257, and HYP-3258 packets.
2. Compute row collisions under candidate column subsets; every collision must
   be a genuine same-status pair or a named missing sidecar.
3. Build the first finite Morse matching on HYP-3202 trap-neighborhood rows,
   using Toeplitz slack, Green resistance, lag transport, index-equality
   defect, Gauss-sum word depth, Borsuk-Ulam forcing-gap status, and
   odd-negative debt as coordinates.
4. Audit whether S289's Roth-Halasz and Hensel-Krasner columns separate real
   rows or only restate schematic coordination language.
