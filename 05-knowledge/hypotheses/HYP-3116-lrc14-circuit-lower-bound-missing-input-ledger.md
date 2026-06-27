---
id: HYP-3116
title: LRC14 circuit lower-bound and missing-input ledger
status: EVIDENCE / executable missing-input, endpoint-kernel, and proof-carrier synthesis; not a proof
source: codex-2026-06-27-S266
tangent: T1191
technique: LTI-252
tournament_technique: LTT-150
scripts:
  - 04-computation/lrc14_circuit_missing_input_ledger_codex_s266.py
  - 04-computation/lrc14_circuit_complexity_past_work_synthesis_codex_s266.py
results:
  - 05-knowledge/results/lrc14_circuit_missing_input_ledger_codex_s266.out
  - 05-knowledge/results/lrc14_circuit_complexity_past_work_synthesis_codex_s266.out
reflections:
  - lrc14-circuit-missing-input-ledger-codex-s266
  - lrc14-circuit-lower-bound-missing-input-ledger-codex-s266
  - lrc14-circuit-complexity-past-work-synthesis-codex-s266
related:
  - HYP-3119
  - HYP-3117
  - HYP-3115
  - HYP-3114
  - HYP-3113
  - HYP-3112
  - HYP-3111
  - HYP-3108
  - HYP-3107
  - HYP-2997
  - HYP-2991
  - HYP-2989
  - HYP-2981
  - HYP-2974
  - HYP-2791
  - HYP-2790
  - HYP-2744
  - HYP-2112
  - HYP-2108
  - HYP-3098
  - HYP-3083
  - HYP-3082
  - HYP-3077
  - HYP-3074
  - HYP-3054
  - HYP-3023
  - HYP-3016
  - HYP-2109
  - HYP-2963
  - THM-572
  - OPEN-Q-108
---

# HYP-3116: LRC14 Circuit Lower-Bound And Missing-Input Ledger

## Claim

Circuit complexity should be used in the LRC14 proof search as a
missing-input discipline, not as a metaphor for hardness.  A proposed shortcut
is proof-facing only after it names the Boolean/proof circuit, its input basis,
which inputs are essential, which certificate minterms close the route, and
which deleted inputs are reconstructible or paid for by sidecars.

The active synthesis connects three older repo threads:

1. **Data circuits for tournament invariants.**  The staircase/tiling Walsh
   work treats Hamiltonian-path counts as a low-degree Boolean/carry circuit
   over fixed-path tiles.
2. **Proof circuits for LRC14.**  HYP-3111/HYP-3115 treat the current frontier
   as a monotone certificate DAG whose outputs are finite-address,
   observer-gluing, root/ear, lattice, and algebraic fold exits.
3. **Endpoint-cover activation circuits.**  HYP-2108 gives the max gate
   `P(S)>0`, while HYP-2112 gives the sum gate `Phi(C)=G(v)`, exactly equal
   to the LRC gap on the tested circuit.

The live bridge is therefore sharper than a generic lower-bound analogy:
LRC14 circuit complexity is a kernel-exclusion problem for the endpoint-cover
activation circuit.  A low-depth quotient may be useful, but it cannot be
terminal unless it retains or reconstructs the `Phi`/`P` activation data and
feeds the HYP-3107 finite-address/observer-gluing proof interface.

S266 executable integration turns that reservation into a concrete
proof-carrier synthesis.  The useful invariant is:

```text
circuit_certificate_vector =
  (input_packet_schema,
   gate_basis,
   sidecar_closure,
   exact_gap_functional,
   route_purity,
   bridge_safety,
   uniform_family_parameter,
   terminal_exit)
```

This vector treats equinumerosity, equidecomposability, and equidistribution as
distinct circuit input types.  Cardinal/fiber shadows, scissors/branch
decompositions, and measure/root/Phi shadows may interact, but a proof circuit
may forget one only after reconstruction, dual annihilation, or named residual
debt is recorded.

## Executable Ledger

Artifact:

```text
04-computation/lrc14_circuit_missing_input_ledger_codex_s266.py
05-knowledge/results/lrc14_circuit_missing_input_ledger_codex_s266.out
```

The script records `13` proof gates / sidecar ledgers:

```text
endpoint_phi_sum_gap
endpoint_P_max_activation
lean_proof_frontier_dag
finite_address_observer_gluing
route_state_median_closure
boolean_mobius_low_depth_cut
full_boolean_mobius_hierarchy
fourier_toeplitz_psd_dual
haar_zipper_cocycle
lee_yang_root_ear_payload
minkowski_circuit_ising_bridge
endpoint_period_warning
finite_bank_apex7_literal
```

The top proof-payload carrier is HYP-2112's `endpoint_phi_sum_gap`:

```text
score = 149
vector = (exact_gap=5, kernel_exclusion=5, frontier_adjacency=4,
          sidecar_retention=4, quotient_leakage_control=5,
          formalizable=4, global_lift=4, uniformity=5,
          finite_checkability=4)
missing =
  finite_address_packet,
  observer_gluing_certificate,
  residual_kernel_exclusion_certificate,
  endpoint_period_numerator_sidecar
```

That ordering is the main new information from the search.  The old
endpoint-cover work is not merely related to circuit complexity; it is the
repo's concrete proof circuit.  `P(S)` is the max activation gate and
`Phi(C)` is the sum activation gate equal to the exact gap.  HYP-2791's
low-depth signed Boolean cut is a useful threshold subcircuit, but HYP-2790
warns that it leaks endpoint-period numerator and owner/phase coordinates
unless those are explicitly retained.

## Missing-Input Vector

The most frequent missing inputs across the ledger are:

```text
endpoint_cover_activation_vector: 6 gates
finite_address_packet:           5 gates
observer_gluing_certificate:      4 gates
endpoint_period_numerator_sidecar:3 gates
residual_kernel_exclusion_certificate: 2 gates
speed_owner_sidecar:              2 gates
proof_uniformity_schema:          2 gates
```

This suggests a compact packet extension:

```text
circuit_packet =
  (endpoint_cover_activation_vector,
   phi_gap_sum,
   phi_kernel_status,
   P_max_activation,
   simultaneous_resonance_winding_word,
   boolean_mobius_low_depth_cut,
   endpoint_period_numerator_sidecar,
   proof_circuit_missing_input_vector,
   proof_uniformity_schema,
   finite_address_packet,
   observer_gluing_certificate,
   walsh_degree_support_profile)
```

The lower-bound readout is proof-facing and finite:

1. Any shortcut omitting `endpoint_cover_activation_vector` has not reached
   the HYP-2108/HYP-2112 gap circuit.
2. Any shortcut omitting both `finite_address_packet` and
   `observer_gluing_certificate` does not feed the HYP-3107 proof frontier.
3. Any Boolean/type shortcut omitting `endpoint_period_numerator_sidecar`
   repeats the HYP-2790 scalar-transfer failure mode.
4. A terminal circuit must carry `Phi`/kernel data plus packet/gluing sidecars,
   or explicitly route the first missing input to named residual debt.

## Tournament Analysis

Vertices are proof carriers and sidecar ledgers, not runners, arcs, Boolean
gates, roots, or scalar values.  The pairwise observable is weighted
proof-payload retention: exact gap, kernel exclusion, HYP-3107 frontier
adjacency, sidecar retention, quotient-leakage control, formal readiness,
global lift, uniformity, and finite checkability.  The tie Hamiltonian path is
the script's declared gate order.

Fingerprint:

```text
score_hist = {0:1, 1:1, 2:1, 3:1, 4:1, 5:1, 6:1,
              7:1, 8:1, 9:1, 10:1, 11:1, 12:1}
directed_3cycles = 0
scc_sizes = [1,1,1,1,1,1,1,1,1,1,1,1,1]
hamiltonian_path_count = 1
edge_flips_against_novelty_gauge = 55
priority_path =
  endpoint_phi_sum_gap >
  lean_proof_frontier_dag >
  finite_address_observer_gluing >
  endpoint_P_max_activation >
  route_state_median_closure >
  haar_zipper_cocycle >
  fourier_toeplitz_psd_dual >
  boolean_mobius_low_depth_cut >
  lee_yang_root_ear_payload >
  full_boolean_mobius_hierarchy >
  endpoint_period_warning >
  minkowski_circuit_ising_bridge >
  finite_bank_apex7_literal
```

The `55` edge flips against a novelty-first gauge are the warning: fashionable
outside carriers and bounded-bank literals look attractive, but the proof
gauge pushes the exact endpoint activation circuit and Lean frontier above
them.

## Proof Route Shape

The proposed HYP-3116 route is:

```text
residual row
  -> HYP-2963/HYP-3107 packet key
  -> finite_address_packet + observer_gluing_certificate
  -> endpoint_cover_activation_vector
  -> P=max activation and Phi=sum activation
  -> either Phi>0,
     or a low-depth Boolean/Toeplitz/Haar certificate lower-bounds Phi,
     or the first missing input becomes named residual debt.
```

This is compatible with HYP-3111/HYP-3115.  Their
`proof_circuit_missing_input_vector` should now be read specifically as the
vector above, not just as a generic request for size/depth/uniformity data.

## Executable Synthesis

Artifact:

```text
04-computation/lrc14_circuit_complexity_past_work_synthesis_codex_s266.py
05-knowledge/results/lrc14_circuit_complexity_past_work_synthesis_codex_s266.out
```

The script records audited facts from local past work and gives each carrier a
common vector:

```text
(retained LRC predicate,
 exact gap readout,
 uniform family clarity,
 route purity,
 sidecar closure,
 bridge safety,
 proof readiness,
 destroyed-coordinate discharge,
 finite-bank warning)
```

Gate basis extracted from older work:

```text
HYP-2112: Phi(C)=mu(safe set), verified exactly 900/900 for n=6..14.
HYP-2108: P(S)>0 iff loose for endpoint-cover circuits.
HYP-2109: L/M/R automaton proposes M as the wall/tie corridor.
HYP-3023: magnitude_cocycle has 100.0% route purity after automatic fibers mix.
HYP-3077: 41 features, 34 Horn rules, 29791 triples, 0 illegal centers.
HYP-3082: raw scalar star has 5 naked bridges; protected branch graph has 0.
HYP-3111: monotone proof circuit has 10 inputs, size 8, depth 4.
HYP-3115: apex7_error<=5 isolates the max row only as a finite-bank warning.
```

## Tournament Analysis

Vertices are proof carriers and sidecar gates, not runners, arcs, Boolean
literals, Ising spins, or scalar values.

Pairwise observable: which carrier better retains the LRC predicate, gives an
exact gap readout, scales as a uniform family, keeps route purity, enforces
sidecar closure, protects against naked bridges, is proof-ready, discharges
destroyed coordinates, and avoids finite-bank fitting.

Switch: majority comparison of the observable vector.  Ties use weighted
payload and then the declared tie Hamiltonian path:

```text
endpoint_phi_relu_gap
endpoint_cover_circuit_positivity
automatic_magnitude_zipper
route_state_horn_median_hull
protected_branch_graph_no_naked_bridge
pde_weak_form_endpoint_compiler
mciq_monotone_proof_frontier
three_state_middle_automaton
finite_bank_apex_threshold_warning
raw_scalar_p0_shadow
```

Fingerprint:

```text
score_hist={0:1, 1:1, 2:1, 3:1, 4:1, 5:1, 6:1, 7:1, 8:1, 9:1}
directed_3cycles=0
scc_sizes=[1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
hamiltonian_path_count=1
edge_flips_against_smallest_circuit_first=38
```

Proof-payload order:

```text
protected_branch_graph_no_naked_bridge
> route_state_horn_median_hull
> endpoint_phi_relu_gap
> automatic_magnitude_zipper
> pde_weak_form_endpoint_compiler
> endpoint_cover_circuit_positivity
> mciq_monotone_proof_frontier
> three_state_middle_automaton
> finite_bank_apex_threshold_warning
> raw_scalar_p0_shadow
```

The `38` edge flips against the "smallest circuit first" gauge are the main
readout.  Small circuits are attractive, but the proof payload prefers larger
gates when they retain exact gap data, legal sidecar closure, or
no-naked-bridge protection.

## Proof Route

The next theorem-facing experiment is a row ledger over HYP-2963/HYP-3107:

```text
row_id
input_packet_schema
Phi_gap
P_sign
endpoint_owner_word
LMR_terminal_state
magnitude_cocycle
automatic_word
root_or_ear_payload
Horn_sidecar_closure
protected_branch_status
proof_depth_stage
finite_bank_literal_alarm
uniform_family_parameter
terminal_exit_or_named_debt
```

Candidate closure theorem:

> Every primitive residual row either passes through the exact Phi/endpoint
> gap gate, is split by the magnitude/route-purity gate, is legally closed by
> the Horn sidecar gate, reaches the protected no-naked-bridge terminal graph,
> or emits a named THM-572/F7 residual. Any row whose proof uses only a finite
> threshold or raw scalar must be rejected as an unsafe quotient.

## Guardrail

HYP-3115's one-literal finite-bank rule `apex7_error <= 5` is a signal, not a
uniform proof.  This hypothesis exists to prevent that kind of finite fitted
classifier from masquerading as a theorem route.  The proof-carrier scout
sharpens this: `apex7_error<=5` is useful only as a missing-input alarm until
the uniform family parameter and destroyed-coordinate discharge are part of the
circuit schema.  The endpoint-kernel script also ranks that literal last
despite its finite checkability:

```text
finite_bank_apex7_literal:
  score = 31
  missing =
    proof_uniformity_schema,
    endpoint_cover_activation_vector,
    finite_address_packet,
    observer_gluing_certificate,
    train_test_family_split
```

Incoming HYP-3117/mainline S266 also ran a broader monotone proof-DAG audit.
It models `LRC14Statement` as:

```text
direct_witness
OR ap_gw_boundary
OR (finite_address AND observer_gluing AND endpoint_owner AND uniformity
    AND one_of(root_ear, relation_lattice, component_bound,
               cocycle_exactness, state_lift, pde_weak_form))
```

That audit found all `12` inputs essential, tested ten tempting shortcuts,
closed `0/10`, and found broad missing-input frequencies:

```text
finite_address: 10
observer_gluing: 8
endpoint_owner: 7
uniformity: 5
root_ear_sidecar: 3
cocycle_exactness: 2
state_lift: 1
```

The endpoint-cover augmentation refines the same lesson: the broad proof-DAG
inputs must be paired with the exact `Phi`/`P` activation circuit before the
route touches a positive LRC gap.

## Next Work

- Attach the `circuit_packet` fields to a HYP-2963/HYP-3107 residual sample.
- Re-express HYP-2791's low-depth Boolean cut as a lower bound on `Phi`, not
  as a standalone Boolean scalar.
- Use HYP-2974/HYP-2981 and HYP-2991/HYP-2989 as optional dual/cocycle gates
  only when they emit endpoint owner and packet-route keys.
- Formalize a small Lean-facing record whose fields are
  `endpoint_cover_activation_vector`, `phi_gap_sum`, `phi_kernel_status`,
  `finite_address_packet`, and `observer_gluing_certificate`.
- Add `proof_circuit_missing_input_vector`, `input_basis`,
  `essential_input_set`, `minimal_certificate_minterm`, `repair_cover`,
  `reconstructible_coordinate_certificate`, and `terminal_exit_or_named_debt`
  to HYP-2963/HYP-3107/HYP-3098 packet rows.
- Prove a gate-elimination step: a residual packet either hits a known minterm,
  or its missing-input vector strictly decreases after a legal sidecar repair.
