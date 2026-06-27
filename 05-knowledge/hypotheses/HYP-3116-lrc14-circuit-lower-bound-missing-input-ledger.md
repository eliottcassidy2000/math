---
id: HYP-3116
title: LRC14 circuit lower-bound and missing-input ledger
status: EVIDENCE / executable missing-input and endpoint-kernel ledger; not a proof
source: codex-2026-06-27-S266
script: 04-computation/lrc14_circuit_missing_input_ledger_codex_s266.py
result: 05-knowledge/results/lrc14_circuit_missing_input_ledger_codex_s266.out
reflections:
  - lrc14-circuit-missing-input-ledger-codex-s266
  - lrc14-circuit-lower-bound-missing-input-ledger-codex-s266
tangent: T1191
technique: LTI-252
tournament_technique: LTT-150
related:
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
  - HYP-3074
  - HYP-3054
  - HYP-2963
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

## Guardrail

HYP-3115's one-literal finite-bank rule `apex7_error <= 5` is a signal, not a
uniform proof.  This hypothesis exists to prevent that kind of finite fitted
classifier from masquerading as a theorem route.  The endpoint-kernel script ranks that
literal last despite its finite checkability:

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
