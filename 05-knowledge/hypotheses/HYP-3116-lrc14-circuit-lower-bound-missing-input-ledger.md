---
id: HYP-3116
title: LRC14 circuit lower-bound and missing-input ledger
status: EVIDENCE / executable missing-input audit; not a proof
source: codex-2026-06-27-S266
tangent: T1191
technique: LTI-252
tournament_technique: LTT-150
related:
  - HYP-3115
  - HYP-3114
  - HYP-3113
  - HYP-3112
  - HYP-3111
  - HYP-3108
  - HYP-3107
  - HYP-3098
  - HYP-3083
  - HYP-3074
  - HYP-3054
  - HYP-2997
  - HYP-2991
  - HYP-2963
  - OPEN-Q-108
script: 04-computation/lrc14_circuit_missing_input_ledger_codex_s266.py
result: 05-knowledge/results/lrc14_circuit_missing_input_ledger_codex_s266.out
reflection: 07-reflections/lrc14-circuit-lower-bound-missing-input-ledger-codex-s266.md
---

# HYP-3116: LRC14 Circuit Lower-Bound And Missing-Input Ledger

## Claim

Circuit complexity should be used in the LRC14 proof search as a
missing-input discipline, not as a metaphor for hardness.  A proposed shortcut
is proof-facing only after it names the Boolean/proof circuit, its input basis,
which inputs are essential, which certificate minterms close the route, and
which deleted inputs are reconstructible or paid for by sidecars.

The active synthesis connects two older repo threads:

1. **Data circuits for tournament invariants.**  The staircase/tiling Walsh
   work treats Hamiltonian-path counts as a low-degree Boolean/carry circuit
   over fixed-path tiles.
2. **Proof circuits for LRC14.**  HYP-3111/HYP-3115 treat the current frontier
   as a monotone certificate DAG whose outputs are finite-address,
   observer-gluing, root/ear, lattice, and algebraic fold exits.

The live conjectural bridge is that failed scalar proof routes are exactly
low-depth quotient circuits with a nonzero missing-input vector.

## Evidence To Gather

- Mine prior Boolean-circuit, Walsh, cocycle, finite-address, observer-gluing,
  Savitch, and route-center work for reusable gate types.
- Build a small executable ledger of known gates, shortcut attempts, essential
  inputs, and missing-input vectors.
- Run Tournament Analysis on proof-gate families rather than runners, Boolean
  gates, or scalar metrics.
- Use the result to update OPEN-Q-108: any new LRC14 shortcut must either
  provide a certificate minterm or state the first missing input.

## Guardrail

HYP-3115's one-literal finite-bank rule `apex7_error <= 5` is a signal, not a
uniform proof.  This hypothesis exists to prevent that kind of finite fitted
classifier from masquerading as a theorem route.

## S266 Executable Audit

Artifacts:

```text
04-computation/lrc14_circuit_missing_input_ledger_codex_s266.py
05-knowledge/results/lrc14_circuit_missing_input_ledger_codex_s266.out
07-reflections/lrc14-circuit-lower-bound-missing-input-ledger-codex-s266.md
```

The scout builds a deliberately small monotone proof circuit:

```text
LRC14Statement =
  direct_witness
  OR ap_gw_boundary
  OR (finite_address AND observer_gluing AND endpoint_owner AND uniformity
      AND one_of(root_ear, relation_lattice, component_bound,
                 cocycle_exactness, state_lift, pde_weak_form))
```

Readout:

- `input_count=12`, `gate_count=3`, `depth=3`.
- All `12` inputs are essential in the Boolean sense: each can flip the output
  under some assignment.
- Minimal certificate minterms are the two terminal exits
  `direct_witness`, `ap_gw_boundary`, plus six five-input observer exits:
  `finite_address AND observer_gluing AND endpoint_owner AND uniformity AND X`
  for `X` in the six sidecar families.
- The audit tests ten tempting shortcuts and closes `0/10`.  Total missing
  inputs across the shortcut list: `36`.
- Missing-input frequency:

```text
finite_address: 10
observer_gluing: 8
endpoint_owner: 7
uniformity: 5
root_ear_sidecar: 3
cocycle_exactness: 2
state_lift: 1
```

This is the strongest proof-route signal from the session: every audited
shortcut lacks finite-address data, and most also lack observer-gluing and
endpoint-owner data.  Therefore circuit complexity is acting as a lower-bound
language on proof compression: a route that deletes those coordinates must
prove reconstructibility, not merely show a scalar correlation.

## Mined Prior Work

The executable ledger mined these older carriers as proof gates:

- `staircase_walsh_carry_circuit`: the data-circuit lesson from fixed-path
  tournament tilings.  Useful for low-degree/carry intuition, but insufficient
  as proof without LRC packet inputs.
- `haar_zipper_cocycle` and `cocycle_normal_form`: quotient repair via local
  cocycles, coboundaries, endpoint owners, and state-lift residual naming.
- `observer_extension_cut_payload`, `observer_gluing_ledger`, and
  `finite_address_branch_closure`: the dominant proof-circuit gates.
- `route_state_closure_median`: prevents raw route labels from acting like
  legal centers before sidecars close.
- `lee_yang_savitch_ear_lattice`, `lee_yang_ear_payload`,
  `minkowski_circuit_ising_quintic`, and `pde_weak_form_packet`: root/ear,
  relation-lattice, uniformity, and operator sidecars.

## Shortcut Audit

The ten audited shortcuts are:

```text
raw_p0_scalar
one_literal_apex7_threshold
raw_H_gap_transfer
raw_minkowski_volume
raw_ising_energy
raw_demoivre_residual
raw_walsh_low_degree
raw_direct_component_count
raw_pair_pascal_shadow
raw_automaton_language
```

None is closed as stated.  Each receives a repair cover, for example:

- `one_literal_apex7_threshold` needs endpoint owner, finite address,
  observer gluing, and uniformity; it repairs through
  `observer_extension_cut_payload + finite_address_branch_closure`.
- `raw_minkowski_volume` needs finite address, observer gluing, root/ear data,
  and uniformity; it repairs through
  `observer_gluing_ledger + lee_yang_savitch_ear_lattice`.
- `raw_walsh_low_degree` needs cocycle exactness, endpoint owner, and finite
  address; it repairs through `haar_zipper_cocycle + finite_address_branch_closure`.

## Tournament Analysis

Vertices are proof-gate families:

```text
finite_address_branch_closure
observer_gluing_certificate
root_ear_payload_ledger
cocycle_exactness_ledger
route_center_median_ledger
pde_weak_form_operator
relation_lattice_minkowski
staircase_walsh_carry_data_circuit
finite_bank_threshold_signal
raw_scalar_p0
```

Pairwise observable: retained proof inputs, proof readiness, uniformity,
reconstructibility, prior-source support, and whether the carrier is only a
data circuit.

Fingerprint:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1,9:1}
directed_3cycles=0
scc_sizes=[1,1,1,1,1,1,1,1,1,1]
hamiltonian_path_count=1
priority_path=
  finite_address_branch_closure
  -> observer_gluing_certificate
  -> root_ear_payload_ledger
  -> cocycle_exactness_ledger
  -> route_center_median_ledger
  -> pde_weak_form_operator
  -> relation_lattice_minkowski
  -> staircase_walsh_carry_data_circuit
  -> finite_bank_threshold_signal
  -> raw_scalar_p0
```

## Next Work

Add `proof_circuit_missing_input_vector` to HYP-2963/HYP-3107/HYP-3098 packet
rows.  For each proposed proof shortcut, record:

```text
input_basis
essential_input_set
minimal_certificate_minterm
missing_input_vector
repair_cover
reconstructible_coordinate_certificate
required_sidecar_or_exit
terminal_exit_or_named_debt
```

The hoped-for closure route is not a larger scalar search.  It is a gate
elimination proof: show that a residual packet either hits a known minterm or
that the missing-input vector strictly decreases after a legal sidecar repair.
