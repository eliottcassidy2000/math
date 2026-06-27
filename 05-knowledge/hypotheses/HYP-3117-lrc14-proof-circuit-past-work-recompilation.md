---
id: HYP-3117
title: LRC14 proof-circuit past-work recompilation atlas
status: SYNTHESIS / executable scout and proof-compiler proposal; not a proof
source: codex-2026-06-27
script: 04-computation/lrc14_proof_circuit_past_work_scout_codex_20260627.py
result: 05-knowledge/results/lrc14_proof_circuit_past_work_scout_codex_20260627.out
reflection: 07-reflections/lrc14-proof-circuit-past-work-synthesis-codex-20260627.md
related:
  - HYP-3116
  - HYP-3115
  - HYP-3111
  - HYP-3112
  - HYP-3108
  - HYP-3107
  - HYP-3102
  - HYP-3016
  - HYP-2963
  - HYP-2962
  - HYP-2961
  - HYP-2115
  - HYP-2114
  - HYP-2112
  - HYP-2108
  - THM-573
  - OPEN-Q-108
---

# HYP-3117: LRC14 Proof-Circuit Past-Work Recompilation Atlas

## Claim

The useful circuit-complexity object in the LRC14 proof search is a uniform
proof-state compiler, not a finite fitted maximizer rule.  The old
endpoint-cover and additive-circuit work already supplies real proof gates,
while the newer packet, automaton, obstruction, Lee-Yang, and Lean frontiers
supply sidecars and legality checks.

This extends the upstream HYP-3116 missing-input ledger by supplying the
executable past-work gate map and a first residual compiler sketch.

Post-rebase integration: HYP-3116 now supplies the base monotone
missing-input circuit

```text
LRC14Statement =
  direct_witness
  OR ap_gw_boundary
  OR (finite_address AND observer_gluing AND endpoint_owner AND uniformity
      AND one_of(root_ear, relation_lattice, component_bound,
                 cocycle_exactness, state_lift, pde_weak_form))
```

with `input_count=12`, `gate_count=3`, `depth=3`, all inputs essential, and
`0/10` tempting shortcuts closed.  HYP-3117 supplies the older proof gates
that should implement or repair those missing inputs.  In HYP-3116's shortcut
audit, the dominant missing coordinates were `finite_address` (`10` hits),
`observer_gluing` (`8`), `endpoint_owner` (`7`), and `uniformity` (`5`); the
past-work map says the first places to seek repairs are `Phi/P`, hidden fold
nodes, magnitude cocycles, first obstruction cochains, and root/ear payloads.

Recompiled this way, a residual row should emit a circuit packet:

```text
proof_circuit_packet(S) =
  (input_basis_id,
   missing_input_vector,
   endpoint_cover_P_gate,
   Phi_gap_output_wire,
   fold_gate_depth,
   hidden_virtual_sum_nodes,
   automaton_magnitude_cocycle,
   Lee_Yang_ear_payload,
   relation_wall_class,
   first_obstruction_class,
   finite_address_or_observer_gluing_exit)
```

The proof target is not a small circuit by itself.  It is a fixed, uniform
circuit over legal inputs that either computes a positive LRC gap, routes the
row to a named finite-address or observer-gluing certificate, or exposes the
exact missing sidecar.

## Map 1: Past Work As Gates

The past-work search found the following reusable gates.

1. HYP-2108 gives the endpoint-cover circuit gate `P(S)`.  Tightness requires
   one integer `v` to resonate with all component midpoints at once, so the
   residual theorem is a simultaneous-resonance infeasibility statement.
2. HYP-2112 gives the exact output wire `G(v)=Phi(C)`.  `Phi>0` is an exact
   positive gap, while `ker Phi` is the tight/worry set.
3. HYP-2114/HYP-2115 say the hard additive boundary is the 3-term fold gate
   plus hidden virtual-sum nodes.  Raw 4-term energy is only a shadow unless
   the hidden sum node is retained.
4. HYP-2961/HYP-2962/HYP-2963 give a deterministic labelled-packet decision
   tree.  This is already circuit-shaped: every branch states the predicate it
   preserves and the data it destroys.
5. HYP-3016 says finite automaton or branching-program shadows are legal only
   with a magnitude/Farey cocycle, because coarse residue-language fibers mix
   boundary and strictly open rows.
6. HYP-3102 turns illegal forgetting into a first obstruction cochain.  This
   is the missing-input detector for the whole compiler.
7. HYP-3107 gives the Lean-facing terminal bus: the circuit must ultimately
   feed finite-address packets, observer-gluing certificates, or still-open
   theorem fields.
8. HYP-3109/HYP-3112 give the whole PGF root curve and one-runner ear payload.
   They should be treated as sidecar inputs, not as scalar `p0` replacements.
9. HYP-3111/HYP-3115 give the circuit-uniformity guardrail.  The finite-bank
   rule `apex7_error<=5` is a signal to test, not a theorem-grade circuit.

## Map 2: Residual Compiler

The proposed compiler route is:

```text
primitive_residual_row
  -> labelled_packet_decision_tree
  -> endpoint_cover_circuit_gate
  -> Phi_gap_output_wire
  -> fold_gate_virtual_sum
  -> automaton_branch_program_guard
  -> Lee_Yang_ear_payload_gate
  -> Minkowski_relation_wall_gate
  -> first_obstruction_cocycle_gate
  -> Lean_frontier_obligation_bus
  -> LRC14Statement or named residual debt
```

This route is not meant to be strictly linear in a future implementation.  It
is a circuit-basis ordering: each stage names one input, gate, or missing-input
check that a shortcut must preserve or discharge.

## Scout Results

Executable scout:

```text
04-computation/lrc14_proof_circuit_past_work_scout_codex_20260627.py
05-knowledge/results/lrc14_proof_circuit_past_work_scout_codex_20260627.out
```

The scout ranks eleven proof carriers by proof-readiness features: LRC
predicate retention, uniformity, exact measure, sidecar retention, terminal
exit, Lean-facing status, route purity, residual separation, missing-input
detection, counterexample filtering, and a penalty for raw scalar rules.

Tournament fingerprint:

```text
vertices = proof gates/sidecars, not runners or arcs
score_hist = {-3:1,18:1,19:2,20:1,21:1,22:1,24:1,25:1,28:1,31:1}
directed_3cycles = 0
scc_sizes = [1,1,1,1,1,1,1,1,1,1,1]
hamiltonian_path_count = 1
priority_path =
  first_obstruction_cocycle_gate
  -> lean_frontier_obligation_bus
  -> labelled_packet_decision_tree
  -> Phi_gap_output_wire
  -> automaton_branch_program_guard
  -> endpoint_cover_circuit_gate
  -> Minkowski_relation_wall_gate
  -> Lee_Yang_ear_payload_gate
  -> mciq_uniformity_guard
  -> fold_gate_virtual_sum
  -> raw_scalar_maximizer_rule
```

The top score for `first_obstruction_cocycle_gate` is the main lesson.  Circuit
complexity should audit illegal forgetting before optimizing for small
size/depth.  A shallow proof DAG is useful only after every input wire is a
legal proof sidecar.

## New Signals

Add these fields to the next HYP-2963/HYP-3107 residual-row ledger:

```text
proof_circuit_input_basis_id
proof_circuit_missing_input_vector
uniformity_guard_status
endpoint_cover_P_gate
Phi_gap_output_wire
fold_gate_depth
hidden_virtual_sum_count
automaton_fiber_mixing_bit
magnitude_cocycle_height
first_obstruction_class
Lee_Yang_ear_payload_mean_level
root_motion_reconstruction_status
relation_wall_class
sidecar_fanin_profile
minimal_certificate_depth
gate_route_purity
terminal_exit_kind
```

The two most important derived signals are:

```text
missing_input_vector =
  which named sidecars are absent after a proposed proof shortcut;

minimal_certificate_depth =
  shortest legal route from a residual packet to positive Phi/P,
  finite-address packet, observer-gluing certificate, or named state-lift debt.
```

## Assumption Challenge

Candidate tournament vertices considered:

```text
runners
gaps
fixed circle sections
section boundaries
wall-crossing events
residues
cover arcs
Fourier modes
matroid circuits
proof obligations
Boolean gates
automaton states
cocycle basis atoms
Lee-Yang root-motion ears
```

The selected vertices are proof gates and sidecar obligations.  This quotient
preserves residual-packet progress toward `LRC14Statement` through positive
`Phi/P`, finite-address packets, observer-gluing certificates, or named
state-lift debt.  It destroys raw row identity, runner order, raw time, some
owner history, and scalar maximizer labels unless a sidecar retains them.

## Bold Hypotheses

1. The remaining LRC14 proof obstruction is a missing-input vector, not a new
   scalar invariant.  A live residual should fail by naming exactly one absent
   sidecar: fold node, magnitude cocycle, first obstruction class, root-motion
   ear payload, relation wall, or observer-gluing exit.
2. `Phi` and `P(S)` are the old gates that should be closest to the terminal
   positive-gap exit.  The Lee-Yang/root and Minkowski/relation lanes should
   explain when the compiler reaches those gates, not replace them.
3. The finite-bank circuit rule from HYP-3115 becomes useful only after it is
   forced through a train/test residual bank with the input basis fixed in
   advance.  Any threshold that changes its basis is a sidecar-discovery
   heuristic, not a proof.
4. Hidden virtual folds are likely the recurring structure behind several
   current signals: additive 3-folds, ear payloads, relation walls, and
   first-obstruction cocycles may be different projections of the same
   "retained sum node" requirement.

## Next Work

1. Attach the new proof-circuit fields to a sample of HYP-2963 residual rows
   and HYP-3107 frontier nodes.
2. For each row, emit `missing_input_vector` and `minimal_certificate_depth`.
3. Reject any circuit classifier whose input basis changes between the
   anchored bank and the residual bank.
4. Test whether every live row reaches one of four exits:

```text
positive Phi/P gap
finite-address packet
observer-gluing certificate
named K33/F7/THM-572 or first-obstruction debt
```

If a row reaches none of these exits, HYP-3117 predicts that the missing field
is itself the next theorem target.
