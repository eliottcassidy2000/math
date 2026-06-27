# LRC14 Branch-Kernel Orientation Audit

This post adds HYP-3082 / LTI-231 / LTT-129.

The S250 script runs the S249/HYP-3081 Robbins no-bridge orientation idea,
downstream of the S247/HYP-3079 Lean q-cusp finite-tail ledger, S246/HYP-3078
q-cusp principal-part gate, S245/HYP-3077 median-hull scheduling, and
S240/HYP-3074 route-state closure, on the existing HYP-2963
packet bank.  It audits `21913` packets, including `7235` hard non-AP/GW
packets, and finds `0` candidate F7 rows in this bounded bank.

The raw scalar graph is unsafe:

```text
nodes=6 edges=5 bridges=5 naked_bridges=5
```

After attaching the known branch sidecars, the protected graph is safe in the
bounded audit:

```text
nodes=80 edges=83 bridges=69 naked_bridges=0
contracted_core_nodes=1
contracted_core_bridges=0
strong_orientation_exists=True
```

The protected exits are not one scalar reason.  They include direct
q-witnesses, AP/GW boundary stops, C27 owner-strip discharge, K33 state-lift
debt, section exits, Haar grid exits, positive open/nested-refinement exits,
power-lift sidecars, route kernels, Desargues/Beal finalizer gates, and named
residual debt.

Branch ledger:

```text
Q-WITNESS                rows=14676 exit=direct-q-witness
COVERING-MOMENT          rows= 6040 exit=strict-safe-component
COVERING-MOMENT          rows= 1188 exit=strict-safe-component
BOUNDARY-PETAL-SPORADIC  rows=    4 exit=petal-discharge
K33-STATE-LIFT           rows=    3 exit=state-lift-debt
BOUNDARY-AP-GW           rows=    2 exit=boundary-equality
```

Tournament Analysis vertices are branch carriers, section exits, gates, and
proof debts, not runners.  The retained-payload tournament is transitive with
one Hamiltonian path:

```text
observer_cut_payload_orbit
> Robbins_no_bridge_assembly
> labelled_packet_bank
> residual_section_exit
> Haar_grid_exit
> endpoint_owner_strip
> C27_petal_branch
> K33_state_lift_branch
> covering_moment_branch
> Fermat_Catalan_no_lift_guard
> Desargues_Beal_finalizer_gate
> named_residual_debt
> raw_scalar_shadow
```

What this proves: in the audited HYP-2963 bank, every hard non-AP/GW packet has
a declared section/grid exit or named state-lift debt before the Robbins bridge
test, and the protected branch graph has no naked bridge.

What this does not prove: the global reduction to the HYP-2963 bank, THM-572,
K33 discharge, or the family theorem for all covering tails.  The next proof
obligation is now precise: prove every primitive residual reaches this branch
graph, then discharge K33/THM-572 and covering-family exits without creating a
new naked bridge.
