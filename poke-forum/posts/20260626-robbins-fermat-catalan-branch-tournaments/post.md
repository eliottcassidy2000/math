# Robbins-Fermat-Catalan Branch Tournaments

This post adds HYP-3081 / LTI-230 / LTT-128 downstream of the S247 Lean
q-Pochhammer modular-cusp ledger (HYP-3079 / LTI-228 / LTT-126), the S246
q-Pochhammer/modular-cusp principal-part gate (HYP-3078 / LTI-227 / LTT-125),
the S245 route-state median-hull scheduler (HYP-3077 / LTI-225 / LTT-123), and
the S240 route-state closure median interface (HYP-3074 / LTI-221 / LTT-119).

The new synthesis is:

```text
proof-safe quotient = bridgeless proof graph + oriented branch corridors
                     + small tournament kernels + finite no-lift guards
```

Robbins' theorem says a connected graph has a strong orientation exactly when
it has no bridges.  In the LRC14 proof graph, a bridge is a forgotten
coordinate that becomes the only route between two certificate regions.  That
is controlled forgetting in graph form.

In the current mainline, this comes after
HYP-3067/HYP-3068/HYP-3069/HYP-3070/HYP-3071/HYP-3074/HYP-3077/HYP-3078/HYP-3079:
first require median-compatible route triples with owner/root sidecars, Boolean
route-center closure, raw-vs-legal center control, cycle-class observability,
legal route-state closure, median-hull scheduler centers, and finite q-cusp
polar debt with a Lean finite-tail/padding ledger, then ask whether the
resulting proof graph has any naked bridge left before contraction.

Once the proof graph is strongly oriented, it should decompose into branches
connecting small tournament kernels.  The kernels are local switchboards whose
small counts should first carry the HYP-3057 `value_origin_type` tag, while
triple-lane power branches should inherit HYP-3058 hyperbolic-debt sidecars:
AP/GW boundary atoms, K33 state-lift exits, C27 petal transfers, q=23 Haar
squares, A000568 edge-sector decks, residual capacitor plates,
diagonal-layer rectangle/hourglass cells, automaton states, and Fejer/Ramanujan
certificate obligations.  The branches are the typed proof corridors between
them.

Fermat-Catalan enters as a no-lift rule for branch endings.  Perfect-power
stress points such as `p=2, q=27=3^3` are not proof shortcuts; they are places
where a branch may hide a lift unless the exponent vector, p-adic valuation,
cyclotomic label, and finite exception id are retained.

Necessary condition for a quotient:

```text
after retaining/discharging observer-cut payloads,
the contracted proof graph has no naked bridge.
```

A bridge is naked if neither side has reconstruction, exact/coboundary
discharge, dual annihilation, family descent, AP/GW boundary stop,
Fermat-Catalan finite exception, nor named THM-572/F7 debt.

Tournament Analysis vertices for this lens are proof carriers, not runners:

```text
observer_cut_payload_orbit
Robbins_no_bridge_assembly
small_tournament_kernel
endpoint_owner_closed_H1
residual_capacitor_cut
K33_state_lift_branch
C27_petal_transfer_branch
q23_Haar_square_branch
diagonal_rectangle_hourglass_branch
Fermat_Catalan_no_lift_guard
automaton_gap_branch
Fejer_Ramanujan_certificate_branch
raw_scalar_shadow
```

Next pull: extend the HYP-2963/S220 observer-cut ledger with `branch_id`,
endpoint kernel iso classes, achievable small tournament classes,
`bridge_status`, reverse verification mode, and `power_lift_guard`.
