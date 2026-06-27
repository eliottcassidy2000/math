# LRC14 Circuit Missing-Input Ledger - Codex S266

The useful circuit-complexity connection was hidden in older LRC work, not in
generic lower-bound language.  HYP-2108 and HYP-2112 already define the
proof-facing circuit: endpoint-cover ramp activations.  `P(S)` is the max gate
that detects a positive activation, and `Phi(C)` is the sum gate equal to the
exact LRC gap.

That changes the shape of the remaining proof.  The central question is not
"is LRC14 hard for small circuits?" but "can a residual packet make every
endpoint-cover activation vanish after the legal packet sidecars are
retained?"  In other words, the proof has become kernel exclusion for a
finite activation circuit.

## What The Ledger Found

The S266 script ranks `endpoint_phi_sum_gap` first, above the Lean proof
frontier and above the attractive newer Minkowski/circuit/Ising/De Moivre
source bridge.  This is the right correction: outside carriers help only when
they deliver the inputs needed by the exact gap circuit.

The repeated missing inputs are also precise:

```text
endpoint_cover_activation_vector
finite_address_packet
observer_gluing_certificate
endpoint_period_numerator_sidecar
residual_kernel_exclusion_certificate
speed_owner_sidecar
proof_uniformity_schema
```

The first three are the structural core.  The endpoint activation vector
connects HYP-2108/HYP-2112 to a positive gap.  The finite-address and
observer-gluing packets connect that gap to the HYP-3107 proof frontier.  The
endpoint-period numerator sidecar is the HYP-2790 warning: Boolean/type
shortcuts leak phase unless this coordinate travels with them.

## How This Augments The Proof Attempts

HYP-2791's low-depth signed Boolean cut should be reinterpreted as a possible
lower bound on `Phi`, not as a terminal Boolean proof.  Its atom cut is small
and formal, but it destroys height, owner, and phase data.  The repair is to
pair it with endpoint-period numerator and owner sidecars, then ask whether it
forces at least one positive endpoint activation or a positive `Phi` sum on
the residual family.

HYP-2974/HYP-2981 and HYP-2991/HYP-2989 fit the same pattern.  Toeplitz/Fejer
and Haar zipper certificates are dual or cocycle gates; they become proof
currency only when packet keys and endpoint-owner decodings are retained.

HYP-3111/HYP-3115's `proof_circuit_missing_input_vector` now has a concrete
payload.  It should include `endpoint_cover_activation_vector`, `phi_gap_sum`,
`phi_kernel_status`, `P_max_activation`, `endpoint_period_numerator_sidecar`,
`finite_address_packet`, and `observer_gluing_certificate`, not only size,
depth, and uniformity.

## Honest Status

This does not close LRC14.  It sharpens what remains: attach the endpoint
activation circuit to HYP-2963/HYP-3107 residual packets and prove that the
legal residual kernel is empty, or name the first missing input as residual
debt.

The finite-bank literal `apex7_error <= 5` is useful as a stress signal and
bad as a proof.  Its missing vector contains proof uniformity, endpoint
activation data, finite address, observer gluing, and family split data.  It
stays last in the proof-payload tournament despite being easy to check.
