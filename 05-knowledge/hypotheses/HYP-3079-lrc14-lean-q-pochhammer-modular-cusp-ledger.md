---
id: HYP-3079
status: SYNTHESIS / Lean q-series sidecar ledger; not a proof of LRC14
source: codex-2026-06-27-S247
tangent: T1163
tags: [lrc14, lean, formalization, q-series, q-pochhammer, modular-functions, eta, hurwitz, cusp, controlled-forgetting, tournament-analysis]
related:
  - HYP-3078
  - HYP-3075
  - HYP-3077
  - HYP-3076
  - HYP-3074
  - HYP-3073
  - HYP-3072
  - HYP-3071
  - HYP-3070
  - HYP-3066
  - HYP-3064
  - HYP-3062
  - HYP-3039
  - HYP-3026
  - HYP-2645
  - HYP-2627
  - HYP-2428
  - THM-572
  - OPEN-Q-108
---

# HYP-3079: Lean q-Pochhammer Modular Cusp Ledger for LRC14

This is the Lean-facing companion to HYP-3078's q-Pochhammer modular-cusp
principal-part gate.  It keeps the q-Pochhammer product as a modular-cusp
sidecar, not as a raw coefficient oracle, and records exactly which finite-tail
facts are currently formalized.

The basic dictionary is:

```text
(q;q)_inf = prod_{n>=1} (1-q^n)
eta(tau) = q^(1/24) (q;q)_inf,       q = exp(2*pi*i*tau)
Delta = eta^24 = q (q;q)_inf^24
j = E4^3 / Delta
```

The correction is that `eta` and `(q;q)_inf` are not full modular functions by
themselves.  They carry weight, multiplier, and product-tail data.  A modular
function for the full modular group must also be invariant under the full
`SL2(Z)` action and meromorphic at the cusp `infinity`.  In q-language that
means its Laurent expansion has only finitely many negative-power terms.

## Exact scout

Script:

- `04-computation/lrc14_q_pochhammer_modular_cusp_s246.py`
- output: `05-knowledge/results/lrc14_q_pochhammer_modular_cusp_s246.out`

HYP-3078's exact scout verifies through `q^30`:

```text
(q;q)_inf = 1 - q - q^2 + q^5 + q^7 - q^12 - q^15 + q^22 + q^26
Delta = q (q;q)_inf^24
1/Delta = q^-1 + 24 + 324 q + 3200 q^2 + ...
j = q^-1 + 744 + 196884 q + 21493760 q^2 + ...
```

Both `1/Delta` and `j` have principal part supported only at exponent `-1` in
this scout.  The important proof readout is not the displayed coefficients; it
is the sidecar split:

```text
q_pochhammer_tail
eta_multiplier_weight_balance
full_modular_group_ST_law
cusp_principal_part_finite
zero_or_pole_divisor_ledger
j_rational_function_exit
```

## Hurwitz role

Here "Hurwitz" should first be read as the complex-analysis zero-persistence
theorem: a locally uniform limit of nonzero holomorphic functions is either
nonzero or identically zero.  Since finite q-Pochhammer products converge
locally uniformly on `|q|<1` to a nonzero product, any proof that passes to a
limit must retain zero/pole divisor data and cannot silently create hidden
interior zeros.

This is adjacent to, but distinct from, the older Markov-Hurwitz and
`(2,3,7)` geometry lanes.  The full modular group has orbifold signature
`(2,3,infinity)`, so the cusp `infinity` changes the proof obligation from a
finite triangle symmetry bound into a finite principal-part ledger.

## LRC14 proof use

For the current LRC14 formalization frontier, this gives a new packet-sidecar
rule:

```text
raw q-series tail is legal only after transformation law,
finite cusp principal part, multiplier balance, and zero/pole persistence
are retained, reconstructed, dual-annihilated, descended, boundary-stopped,
or emitted as named THM-572/F7 debt.
```

The Lean implication is concrete.  The existing center-control packet shell
should not receive a q-series field alone.  It should receive a proof-bearing
record with at least:

```text
q_pochhammer_tail_id
eta_multiplier_balance_status
sl2z_transformation_status
cusp_principal_part_order
finite_negative_tail_proof
hurwitz_zero_persistence_status
j_rational_exit_status
formal_series_truncation_bound
```

Only after those fields are populated should a modular/q-series certificate be
allowed to serve as a route center or residual discharge.

## Lean ledger

The first Lean-facing layer is now:

- `04-computation/lean/TournamentH7/TournamentH7/LRCModularCuspLedger.lean`
- aggregate import: `TournamentH7.LRCModularCuspLedger`

This file does not assert the analytic modular-forms theorem.  It separates the
proved finite-tail glue from the open analysis:

```text
HasOnlyFiniteNegativePowers
FullModularCuspExpansionObligation
HurwitzQExpansionGate
qPochhammerEulerCoeffTo12_finite_negative_tail
jPrincipalCoeffS247_finite_negative_tail
invDeltaPrincipalCoeffS247_finite_negative_tail
```

The checked Lean statement for the full modular-function quote is:

```text
full modular invariance + meromorphic at the cusp
  -> finite negative q-tail
```

but this is represented as the named obligation
`FullModularCuspExpansionObligation`, not assumed as a theorem.

## Sixth-power collision merge

HYP-3076 supplies the separate sixth-power sidecar for

```text
a^6+b^6+c^6=d^6+e^6+f^6
a^6+b^6=d^6+e^6
```

S247 merges that idea into the Lean ledger by proving the safe face map:

```text
SixthPowerTwoCollision
  -> twoCollisionToThreeWithTail
  -> SixthPowerThreeCollision
```

The map adds the same sixth-power tail to both sides.  Thus a two-term
collision can enter the three-term ledger only as explicitly padded data; it is
not counted as a native six-slot support obstruction unless another packet
proves the padding legal.

## Tournament Analysis

The scout uses proof sidecars as vertices, not runners, bases, or raw
q-coefficients.

Pairwise observable:

```text
retained LRC status
SL2Z law
finite principal part
eta multiplier
zero divisor / Hurwitz persistence
Lean payload readiness
penalty for raw tail only
```

Fingerprint:

```text
score_hist = {0:1,1:1,2:1,3:1,4:1,5:1,6:1}
directed_3cycles = 0
scc_sizes = [1,1,1,1,1,1,1]
hamiltonian_path_count = 1
tie_hamiltonian_path =
  full_modular_function_packet
  > j_rational_exit
  > finite_principal_part_ledger
  > eta_multiplier_balance
  > hurwitz_zero_persistence_gate
  > q_pochhammer_tail
  > raw_q_coefficients
```

Assumption challenge: considered runners, gaps, residues, Fourier modes,
q-series coefficients, cusp principal parts, zero divisors, and proof
obligations.  The proof-obligation sidecars are the right vertices because raw
q-coefficients destroy transformation law, pole order, multiplier balance, and
zero/pole persistence data.

## Status

Open / synthesis.  HYP-3079 does not prove LRC14.  It makes the modular/q-series
import legally shaped for the current proof stack: q-Pochhammer is a product
tail, eta is a weighted/multiplier object, and a full-level modular-function
exit needs finite principal part plus transformation data.

Next pull: instantiate one non-tautological modular sidecar packet in the Lean
center-control interface from `TournamentH7.LRCModularCuspLedger`, probably as
a record that can hold `j`-rational exit data and a finite principal-part proof
without invoking a full modular-forms library.
