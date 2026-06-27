# LRC14 Lean q-Pochhammer Modular Cusp Ledger

**codex-2026-06-27-S247.**  User asked to merge q-Pochhammer series, the full
modular-group q-expansion rule, and a Hurwitz-theorem lens into the LRC14 proof
work.

The useful correction is a type separation:

```text
(q;q)_inf             product tail
eta                   product tail plus weight/multiplier
full modular function transformation law plus finite cusp principal part
j                     clean full-level modular-function exit
```

The upstream HYP-3078/S246 exact scout verifies the q-series dictionary:

```text
(q;q)_inf = 1 - q - q^2 + q^5 + q^7 - q^12 - q^15 + q^22 + q^26 + ...
Delta = q (q;q)_inf^24
1/Delta = q^-1 + 24 + 324 q + 3200 q^2 + ...
j = q^-1 + 744 + 196884 q + 21493760 q^2 + ...
```

The principal-part point is the live one.  For a full-level modular function,
the q-expansion is allowed an infinite positive tail but only finitely many
negative powers.  That gives the exact controlled-forgetting sidecar:

```text
q_pochhammer_tail
eta_multiplier_weight_balance
full_modular_group_ST_law
cusp_principal_part_finite
hurwitz_zero_persistence_status
j_rational_exit_status
```

Hurwitz should be used carefully.  The complex-analysis Hurwitz theorem gives a
noncollapse guard for the q-product limit: finite products converge locally
uniformly on `|q|<1` to a nonzero product, so a proof cannot introduce or delete
interior zeros without a divisor ledger.  This is related to the older
`(2,3,7)`/Markov-Hurwitz lanes only by analogy; the full modular group has
signature `(2,3,infinity)`, and the cusp changes the obligation into a finite
principal-part record.

For the Lean frontier, this means a future packet should not contain a bare
q-series field.  It should contain transformation status, principal-part order,
finite-negative-tail proof, eta multiplier status, and zero/pole persistence
status.  Then it can be linked to `CenterControlPacket` as a non-tautological
certificate source.

The first Lean ledger now exists at
`04-computation/lean/TournamentH7/TournamentH7/LRCModularCuspLedger.lean`.
It checks only the finite-tail glue and sidecar typing:
`HasOnlyFiniteNegativePowers`, `FullModularCuspExpansionObligation`,
`HurwitzQExpansionGate`, the S247 q-Pochhammer/`j`/`1/Delta` principal-part
readouts, and the sidecar-vertex count.  It also merges the sixth-power prompt
equations by proving that a `a^6+b^6=d^6+e^6` collision maps into the
`a^6+b^6+c^6=d^6+e^6+f^6` ledger only by adding the same sixth-power tail to
both sides.  That keeps HYP-3076's native-vs-padded distinction intact.

Assumption challenge: I considered runners, gaps, fixed sections, residues,
Fourier modes, q-coefficients, cusp principal parts, zero divisors, and proof
obligations.  Raw q-coefficients preserve too little.  The modular proof object
is the proof-obligation sidecar that says exactly why the infinite tail is
legal and why the negative tail is finite.
