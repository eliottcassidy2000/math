---
id: HYP-3238
title: LRC14 even/odd and positive/negative duality bridge
status: RESERVED / synthesis bridge; exact scout planned; not an LRC14 proof
source: codex-2026-06-28
tangent: T1338
technique: LTI-338
tournament_technique: LTT-238
script: 04-computation/lrc14_even_odd_positive_negative_bridge_codex_20260628.py
result: 05-knowledge/results/lrc14_even_odd_positive_negative_bridge_codex_20260628.out
reflection: 07-reflections/lrc14-even-odd-positive-negative-duality-bridge-codex-20260628.md
related:
  - HYP-3237
  - HYP-3236
  - HYP-3235
  - HYP-3234
  - HYP-3233
  - HYP-3232
  - HYP-3231
  - HYP-3230
  - HYP-3228
  - HYP-3227
  - HYP-3225
  - HYP-3224
  - HYP-3223
  - HYP-3222
  - HYP-3221
  - HYP-3219
  - HYP-3218
  - HYP-3217
  - HYP-3216
  - HYP-3214
  - HYP-3205
  - HYP-3204
  - HYP-3202
  - HYP-3201
  - HYP-3200
  - HYP-3163
  - HYP-3153
  - OPEN-Q-108
---

# HYP-3238: LRC14 Even/Odd And Positive/Negative Duality Bridge

## Claim

The current LRC14 proof frontier should be organized as a crossed duality,
not as a search for one final scalar:

```text
even / positive side:
  Fejer square, cyclotomic SOS, pair-Pascal cap mass,
  covariance layers, positive Green conductance, Perron coherent mode,
  bulk positive-measure equidistribution

odd / negative side:
  Worpitzky associator, Brouwer trace sign, Hermite-Biehler interlacing,
  negative covariance leakage, chart-change cancellation debt,
  measure-zero cyclotomic core witnesses
```

The conjectural proof packet is therefore:

```text
row E
  -> even_positive_certificate(E)
  -> odd_negative_sidecar(E)
  -> bulk_core_transfer(E)
  -> terminal AP witness or named residual debt
```

The key compression rule is the HYP-3201 no-free-slider rule in this specific
duality language:

```text
even/positive compression is legal only when the odd/negative coordinate is
zero, reconstructible, dual-annihilated, or explicitly retained as sidecar.
```

Thus HYP-3236's positive-part Green graph is a powerful certificate face, but
not a proof object by itself; HYP-3219's Brouwer sign and HYP-3222's
Hermite-Biehler even/odd interlacing are the complementary data that explain
which clipped coordinates may be discharged.

## Why This Extends The Previous Session

HYP-3236 produced a strong positive/Green coordinate:

```text
C_E -> G_+(E) -> L_E -> L_E^+ -> lambda2 / resistance / current entropy.
```

AP/consecutive and doubled AP are the only all-bank extremals for this
positive conductance face.  The dangerous arrow is still
`C_E -> G_+(E)`, because it clips negative covariance.  HYP-3238 names the
missing bridge:

```text
positive Green slack must either factor through the even Fejer/SOS face
or carry the negative/odd leakage that was clipped.
```

Incoming HYP-3237 supplies the bulk/core transfer law.  Bulk equidistribution
has positive measure information; AP is a measure-zero core where measure is
blind and cyclotomic arithmetic information turns on.  That is the same
information pattern as positive/negative covariance clipping: the coarse
positive channel becomes useful only when the missing core coordinate is
reintroduced.

Incoming HYP-3219 supplies the topological sign law.  The de Moivre cubic
obstruction factors as:

```text
non-SOS cubic obstruction = Brouwer/trace sign * SOS magnitude.
```

The magnitude is even/positive; the sign is odd/topological.  Therefore the
odd half should not be forced into an SOS certificate.  It should be carried
as a sign, degree, or interlacing sidecar until a dual certificate annihilates
it.

## Compression Failures Beyond Commutativity

The useful abstraction is a function audit.  A quotient
`q:X -> Y` may support a theorem about `f:X -> Z` only when `f` factors
through `q` or when the residual `H(f | q)` is named and discharged.

For LRC14, the relevant failures are:

```text
commutativity failure:
  unordered pair compression forgets ordered Worpitzky / exponent sidecars.

associativity failure:
  Schur complements, star-mesh reductions, and exchange moves compose only
  with boundary terminals and eliminated coordinates retained.

positivity failure:
  C_E -> max(C_E,0) gives a conductance graph but forgets negative covariance.

evenness failure:
  Fejer/SOS and covariance quadratic shadows forget the odd Worpitzky sign.

measure failure:
  bulk measure sees open lonely mass, but AP core witnesses are measure-zero.

scalarization failure:
  lambda2, total conductance, Perron alignment, or L_y can be AP-tight yet
  still forget chart, sign, core, or sidecar data needed for proof transport.
```

The proof object must therefore be a packet of functions, not a single
function value.

## Candidate Proof Packet

The bridge suggests three named components:

```text
even_positive_certificate:
  Fejer/de Moivre square, pair-Pascal cap mass, AP-support normal,
  Toeplitz moment cone, covariance distance layers,
  positive Green conductance and Rayleigh/Thomson slack.

odd_negative_sidecar:
  Worpitzky q3 center debt, associator cocycle, Brouwer trace sign,
  Hermite-Biehler odd polynomial, negative covariance leakage,
  signed-address chart-change cancellation debt.

bulk_core_transfer:
  Vitali wall, positive-measure equidistribution,
  Phi_14 closed core witnesses, Phi_7 equioscillation,
  cyclotomic arithmetic information replacing vanished measure information.
```

The theorem-facing target is:

```text
For every primitive non-AP k=8 row, either the even_positive_certificate has
strict AP slack, or the odd_negative_sidecar is nonzero and is discharged by
Brouwer/Hermite-Biehler/ordered-tail exchange, or the row belongs to a named
bulk/core residual packet.  Equality forces AP; all-bank equality also permits
the doubled AP dilation.
```

This is not a proof yet.  The planned exact scout should measure which rows
look falsely terminal under positive/even compression, especially rows with no
negative covariance leakage and rows that are local exchange traps.

## Relation To Earlier Packets

HYP-3222 gives the exact local even/odd polynomial clue:

```text
E(x) = x^2 + 5x + 4
O(x) = x^2 + 4x + 1
```

with strict interlacing and positive Wronskian.  This is the algebraic shape
expected of a legal even-positive / odd-negative gluing certificate.

HYP-3204 prices central `q3` debt by the `q0+q6` bimodality loss.  In the
bridge language, that is an exchange-rate map from odd center mass into an
even endpoint-cap deficit.

HYP-3235 and HYP-3218 make the even side more rigid: the cap lies in
`Q(cos(2*pi/7))` and the Fejer magic is a totally positive cyclotomic square.
HYP-3219 says the remaining obstruction is not lack of magnitude, but the
sign/degree factor.

HYP-3237 supplies the information-theoretic wall:

```text
measure information -> 0
cyclotomic arithmetic information -> on
```

The bridge is the statement that Green conductance, Fejer positivity, and
pair-Pascal mass live on the positive/even face, while Brouwer sign,
negative covariance, and core witnesses are the sidecar that prevents illegal
compression.

## Tournament Analysis

Vertices are proof obligations and retained information channels, not runners:

```text
even_positive_fejer_square
green_dirichlet_positive_face
bulk_vitali_positive_measure
pair_pascal_cap_mass
toeplitz_normal_fan_slack
odd_negative_brouwer_sign
hermite_biehler_interlacing
negative_covariance_leakage
signed_chart_change_debt
raw_lambda2_scalar
raw_positive_association
```

Pairwise observable: which carrier preserves the LRC14 predicate while losing
less of the crossed even/odd and positive/negative packet.

Switch/gauge: orient `A -> B` when `A` can discharge a destroyed coordinate of
`B` or retains both the positive/even certificate and the odd/negative sidecar.
If both preserve the same predicate, prefer the one with lower conditional
residual `H(lost_duality_payload | carrier)`.

Reserved fingerprint for the planned scout:

```text
score_hist = to be computed
directed_3cycles = to be computed
sccs = to be computed
tie_hamiltonian_path =
  even_positive_fejer_square
  -> green_dirichlet_positive_face
  -> bulk_vitali_positive_measure
  -> pair_pascal_cap_mass
  -> toeplitz_normal_fan_slack
  -> odd_negative_brouwer_sign
  -> hermite_biehler_interlacing
  -> negative_covariance_leakage
  -> signed_chart_change_debt
  -> raw_lambda2_scalar
  -> raw_positive_association
```

## Assumption Challenge

Alternate vertex sets considered: runners, gaps, fixed circle sections,
section boundaries, wall-crossing events, residues, cover arcs, Fourier modes,
Fiedler modes, negative covariance edges, Brouwer signs, Fejer coefficients,
matroid circuits, and proof obligations.

Chosen vertices: proof obligations / information channels.  This choice
preserves the crossed duality structure and directly records which quotient
was legal.  It destroys row identity, but row identity is recoverable in the
planned scout output by keeping exemplar rows and trap sidecars.

Challenged assumption: the positive/even side is the proof and the
odd/negative side is error.  The stronger reading is that the proof needs a
two-factor certificate: positive/even magnitude plus odd/negative sign or
core data.

## Next Work

1. Build the exact HYP-3238 scout on the same `3432` anchored rows as
   HYP-3236.
2. Count false terminals for positive/even compression: rows with no negative
   covariance leakage, rows with high Green connectivity, and HYP-3202 traps.
3. Measure exchange-rate links between odd central debt (`q3`), endpoint
   bimodality (`q0+q6`), Green slack, and negative leakage.
4. Update this hypothesis from reserved synthesis to evidence after the scout
   produces a bounded-bank table.
