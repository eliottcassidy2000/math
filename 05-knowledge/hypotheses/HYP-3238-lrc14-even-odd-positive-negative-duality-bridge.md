---
id: HYP-3238
title: LRC14 even/odd and positive/negative duality bridge
status: EVIDENCE / exact bounded-bank scout plus synthesis; not an LRC14 proof
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
  - HYP-3220
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

This is not a proof yet.  The exact scout below measures which rows look
falsely terminal under positive/even compression, especially rows with no
negative covariance leakage and rows that are local exchange traps.

## Exact Scout

The script recomputes the same anchored bounded bank as HYP-3236:

```text
E = {0} union A, A subset {1,...,14}, |A|=7
rows_all = 3432
rows_primitive = 3431
```

It adds the duality-facing coordinates:

```text
endpoint_bimodality = q0 + q6
odd_center_debt = max(0, q3(E) - q3(AP))
L_y = q0 + q6 + q3/10
positive Green slack = lambda2 / Kirchhoff / maxR deficits
negative leakage = negative_edges / negative_mass
```

The AP profile is:

```text
AP_q0 = 481/1470 = 0.327210884354
AP_q3 = 26/245 = 0.106122448980
AP_q6 = 1/49 = 0.020408163265
AP_q0_plus_q6 = 73/210 = 0.347619047619
AP_L_y = 2633/7350 = 0.358231292517
AP_lambda2 = 0.192033074001
AP_kirchhoff = 108.654718079151
AP_negative_edges = 0
```

HYP-3220 is not merely adjacent.  The scout records its parity certificate:

```text
de_moivre_power_sums p1..p8 = -1,5,-4,13,-16,38,-57,117
sign pattern = (-1)^k
dominant negative Perron period = -2cos(pi/7) = -1.801937735805
complement pairs = (1,6),(2,5),(3,4)
```

So even/odd and positive/negative are the same `Z/2` complement/parity
operator in the de Moivre chart.  HYP-3238 should therefore be read as a
single-duality packet with two coordinate languages, not as two parallel
analogies.

AP is uniquely primitive-tight for the main positive/even coordinates:

```text
L_y_MAX:
  all_beaters=0 primitive_beaters=0 all_ties=2 primitive_ties=1

endpoint_bimodality_q0_plus_q6_MAX:
  all_beaters=0 primitive_beaters=0 all_ties=2 primitive_ties=1

lambda2_MAX:
  all_beaters=0 primitive_beaters=0 all_ties=2 primitive_ties=1
```

The false-terminal warning is equally concrete:

```text
negative_edges_MIN:
  all_ties=20 primitive_ties=19 primitive_false_ties=18

negative_mass_MIN:
  all_ties=20 primitive_ties=19 primitive_false_ties=18

primitive_zero_negative_edges_rows = 19
primitive_zero_negative_edges_nonAP_false_terminals = 18
primitive_connected_positive_graph_nonAP_rows = 2754
```

Thus "no negative leakage" and "connected positive graph" are not proof
conditions.  They are sidecar diagnostics.  AP is still separated by Green
slack and `L_y`, but the negative/odd coordinate cannot be erased.

The ordered-tail exchange proxy is clean on the bank:

```text
primitive_rows_with_positive_q3_debt = 2879
q3_exchange_margin_violations_for_debt_rows = 0
Ly_AP_max_violations = 0
worst_endpoint_bimodality_gap_per_q3_debt =
  17161/12882 = 1.332168917870
  at E=(0,1,4,5,9,10,13,14)
```

Equivalently, the `L_y` dual prices all positive `q3` central debt by endpoint
bimodality loss.  This is the finite table where HYP-3204's odd-center
exchange rate meets HYP-3222/HYP-3220 sign/interlacing data.

The HYP-3202 trap boundary splits through the bridge as:

```text
non_AP_trap_bridge_classes =
  negative_leakage_plus_odd_q3_debt: 8
  odd_q3_debt_without_negative_covariance: 3
```

The `3` traps with no negative covariance leakage are important: even on the
finite trap boundary, negative leakage alone is not the odd sidecar.  Odd
`q3` debt and Green resistance excess remain live.

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

HYP-3220 sharpens this into an identity of dualities.  The de Moivre power
sums have sign exactly `(-1)^k` because the dominant period is the negative
Perron root `-2cos(pi/7)`.  The complement involution `x -> -x` swaps sector
pairs `(1,6),(2,5),(3,4)` and is simultaneously the positive/negative fold
and the even/odd parity operator.  For the LRC(2p) family, the sign obstruction
is the `p mod 4` / imaginary-quadratic wall.

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

Exact scout fingerprint:

```text
score_hist = {0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1,9:1,10:1,11:1}
directed_3cycles = 0
sccs = singletons
hamiltonian_path_count = 1
priority_path =
  even_positive_fejer_square
  -> green_dirichlet_positive_face
  -> odd_negative_brouwer_sign
  -> hermite_biehler_interlacing
  -> bulk_core_vitali_transfer
  -> pair_pascal_cap_mass
  -> toeplitz_normal_fan_slack
  -> negative_covariance_leakage
  -> signed_chart_change_debt
  -> raw_lambda2_scalar
  -> raw_connected_positive_graph
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
scout output by keeping exemplar rows and trap sidecars.

Challenged assumption: the positive/even side is the proof and the
odd/negative side is error.  The stronger reading is that the proof needs a
two-factor certificate: positive/even magnitude plus odd/negative sign or
core data.

## Next Work

1. Prove the `q3` exchange-rate inequality symbolically from the shell
   Delsarte dual or HYP-3204 ordered-tail exchange, not only by bounded-bank
   audit.
2. Attach HYP-3222 Hermite-Biehler interlacing directly to the three
   zero-negative-leakage exchange traps with positive `q3` debt.
3. Compare the HYP-3220 negative Perron sign with Fiedler vectors and Green
   bottleneck currents from HYP-3236.
4. Lift the bridge beyond the bounded bank: determine which pieces are
   exact finite k=8 facts and which persist as the LRC(2p) `p mod 4`
   family-law sidecar.
