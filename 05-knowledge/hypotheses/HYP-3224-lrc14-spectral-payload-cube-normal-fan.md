---
id: HYP-3224
title: LRC14 spectral payload cube and normal-fan support certificate
status: EVIDENCE / exact bounded-bank synthesis; not an LRC14 proof
source: codex-2026-06-28
tangent: T1307
technique: LTI-307
tournament_technique: LTT-207
script: 04-computation/lrc14_spectral_payload_cube_codex_20260628.py
result: 05-knowledge/results/lrc14_spectral_payload_cube_codex_20260628.out
reflection: 07-reflections/lrc14-spectral-payload-cube-normal-fan-codex-20260628.md
related:
  - HYP-3236
  - HYP-3235
  - HYP-3234
  - HYP-3233
  - HYP-3232
  - HYP-3231
  - HYP-3230
  - HYP-3229
  - HYP-3228
  - HYP-3227
  - HYP-3226
  - HYP-3225
  - HYP-3223
  - HYP-3222
  - HYP-3221
  - HYP-3218
  - HYP-3217
  - HYP-3214
  - HYP-3216
  - HYP-3213
  - HYP-3212
  - HYP-3211
  - HYP-3210
  - HYP-3205
  - HYP-3204
  - HYP-3203
  - HYP-3202
  - HYP-3201
  - HYP-3200
  - HYP-3163
  - HYP-3162
  - HYP-3161
  - HYP-3160
  - HYP-3154
  - HYP-3153
  - HYP-3152
  - HYP-3151
  - HYP-3150
  - HYP-3138
  - HYP-3132
  - OPEN-Q-108
---

# HYP-3224: LRC14 Spectral Payload Cube And Normal-Fan Support Certificate

## Claim

The strongest current k=8 signals should be read as projections of one
normal-fan certificate, not as unrelated scalar maxima:

```text
AP support hyperplane
Caratheodory-Toeplitz lambda_min moment-cone interior
cyclic-distance covariance layers D1,D2,D3
total covariance Sigma_kappa2
finite exchange-trap discharge
ordered-tail exchange-rate face
HYP-3222 Perron/Joukowski/Hermite-Biehler exact-certificate sidecar
Chebyshev/Cohn-Elkies magic-function dual
```

The proof target is to construct the dual certificate whose visible faces are
these signals.  In this view, AP/consecutive is an exposed point of a payload
polytope or cone, and the doubled AP row is the known nonprimitive equality
artifact.

## Exact Scout Readout

The scout recomputes every signal over the anchored bounded k=8 bank

```text
E = {0} union A, A subset {1,...,14}, |A|=7
rows_all = 3432
rows_primitive = 3431
```

For consecutive speeds,

```text
q = (481/1470, 359/1470, 25/147, 26/245, 17/210, 5/98, 1/49)
AP_support = 39766/540225
Toeplitz_lambda_min = 0.042304730706
Sigma_kappa2 = 6237419/8643600
D1,D2,D3 =
  308509/1080450,
  547577/2160900,
  225577/1234800
Perron_alignment = 0.996784201609
bimod = q0+q6 = 73/210
bimod_plus_q3 = 667/1470
```

Consecutive and doubled AP tie for maximum coverage, `L_y`,
`q0+q6`, `q0+q6+q3`, `Sigma_kappa2`, all three covariance layers, AP support
projection, and Toeplitz `lambda_min`.  Raw cyclotomic energy is again a
guardrail: consecutive has rank `19` for minimum energy, so "closest to
uniform" is not the certificate.

## Main New Evidence

The payload cube uses the six max-facing metrics

```text
AP_support
Toeplitz_lambda_min
D1
D2
D3
Sigma_kappa2
```

and computes the all-bank Pareto skyline.  The result is:

```text
pareto_skyline_size = 2
pareto_skyline_rows =
  (0,1,2,3,4,5,6,7)
  (0,2,4,6,8,10,12,14)
```

This is the cleanest synthesis so far.  The AP support functional, Toeplitz
moment-cone margin, distance-layer inequalities, and total covariance all
select exactly the same visible face.  The conjectural theorem should be
phrased as a normal-cone or Delsarte/Farkas support certificate:

```text
for every legal row E,
  payload(E) lies below the AP normal cone;
equality forces AP or named dilation.
```

## Trap Discharge

HYP-3202 found `11` primitive arbitrary-exchange local maxima besides
consecutive.  HYP-3224 tests all of them against the full payload cube.

Readout: every nonconsecutive exchange trap has a strict Toeplitz
`lambda_min` deficit.  In the scout's normalized gap classifier, all `11`
traps discharge through the Caratheodory-Toeplitz moment-interior signal.

This suggests a useful division of labor:

```text
exchange-gradient route:
  prove bulk improvement away from the finite trap manifold

moment-cone route:
  discharge the finite traps by Toeplitz lambda_min / Schur / Verblunsky
```

So the trap manifold is not merely an exception list.  It is a boundary chart
where local covariance compression fails but the moment-cone interior sees the
missing curvature.

## Ordered-Tail Exchange Face

Post-rebase incoming work added HYP-3204, which tests a coefficient-level
route to the same `L_y=q0+q6+q3/10` target.  HYP-3224 now reproduces that
check inside the spectral payload scout:

```text
(q3 - q3_consec)_+ <= (q0+q6)_consec - (q0+q6)
primitive rows with q3 gain = 2879
exact violations = 0
worst ratio = 12882/17161
worst row = (0,1,4,5,9,10,13,14)
```

This makes HYP-3204 the coefficient face of the same normal-fan cube.  The
central odd/Worpitzky mass `q3` may increase, but it is priced by a larger
loss in ordered-state bimodality `q0+q6`.  Therefore the payload diagram now
has a concrete coefficient chart:

```text
q0+q6 bimodality atom
+ central exchange-rate lemma
+ upper ordered-tail barrier
=> L_y extremality
```

This chart should join the Toeplitz trap-discharge chart rather than replace
it.

## Incoming Arithmetic Synthesis

The same fetches that introduced HYP-3204 also added HYP-3205, HYP-3211,
HYP-3212, HYP-3213, HYP-3221, and the rebase-added HYP-3222 exact-certificate
packet.  They sharpen the interpretation of the dual certificate:

- HYP-3205 is the spectral dictionary compatibility layer: AP and doubled AP
  are the only simultaneous maximizers across the dictionary, and primitive
  normal form leaves AP unique.  HYP-3224 treats that dictionary intersection
  as one face of a larger normal-fan/magic-function certificate.
- HYP-3222 supplies the exact arithmetic gluing layer: the minimal even fold
  `v^2-5v+4` becomes `E(x)=x^2+5x+4`, the odd Worpitzky leg is
  `O(x)=A_3(x)=x^2+4x+1`, the roots strictly interlace, and
  `E O'-E' O=(x+3)^2+2>0`.  Its ideal C6 Perron quotient also turns the
  HYP-3202 layer sums into `lambda0=(1^T C 1)/6=lambda_max`.  HYP-3224
  treats those as exact faces that must be glued to AP support, Toeplitz
  curvature, and ordered-tail exchange rather than as standalone scalars.
- HYP-3212 proves `V_7(u)-2=(u-2)(u^3+u^2-2u-1)^2`, so the de Moivre cubic
  appears as a Chebyshev equioscillation double root.  This is exactly the
  kind of extremal/magic-function dual HYP-3224 is looking for.
- HYP-3213 packages the cap/dip as arithmetic in `Q(cos 2pi/7)`: cyclotomic
  ideal, Chebyshev equioscillation, rationality, discriminant `7^2`,
  Joukowski trace, and Delsarte/Cohn-Elkies magic function are one field.
- HYP-3221 warns that config-blind algebraic certificates stall at the
  apex-7 Lee-Yang obstruction and that the final closure must carry
  config-specific analytic/equidistribution information.
- HYP-3223 supplies two downstream certificate languages for this packet's
  finite trap-discharge face: a Green-current/effective-resistance reading of
  the covariance matrix, and a Lorentzian/valuated-exchange reading of
  co-emptiness finite differences.  HYP-3224 should treat those as candidate
  additional coordinates of the same normal fan, especially for classifying
  the `11` non-AP exchange traps.
- HYP-3236 executes the Green-current/effective-resistance candidate.  The
  positive covariance conductance graph adds an algebraic-connectivity face to
  the payload cube: AP and doubled AP are the only rows maximizing `lambda2`
  and total positive conductance and the only rows minimizing Green resistance
  profiles.  All `11` non-AP exchange traps discharge by Green resistance
  excess, split between Kirchhoff excess and max-resistance bottlenecks.  This
  new face is legal only with negative covariance leakage and odd
  Worpitzky/Hermite-Biehler sidecars retained.
- HYP-3225 supplies the trap-local refinement of the same face.  Toeplitz
  `lambda_min` remains the universal first discharge on the `12` local maxima,
  while residual mechanisms split into rank-2 Plucker, Green low-connectivity,
  AFM/Rayleigh, and mixed sidecars.  HYP-3236 should be read as the global
  conductance face, HYP-3225 as its finite-boundary taxonomy, and HYP-3226 as
  the motif atlas that should ingest both.
- HYP-3214 identifies the 7-sector cyclotomic magic function as the Fejer
  kernel `F_7=(de Moivre cubic)^2`, positive-definite and equal to AP
  autocorrelation, while honestly separating it from the 14-clock Johnson
  pair-Pascal cap.  This gives the payload cube a harmonic face adjacent to
  HYP-3236's Green Dirichlet face: Fejer controls autocorrelation/coverage,
  Green controls finite covariance conductance, and the cap remains a distinct
  pair-normalized coordinate until a dual certificate glues them.
- HYP-3211 separates the additive cyclotomic face from the multiplicative
  octonion/G2 face; for HYP-3224 this says the relevant normal cone is
  additive/cyclotomic, while octonion numerology is a productive negative
  unless it preserves the q-vector payload.

Thus the best current target is more specific than "find a Toeplitz dual":
construct the level-7 Chebyshev/Cohn-Elkies/Delsarte magic function whose
finite shadows are AP support, Toeplitz lambda-min, ordered-tail exchange,
covariance layers, and Hermite-Biehler interlacing.

## Correlation And Non-Equivalence

The signals are aligned but not identical:

```text
AP_support vs Toeplitz_lambda_min: pearson 0.605024, spearman 0.576158
AP_support vs Sigma_kappa2:        pearson 0.777602, spearman 0.731513
Toeplitz_lambda_min vs Sigma_k2:   pearson 0.858314, spearman 0.866246
Sigma_kappa2 vs Perron_alignment:  pearson 0.788154, spearman 0.669266
AP_support vs cyclotomic_energy:   pearson -0.323500, spearman -0.251094
```

This matters.  The route should not collapse the proof to one scalar.  The
right object is a commuting diagram with different faces:

```text
root-locus face      -> AP support
moment-cone face     -> Toeplitz lambda_min
coefficient face     -> q0+q6 exchange-rate pricing for q3
ferromagnetic face   -> D1,D2,D3 and Sigma_kappa2
spectral face        -> Perron alignment / HB interlacing
local-move face      -> exchange traps and trap discharge
```

## Creative Proof Angles

### Angle A: Moment-cone normal fan

Try to prove that the AP support hyperplane is the normal cone of the
Caratheodory-Toeplitz moment cone at the consecutive row.  In practical terms,
seek a PSD Toeplitz dual, Christoffel function, Fejer-Riesz factor, Schur
parameter, or Verblunsky certificate whose slack is exactly the AP support
gap plus named dilation equality.

This is stronger than "Toeplitz lambda_min is maximal."  It asks for the dual
functional behind the maximum.

### Angle B: Finite trap curvature

Use HYP-3202's exchange graph only until it reaches the `11` traps, then stop
using exchange monotonicity.  At the traps, switch to Toeplitz curvature.
The traps should become the finite boundary of a proof by charts:

```text
bulk chart: exchange/covariance improvement
boundary chart: Toeplitz lambda_min deficit
odd chart: Worpitzky/Hermite-Biehler interlacing debt
```

### Angle C: Chebyshev/Cohn-Elkies magic function

HYP-3212/HYP-3213 turn the normal-cone hunt into a named classical object:
the de Moivre cubic is the double-root factor of the seventh Chebyshev
equioscillation polynomial, and the cover bound should be a one-dimensional
level-7 Delsarte/Cohn-Elkies magic-function certificate.  This is the most
concrete candidate dual for the AP support hyperplane.

### Angle D: HYP-3222 Hermite-Biehler / Perron gluing

HYP-3210 says the LRC cover bound and tournament TRRT are one
Hermite-Biehler theorem under Joukowski; HYP-3222 makes the minimal local
certificate exact by giving the interlacing `E/O` legs and the C6 Perron
quotient.  HYP-3224 gives the finite spectral payload that the gluing theorem
should consume: AP support, Toeplitz interiority, covariance layers, trap
curvature, ordered-tail exchange, and Perron alignment are the measurable
shadows of the even/odd interlacing slack.

## Tournament Analysis

Vertices are proof currencies, not runners, arcs, or sectors:

```text
normal_fan_support_certificate
moment_cone_toeplitz_lambda_min
covariance_distance_layer_certificate
spectral_trap_discharge
joukowski_hermite_biehler_transport
perron_uniform_mode_alignment
exchange_gradient_bulk
raw_total_covariance
raw_cyclotomic_norm
raw_left_compression
```

Pairwise observable: retained LRC proof payload under quotienting.

Switch/gauge: orient toward the currency that keeps more support, moment,
layer, spectral, and sidecar payload.

Fingerprint:

```text
directed_3cycles = 0
SCC_sizes = [1,1,1,1,1,1,1,1,1,1]
Hamiltonian_path_count = 1
priority_path =
normal_fan_support_certificate
-> moment_cone_toeplitz_lambda_min
-> covariance_distance_layer_certificate
-> spectral_trap_discharge
-> joukowski_hermite_biehler_transport
-> perron_uniform_mode_alignment
-> exchange_gradient_bulk
-> raw_total_covariance
-> raw_cyclotomic_norm
-> raw_left_compression
```

## Assumption Challenge

Alternate vertices considered: proof currencies, moment-cone faces,
covariance layers, exchange traps, sidecar obligations, and spectral modes.

Preserved predicate: k=8 AP/consecutive coverage and covariance extremality.

Destroyed information: raw q-vector order, trap identity, dilation status,
and odd Worpitzky/Hermite-Biehler interlacing unless explicit sidecars are
kept.

Challenged assumption: one scalar, one compression move, or one tournament
vertex model is the proof object.  The stronger hypothesis is that LRC14's
k=8 node is a small normal-fan diagram, and the proof must retain enough
coordinates to see that diagram.

## Green-Current Addendum: HYP-3227

HYP-3225 classifies the local Green/Lorentzian fingerprints of the trap
manifold; HYP-3227 adds full-bank Green-current conductance coordinates to
this normal-fan story.  The AP/doubled-AP equality face survives several
electrical metrics:

```text
all_ones_green_energy
cov_positive_lambda2
cov_positive_kirchhoff
precision_lambda2
precision_kirchhoff
precision_killing_abs
```

But `precision_mmatrix_defect` is a guardrail with `181` primitive beaters,
so the Green precision graph is a sidecar with named leakage, not a new
terminal scalar.

The trap result strengthens the two-chart theorem below.  The finite
trap/certificate graph remains connected without Toeplitz
(`lambda2=2.537866286`) and remains connected with only Green coordinates
(`lambda2=1.208613477`).  Thus the moment-cone trap chart is supported by an
independent Green-current discharge chart.  The next normal-fan rerun should
include AP-tight Green coordinates while excluding or sidecar-marking the
M-matrix defect.

## Next Tests

1. Search for a Toeplitz/Fejer-Riesz/Verblunsky dual whose slack equals or
   dominates the AP support gap, now with the Chebyshev/Cohn-Elkies
   magic-function interpretation from HYP-3212/HYP-3213.
2. Add HYP-3204's ordered-tail exchange-rate face to any future normal-cone
   certificate, especially `q0+q6` and the central `q3` pricing lemma.
3. Rerun the payload skyline after adding higher odd/Worpitzky or
   Hermite-Biehler interlacing shadows; check whether AP/dilation remains
   the only exposed face.
4. Prove a two-chart theorem: exchange-gradient improvement off the finite
   trap manifold, and Toeplitz moment-cone discharge on the traps.
5. Test whether the same payload cube survives beyond the bounded bank or
   after quotienting to primitive normal forms.
