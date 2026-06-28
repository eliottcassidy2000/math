---
id: HYP-3205
title: Spectral dictionary compatibility for the k=8 LRC14 frontier
status: EVIDENCE / exact bounded-bank scout plus synthesis; not an LRC14 proof
source: codex-2026-06-28
tangent: T1305
technique: LTI-305
tournament_technique: LTT-205
script: 04-computation/lrc_k8_spectral_dictionary_compatibility_codex_20260628.py
result: 05-knowledge/results/lrc_k8_spectral_dictionary_compatibility_codex_20260628.out
reflection: 07-reflections/spectral-dictionary-compatibility-codex-20260628.md
related:
  - HYP-3204
  - HYP-3203
  - HYP-3202
  - HYP-3210
  - HYP-3201
  - HYP-3200
  - HYP-3163
  - HYP-3162
  - HYP-3161
  - HYP-3160
  - HYP-3154
  - HYP-3153
  - HYP-3152
  - HYP-3150
  - HYP-3138
  - HYP-3132
  - HYP-3117
  - HYP-3116
  - HYP-3115
  - OPEN-Q-108
---

# HYP-3205: Spectral Dictionary Compatibility

## Claim

The next theorem-facing object should be a compatibility vector, not a single
scalar.  On the k=8 bounded bank, AP/consecutive and its doubled dilation are
the only rows simultaneously coherent across the current certificate
languages:

```text
coverage q0
L_y = q0 + q6 + q3/10
total covariance Sigma kappa_2
three cyclic-distance covariance layers D1,D2,D3
AP-polarized cyclotomic support
Caratheodory-Toeplitz lambda_min margin
```

The proposed proof shape is therefore a certificate-Helly theorem:

```text
Every primitive non-AP row violates at least one certificate coordinate by a
usable margin; the finite exchange/compression traps are discharged by the
first certificate coordinate where they visibly fail.
```

This merges HYP-3202's covariance layer/trap program, HYP-3203's polarized
cyclotomic support, the S31al/S73d Toeplitz/moment route, and HYP-3210's
Joukowski/Hermite-Biehler/Perron bridge.

Post-rebase integration: incoming HYP-3204 ordered-tail exchange is not a
collision after renumbering; it is a coefficient chart inside this dictionary.
Its central exchange-rate lemma prices `q3` gains against `q0+q6` loss, while
HYP-3205 says that any such coefficient chart should remain compatible with
the layer/support/Toeplitz certificate vector.  In this reading, HYP-3204 is a
strong candidate discharge route for the `L_y` coordinate, and HYP-3205 is the
cross-chart compatibility audit that keeps it from becoming a raw `q3` or
entropy shortcut.

## Exact Scout Readout

The executable scout scans the anchored bounded bank

```text
E = {0} union A,  A subset {1,...,14}, |A|=7
rows_all = 3432
rows_primitive = 3431
```

It recomputes the AP miss vector

```text
q_AP =
(481/1470, 359/1470, 25/147, 26/245, 17/210, 5/98, 1/49).
```

Against AP, the channel ranks are:

```text
coverage_q0:             0 beaters, 2 all-bank ties, 1 primitive tie
L_y:                     0 beaters, 2 all-bank ties, 1 primitive tie
Sigma_kappa2:            0 beaters, 2 all-bank ties, 1 primitive tie
distance_1:              0 beaters, 2 all-bank ties, 1 primitive tie
distance_2:              0 beaters, 2 all-bank ties, 1 primitive tie
distance_3:              0 beaters, 2 all-bank ties, 1 primitive tie
AP_residual_projection:  0 beaters, 2 all-bank ties, 1 primitive tie
Toeplitz_lambda_min:     0 beaters, 2 all-bank ties, 1 primitive tie
```

The two all-bank ties are exactly

```text
(0,1,2,3,4,5,6,7)
(0,2,4,6,8,10,12,14).
```

The simultaneous exact maximizers across coverage, `L_y`, total covariance,
all three distance layers, and AP support are the same two rows; adding the
Toeplitz margin leaves the same two rows.  In primitive normal form, AP is the
unique simultaneous maximizer across the dictionary.

Two warnings remain important:

```text
Perron_uniform_alignment: AP has rank 2, so this is diagnostic, not terminal.
raw_cyclotomic_energy_MIN: AP has rank 19, so raw uniformity is false.
```

This agrees with the incoming spectral synthesis: Perron alignment explains
the ferromagnetic mode, but the theorem-facing coordinate is still a
payload-bearing certificate such as a layer inequality, support hyperplane, or
Toeplitz slack.

## Nearest Decoys

The scout ranks non-AP rows by smallest maximum normalized deficit across

```text
q0, L_y, Sigma kappa_2, D1, D2, D3, AP support, Toeplitz margin.
```

The closest primitive decoy is

```text
E=(0,2,3,4,5,6,7,8)
max_norm=0.206062
def_Sigma_kappa2=0.1434055544
def_AP_support=0.0057647277
def_Toeplitz=0.0844938805
layer_defs=(0.0817898503,0.0354682945,0.0261474096)
```

No non-AP decoy has a zero deficit in any of the eight dictionary coordinates.
This is the key new evidence: the strongest-looking decoys are close in
different ways, but none is invisible to the whole dictionary.

## Trap Sheaf

HYP-3202 left `11` primitive arbitrary-exchange traps plus AP.  HYP-3205
adds a discharge table for those traps.  Each trap has a visible deficit in
AP support and Toeplitz margin, and at least one visible distance-layer
deficit.  For example:

```text
(0,2,4,6,7,8,10,12):
  def_AP_support = 0.0032783794
  def_Toeplitz   = 0.1458043926
  layer_defs     = (0.1364451155, 0.0772303207, 0.0491342438)

(0,1,3,5,7,9,11,13):
  def_AP_support = 0.0118544235
  def_Toeplitz   = 0.1139353483
  layer_defs     = (0.1499399102, 0.0909587958, 0.0671140852)
```

This suggests replacing "one monotone exchange path" by a sheaf-style proof:
local charts are exchange or compression neighborhoods, and the gluing data
is the certificate vector.  A row can be locally maximal in the exchange
chart, but then the Toeplitz chart, AP-support chart, or a distance-layer
chart exposes the obstruction.

Primary-discharge histograms from the normalized certificate deficits:

```text
exchange_trap_primary_discharge =
  {'ap_support': 5, 'd1': 4, 'd2': 2}

left_compression_non_AP_primary_discharge =
  {'Ly': 1, 'ap_support': 8, 'd1': 4, 'd2': 2,
   'q0': 2, 'sk2': 1, 'toeplitz': 1}
```

So the exchange-trap manifold is mostly a support/layer problem, while
left-compression traps need a wider dictionary.  That is a concrete reason to
separate exchange-gradient proof work from orbit-aware compression proof work.

## Two Maps

### Map A: Certificate Dictionary

```text
row E
  -> miss vector q
  -> coverage/L_y
  -> covariance matrix C and distance layers D1,D2,D3
  -> AP residual support <q-1/7,q_AP-1/7>
  -> Toeplitz moment matrix [q_|i-j|]
  -> Perron uniform-mode diagnostic
  -> trap flags and orbit sidecars
```

The proof target is to show that only AP/dilation can be maximal in every
coordinate.  This is a finite Helly/separating-hyperplane program: every
non-AP row should be separated by one member of a small certificate family.

### Map B: Adjoint Transport

```text
Toeplitz moment cone
  <-> Caratheodory/Fejer/Szego/Verblunsky slack
  <-> Lee-Yang circle and Joukowski image
  <-> Hermite-Biehler even/odd interlacing
  <-> covariance/Perron ferromagnetic mode
  <-> tournament proof-channel order
```

This is where HYP-3210 matters.  The Joukowski bridge should not be used as a
raw root-distance slogan.  It is a transport map: Toeplitz slack and AP
support are the circle/moment side; Hermite-Biehler interlacing is the
real-axis/tournament side; covariance layers are the finite bounded-bank
coordinates where the transport can be tested.

## Creative Extensions Worth Testing

1. **Certificate-Helly lemma.**  Prove that the intersection of the AP-tight
   faces of the layer, support, and Toeplitz certificates is exactly the
   AP/dilation orbit.  Then primitive uniqueness follows by removing dilation.

2. **Trap-sheaf discharge.**  Treat the `11` exchange traps and `19` non-AP
   left-compression traps as stalks.  The local chart may not improve
   covariance, but a transition function to Toeplitz/AP-support/layer
   coordinates should discharge it.

3. **Random-current compatibility.**  The Griffiths/FKG analogy is too weak
   on the speed lattice, but a random-current representation might define the
   correct coupling coordinates.  HYP-3205 says the current must preserve the
   same certificate vector, not only positive pair covariances.

4. **Circuit lower-bound guard.**  In the HYP-3116/HYP-3117 language, a
   scalar proof is a too-small circuit.  The input basis should be
   `certificate_vector=(D1,D2,D3,AP_support,Toeplitz_slack,trap_orbit,odd_debt)`.
   A shortcut is legal only if its circuit reconstructs those inputs or names
   the destroyed sidecar.

5. **Minkowski/Bravais relation wall.**  Read AP as the central lattice point
   in a low-dimensional relation body for the q-vector.  A Minkowski-style
   proof would show every non-AP relation vector exits at least one certificate
   slab.  The current scout supplies the slab coordinates.

6. **Irrational approximation sidecar.**  The 7th-cyclotomic ideal lives in
   the cubic real subfield generated by `2cos(2pi/7)`, while all bounded-bank
   q-data are rational.  The AP support functional may be the rational
   separating hyperplane for the cubic Joukowski ideal.

7. **Descent `14 -> 7 -> 2`.**  The composite `k+1=14` obstruction mirrors
   the paper's composite fallback: lifts at `c=7` and `c=2`.  HYP-3205 says
   the descent should carry a certificate vector through the lift, not just a
   witness time.  This is also the honest relation to prime-gap/sexy-prime
   analogies: residue-pair survival is useful as a sidecar language, but it
   does not by itself prove an LRC extremal inequality.

## Tournament Analysis

Vertices are proof channels/certificates:

```text
cyclic_distance_layer_bundle
AP_residual_support_hyperplane
Toeplitz_moment_margin
total_covariance_scalar
exchange_trap_sheaf
orbit_aware_left_compression
Perron_uniform_mode_diagnostic
entropy_or_min_description_scalar
one_seventh_associator_scalar
raw_cyclotomic_norm
```

Pairwise observable: retained certificate payload for k=8
coverage/covariance extremality.

Switch/gauge: orient toward the channel with fewer AP beaters, fewer
nontrivial ties, better obstruction localization, and better sidecar
retention.

Fingerprint:

```text
score_hist={-100:1,-19:1,-11:1,-5:1,19:1,37:1,69:1,75:1,76:1,79:1}
directed_3cycles=0
sccs=singletons
hamiltonian_path_count=1
priority_path:
cyclic_distance_layer_bundle
-> AP_residual_support_hyperplane
-> Toeplitz_moment_margin
-> total_covariance_scalar
-> exchange_trap_sheaf
-> orbit_aware_left_compression
-> Perron_uniform_mode_diagnostic
-> entropy_or_min_description_scalar
-> one_seventh_associator_scalar
-> raw_cyclotomic_norm
```

## Assumption Challenge

Alternate vertex sets considered:

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
certificate channels
```

Chosen vertices: certificate channels, with bounded-bank rows as test
objects.  This preserves the LRC coverage/covariance extremality predicate,
the dilation exception, trap identities, and root/moment/covariance sidecar
coordinates.  It destroys literal runner identity and raw arc orientation on
purpose.

The challenged assumption is that the obstruction is a single scalar maximum.
The better model is a compatibility failure across adjoint certificate
languages.

## Next Tasks

1. Prove the certificate-Helly statement on the exact bounded bank, then try
   to lift it to the LRC14 proof state.
2. Classify the `11` exchange traps and `19` non-AP left-compression traps by
   their first failing dictionary coordinate.
3. Try to derive the AP support inequality from Toeplitz moment-cone duality
   or signed SPEC, and then transport it through the Joukowski/Hermite-Biehler
   bridge.
4. Keep the odd Worpitzky/associator channel as named debt; do not fold it
   into `1/7`, entropy, or raw cyclotomic norm.
