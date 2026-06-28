# HYP-3202: k=8 Covariance Attack Angles

**Status:** EVIDENCE / exact bounded-bank scout; not an LRC14 proof.

## Claim Tested

HYP-3200 leaves the clean theorem-facing target:

```text
Primitive k=8 bounded-core rows maximize Sigma kappa_2 at
E=(0,1,2,3,4,5,6,7).
```

HYP-3202 asks how to attack that inequality without falling back to entropy,
the false exact `1/7` associator law, or plain positive association.  It tests
two routes on the same anchored bank

```text
E = {0} union A,  A subset {1,...,14}, |A|=7
```

with `3432` rows and `3431` primitive rows.

## Incoming Mainline Integration

Push-time fetch revealed that mainline had already claimed `HYP-3201` for a
law-defect entropy packet, and also added k=8 Perron, Toeplitz/Szego, and
covariance-Laplacian proof angles.  This packet is therefore renumbered to
HYP-3202 and read as a complementary finite localization:

- HYP-3163 asks for covariance-kernel, Monge, conditionally-positive, or PSD
  certificates.  The distance layers below are concrete candidate coordinates
  for its `covariance_kernel_distance_profile` and `monge_four_point_defect`
  fields.
- The S73c Perron route explains why the all-ones direction is the right
  ferromagnetic covariance mode.  The layer split tests whether that Perron
  alignment can be sharpened into three distance-wise inequalities.
- The S31al Caratheodory-Toeplitz route gives a global moment/spectral
  certificate.  The exchange-trap manifold below is a finite obstruction list
  for any local move or random-current coupling proof.
- HYP-3201 law-defect entropy remains the quotient-legality meter: any
  compression of the layer profile, exchange graph, or trap identities must
  have zero target-function residual or retain a named sidecar.

## Angle B: Cyclic-Distance Covariance Layers

Split the `15` pair covariances between the six inner sectors by cyclic
distance `d=1,2,3` on the heptagon.  For consecutive speeds:

```text
Sigma kappa_2 = 6237419/8643600 = 0.721622819196
distance 1    = 308509/1080450
distance 2    = 547577/2160900
distance 3    = 225577/1234800
```

Exact result: in the all-bank, each distance layer has `0` beaters and `2`
ties, where the second tie is the nonprimitive even-AP dilation
`(0,2,4,6,8,10,12,14)`.  In primitive normal form, each distance layer has
`0` beaters and `1` tie: consecutive speeds alone.

This suggests a sharper proof program than one total-covariance scalar:

```text
prove D1(E) <= D1(consec)
prove D2(E) <= D2(consec)
prove D3(E) <= D3(consec)
sum the three inequalities to get Sigma kappa_2 extremality.
```

The layer split matches the earlier reflection/Perron route, but it adds a
new local target: each distance class should have its own majorization or
transport proof, with dilation identified as the only all-bank equality
artifact.

## Angle A: Exchange Gradient and Critical Traps

Treat `Sigma kappa_2` as an energy on the primitive bounded bank and look for
improving moves toward consecutive speeds.

```text
adjacent +/-1 moves:       3066/3430 non-target rows improve; 364 stuck
gap-fill moves:            3411/3430 non-target rows improve; 19 stuck
arbitrary one-point swap:  3419/3430 non-target rows improve; 11 stuck
```

The arbitrary-swap local maxima are exactly `12`: consecutive plus the `11`
trap rows

```text
(0,2,4,6,7,8,10,12)
(0,1,3,5,7,9,11,13)
(0,1,2,3,11,12,13,14)
(0,3,6,8,9,11,12,14)
(0,2,3,5,6,8,11,14)
(0,6,7,8,9,10,11,12)
(0,4,5,8,9,10,13,14)
(0,8,9,10,11,12,13,14)
(0,1,2,3,7,8,9,10)
(0,2,5,7,9,10,12,14)
(0,1,4,5,7,8,11,12)
```

The best-improving exchange climb has `12` endpoints and maximum depth `7`.
The largest basin is consecutive with size `1853`; the largest decoy basin is
`(0,2,4,6,7,8,10,12)` with size `475`.

The proof shape is therefore not "find one monotone exchange that always
works."  The more realistic target is:

```text
bulk exchange lemma:
  outside the finite critical-trap manifold, there is an improving exchange

trap discharge lemma:
  every trap is below consecutive, or routes to one distance-layer deficit
```

## Ferromagnetic Guardrail

Consecutive speeds have all `15` pair covariances nonnegative and
`min_cov=1873/180075`.  But positive association is not enough: `19`
primitive rows have all `15` pair covariances nonnegative, including
structured decoys such as `(0,1,3,5,7,9,11,13)` and
`(0,6,7,8,9,10,11,12)`.

Thus the useful ferromagnetic statement is not raw FKG.  It needs either
cyclic-distance layer dominance or a trap-discharge mechanism that sees where
the positive covariances are located.

## Tournament Analysis

Vertices are proof moves and retained signals, not runners, sectors, or arcs.
The pairwise observable is which signal preserves a usable route to
`Sigma kappa_2` extremality.  The scout's transitive priority path is:

```text
cyclic_distance_covariance_layers
-> exchange_gradient_bulk
-> finite_critical_trap_manifold
-> ferromagnetic_positive_pair_sidecar
-> raw_total_covariance_scalar
-> plain_FKG_monotonicity
-> entropy_min_description
-> one_seventh_associator_scalar
```

This challenges the default mapping explicitly: for this proof obligation,
the useful tournament is over proof channels and route quality.  Runner- or
arc-level tournaments may still be sidecars, but they are not the primary
vertices for the covariance inequality.

## Proof-Frontier Use

HYP-3202 gives two concrete next tasks:

1. Prove the three cyclic-distance covariance inequalities separately.
2. Prove an exchange-gradient lemma except for the `11` named traps, then
   discharge those traps using distance-layer deficits, reflection folds, or
   finite-resolvent sidecars.

The odd Worpitzky/associator channel remains separate.  Do not scalarize it
to `1/7`; keep `S3`, minority-edge gates, PGF/root data, and terminal debt as
sidecars.

## Post-HYP-3224 Integration

HYP-3224 reruns the `11` arbitrary-exchange traps against a larger spectral
payload cube.  Every nonconsecutive trap has a strict Toeplitz `lambda_min`
deficit, and the combined skyline of

```text
AP_support, Toeplitz_lambda_min, D1, D2, D3, Sigma_kappa2
```

is exactly consecutive plus the nonprimitive doubled AP dilation.  This
sharpens the route: use exchange/covariance as the bulk chart, then use the
Caratheodory-Toeplitz moment-cone chart to discharge the trap boundary.

## Artifacts

- Script: `04-computation/lrc_k8_covariance_attack_angles_codex_20260628.py`
- Result: `05-knowledge/results/lrc_k8_covariance_attack_angles_codex_20260628.out`

## Links

HYP-3202 extends HYP-3200 and should be read with HYP-3161, HYP-3160,
HYP-3163, HYP-3201, HYP-3154, HYP-3153, HYP-3152, HYP-3151, HYP-3150,
HYP-3147, HYP-3144, HYP-3142, HYP-3139, HYP-3138, HYP-3132, HYP-3122,
THM-577, LTI-302, LTT-202, T1302, and OPEN-Q-108.
