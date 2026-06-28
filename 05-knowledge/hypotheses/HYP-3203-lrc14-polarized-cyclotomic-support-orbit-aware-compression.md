---
id: HYP-3203
title: Polarized cyclotomic support and orbit-aware compression are two new k=8 proof angles
status: EVIDENCE / exact bounded-bank scout; not an LRC14 proof
source: codex-2026-06-28
tangent: T1303
technique: LTI-303
tournament_technique: LTT-203
script: 04-computation/lrc14_new_angles_polarized_cyclotomic_compression_codex_20260628.py
result: 05-knowledge/results/lrc14_new_angles_polarized_cyclotomic_compression_codex_20260628.out
reflection: 07-reflections/polarized-cyclotomic-support-and-orbit-aware-compression-codex-20260628.md
related:
  - HYP-3162
  - HYP-3163
  - HYP-3200
  - HYP-3201
  - HYP-3202
  - HYP-3154
  - HYP-3153
  - HYP-3161
  - HYP-3160
  - HYP-3152
  - HYP-3150
  - HYP-3138
  - HYP-3132
  - HYP-3108
  - OPEN-Q-108
---

# HYP-3203: Polarized Cyclotomic Support And Orbit-Aware Compression

## Claim

Two new proof angles should be tried on the remaining LRC14 k=8 targets:

```text
Angle A: prove a polarized 7-cyclotomic support inequality
Angle B: replace naive left-compression by orbit-aware compression
```

The first attacks the Joukowski/Lee-Yang root-locus target from HYP-3154.  The
second attacks the ferromagnetic/covariance extremality target from
HYP-3160/HYP-3161.

## Exact Scout Readout

The bounded anchored k=8 bank has `3432` rows.  The AP row has miss-vector

```text
q_AP = (481/1470, 359/1470, 25/147, 26/245, 17/210, 5/98, 1/49).
```

Raw closeness to the uniform 7-fold/cyclotomic profile is the wrong target:
AP has rank `19` for minimum nontrivial cyclotomic energy.  Split-block rows
can be more uniform while having worse coverage and covariance.

The corrected support target is:

```text
<q(E)-1/7, q_AP-1/7> <= ||q_AP-1/7||^2.
```

On the bounded bank this is maximized exactly by `{0,1,2,3,4,5,6,7}` and its
doubled dilation `{0,2,4,6,8,10,12,14}`, with value

```text
39766/540225.
```

The same two rows maximize coverage `q0`, `L_y=q0+q6+q3/10`, and total
covariance `Sigma-kappa_2`.

## Angle A: Polarized Cyclotomic Support

HYP-3154 says the on-circle ideal is the seventh cyclotomic profile.  The
scout shows that "closest to uniform/cyclotomic" is not the LRC extremality.
The AP is not the minimum-norm residual; it is the exposed point in the AP
residual direction.

This turns the root-locus problem into a linear inequality over the whole
miss-count vector.  If this support inequality can be proven globally or on
the primitive bounded-core bank, then the Joukowski sidecar gets a theorem
facing shape:

```text
circle-root data -> AP residual support halfspace -> coverage/covariance cap.
```

This is a sharper instruction than "minimize off-circle error."  The
off-circle defect is directional.

## Angle B: Orbit-Aware Compression

A naive local proof route is false.  Replacing a larger runner by a smaller
missing runner does not always increase total covariance.  In the bounded k=8
bank:

```text
local no-improvement traps = 19
greedy left-compression stuck states = 919
```

The doubled AP dilation is one of the local traps, so any compression proof
that does not first quotient dilation/orbit information is aiming at the wrong
predicate.  Other traps are two-block or tail-cluster rows.  The repair is to
use a compression lemma with sidecars:

```text
dilation_orbit_sidecar
mirror_orbit_sidecar
two_block_trap_sidecar
AP_residual_projection_Lyapunov
```

That makes compression a legal controlled-forgetting move instead of a raw
coordinate-left-shift.

## Tournament Analysis

Vertices are proof moves, not runners:

```text
AP_polarized_cyclotomic_support
ferromagnetic_covariance
orbit_aware_compression
dilation_orbit_sidecar
two_block_trap_sidecar
raw_cyclotomic_norm
raw_left_compression
raw_scalar_p0
```

Pairwise observable: which carrier preserves the LRC14 coverage/dip predicate
while keeping the needed root-locus and compression sidecars.

Switch/gauge: orient toward the carrier with stronger predicate preservation
and less destroyed-coordinate debt.

Fingerprint from the scout:

```text
score_hist = {0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1}
directed_3cycles = 0
hamiltonian_path_count = 1
selected path:
AP_polarized_cyclotomic_support
-> ferromagnetic_covariance
-> orbit_aware_compression
-> dilation_orbit_sidecar
-> two_block_trap_sidecar
-> raw_cyclotomic_norm
-> raw_left_compression
-> raw_scalar_p0
```

## Assumption Challenge

The challenged assumption is that "more cyclotomic/uniform" or "more
left-compressed" should be the proof invariant.  Both fail on the exact bank.
The preserved predicate is not raw root distance or coordinate sum; it is
AP-directed coverage/covariance support.  The destroyed information is orbit
type: dilation, mirror, and two-block trap status.

## Next Tests

1. Try to prove the support inequality as a finite signed-SPEC or
   moment-cone inequality.
2. Classify the `19` local compression traps into dilation, mirror, and
   two-block sidecar families.
3. Test whether AP residual projection remains sharp on larger bounded banks
   and primitive representatives after quotienting dilation.

## Post-HYP-3224 Integration

HYP-3224 confirms that the AP support functional is not isolated: when paired
with Toeplitz `lambda_min`, covariance layers `D1,D2,D3`, and total
covariance, the bounded-bank Pareto skyline is still exactly AP plus doubled
AP.  It also shows that the `11` strongest primitive exchange traps are all
discharged by strict Toeplitz moment-cone deficits.  Therefore the next
support proof should look for a normal-cone dual certificate rather than only
a linear inequality in the q-vector.
