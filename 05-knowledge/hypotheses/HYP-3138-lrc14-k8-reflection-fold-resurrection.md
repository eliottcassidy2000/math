---
id: HYP-3138
title: LRC14 k=8 reflection fold admits finite coordinate resurrection
status: EVIDENCE / executable bounded-bank quotient audit; not a proof
source: codex-2026-06-27-k8-reflection-fold
tangent: T1203
technique: LTI-264
tournament_technique: LTT-162
script: 04-computation/lrc14_k8_reflection_fold_resurrection_codex_20260627.py
result: 05-knowledge/results/lrc14_k8_reflection_fold_resurrection_codex_20260627.out
related:
  - HYP-3134
  - HYP-3133
  - HYP-3132
  - HYP-3122
  - HYP-3118
  - HYP-3116
  - HYP-3110
  - HYP-3085
  - THM-577
  - OPEN-Q-108
---

# HYP-3138: LRC14 k=8 Reflection-Fold Coordinate Resurrection

## Claim

The k=8 De Moivre/phi4 hard-row reduction should be treated as a legal
reflection fold plus coordinate-resurrection table, not as a claim that odd
coordinates vanish.

Fold the miss-count distribution `q_0,...,q_6` by the reflection `t <-> 6-t`:

```text
even_fold(q) = (q0+q6, q1+q5, q2+q4, q3).
```

The gK8 Delsarte functional

```text
L_yK8 = 10*q0 + q3 + 10*q6
```

depends only on the even fold.  However endpoint `Phi/P`, observer gluing,
finite-address exits, and packet legality may still need the destroyed odd
coordinates:

```text
odd_leakage = (q0-q6, q1-q5, q2-q4).
```

HYP-3138 tests whether the even fold has a finite adjoint/lookup on bounded
k=8 banks, in the sense of HYP-3118 coordinate resurrection.

## Evidence

The executable scout uses an exact integer breakpoint grid for primitive rows
`E={0}+7 speeds` in spans `14`, `15`, and `16`.

```text
span<=14: primitive_rows=3431, folded_signatures=3431, collision_fibers=0
span<=15: primitive_rows=6434, folded_signatures=6434, collision_fibers=0
span<=16: primitive_rows=11432, folded_signatures=11432, collision_fibers=0
```

Thus the even fold is injective on every tested bounded bank.  In this range,
the fold behaves like a finite coordinate-resurrection lookup, not a lossy
scalar quotient.

The best row remains the consecutive row:

```text
best_row=(0,1,2,3,4,5,6,7)
L_yK8=2633/735
10*cap8-L_yK8=683/2940
even_fold=(73/210,31/105,123/490,26/245)
odd_leakage=(451/1470,142/735,131/1470)
```

The nonzero `odd_leakage` is the important warning.  HYP-3132's biquadratic
symmetry gives the legal even fold:

```text
(t-1)(t-2)(t-4)(t-5), t=u+3
  = u^4 - 5u^2 + 4
  = (u^2-4)(u^2-1).
```

It does not erase odd proof coordinates.  It says which quotient is natural
for the Delsarte/phi4 part, while HYP-3116/HYP-3118 must supply the endpoint
activation and coordinate-repair sidecars before final quotienting.

## Candidate Invariant

```text
k8_reflection_fold_certificate =
  even_folded_miss_distribution
  + odd_coordinate_resurrection_table
  + de_moivre_biquadratic_resolvent_word
  + endpoint_phi_activation_status
  + observer_gluing_or_named_debt
```

## Niche-Work Transfer

- HYP-3110's Jacobi theta channels become the even/odd signed residue tails,
  not a standalone analytic shortcut.
- The `17` wallpaper groups and `230` space groups become finite quotient
  audits: a symmetry quotient is proof-legal only after stabilizers and
  destroyed coordinates are named.
- HYP-3116's circuit work supplies the exact endpoint `Phi/P` proof circuit.
- HYP-3118 supplies the rule that a quotient is legal only with a right
  adjoint/repair section for the coordinate demanded next.
- HYP-3134 supplies the A000568 global-consistency warning: do not quotient
  local data until the global gluing rule is named.

## Tournament Analysis

Tournament vertices are proof carriers and quotient operators, not runners or
raw roots.  The pairwise observable compares k=8 predicate retention,
destroyed-coordinate control, quotient legality, formal next step, exactness,
and niche bridge value.

Fingerprint:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1}
directed_3cycles=0
scc_sizes=[1,1,1,1,1,1,1,1,1]
hamiltonian_path_count=1
selected_path =
  endpoint_phi_activation_circuit
  -> coordinate_resurrection_adjoint
  -> k8_even_reflection_fold_table
  -> gK8_even_delsarte_functional
  -> A000568_global_consistency_quotient
  -> jacobi_theta_even_odd_channels
  -> de_moivre_biquadratic_resolvent
  -> raw_LyK8_scalar
  -> wallpaper17_space230_orbifold_audit
```

## Next Theorem Target

Prove a finite k=8 fold-adjoint lemma: on the bounded-core bank, the
`even_fold` determines the odd leakage or routes it to named finite-address /
observer-gluing debt.  Then the gK8/phi4 dip bound can use the even
biquadratic coordinate, while HYP-3116/HYP-3118 supply the legal repair for
endpoint `Phi/P` coordinates before final quotienting.
