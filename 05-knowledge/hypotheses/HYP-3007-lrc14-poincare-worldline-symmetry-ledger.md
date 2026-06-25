---
id: HYP-3007
title: LRC14 Poincare worldline symmetry ledger
status: SYNTHESIS / symmetry guardrail and proof-interface carrier; not a proof
source: codex-2026-06-25-S172
script: 04-computation/lrc14_poincare_symmetry_ledger_codex_s172.py
result: 05-knowledge/results/lrc14_poincare_symmetry_ledger_codex_s172.out
related:
  - HYP-3006
  - HYP-3005
  - HYP-3004
  - HYP-3002
  - HYP-2997
  - HYP-2996
  - HYP-2995
  - HYP-2990
  - HYP-2963
  - HYP-2953
  - HYP-2924
  - HYP-2913
  - HYP-2486
  - HYP-2291
  - THM-381
  - THM-385
  - THM-572
  - OPEN-Q-108
---

# HYP-3007: LRC14 Poincare Worldline Symmetry Ledger

## Claim

The useful part of the Poincare-group analogy is not a new scalar invariant.  It
is a symmetry ledger for the time/phase cylinder:

```text
runner worldline:      x_i(t) = v_i t mod 1
observer worldline:    x_0(t) = u t mod 1
danger tube:           ||x_i(t)-x_0(t)|| < 1/14
safe time:             one horizontal time slice misses all danger tubes
```

The standard LRC fixes `u=0`, so it is observer-anchored.  A Poincare or
Galilean "boost" is legal only in the enlarged observer-coupled worldline
groupoid where `u` is retained and all verdicts are recentered to relative
speeds `v_i-u`.  If the observer velocity is forgotten, the same boost becomes
ordinary speed translation `v_i -> v_i+c`, and it is not an LRC automorphism.

So the proof lesson is:

```text
true standard-LRC symmetries:
  runner permutation,
  independent sign flips v_i -> +/- v_i,
  global reflection / time reversal,
  integer time dilation plus primitive gcd normalization.

observer-coupled symmetries:
  boosts of all worldlines, legal only while observer velocity is retained.

useful non-symmetries:
  stationary speed translation and Lorentz-like velocity addition,
  because their failures identify lost observer, metric, or lattice labels.
```

## Computation

Script:

```text
04-computation/lrc14_poincare_symmetry_ledger_codex_s172.py
```

Output:

```text
05-knowledge/results/lrc14_poincare_symmetry_ledger_codex_s172.out
```

The script computes exact safe measure for named LRC14 rows under six
transformations.  The exact stress test is deliberately simple:

```text
AP13 base_safe=0
  sign flips:                 0, unchanged
  scale by 5:                 0, unchanged
  stationary +5 translation:  41299/299880, changed
  observer-coupled +5 boost:  0, unchanged

GW13 base_safe=0
  stationary +5 translation:  218809/1662570, changed

K33_12_to_36 base_safe=1/1260
  stationary +5 translation:  534439/3760848, changed

q27_two_block_probe base_safe=2556/245245
  stationary +5 translation:  1935391862213/16530724395900, changed
```

The key point is not the particular values.  It is the invariant/failure split:
signs, reflection, and scale preserve exact safe measure; speed translation
does not, unless it is recorded as an observer-coupled boost and recentered.

## Tournament Analysis

Vertices are symmetry/proof carriers, not runners:

```text
observer_coupled_worldline_tube_groupoid
individual_sign_flip_parity_kernel
integer_time_dilation_primitive_scale
anchored_metric_winding_tournament
stationary_velocity_translation
bare_winding_iso_class
lorentz_velocity_addition_shadow
raw_speed_scalar
```

Pairwise observable:

```text
predicate retention
observer anchor retention
integer lattice retention
exact scale retention
packet label retention
boost compatibility
anti-scalar guardrail
proof maturity
```

The conservative gauge is transitive:

```text
score_hist=[(0,1),(1,1),(2,1),(3,1),(4,1),(5,1),(6,1),(7,1)]
directed_3cycles=0
Hamiltonian path:
  individual_sign_flip_parity_kernel
  > observer_coupled_worldline_tube_groupoid
  > integer_time_dilation_primitive_scale
  > anchored_metric_winding_tournament
  > stationary_velocity_translation
  > bare_winding_iso_class
  > lorentz_velocity_addition_shadow
  > raw_speed_scalar
```

The high rank of `individual_sign_flip_parity_kernel` is intentional.  It is an
exact kernel of the anchored LRC predicate but not of winding order.  Therefore
any proof route that uses pairwise order, signed tournaments, or directed
worldline crossings must either quotient by this parity kernel or explicitly
retain the sign/orientation debt.

## Relation To The Existing Stack

This sharpens the three observer categories:

```text
category 1: observer-relative tiling / LRC coverage
category 2: metric winding order plus gap widths
category 3: observer-blind tournament or affine scalar
```

The Poincare ledger says to add a fourth line to the checklist:

```text
worldline-frame data: observer velocity, tube metric, and lattice embedding.
```

HYP-3002's curried packet evaluator becomes:

```text
E(S)(packet)(observer_frame)(time_phase_tube)(lane)(certificate)(verdict).
```

HYP-2997/HYP-3006's cocycle language becomes:

```text
forgotten observer velocity = boost cocycle,
forgotten tube metric       = cone/tube deformation cocycle,
forgotten sign orientation  = parity-kernel cocycle for winding carriers.
```

This also explains why negative speeds are a genuine trick but not a complete
proof.  Individual direction reversal is an exact LRC symmetry, so it can
simplify observer-distance predicates.  But it destroys pairwise orientation
data, so tournament or winding routes must record what was lost.

## Proof Targets

1. **Sign-kernel lemma.**  The standard anchored LRC predicate factors through
   the absolute relative speed multiset.  Any winding/tournament route that
   distinguishes sign choices must discharge the sign debt by parity-kernel
   exactness or reattach orientation labels.

2. **Boost admissibility lemma.**  A common velocity translation is a legal
   quotient only for observer-coupled packets carrying observer velocity `u`;
   after recentering by `v_i-u`, the standard predicate is recovered.  Without
   `u`, stationary speed translation changes safe measure and is forbidden as
   a proof quotient.

3. **Tube-cone deformation lemma.**  Lorentz/Poincare-style velocity addition is
   not an automorphism of the integer-speed LRC lattice with fixed Euclidean
   circle metric.  It can be used only if the deformed tube metric and lattice
   embedding are carried as named cochains or routed to F7/THM-572 residual
   debt.

4. **Packet-schema addition.**  Add these fields to HYP-2963-style records:

```text
observer_velocity_label
relative_speed_normal_form
sign_kernel_status
primitive_scale_gcd
tube_metric_label
worldline_frame_label
boost_cocycle_status
orientation_debt_for_winding
```

The intended LRC14 use is narrow but important: it prevents a proof from using
a Poincare/boost intuition while silently jumping from the observer-coupled
worldline groupoid back to the stationary-observer problem.
