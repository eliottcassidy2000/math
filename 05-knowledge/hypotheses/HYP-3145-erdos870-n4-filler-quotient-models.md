---
id: HYP-3145
title: Erdos-870 filler interfaces and n=4 tournament quotient models
status: EVIDENCE / executable finite interface scout; not a proof
source: codex-2026-06-27; prompted by davidturturean/erdos-870 and the two n=4 tables
tangent: T1210
technique: LTI-271
tournament_technique: LTT-169
script: 04-computation/lrc14_erdos870_n4_filler_models_codex_20260627.py
result: 05-knowledge/results/lrc14_erdos870_n4_filler_models_codex_20260627.out
external:
  - https://github.com/davidturturean/erdos-870
related:
  - HYP-3144
  - HYP-3147
  - HYP-3143
  - HYP-3136
  - HYP-3135
  - HYP-3134
  - HYP-3133
  - HYP-3132
  - HYP-3129
  - HYP-3125
  - HYP-3124
  - HYP-3118
  - HYP-3107
  - HYP-3054
  - HYP-3053
  - HYP-3049
  - OPEN-Q-108
---

# HYP-3145: Erdos-870 Filler Interfaces And n=4 Quotient Models

## Claim

The useful import from the Erdos-870 formalization is architectural:
separate a small retained core from deterministic finite fillers, prove the
interface lemmas at that boundary, and only then transfer to the terminal
theorem.  The n=4 tournament tables give a miniature version of the same
choice.

Incoming HYP-3143 already captures the exact-order packet-subbasis audit:
the Hamiltonian-path cube has lower-order leakage, while a two-bit matching
core has one squarefree representation per class.  This card extends that
audit toward the LRC14 floor closure by naming the filler/core boundary fields
that have to survive before a factorized estimate can be used as proof data.
HYP-3144/HYP-3147 supply the smaller K3 preflight: a quotient may preserve the
class kernel while losing the minority edge, Worpitzky descent word, and
ordered-function channel.  The n=4 filler/core interface below should therefore
be read as the first place where that local sidecar rule becomes a scaffold
choice.

In the fixed-Hamiltonian-path model, three off-path arcs `a,b,c` remain.  With
the transitive state `E` and class names

```text
T=(0,1,2,3), S=(1,1,2,2), +=(0,2,2,2), -=(1,1,1,3),
```

the representative table is

```text
* | E  a  b  c
--+------------
E | T  +  -  S
a | +  T  S  S
b | -  S  T  S
c | S  S  S  T
```

This is a good tiling atlas but not a class algebra.  The full fixed-path
fiber has

```text
T: E
+: a
-: b
S: c, ab, ac, bc, abc
```

so multiplying from an `S` representative is ambiguous.

In the second model, four arcs are fixed so their partial score profile is
`(0,1,1,2)`.  One exact realization fixes

```text
(0,2) as 2->0
(0,3) as 0->3
(1,2) as 1->2
(1,3) as 1->3
```

and leaves `x=(0,1)` and `y=(2,3)` as the two free arcs.  The user-facing
table is

```text
* | E  x  y
--+---------
E | T  +  -
x | +  T  S
y | -  S  T
```

and the closed table with `xy=S` is a Klein-square interface.  Unlike the
tiling table, this two-free-arc presentation is congruent to the four
isomorphism classes.

## Erdos-870 Transfer

The referenced Erdos-870 repository formalizes a negative answer to the
minimal-subbasis question by packaging a sparse order-two source and
deterministic endings.  For `k>=4` it uses a finite residue-class filler
gadget; for `k=3` it uses a parity reduction and clustered construction.  The
local lesson for LRC14 is not the additive-basis theorem itself, but the
interface discipline:

- pick fillers that force the residue/score scaffold;
- retain the small core where representations or signs actually vary;
- prove deletion/nonminimality or quotient failure at the interface;
- expose the proof boundary as a formal target before using it downstream.

Before scalarizing a filler/core table, run the HYP-3144 pair-function test:
sum/product-like observables may factor through the quotient, while
exponent-like observables need an ordered or deletion sidecar.  HYP-3147's
minority-edge gate is the local witness for exactly that missing coordinate.

In LRC14 language, HYP-3136 should not be pushed as a scalar
`Rprime*meas(R-safe)*meas(Q-lonely)` product without such a boundary.  The
finite filler side should retain normalized residue, global-consistency,
endpoint-owner, and tail/tip child data.  The small core should retain the
signed SPEC low modes and the HYP-3132 k=8 De Moivre/phi4 resolvent packet.
If a fixed-path tiling coordinate has a many-to-one `S` fiber, it is a
nonminimal/deletable-data alarm rather than a proof algebra.

## Evidence

The executable scout verifies the two user tables and finds the explicit
partial-score realization above.  Its tournament-analysis vertices are proof
carriers rather than runners:

```text
erdos870_filler_interface
-> quotient_congruence_audit
-> partial_score_two_arc_core
-> edge_witness_SPEC_packet
-> nonminimal_deletable_fiber_alarm
-> fixed_path_tiling_cube
-> raw_n4_class_table
```

Fingerprint:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1}
directed_3cycles=0
```

## Next Use

Add a `filler_core_interface` field to HYP-3125/HYP-3129 edge-floor rows:

```text
finite_filler_scaffold
partial_score_or_residue_profile
core_variable_pair
quotient_congruence_status
nonminimal_fiber_alarm
deletion_or_forgetting_exit
formal_interface_target
```

The next concrete test is to rewrite one multi-far covering row in this form:
deterministic Q/apex and residue fillers first, then a small signed core whose
two-variable table is SPEC/De Moivre-facing rather than a raw fixed-path
tiling cube.
