---
id: HYP-3149
title: Tournament-4 canary/filler tables import the Erdos-870 order-two source interface
status: EVIDENCE / exact n=4 table scout and proof-interface transfer; not a proof
source: codex-2026-06-28
tangent: T1214
technique: LTI-275
tournament_technique: LTT-173
script: 04-computation/lrc14_erdos870_tournament4_canary_filler_codex_20260628.py
results:
  - 05-knowledge/results/lrc14_erdos870_tournament4_canary_filler_codex_20260628.out
related:
  - HYP-3143
  - HYP-3144
  - HYP-3145
  - HYP-3146
  - HYP-3147
  - HYP-3148
  - HYP-3142
  - HYP-3141
  - HYP-3140
  - HYP-3138
  - HYP-3137
  - HYP-3134
  - HYP-3133
  - HYP-3124
  - HYP-3118
  - HYP-3116
  - HYP-3093
  - HYP-3097
  - OPEN-Q-108
external_sources:
  - https://github.com/davidturturean/erdos-870
  - https://github.com/davidturturean/erdos-870/blob/main/problem_statement.md
  - https://github.com/davidturturean/erdos-870/blob/main/Erdos870/MainTheorem.lean
---

# HYP-3149: Tournament-4 Canary/Filler Tables

## Claim

The two prompt tables for tournaments on four vertices are not two unrelated
ways to name `S,T,+,-`.  They are a small exact model of the same proof
interface used in the formalized negative answer to Erdos Problem #870:

```text
order-two source + deterministic filler/canary + finite-shift deletion rule.
```

This refines the adjacent S276/HYP-3143 n=4 tournament subbasis packet.
HYP-3143 records the exact-order representation lesson: the four-fixed-arc
model is a minimal two-bit class basis, while the Hamiltonian-path cube has
lower-order leakage into `S`.  Upstream HYP-3146 records the broader
cover-versus-scaffold shift-package policy.  HYP-3149 names the narrower
canary/filler mechanism: in the fixed-Hamiltonian-path scaffold itself, the
extra coordinate `c` is what turns the two-bit table from a collision-prone
quotient into a proof-facing slice.

It also sits after S274/HYP-3144, HYP-3145, HYP-3146, upstream HYP-3147, and
upstream HYP-3148.  HYP-3144 asks when pair-function and three-edge quotients may be
scalarized without losing order-sensitive data; HYP-3145 names the broader
filler/core boundary; HYP-3146 packages the redundant fixed-path cover against
finite scaffolds; and HYP-3147 supplies the normalized n=3 edge-flip /
Worpitzky kernel.  HYP-3148 adds the live-core/deletable-coordinate audit;
here the order-sensitive datum is the exact fixed-path canary-slice collapse.

In the fixed-Hamiltonian-path tiling model, take the path

```text
0 -> 1 -> 2 -> 3
```

and the remaining arcs

```text
a=(0,2), b=(1,3), c=(0,3).
```

Starting from the transitive orientation, the quotient map from the flip cube
`F_2^3` to isomorphism classes has fibers

```text
T: {E}
+: {a}
-: {b}
S: {c, ab, ac, bc, abc}.
```

Thus the table

```text
* E a b c
E T + - S
a + T S S
b - S T S
c S S S T
```

is not a group law on tournament classes.  It is the pairwise visible part of
a quotient map with a large `S` collision fiber.

The second model fixes the extra arc `c=(0,3)` in its unflipped orientation,
so the fixed arcs have partial score sequence

```text
(0,1,1,2).
```

The two remaining coordinates `x=a`, `y=b` give the exact slice

```text
* E x y
E T + -
x + T S
y - S T
```

This `c=0` slice is an exact four-class transversal.  The complementary
`c=1` slice collapses every `x,y` completion to `S`.  That makes `c` a
canary/filler coordinate: proof-legal when fixed and audited, proof-illegal
when erased as raw symmetry.

## Evidence

The executable scout verifies the prompt tables, then checks all `64` labelled
four-vertex tournaments.  Labeled class counts and Hamiltonian-path counts are

```text
class labeled_count H_per_labeled_tournament fixed_path_conditioned_count
T            24                        1                            1
+             8                        3                            1
-             8                        3                            1
S            24                        5                            5
```

Conditioning on the fixed Hamiltonian path reproduces the three-bit tiling
fiber sizes `T=1,+ =1,-=1,S=5`.  Fixing `c` removes the path-conditioning
bulk and leaves the order-two `x,y` source.

The scout also records Tournament Analysis over proof carriers rather than
runners or raw arcs.  The proof-payload gauge is transitive:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1}
directed_3cycles=0
hamiltonian_path_count=1
edge_flips_vs_raw_symmetry_gauge=20
selected_path =
  fixed_c_xy_transversal
  -> erdos870_order2_plus_filler_interface
  -> tiling_hamiltonian_path_cube
  -> edge_tip_tail_information_packet
  -> S_bulk_collision_fiber
  -> fiber_pgf_or_distribution_sidecar
  -> raw_score_sequence_scalar
  -> raw_einheit_group_table_numerology
```

The twenty edge flips against the raw-symmetry gauge are the warning: the
prettiest group-like table ranks the wrong object first.

## Erdos-870 Transfer

The external repository `davidturturean/erdos-870` formalizes a negative
answer to Erdos Problem #870.  Its public theorem splits into `k=3` and
`k>=4`; the `k>=4` route consumes an order-two Larsen-Larsen source and adds
finite fillers, while the `k=3` route uses a clustered finite-shift canary
construction.  The useful import is not the additive-basis theorem itself,
but this interface discipline:

```text
low-order source
+ finite filler slots
+ canary exactness / finite-shift uniformity
+ deletion property
+ formal audit boundary.
```

HYP-3143 imports the exact-order / no-lower-order language, HYP-3144 adds the
pair-function scalarization alarm, and HYP-3145 names the broader filler-core
interface.  This note imports the deletion audit: if the filler coordinate is
not explicitly fixed or carried as a sidecar, the quotient is not legal proof
evidence.

The tournament-4 dictionary is:

```text
order-two source        = free coordinates x,y after c is fixed
deterministic filler    = fixed canary arc c plus the Hamiltonian path boundary
canary exactness        = c=1 makes every x/y completion land in S
finite-shift uniformity = the same x/y table must survive relabelled edge packets
deletion property       = if one core coordinate is removed, name restoration debt
```

## LRC14 Use

The next LRC14 packet should add

```text
tournament4_canary_filler_certificate =
  fixed_path_word
  + c_canary_status
  + xy_completion_table
  + S_bulk_fiber_words
  + deletion/restoration_sidecar
  + edge_tip_tail_exit_or_named_debt.
```

This connects directly to HYP-3141: a directed edge is not proof evidence
until its tail payload, tip payload, commutator defect, observer orbit, and
terminal exit are named.  The n=4 canary table is the smallest exact warning
that a raw edge/tile quotient can hide a whole `S` fiber unless one extra
coordinate is fixed as a filler.

The theorem target is a finite interface lemma:

```text
For every legal local edge packet, either the order-two completion table is an
exact transversal after fixing the canary/filler coordinate, or the packet emits
the S-bulk collision fiber as named restoration / observer-gluing debt.
```

This is a proof-interface improvement, not yet a covering proof.  It says how
to make a two-coordinate source legal before using it inside HYP-3141 edge
packets, HYP-3140 fiber-PGF rows, or HYP-3138/HYP-3139 k=8 reflection folds.

## Guardrail

Do not promote `E,a,b,c` or `E,x,y` to a literal group of isomorphism classes.
The table entries are values of a quotient map from a flip cube.  The relevant
invariant is not the multiplication table; it is the missing-coordinate fiber:

```text
visible quotient + canary status + fiber/restoration debt.
```
