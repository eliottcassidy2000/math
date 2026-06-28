---
id: HYP-3143
title: n=4 tournament class quotients should be audited as packet subbases
status: EVIDENCE / exact n=4 enumeration and proof-packet proposal; not an LRC proof
source: codex-2026-06-28-S276
tangent: T1208
technique: LTI-269
tournament_technique: LTT-167
script: 04-computation/tournament_n4_dual_basis_erdos870_codex_s276.py
result: 05-knowledge/results/tournament_n4_dual_basis_erdos870_codex_s276.out
reflection: 07-reflections/n4-tournament-subbasis-erdos870-codex-s276.md
related:
  - HYP-3142
  - HYP-3141
  - HYP-3140
  - HYP-3139
  - HYP-3138
  - HYP-3137
  - HYP-3136
  - HYP-3135
  - HYP-3134
  - HYP-3133
  - HYP-3132
  - HYP-3129
  - HYP-3124
  - HYP-3118
  - HYP-3106
  - HYP-3002
  - HYP-2998
  - OPEN-Q-108
links:
  - https://github.com/davidturturean/erdos-870
---

# HYP-3143: n=4 Tournament Subbasis Packet

## Claim

The two n=4 tournament tables in the prompt are not two cosmetic
presentations of the same data.  They are two quotient policies.

The Hamiltonian-path tiling model fixes the path

```text
0 -> 1 -> 2 -> 3
```

and leaves three chords `a=02`, `b=13`, `c=03`.  It is faithful to the raw
tiling geometry, but the induced class map from the 3-bit chord cube is not a
unique representation system:

```text
E  -> T
a  -> +
b  -> -
c  -> S
ab -> S
ac -> S
bc -> S
abc -> S
```

So the class `S` appears at squarefree flip orders `1`, `2`, and `3`.  This
is a lower-order leakage failure.

The second model fixes four arcs whose partial score sequence is `0,1,1,2`
and leaves two arcs `x,y`.  The exact finite search finds `12` witnesses, all
with free arcs a perfect matching.  A canonical representative fixes the
cycle-like filler

```text
01, 03, 12, 23
```

all forward and leaves `x=02`, `y=13`.  Then

```text
E  -> T
x  -> +
y  -> -
xy -> S
```

is a minimal two-bit class basis: each n=4 isomorphism class is represented
once, and `S` first appears at squarefree order `2`.

## Erdős-870 Transfer

The `erdos-870` proof package solves Erdős problem 870 by building minimal
subbases: every sufficiently large integer has a unique representation of the
declared order `h`, and no representation of smaller order.  The n=4 table
above gives a finite tournament analogue:

```text
packet_subbasis =
  fixed filler
  + free obstruction basis
  + exact-order / no-lower-order certificate
  + collision sidecar or named debt
```

Scheme B is the order-2 subbasis model.  Scheme A is the useful tiling
coordinate model, but it is not theorem-facing until its `S` collision fiber
is repaired or explicitly carried.

## LRC14 Use

This adds a new guardrail for the labelled-packet LRC14 proof route:

```text
A quotient is allowed to forget only if every terminal obstruction class
still has a unique squarefree packet word of the declared order, or if the
collision fiber is named and repaired.
```

That rule lines up with the current frontier:

- HYP-3141 says directed edges must keep tail/tip information before becoming
  proof witnesses.
- HYP-3142 says the terminal k=8 bounded-core sidecar must be reached with the
  exact moment/resolvent packet, not by a lower-order scalar shortcut.
- HYP-3138 and HYP-3139 already found concrete versions of collision debt:
  odd-coordinate leakage and center/boundary leakage.
- HYP-3133/HYP-3134 make A000568 a global-consistency quotient, not a raw
  class-count eraser.

## Assumption Challenge

This session explicitly did not take tournament vertices to be LRC runners.
The tested vertex and carrier choices were:

- labelled arc-flip states,
- fixed filler arcs,
- free matching bits,
- class fibers,
- squarefree representation words,
- proof-carrier sidecars,
- A000568 global class shadows,
- bounded-core LRC exit packets.

The quotient preserves the n=4 isomorphism class and, in Scheme B, the first
representation order.  It destroys raw labels, endpoint roles, and raw chord
geometry unless those are retained as filler or collision sidecars.  The
challenged assumption is that a fixed Hamiltonian path plus all remaining
chords is automatically the right tournament coordinate system; the exact
enumeration says it is visually faithful but not representation-faithful.

## Next Tests

1. Enumerate n=5 and n=6 partial assignments and search for minimal
   class-bases whose fibers separate A000568 classes with lower-order
   exclusion.
2. Add a `packet_order` column to HYP-3141 edge witnesses: at what order can
   each proof class first be represented after fixed filler?
3. Audit HYP-3142 bounded-core U4 sidecars by packet order: the terminal k=8
   obstruction should not be reachable by a lower-order quotient word.
4. Treat q=3 unital four-point blocks as `C4 + matching` packets: boundary
   cycle as deterministic filler, matching diagonals as basis bits, and
   AP/Goddyn-Wong labels as sidecar colors.
