---
id: HYP-3147
title: n=3 edge flips form a two-class Markov kernel with Worpitzky/function sidecars
status: EVIDENCE / exact n=3 enumeration and proof-packet proposal; not an LRC proof
source: codex-2026-06-28-S277
tangent: T1212
technique: LTI-273
tournament_technique: LTT-171
script: 04-computation/tournament_n3_edge_flip_worpitzky_codex_s277.py
result: 05-knowledge/results/tournament_n3_edge_flip_worpitzky_codex_s277.out
reflection: 07-reflections/n3-edge-flip-worpitzky-function-kernel-codex-s277.md
related:
  - HYP-3144
  - HYP-3146
  - HYP-3145
  - HYP-3143
  - HYP-3142
  - HYP-3141
  - HYP-3139
  - HYP-3138
  - HYP-3134
  - HYP-3133
  - HYP-3129
  - HYP-3124
  - HYP-3106
  - THM-084
  - OPEN-Q-108
---

# HYP-3147: n=3 Edge-Flip Worpitzky/Function Kernel

## Claim

The tournament on three vertices, viewed from its three edges, is the local
transition kernel that HYP-3143 needs beneath its n=4 packet-basis audit.

Use a fixed cyclic reference

```text
0 -> 1 -> 2 -> 0
```

and record each edge as a coin flip against that reference.  Then:

- `HHH` and `TTT` are the two cyclic tournaments, score class `C=(1,1,1)`.
- The six `2:1` coin mixes are transitive tournaments, score class
  `T=(0,1,2)`.

Flipping one edge gives the exact collapsed two-class kernel:

```text
P(C -> T) = 1
P(T -> C) = 1/3
P(T -> T) = 2/3
```

with stationary distribution

```text
pi(C,T) = (1/4, 3/4),
```

matching the straight/mixed coin split `2/8,6/8`.  The nontrivial eigenvalue
is `-1/3`, the local signed oscillation left after collapsing the cube to the
two tournament classes.

The edge-level rule is sharper:

```text
cyclic state: every edge flip is a C -> T gate;
transitive state: exactly one minority edge is a T -> C gate, and two majority
edges are T -> T self-class flips.
```

## Worpitzky Sidecar

Inside the transitive fiber, the source-to-sink linear order has descent
counts

```text
1, 4, 1.
```

This is the Eulerian row `A(3,k)`, and it verifies the n=3 Worpitzky identity

```text
x^3 = binom(x+2,3) + 4 binom(x+1,3) + binom(x,3).
```

Thus Worpitzky does not replace the `C/T` edge kernel.  It refines the
transitive fiber by ordered path descents.  This is a useful correction to a
too-coarse quotient: after a packet knows it is in `T`, it still needs a
descent/linear-order sidecar if later proof moves depend on ordered paths.

## Function Quartet

For an edge with endpoint values `(a,b)`, the four functions in the prompt
split into two quotient types:

```text
a+b, a*b        symmetric / unordered edge shadow
a^b, b^a        ordered / orientation-sensitive channels
```

Sum and product cannot see an arc flip.  The ordered exponential pair can see
which endpoint is base and which is exponent.  The LRC guardrail is therefore:

```text
symmetric pair functions alone cannot certify an edge-flip quotient;
an ordered channel or a named orientation sidecar is required.
```

## LRC14 Use

Candidate local packet:

```text
n3_edge_kernel_packet =
  C/T class
  + minority_edge_gate
  + Worpitzky_descent_word
  + ordered_function_channel
```

Transfer to the current proof route:

- HYP-3141 edge witnesses should add `edge_flip_class_kernel`,
  `minority_edge_gate`, `worpitzky_descent_word`, and
  `ordered_function_payload`.
- HYP-3142 should treat the k=8 antisymmetric nonmax block as
  order-sensitive payload, not disposable noise.
- HYP-3143 should compute `packet_order` only after straight/mixed separation
  and Worpitzky refinement of transitive subtriangles.
- HYP-3129's signed SPEC low modes should be tested against the local `-1/3`
  edge-kernel eigenmode.

## Assumption Challenge

The vertices here are not runners and not even tournament classes alone.  The
live objects are edge flips, coin words, minority/majority edge roles,
source-to-sink orders, descent words, and functions on ordered endpoint pairs.

The quotient preserves only the `C/T` score class and the one-step edge-flip
transition probabilities.  It destroys the linear order inside `T`, the
minority edge, and the ordered endpoint function channel unless these are
attached as sidecars.

## Next Tests

1. Add these fields to the n=4 packet-basis scout: minority/majority flip role,
   Worpitzky descent word on each transitive subtriangle, and ordered-function
   payload.
2. For HYP-3141 edge witnesses, compute the two-state edge kernel on every
   3-vertex shadow around a directed proof edge.
3. For HYP-3142, test whether the k=8 antisymmetric shell can be bounded by
   summing n=3 ordered-function oscillations with Worpitzky descent weights.
4. For LRC `Rprime` rows, test whether the `-1/3` eigenmode aligns with the
   HYP-3129 signed SPEC low modes.
