---
id: HYP-3150
title: LRC14 function-compression packets may keep the hard core below the Abel-Ruffini wall
status: RESERVED / executable factor-through scout pending; not a proof
source: codex-2026-06-28
tangent: T1215
technique: LTI-276
tournament_technique: LTT-174
script: 04-computation/lrc14_function_compression_resolvent_wall_codex_20260628.py
result: 05-knowledge/results/lrc14_function_compression_resolvent_wall_codex_20260628.out
related:
  - HYP-3149
  - HYP-3148
  - HYP-3147
  - HYP-3146
  - HYP-3145
  - HYP-3144
  - HYP-3143
  - HYP-3142
  - HYP-3141
  - HYP-3140
  - HYP-3139
  - HYP-3138
  - HYP-3137
  - HYP-3135
  - HYP-3134
  - HYP-3133
  - HYP-3132
  - HYP-3129
  - HYP-3122
  - THM-084
  - THM-577
  - OPEN-Q-108
external_sources:
  - https://github.com/davidturturean/erdos-870
  - https://www.erdosproblems.com/870
---

# HYP-3150: Function-Compression Below the Resolvent-Degree Wall

## Reserved Claim

The next scout will test the user's suggested synthesis as a single legality
calculus:

```text
compression q: X -> Y
observable  f: X -> Z
legal iff   f factors through q,
or the missing part is carried by a named sidecar.
```

This is meant to unify the nearby exact packets:

- HYP-3147/HYP-3144: K3 edge flips collapse to a two-class kernel, but the
  minority edge, Worpitzky descent word, state-level PGF curve, and ordered
  pair functions may not factor through the score-class quotient.
- HYP-3149/HYP-3146/HYP-3148: the K4 fixed-path cube factors to a legal
  two-bit class table only through a non-linear compression
  `x = a OR c`, `y = b OR c`, or by fixing the canary/filler coordinate.
- HYP-3142/HYP-3132/HYP-3139: the k=8 bounded-core dual is a quartic
  resolvent, but after the reflection coordinate `u=t-3` it is even:
  `u^4 - 5u^2 + 4 = (u^2-1)(u^2-4)`, so the proof-facing variable is
  the quadratic `v=u^2` together with resurrection sidecars for odd data.

The speculative meta-point to test is:

```text
LRC14 stays solvable because every hard compression seen so far has
effective degree <= 4, and the deepest k=8 node is a biquadratic below the
generic quintic / Abel-Ruffini wall.
```

This is not a theorem claim.  The executable scout should separate exact
factor-through facts from creative proof-route hypotheses.

## Planned Measurements

1. Audit pair functions on ordered pairs:
   `a+b` and `a*b` should factor through unordered-pair quotienting, while
   `a^b` and `b^a` should fail except at accidental equalities.
2. Recompute the K3 two-class edge-flip graph and record which observables
   factor through the `C/T` quotient, which require the minority-edge or
   ordered-function sidecar, and which only agree after aggregation.
3. Recompute the K4 fixed-path cube and the two-bit canary/filler scaffold.
   Test `x=a OR c`, `y=b OR c` as a class-preserving compression and record
   which information it destroys: flip word, fiber PGF, deletion stability,
   and canary status.
4. Recompute the k=8 quartic fold:
   `g(t)=(t-1)(t-2)(t-4)(t-5)`, `u=t-3`,
   `g(u+3)=u^4-5u^2+4`, `v=u^2`, `v^2-5v+4`.
   Record the degree drop from quartic to quadratic and the sidecar debt for
   odd coordinates.
5. Run Tournament Analysis over proof carriers rather than runners or raw
   arcs.  Candidate vertices: factor-through audit, K3 kernel, K4 OR
   compression, canary sidecar, fiber-PGF curve, Worpitzky/order sidecar,
   resolvent-even fold, raw scalar value, and quintic-wall alarm.

## Assumption Challenge

The vertex set for the scout is not fixed in advance.  It should explicitly
compare at least these carriers:

```text
functions, quotient fibers, edge roles, fixed-path arc words,
canary/filler coordinates, generating-function curves, resolvent variables,
Fourier/SPEC modes, and proof obligations.
```

The quotient preserves the LRC predicate only when the target observable is
constant on visible fibers or reconstructible from sidecars.  It destroys
orientation, order, deletion robustness, odd coordinates, or full PGF curves
unless those are named.

## Guardrail

Do not promote the Abel-Ruffini wall phrase into proof.  The exact content to
verify is smaller and stronger: every compression used in the current LRC14
endgame must either be a finite factor-through map of bounded degree or carry
the missing coordinate as an explicit sidecar.
