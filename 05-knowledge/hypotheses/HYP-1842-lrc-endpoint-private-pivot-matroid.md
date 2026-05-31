---
id: HYP-1842
status: EXPLORATORY
source: codex-2026-05-31-S385
related:
  - HYP-1811
  - HYP-1815
  - HYP-1836
  - HYP-1843
  - THM-359
---

# HYP-1842: LRC endpoint protection is controlled by a private-pivot matroid

## Statement

For a Lonely Runner speed set, form the binary endpoint-protection matrix over
`GF(2)` whose rows are forbidden endpoints and whose columns are pulled-back
forbidden intervals.  A counterexample should require an unusually
self-supporting representation of this matrix: no unprotected row, no peelable
leaf, and no private pivot after quotienting the scalar equality spine.

The working hypothesis is stronger:

```text
Every integer-realizable LRC endpoint-protection matrix has a private pivot or
peels to the scalar boundary skeleton.
```

In this framing, the obstruction is matroidal before it is metric.  The metric
data chooses a realizable circular-arc representation; the proof should show
that its protection matroid cannot be fully self-supported.

## Evidence

S385 computed endpoint-incidence fingerprints for tight examples, fourteen-
runner near-disproofs, and fifteen-runner analogues.  Every sampled system
peeled to terminal `coreE=0`, even when the visible gap was extremely small:

```text
initial n=14:       coreE=0, rank/null=79/91,    privateE=128
n14 seven-ladder:  coreE=0, rank/null=514/552,  privateE=710
n14 14-ladder:     coreE=0, rank/null=1025/1105, privateE=1420
n15 3x5 ladder:    coreE=0, rank/null=229/221,  privateE=204
n15 mixed gates:   coreE=0, rank/null=224/298,  privateE=232
```

The strongest rank pressure in the sample was the sporadic tight `n=8`
example, with `rankP=0.511905`, but it still had exposed endpoint rows and
`coreE=0`.  High rank pressure alone is not enough; the endpoint system still
finds private pivots.

## Predictions

1. Exhaustive primitive boxes should continue to show `coreE=0` after
   endpoint peeling, even in tiny-gap quotient-ladder families.
2. Abstract circular-arc endpoint systems with nonzero all-protected cores may
   exist, but most should fail integer-speed realizability.
3. A proof should be possible by finding a canonical private pivot in each
   quotient layer, then descending through THM-359 peeling.
4. The S367/S371 eight alpha stencils should be describable as private-pivot
   types in the normalized fourteen-runner endpoint matroid.

## Test Plan

1. Search abstract circular-arc endpoint systems for nonzero protection cores.
2. Add integer-realizability constraints and see which abstract cores survive.
3. Compute Smith or Tutte-style invariants of the endpoint-protection matrix
   for the S373-S385 near-disproof families.
4. Try to prove that any primitive integer realization with full forbidden
   measure has a private pivot unless it is the scalar unit skeleton.

## Sources

- `04-computation/lonely_runner_hypothesis_noise_s385.py`
- `05-knowledge/results/lonely_runner_hypothesis_noise_s385.out`
- `07-reflections/lonely-runner-hypothesis-noise-s385.md`
- HYP-1815
- HYP-1836
- THM-359
