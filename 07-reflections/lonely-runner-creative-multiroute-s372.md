---
source: codex-2026-05-31-S372
status: research sprint
tags:
  - lonely-runner
  - fourteen-runners
  - micro-staircase
  - scalar-ramp
  - endpoint-descent
---

# Lonely Runner Creative Multiroute Sprint

The useful move was to stop treating the S364 near-blocker as a mysterious
non-scalar vector.  It is a punctured scalar ramp.

S364's best `n=14` near-blocker was

```text
(8,2,10,4,12,13,0,8,2,10,4,12,6).
```

The scalar ramp with multiplier `8` is

```text
(8,2,10,4,12,6,0,8,2,10,4,12,6).
```

So the whole obstruction is one coordinate change:

```text
i=6: 6 -> 13 = 6+7 mod 14.
```

That matters because the scalar ramps are the Dirichlet equality spine.  If
their nearest non-scalar shadows already leak many cells, then the proof can
surround the spine with a moat instead of trying to classify all residue
vectors uniformly.

## Routes Tried

### 1. Scalar-puncture atlas

I computed the best one-coordinate scalar-ramp puncture for every `n=4..22`.
For the frontier values:

```text
n=14: best_h1_missed=56,  positions=(6,),    deltas=(7,)
n=15: best_h1_missed=120, positions=(6,14),  deltas=(5,10)
```

The positions/deltas are divisor-layer data, not random residue noise.  For
`14=2*7`, the best leak shifts by the half-modulus `7`; for `15=3*5`, the best
leaks shift by the `5`-layer.

### 2. Exact radius 1 and 2

For `n=14`, the script exhausts all punctures around every scalar ramp:

```text
radius 1 checked 2366 vectors,   best missed 56
radius 2 checked 184548 vectors, best missed 112
```

For `n=15`:

```text
radius 1 checked 2940 vectors,   best missed 120
radius 2 checked 267540 vectors, best missed 220
```

So at least through radius 2, moving away from the scalar spine opens more
witness cells, not fewer.

### 3. Greedy search

Fresh greedy searches from random non-scalar starts fell back to distance-one
punctures:

```text
n=14: missed=56,  distance_to_scalar=1
n=15: missed=120, distance_to_scalar=1
```

That says the local search landscape is being pulled toward the scalar moat.
The best near-blockers are not evidence for hidden full blockers nearby; they
are evidence that the scalar line is locally isolated.

### 4. Missed-cell anatomy

The zero-ramp `n=14` puncture

```text
(0,0,0,0,0,7,0,0,0,0,0,0,0)
```

has exactly `56` witness cells, organized as:

```text
7 odd s-layers * 8 alpha-cells.
```

The total alpha width of the eight-cell package is `7/858`.  The transported
S364 vector has the same `56` cells spread over transported alpha packets.

### 5. Gate pressure

Hand-built and random speed-set constructions still failed to make an open
cover.  The closest hand-built `n=14` route in this run was

```text
(1,2,3,4,5,7,8,9,10,11,12,13,14)
```

with a positive gap `5/1848`, only eight boundary witnesses, peel depth `27`,
and empty endpoint core.  It is a nice pressure case, but it still leaks.

## New Proof Target

The next attempt should prove the scalar-puncture moat directly:

```text
Any n=14 scalar ramp with one coordinate changed opens at least 56
micro-staircase witness cells.
```

Then prove the radius-2 moat, and use scalar-distance as a ranking invariant.
If a putative non-scalar blocker exists, it must live far from every scalar
ramp.  That is exactly where endpoint-protection descent should be strongest.

The slogan for the next session:

```text
scalar spine, puncture moat, far-field endpoint descent
```

## Artifacts

- `04-computation/lonely_runner_creative_multiroute_s372.py`
- `05-knowledge/results/lonely_runner_creative_multiroute_s372.out`
- `05-knowledge/hypotheses/HYP-1827-lrc-scalar-puncture-moat.md`
