---
id: HYP-1827
status: EXPLORATORY
source: codex-2026-05-31-S372
related:
  - HYP-1817
  - HYP-1818
  - HYP-1819
  - HYP-1829
  - THM-358
  - THM-360
---

# HYP-1827: LRC scalar-puncture moat

## Statement

In the composite micro-staircase system for the Lonely Runner frontier, the
scalar-ramp line

```text
v_i = m i mod n,  i=1,...,n-1
```

is an exact Dirichlet-equality blocker, but it is surrounded by a moat: small
non-scalar perturbations force many explicit witness cells to open.

For the fourteen-runner case `n=14`, the exact full-cell system has:

```text
all scalar ramps:        0 witness cells
best radius-1 puncture: 56 witness cells
best radius-2 puncture: 112 witness cells
```

The best radius-1 punctures are precisely coordinate `i=6` shifted by `+7`
modulo `14`, transported along the scalar ramp family.  The S364 best
non-scalar near-blocker is one of these transported punctures:

```text
(8,2,10,4,12,13,0,8,2,10,4,12,6).
```

## Evidence

`04-computation/lonely_runner_creative_multiroute_s372.py` builds the exact
full discontinuity arrangement for

```text
floor(n * {i alpha}),  i=1,...,n-1,
```

then exhausts all radius-1 and radius-2 punctures around every scalar ramp for
`n=14` and `n=15`.

For `n=14`, radius 1 checks `2366` vectors and radius 2 checks `184548`
vectors.  The exact minima are:

```text
radius 1: best_missed=56,  best_count=14, positions=(6,),    deltas=(7,)
radius 2: best_missed=112, best_count=14, positions=(8,12),  deltas=(7,)
```

For `n=15`, the analogous minima are:

```text
radius 1: best_missed=120, best_count=60, positions=(6,14), deltas=(5,10)
radius 2: best_missed=220, best_count=60, positions=(6,14), deltas=(5,10)
```

The `n=14` zero-ramp puncture moat has a clean cell anatomy: seven odd
`s`-layers, each with eight alpha-cells, and total alpha width `7/858`.

Speed-set pressure in the same run did not produce open-cover candidates.
Hand-built and random `14`- and `15`-gated families leaked through positive
gaps or boundary witnesses; all endpoint peels ended with empty cores.

## Interpretation

HYP-1818 explained why scalar ramps must be excised.  HYP-1827 strengthens the
picture: the scalar line is not merely an obstruction family; it is locally
isolated.  A proof can try to show that any vector close to the scalar line
automatically exposes a finite package of witness cells.

This turns the micro-staircase problem into two parts:

1. Prove the scalar-puncture moat symbolically.
2. Prove that any far-from-scalar blocker descends through quotient/endpoint
   protection and cannot maintain an all-protected core.

## Test Plan

1. Prove the `n=14` radius-1 moat exactly from the eight alpha-cell intervals.
2. Extend the proof to radius 2, where the minimum doubles to `112` cells.
3. Identify the general divisor-layer rule behind the atlas positions and
   deltas, especially `n=14` shift `+7` and `n=15` shifts `+5,+10`.
4. Add a SAT/backtracking search constrained to vectors at distance at least
   three from every scalar ramp.
5. Combine scalar-distance rank with endpoint peel depth when ranking possible
   counterexample constructions.

## Sources

- `04-computation/lonely_runner_creative_multiroute_s372.py`.
- `05-knowledge/results/lonely_runner_creative_multiroute_s372.out`.
- `07-reflections/lonely-runner-creative-multiroute-s372.md`.
- `05-knowledge/hypotheses/HYP-1818-lrc-scalar-ramp-excision.md`.
- `05-knowledge/hypotheses/HYP-1829-lrc-scalar-distance-endpoint-closure.md`.
