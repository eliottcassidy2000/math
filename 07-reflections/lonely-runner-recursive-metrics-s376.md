# Lonely Runner Recursive Metrics S376

The user asked for the Lonely Runner analogue of the tournament-structure
program: not just "try the next case", but understand what changes as `n`
changes.  S376 is a first atlas.

The best comparison is:

```text
tournaments:  n vertices -> class counts, SC layers, endpoint transfer, residues
LRC:          denominator n -> cells, unit endpoints, scalar moats, endpoint closure
```

## The Metric Menu

The script tracks five layers.

`C(n)`: the number of micro-staircase alpha floor-pattern cells.  This is the
closest analogue of a tournament class-count sequence.  It grows irregularly
and jumps at divisor-rich values:

```text
C(14)=812, C(15)=960, C(16)=1152, C(17)=1360, C(18)=1728.
```

The transition jumps are more informative than the raw counts:

```text
13->14: +214
14->15: +148
15->16: +192
16->17: +208
17->18: +368
```

`phi(n)/(n-1)`: the density of the initial-segment unit skeleton.  Primes
maximize this at `1.0`; composites expose quotient layers immediately.  This
is the LRC version of a symmetry/regularity axis.

Scalar-puncture moat: the best radius-1 deviation from a scalar ramp.  It is
a local curvature invariant around the Dirichlet equality spine.  Composite
`n` often shows nonunit deltas:

```text
n=14: moat=56,  delta_gcd=7
n=15: moat=120, delta_gcd=5
n=16: moat=128, delta_gcd=8
n=18: moat=198, delta_gcd=9
```

Endpoint closure: the initial-segment endpoint system peels to empty for all
checked `n`, but the peel depth itself is a recursive signal.  At `n=14`,
the initial segment has `176` endpoint values and peel depth `15`; at `n=18`,
it has `300` endpoint values and depth `18`.

Gate and quotient-ladder pressure: these are the counterexample-family
metrics.  Quotient ladders start producing tiny global gaps at composite
values:

```text
n=14: gap/th=0.005411, unprotected=84
n=16: gap/th=0.007576, unprotected=140
n=18: gap/th=0.005682, unprotected=176
```

That is the S374 lesson in recursive form: tiny gap does not mean close to a
counterexample unless endpoint closure also improves.

## What Changed Conceptually

Before this run, the LRC thread had strong fixed-case objects: unit boundary
skeletons, scalar ramps, scalar punctures, endpoint cores, and the
fourteen-runner seven-ladder.  S376 turns those into coordinates of one
recursive object.

The likely proof object is not a single scalar invariant.  It is a vector:

```text
(arrangement complexity,
 unit density,
 scalar moat curvature,
 endpoint closure defect,
 quotient-ladder gap pressure).
```

This feels much closer to how the tournament side matured.  Counts were never
enough there either; the useful information came from residues, transfer
matrices, projection defects, and special families.  For LRC, the analogous
warning is: scalar moat data without endpoint incidence is incomplete, and
endpoint incidence without scalar distance is equally incomplete.

## Next Move

The next constructive step is an endpoint-transfer matrix for LRC.  Rows
should be leak layers at denominator `n`, columns should be leak layers at
`n+1`, and entries should count or weight protecting speeds that move one
layer into another.  If that object has rank or peelability behavior like the
tournament endpoint-transfer matrices, it may finally explain why quotient
ladders get so close and still fail.

Artifacts:

- `04-computation/lonely_runner_recursive_metrics_s376.py`
- `05-knowledge/results/lonely_runner_recursive_metrics_s376.out`
- `05-knowledge/hypotheses/HYP-1831-lrc-recursive-metric-tower.md`
