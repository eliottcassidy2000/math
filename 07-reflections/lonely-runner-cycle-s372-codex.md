# Lonely Runner 14/15/Disproof Cycle S372

**Session:** codex-2026-05-31 S372

This run cycled through the requested three routes:

1. attack the fourteen-runner scalar-gauge quotient,
2. force the obstruction shape into the fifteen-runner case,
3. lift the resulting quotient near-blockers into actual speed-set disproof
   pressure tests.

The computation is reproducible at:

- `04-computation/lonely_runner_cycle_s372.py`
- `05-knowledge/results/lonely_runner_cycle_s372.out`

## Fourteen-Runner Route

The `n=14,k=13` quotient system has `812` alpha patterns and `11368`
candidate `(shift,cell)` pairs.  Exact normalized support scans through
support `4`, the complete two-torsion cube, and local search all keep the same
best nonzero quotient vector:

```text
v = (0,0,0,0,0,7,0,0,0,0,0,0,0)
missed = 56
```

New extension beyond S367/S371: the exact support-4 scan covers `14,137,695`
vectors.  Its best vector misses `182` cells:

```text
(0,0,0,7,7,7,0,7,0,0,0,0,0)
```

So adding more small support does not approach a full blocker; it moves away
from the coordinate-6 half-turn extremal.

## Fifteen-Runner Route

For `n=15,k=14`, the scalar-gauge quotient analogue has `960` alpha patterns
and `14400` candidate pairs.  The best nonzero vectors are single-coordinate
order-3 leaks:

```text
coord 14, residue 5 or 10: missed = 120
coord  6, residue 5 or 10: missed = 120
```

Exact support scans found:

```text
support 1: best missed 120
support 2: best missed 220
support 3: best missed 280
full {0,5,10}^{13} subgroup: best missed 120
```

This gives a concrete transfer principle: the `n=14` half-turn leak becomes an
`n=15` order-3 leak at staircase coordinates `6` and `14`.

## Disproof Pressure

The new disproof test was residue-template lifting.  Instead of random gated
sets, it took quotient near-blockers, added scalar ramps, lifted those residue
multisets to integer speeds, and scored exact forbidden-interval coverage.

Best first-pass lift:

```text
n=15 speeds=(3,12,14,15,36,48,54,57,60,66,72,78,81,84)
gap/thresh = 0.023504
```

Best second-pass lift:

```text
n=15 speeds=(2,3,6,9,12,15,18,21,24,27,30,33,36,39)
gap/thresh = 0.025641
```

No lifted template became an open-cover candidate.  The guided search is more
dangerous than unguided gate sampling, but it still leaks a positive complement
gap.

## Next Move

Classify single-coordinate subgroup defects for composite `n`: half-turns for
even `n`, order-3 leaks for multiples of `3`, and the product-sum target
coordinates suggested by the neighboring S372 notes.  The proof target is a
finite stencil certificate for these subgroup leaks before attempting arbitrary
normalized vectors.
