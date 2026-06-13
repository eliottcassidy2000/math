# HYP-1836: LRC is fundamentally an endpoint-incidence problem

**Status:** EXPLORATORY; supported by S379 exact representation comparison.

## Claim

The most fundamental finite object in the Lonely Runner Conjecture is not the
speed set itself.  It is the circular-arc endpoint-incidence system generated
when the forbidden intervals are pulled back to the common time circle.

In this representation:

- a speed is a character `t -> v*t` on `R/Z`;
- a speed set is a generator of a family of pulled-back unsafe arcs;
- an endpoint is a boundary event where one runner is exactly at distance
  `1/(k+1)`;
- a protected endpoint is a boundary event lying strictly inside another
  unsafe arc;
- a tight example is an open cover failure supported only at endpoints;
- a counterexample would be a full open cover whose endpoint graph has no
  unprotected boundary point;
- scalar ramps are the Dirichlet equality spine/gauge orbit;
- torsion and product-sum gates mark quotient layers where protection leaks.

Thus the proof question should be reframed from "where is the lonely time?" to:

```text
Why must the endpoint-protection graph have a leaf after the scalar equality
spine is quotiented?
```

## Evidence

`lonely_runner_shape_questions_s379.py` compares the same examples as speed
sets, interval covers, endpoint-protection graphs, interval-overlap shadows,
and exact max-min loneliness landscapes.

The tight examples in the sample are endpoint failures, not interval failures:

```text
initial n=8:  critical/th=1, boundary_witnesses=4, unit_mod_8
sporadic n=8: critical/th=1, boundary_witnesses=4, unit_mod_8
initial n=14: critical/th=1, boundary_witnesses=6, unit_mod_14
initial n=15: critical/th=1, boundary_witnesses=8, unit_mod_15
```

Tiny-gap near-disproofs are not actually close in the max-min landscape or the
endpoint graph:

```text
n14 seven-ladder: gap/th=0.005411, critical/th=1.217391, unprotected=84
n15 3x5 ladder:  gap/th=0.080000, critical/th=1.363636, unprotected=150
n15 mixed gates: gap/th=0.120000, critical/th=2.142857, unprotected=112
```

The seven-ladder is close only in horizontal gap-width space.  Vertically, its
best lonely radius is still far above the conjectural threshold; structurally,
it has `84` unprotected endpoints in a higher-denominator quotient layer.

## Predictions

1. Apparent near-counterexamples with tiny visible gaps will often have large
   critical-radius surplus or many unprotected endpoints.
2. A successful proof-search invariant should use endpoint leaf/peel structure,
   not only interval gap width.
3. Abstract circular-arc endpoint systems with no peelable leaf may exist but
   fail realizability by integer speed sets; realizability is a separate
   arithmetic constraint.
4. The S367/S371 eight alpha stencils should have a characterization as
   endpoint-incidence leaf types, not only as micro-staircase cells.
5. A useful LRC homology theory, if any, should be built from the protected
   endpoint incidence/repair graph rather than from the interval nerve alone.

## See

- `04-computation/lonely_runner_shape_questions_s379.py`
- `05-knowledge/results/lonely_runner_shape_questions_s379.out`
- `07-reflections/lonely-runner-shape-questions-s379.md`
- HYP-1811, HYP-1813, HYP-1831, HYP-1835
