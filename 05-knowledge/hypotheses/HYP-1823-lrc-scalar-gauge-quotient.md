---
id: HYP-1823
status: EXPLORATORY
source: codex-2026-05-31-S367
related:
  - THM-363
  - THM-358
  - HYP-1817
  - HYP-1818
  - HYP-1819
  - HYP-1824
  - HYP-1825
  - HYP-1832
  - HYP-1833
  - HYP-1828
---

# HYP-1823: Fourteen-runner scalar-gauge quotient lemma

## Statement

In the `n=14`, `k=13` Lonely Runner micro-staircase system, quotient residue
vectors by scalar ramps:

```text
v_i ~ v_i + m i mod 14.
```

Normalize every class by setting `v_1=0`.  Then the zero class is the only full
blocker.  Equivalently, every nonzero normalized vector
`v in (Z/14Z)^13` has some shift `s mod 14` and some cell of the arrangement
of `floor(14 {i alpha})` such that

```text
s v_i + floor(14 {i alpha}) notin {0,13} mod 14
```

for all `i=1,...,13`.

## Evidence

The S367 script builds the exact `n=14` cell system with `812` alpha patterns
and `11368` candidate `(s, cell)` pairs.  It checks that adding a scalar ramp
preserves raw coverage: in `200` random trials, all `14` gauge shifts preserved
the blocked-cell count.

All scalar ramps

```text
v_i = m i mod 14
```

are full blockers, and all normalize to the zero vector.  After quotienting,
the best local nonzero blocker attempt remains

```text
(0,0,0,0,0,7,0,0,0,0,0,0,0),
```

which misses `56` of `11368` candidates.

Exact subfamily scans strengthen the signal:

```text
support 1: scanned 156,    best missed 56
support 2: scanned 11154,  best missed 112
support 3: scanned 483340, best missed 126
2-torsion: scanned 4095,   best missed 56, best_count 1
```

Thus the complete normalized `2`-torsion cube has a unique extremal, the
coordinate-6 half-turn, and it still has explicit safe cells.

For that extremal, all `56` misses occur at odd shifts: eight alpha cells for
each of the seven odd shifts.  The missed cells have positive margin `1` and
widths

```text
1/728, 1/882, 1/1176, 1/1386.
```

The first cell is

```text
alpha in [29/182, 9/56)
bins=(2,4,6,8,11,13,1,3,6,8,10,12,1).
```

S377 extended this quotient pressure in two directions.  The exact `n=14`
two-torsion cube still has the coordinate-6 half-turn as a unique global
minimum, and adding a second defect to it raises the miss count to at least
`126`.  The natural-gate scan also shows that the hardest one-defect
coordinates are product-sum target modes: coordinate `6` for `n=14`, and
coordinates `6` and `14` for `n=15`.

## Why It Matters

HYP-1818 says scalar ramps must be excised before the composite
micro-staircase proof can work.  HYP-1823 is the sharper quotient version: the
scalar ramps are one gauge orbit, and the real finite problem is to show that
no nonzero quotient class can block every cell.

If proved, this would convert the S363/S364 computational obstruction into a
finite lemma suited for the fourteen-runner case: scalar zero is the
Dirichlet-equality spine, and every other class has a micro-staircase witness.

## Test Plan

1. Use THM-363 to treat scalar addition as a genuine gauge quotient rather
   than an exceptional family.
2. Turn the exact `2`-torsion scan into a human-readable interval proof, using
   the coordinate-6 extremal and its eight mirror cells as the hard case.
3. Build a branch-and-bound certificate for all normalized vectors, ordering
   coordinates by distance from the scalar line and pruning once a target
   stencil interval is forced.
4. Split shifts into units, `s=7`, and even nonunits.  Route nonunit behavior
   to quotient descent; use unit shifts to force a positive-margin cell.
5. Compare the resulting quotient classes with the initial bad sets
   `I(13,p,1)` from the public verifier and turn safe cells into cover-search
   pruning rules.

## Sources

- `04-computation/lonely_runner_k13_scalar_gauge_s367.py`
- `05-knowledge/results/lonely_runner_k13_scalar_gauge_s367.out`
- `07-reflections/lonely-runner-fourteen-runner-scalar-gauge-s367.md`
- THM-363
- HYP-1817
- HYP-1818
- HYP-1824
- HYP-1825
- HYP-1828
- HYP-1832
- HYP-1833
