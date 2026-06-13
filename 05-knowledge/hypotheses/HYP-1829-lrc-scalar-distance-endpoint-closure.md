---
id: HYP-1829
status: EXPLORATORY
source: codex-2026-05-31-S374
related:
  - HYP-1823
  - HYP-1827
  - HYP-1828
  - HYP-1831
  - THM-357
  - THM-360
---

# HYP-1829: LRC blockers must optimize scalar distance and endpoint closure together

## Statement

For the fourteen-runner frontier, the two current obstruction languages are
not independent:

1. micro-staircase residue vectors near the scalar-ramp spine hit the
   scalar-puncture moat; and
2. speed-set counterexample constructions that almost cover the circle hit an
   endpoint-closure explosion.

A viable proof or disproof should therefore optimize the pair

```text
(distance from the scalar-ramp quotient, endpoint-protection closure size)
```

rather than either coordinate alone.  In particular, a fourteen-runner
counterexample cannot merely be a tiny-gap quotient ladder; it must also close
the endpoint-protection cycle without falling back into the scalar-puncture
moat.  Conversely, a fourteen-runner proof should try to show that every
attempt to reduce one coordinate forces the other coordinate to grow.

## Evidence

S374 forced a loop through the three active routes: fourteen-runner proof
pressure, fifteen-runner transfer, and direct disproof construction.

On the fourteen-runner cell side, the known best non-scalar near-blocker

```text
(8,2,10,4,12,13,0,8,2,10,4,12,6)
```

is scalar ramp `m=8` with one defect.  The mutation table has exactly one
complete repair: reverting the punctured coordinate `6` from `13` to `6`.
Every other best listed one-coordinate mutation still misses `96` cells, and
the best true second-defect pressure has

```text
missed = 126
count  = 1
```

So local scalar repair is a dead end unless it returns to the scalar equality
spine.

On the fifteen-runner transfer side, the scalar-puncture moat thickens:

```text
n=15 radius-1 moat: 120 missed cells
active coordinates: 6 or 14
active jumps:       5 or 10
```

The preferred positions and jumps are divisor-layer data, not random
coordinates.  This supports the idea that composite scalar moats are organized
by quotient coordinates.

On the disproof side, quotient ladders give very small positive gaps but do
not close endpoints.  The best fourteen-runner seven-ladder remains

```text
(1,7,14,21,28,35,49,56,63,70,77,84,91)
```

with forbidden length `142/143`, max-gap ratio `0.005411`, and `84` exposed
endpoints.  A one-swap neighborhood through speed `252` found no open-cover
repair; the best exact audits preserve the same gap geometry.  Endpoint
protector pressure is also split: speed `42` covers `48` of the seven-ladder
leaks but introduces `84` endpoints, while the best fifteen-runner protector
covers only `34` leaks and introduces at least `150` endpoints.

## Interpretation

HYP-1827 says scalar-near residue blockers are locally isolated.  HYP-1828
says speed-first disproof searches should solve endpoint-protection cycles
before exact interval audits.  S374 ties these together: an endpoint cycle
that uses quotient-ladder speeds may reduce max gap while creating a large
new boundary, and a residue-vector repair that stays near a scalar ramp may
open a fixed package of witness cells.

The missing bridge is a common potential function.  It should measure how a
speed-set endpoint protection move changes the induced residue-vector
distance from scalar ramps, or conversely how a scalar-distance perturbation
forces a finite endpoint leak layer.

## Next Tests

1. Build a finite endpoint graph on the first leak layer
   `{a/14, 9/98, 29/182, 15/182}` and search for directed protection cycles.
2. For every candidate endpoint cycle, compute a scalar-distance signature of
   the realizing speeds before running exact audits.
3. Search the `n=14` residue-vector space at scalar distance at least `3`,
   ranking candidates by induced endpoint closure rather than missed-cell
   count alone.
4. Prove the radius-1 scalar moat as a local lemma, then test whether the
   same proof isolates the endpoint protectors that introduce the most new
   boundary.
5. In the `n=15` analogue, classify why coordinates `6` and `14` with jumps
   `5` and `10` are preferred; use that divisor classification as a guide for
   fourteen-runner quotient layers.

## Sources

- `04-computation/lonely_runner_feedback_cycle_s374.py`
- `05-knowledge/results/lonely_runner_feedback_cycle_s374.out`
- `07-reflections/lonely-runner-feedback-cycle-s374.md`
- HYP-1827
- HYP-1828
- HYP-1831
