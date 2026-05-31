# Lonely Runner Feedback Cycle S374

The session deliberately cycled through three routes:

```text
14-runner attack -> 15-runner transfer -> disproof construction
```

At every dead end, the next route had to inherit a new idea rather than start
fresh.  That constraint was useful.  It prevented the usual split where the
cell-blocking work and the interval-cover work evolve as two separate
languages.

## Cycle Notes

The fourteen-runner route began with the S371/S372 best non-scalar near
blocker

```text
(8,2,10,4,12,13,0,8,2,10,4,12,6).
```

Its anatomy is exactly scalar ramp `m=8` plus one defect at coordinate `6`.
The mutation table was blunt: the only complete repair is to undo the defect.
The best true second-defect attempt misses `126` cells.  So the local proof
route has a solid dead end: scalar-near repair gives the scalar equality
spine back, not a non-scalar full blocker.

The forced fifteen-runner transfer made the moat thicker instead of thinner.
The best scalar-puncture misses `120` cells, and the zero-ramp witnesses are
concentrated at coordinates `6` and `14` with jumps `5` and `10`.  This looks
like divisor-layer behavior: the composite quotient is choosing its leak
coordinates.

The disproof route then tried quotient ladders.  The fourteen-runner
seven-ladder

```text
(1,7,14,21,28,35,49,56,63,70,77,84,91)
```

is still the best tiny-gap near-miss: forbidden length `142/143`, gap ratio
`0.005411`, and first exposed endpoint `9/98`.  But it has `84` exposed
endpoints.  A one-swap repair neighborhood did not find an open cover; the
best swaps mostly reproduce the same geometry.  Endpoint-protector pressure
also split badly: adding speed `42` covers `48` leaks but introduces `84`
new endpoint values.

## Synthesis

The new invariant is not "be far from scalar ramps" and not "close endpoint
cycles" in isolation.  It is the interaction:

```text
scalar distance versus endpoint-closure growth.
```

Scalar-near residue vectors hit the `56`/`120` witness-cell moats.  Speed-set
near-disproofs can make the max gap tiny, but their endpoint closure grows
faster than single protectors can close it.  The missing proof bridge should
show that pushing down one pressure necessarily raises the other.

For a disproof search, this suggests a stricter pipeline:

1. choose a small leak layer, beginning with `a/14`, `9/98`, `29/182`, and
   `15/182`;
2. find a directed endpoint-protection cycle on that finite layer;
3. compute scalar-distance signatures of the realizing speeds;
4. run exact interval audits only after both tests are favorable.

For a proof search, the same data points toward a two-layer lemma: near the
Dirichlet scalar spine, finite cell witnesses are unavoidable; away from it,
endpoint-protection closure fails to stabilize.

## Artifacts

- `04-computation/lonely_runner_feedback_cycle_s374.py`
- `05-knowledge/results/lonely_runner_feedback_cycle_s374.out`
- `05-knowledge/hypotheses/HYP-1829-lrc-scalar-distance-endpoint-closure.md`
