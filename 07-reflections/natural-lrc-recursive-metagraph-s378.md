---
source: codex-2026-05-31-S378
status: research synthesis
tags:
  - natural-operation-graphs
  - product-sum
  - lonely-runner
  - recursive-metagraph
  - scalar-gauge
---

# Natural Modes and the LRC Recursive Metagraph

The old natural-number mode graph has now become a useful test object for the
Lonely Runner work, but only after one correction.

The additive rule

```text
x -> z and y -> z when x+y=z
```

is not really an incomplete tournament once the second input is forgotten.  On
positive naturals, `x -> z` iff `x<z`, so the simple shadow is the complete
transitive order.  The multiplicative analogue keeps structure:

```text
x -> z iff x divides z.
```

That is the divisor DAG, with Hasse edges `x -> xp` for primes `p`.  Addition
is the dense completion.  Multiplication is the sparse skeleton.

The product-sum equations sit exactly at their interface.  If a nonunit core
`C` has

```text
P(C)=product(C),  S(C)=sum(C),  D(C)=P(C)-S(C),
```

then `D(C)` copies of `1` repair the additive fold:

```text
D(C)*1 + S(C) = P(C).
```

This is the natural-number toy model for the LRC obstruction: a dense equality
spine plus a sparse divisibility channel plus a defect that must be repaired
without creating a larger defect elsewhere.

## LRC Translation

For `k` moving speeds the LRC threshold denominator is `n=k+1`.

The initial segment is the Dirichlet equality spine.  In the micro-staircase
residue system, scalar ramps

```text
v_i = m i mod n
```

are the same kind of complete/equality object: they block every cell, but they
are a gauge orbit rather than a real counterexample.  After quotienting by this
scalar spine, the sparse structure reappears.

The S378 metrics make the prime/composite recursion visible:

```text
n=14: best puncture misses 56 cells, delta order 2, delta gcd 7
n=15: best puncture misses 120 cells, delta order 3, delta gcd 5
n=21: best puncture misses 224 cells, delta order 3, delta gcd 7
```

Prime denominators behave differently.  At `n=13,17,19`, the best deltas have
full cyclic order, and the preferred delta layer is broad rather than a proper
torsion subgroup.  This is the LRC version of the additive/multiplicative
split: primes give a unit skeleton; composites give quotient channels.

## Repair Defect

The one-step repair ledger gives a direct analogue of product-sum defect.

For a best scalar puncture, ask for the non-reverting coordinate change that
covers the most old missed cells.  It never repairs for free:

```text
n=14: gain 56 old misses, create 308 new exposed cells, ratio 11/2
n=15: gain 60 old misses, create 280 new exposed cells, ratio 14/3
n=20: gain 180 old misses, create 1220 new exposed cells, ratio 61/9
n=22: gain 220 old misses, create 1386 new exposed cells, ratio 63/10
```

This is the exposed-cell version of the product-sum slack law.  A unit in the
product-sum equation repairs multiplicative excess additively.  A local LRC
coordinate move repairs old exposure only by opening a new exposure package.
The proof target should measure that package, not merely whether old cells are
covered.

## Proposed Metagraph

The tournament work succeeded by replacing one scalar statistic with a
landscape.  LRC seems to need the same move.

For each denominator `n`, define a residue metagraph:

```text
nodes: normalized residue vectors in (Z/nZ)^(n-1) modulo scalar ramps
edges: one-coordinate residue repairs
height: missed-cell count in the micro-staircase system
torsion label: cyclic order/gcd of the changed residue
```

Then couple it to an endpoint-protection graph:

```text
nodes: speed sets
height: max-gap ratio, boundary count, unprotected endpoints, peel depth
operation labels: additive gates, multiplicative gates, divisor edges
```

The recursive object is not a single sequence in `n`.  It is the evolution of a
two-height landscape.  The transition `n -> n+1` changes:

- the micro-staircase arrangement through the denominators `n*i`;
- the unit skeleton size `phi(n)`;
- the divisor channel count `tau(n)`;
- the best torsion order for scalar punctures;
- the product-sum factor-packing type aligned with the runner count.

That explains why the LRC data do not move monotonically with `n`.  The state
is arithmetical, just as the tournament metagraph state was structural rather
than just "larger n means larger graph."

## New Angle

The natural-operation graph says: addition is complete after projection, so
the real data live in fibers and critical pairs.

The LRC analogue says: scalar ramps are complete after gauge projection, so the
real data live in quotient residues and endpoint critical pairs.

This suggests the next proof/search invariant:

```text
LRC state(n, v, D) =
  (missed_cells(v),
   repair_defect(v),
   endpoint_core(D),
   unit_skeleton(n),
   divisor_channel(n),
   operation_closure(D)).
```

That is the first explicit "higher-order concept" for LRC in the same spirit as
the repo's tournament angles of attack: not only whether a blocker exists, but
where it sits in a recursively changing moduli space.

## Artifacts

- `04-computation/lonely_runner_operation_metagraph_s378.py`
- `05-knowledge/results/lonely_runner_operation_metagraph_s378.out`
- `04-computation/natural_lrc_recursive_modes_s378.py`
- `05-knowledge/results/natural_lrc_recursive_modes_s378.out`
- HYP-1835
- HYP-1834 for the sibling operation-shadow feature-vector formulation
