---
id: THM-365
name: lrc-endpoint-protection-cycle-necessity
status: PROVED
date: 2026-05-31
session: codex-2026-05-31-S384
depends_on:
  - THM-357
  - THM-359
results:
  - 05-knowledge/results/lonely_runner_endpoint_cycle_formal_s384.out
---

# THM-365: Endpoint-Protection Cycle Necessity

## Statement

Let `(E,I,B,P)` be a finite endpoint/interval protection system, where:

```text
B(i) subset E       boundary endpoints of interval i,
P(e) subset I       intervals strictly protecting endpoint e.
```

Assume each interval has at least one boundary endpoint.  In the LRC
specialization each interval has exactly two.

Define the endpoint-protection digraph `G` on `E` by drawing an arrow

```text
e -> f
```

when some interval `i in P(e)` has `f in B(i)`.

If `(E,I,B,P)` has a nonempty protection core in the sense of THM-359, then
the induced endpoint-protection digraph on the core contains a directed cycle.

Consequently, by THM-357 and THM-359, every reduced Lonely Runner
counterexample must contain a directed cycle in its labelled arithmetic
endpoint-protection graph.

For a Lonely Runner endpoint

```text
e = (n*m + eps)/(n*u),
eps in {-1,+1},
```

an arrow protected by speed `p` must carry the integer label condition

```text
| p*(n*m + eps) - a*n*u | < u
```

for some integer `a`.

Thus a true counterexample must realize a directed protection cycle whose
arrows satisfy these strict integer inequalities.

## Proof

Let `(E',I')` be a nonempty protection core.  By definition, every endpoint
`e in E'` is protected by at least one interval `i in I'`, so

```text
P(e) cap I' != empty.
```

Also by definition, every interval in `I'` has all of its boundary endpoints
inside `E'`:

```text
B(i) subset E'.
```

Therefore every vertex `e in E'` has at least one outgoing edge in the
endpoint-protection digraph induced by `E'`: choose `i in P(e) cap I'`; since
`B(i)` is nonempty and contained in `E'`, choose any `f in B(i)`.

A finite nonempty directed graph in which every vertex has positive outdegree
contains a directed cycle.  Starting at any vertex and repeatedly following an
outgoing edge, some vertex must repeat; the repeated segment is a directed
cycle.

This proves the abstract statement.

For the Lonely Runner specialization, THM-357 says that a counterexample is a
full-measure forbidden interval system in which every endpoint is protected.
THM-359 says the terminal peeling pair is the largest protection core.  In an
all-protected full system, the whole endpoint system is already a nonempty
core, so the cycle conclusion applies.

Finally, the displayed integer label is exactly the finite protection
criterion from THM-357 applied to owner speed `u` and protector speed `p`.

## Insufficiency of Bare Topology

The cycle condition is necessary but far from sufficient.  Bare circular arcs
have all-protected cores immediately.  On the cyclic quotient `Z/3Z`, the
three open arcs

```text
(0 -> 2), (1 -> 0), (2 -> 1)
```

of length `2` cover every cell, protect every endpoint, and form a nonempty
protection core.  The S384 computation finds such abstract all-protected
systems for `q=3,...,9`, including restricted short-arc variants.

Therefore a proof of Lonely Runner cannot use endpoint-cycle absence as a
pure topological fact.  It must use the arithmetic labels: speed orbits,
endpoint owner labels, strict integer protection inequalities, and the
full-measure balance of the generated intervals.

## Verification Record

`04-computation/lonely_runner_endpoint_cycle_formal_s384.py` verifies:

1. the abstract endpoint-cycle construction on finite cyclic arc systems;
2. the existence of all-protected abstract circular-arc mirages for
   `q=3,...,9`;
3. empty terminal cores, hence no terminal endpoint cycles, for the sampled
   LRC tight and near-disproof systems:
   initial `n=8`, sporadic tight `n=8`, initial `n=14`, the `n=14`
   seven-ladder, the `n=14` single-gate set, initial `n=15`, and the
   `n=15` `3x5` ladder.

The stored output is
`05-knowledge/results/lonely_runner_endpoint_cycle_formal_s384.out`.

## Use

This theorem sharpens the endpoint-protection program:

```text
counterexample
=> nonempty protection core
=> directed endpoint-protection cycle
=> labelled integer cycle with full-measure interval balance.
```

The next disproof search should enumerate labelled cycles first and then ask
whether any speed set realizes them.  The next proof search should find an
arithmetic potential that forces every such labelled cycle to leak an endpoint
or a positive gap.

## Related

- THM-357: Lonely Runner endpoint-protection trichotomy.
- THM-359: endpoint/interval protection core peeling.
- THM-360: unit endpoint divisibility filter.
- HYP-1811: LRC protection peeling.
- HYP-1828: protection-cycle-first disproof search.
- HYP-1836: endpoint incidence is fundamental.
- `04-computation/lonely_runner_endpoint_cycle_formal_s384.py`.
