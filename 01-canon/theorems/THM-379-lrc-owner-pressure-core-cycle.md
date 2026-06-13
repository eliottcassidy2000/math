---
id: THM-379
name: lrc-owner-pressure-core-cycle
status: PROVED
date: 2026-06-01
session: codex-2026-06-01-S503
depends_on:
  - THM-357
  - THM-359
  - THM-365
---

# THM-379: Owner-realized endpoint cores force pressure cycles

## Statement

Let `(E,I,B,P)` be a finite endpoint/interval protection system as in
THM-359:

```text
B(i) subset E       boundary endpoints of interval i,
P(e) subset I       intervals strictly protecting endpoint e.
```

Let `(E',I')` be a nonempty protection core.  Suppose there is a finite owner
set `R` and owner maps

```text
o_E : E' -> R,
o_I : I' -> R.
```

Assume the core is owner-compatible:

1. every interval owner appears as an endpoint owner:

```text
o_I(I') subset o_E(E');
```

2. strict protection is never self-owned:

```text
i in P(e) cap I'  =>  o_I(i) != o_E(e).
```

Define the owner-protection digraph `G_C` on

```text
A = o_E(E')
```

by drawing an edge

```text
o_I(i) -> o_E(e)
```

whenever `i in P(e) cap I'`.

Then `G_C` contains a directed cycle.  In particular, any pressure graph on
`A` that contains all owner-protection edges from the core has a nontrivial
strongly connected component.

## Proof

Since `(E',I')` is a protection core, every endpoint `e in E'` has a strict
protector in the core:

```text
P(e) cap I' != empty.
```

Fix an active owner `r in A`.  By definition of `A`, choose an endpoint
`e_r in E'` with `o_E(e_r)=r`.  Since the core protects every endpoint, choose
some interval

```text
i_r in P(e_r) cap I'.
```

The owner graph therefore has the edge

```text
o_I(i_r) -> r.
```

By the no-self-protection assumption, `o_I(i_r) != r`.  Thus every vertex of
`A` has indegree at least one from another vertex of `A`.  The source
`o_I(i_r)` lies in `A` by owner-compatibility.

It remains only to recall the elementary finite graph fact: a finite directed
graph in which every vertex has indegree at least one contains a directed
cycle.  To prove it, start at any vertex and repeatedly choose an incoming
edge, so `x_{k+1} -> x_k`.  Since the graph is finite, some vertex repeats;
the repeated segment gives a directed cycle by following the displayed arrows
from the later repeat back to the earlier one.

Hence `G_C` contains a directed cycle.  A supergraph containing all of its
edges also contains that cycle, so it has a nontrivial strongly connected
component.

## LRC Specialization

For a Lonely Runner endpoint system, endpoints are owned by the speed whose
forbidden interval has that endpoint, and intervals are owned by their speed.
An interval cannot strictly protect its own endpoint: its own endpoints are
boundary points, not interior points.

Thus a terminal LRC endpoint core gives an owner-protection digraph whenever
we record strict incidences

```text
protector speed -> endpoint-owner speed.
```

If a chosen LRC pressure lift contains these owner-protection edges, then a
nonempty endpoint core forces a pressure cycle.  This is the formal bridge
between endpoint-core certificates and the pressure-SCC search program.

## Verification Record

`04-computation/lrc_endpoint_pressure_formal_s503.py` exhaustively checks the
finite graph lemma for loopless digraphs on up to `4` vertices and checks the
minimal protector-selector model through `8` active owners.  In every case,
the number of acyclic all-owned protector shadows is `0`.

Stored output:

```text
05-knowledge/results/lrc_endpoint_pressure_formal_s503.out
```

## Related

- THM-357: Lonely Runner endpoint-protection trichotomy.
- THM-359: endpoint/interval protection core peeling.
- THM-365: endpoint-protection cycle necessity.
- HYP-1960: LRC pressure searches return DAG certificates.
- HYP-1961: LRC pressure DAG peel layers are endpoint-certificate data.
