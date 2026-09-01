# Endpoint-590 critical 41-obligation response core

**FINITE-EXACT in the fixed labelled rank-8/rank-9 response hypergraph
inherited from the endpoint-590 packet. This is not a physical lonely-runner
entry and does not prove LRC(14).**

## Result

Let `F` be the frozen ordered family of 100 endpoint-590 failure bodies and
let `R` be the 14,368 distinct nonempty response signatures realized by all
active rank-8/rank-9 masks. For a subset `X` of failure ordinals, write
`tau(X)` for the smallest number of members of `R` whose union contains `X`.

The 41-ordinal set

```text
U = {0,5,9,12,14,18,20,22,24,25,29,32,35,37,38,40,43,47,48,
     53,55,57,63,65,68,71,73,76,77,79,80,83,84,88,89,90,94,
     95,96,97,99}
```

is exactly 9-cover-critical:

```text
tau(U) = 9,
tau(U \ {u}) = 8 for every u in U.                       FINITE-EXACT
```

Thus the last unit in the endpoint-590 response optimum is already forced by
41 of the 100 failures, and none of those 41 can be deleted while preserving
that obstruction.

## Exact proof packet

Restriction sends a full response signature `S` to `S intersect U`. This
produces 2,041 distinct nonempty traces. Quotienting those traces by inclusion
leaves 395 maximal traces. This quotient is lossless for cover existence:
replacing a trace by a containing trace never increases the number of chosen
responses and cannot uncover a point of `U`.

The hardened C++ search branches on all maximal traces containing a selected
uncovered pivot. State-local dominated gains are removed only when a containing
gain can replace them. Its sum and integer-dual prunes are necessary
conditions. Independent `-O2` and `-O3` builds agree byte-for-byte:

```text
depth-eight search on U:
  nodes 1,968,373; dead states 1,698,303;
  sum prunes 530,603; dual prunes 1,009,466; result UNSAT.

depth-seven searches on every U \ {u}, aggregated over all 41 u:
  nodes 278,730; dead states 270,531;
  sum prunes 52,383; dual prunes 198,392; 41 UNSAT, 0 SAT.
```

A direct nine-response ledger covers `U`. Separately, the 328-row criticality
ledger gives eight atlas representatives for each of the 41 single deletions;
each union is checked to equal exactly `U \ {u}` on the core coordinates.
These upper witnesses and the exhaustive lower searches prove both displayed
equalities. The MILP used to discover witnesses is not a proof dependency.

## Why eight is almost possible

The existing integer dual has a 21-point support `D` inside `U`, total weight
22, and maximum load 3 on every one of the 14,368 responses. The other 20
points of `U` are essential satellites. Among the 395 maximal core traces, the
dual-load distribution is

```text
load 0: 13,  load 1: 88,  load 2: 177,  load 3: 117.
```

There is an exact saturation identity behind the remaining obstruction.
Suppose eight responses covered `U`. Let

```text
l_j     = dual load of response j,
Delta   = sum_j (3 - l_j),
c_i     = number of chosen responses covering dual point i,
Omega   = sum_{i in D} w_i (c_i - 1).
```

Because every dual point is covered, `c_i >= 1`, and therefore

```text
sum_j l_j = sum_{i in D} w_i c_i = 22 + Omega.
```

But eight responses have total dual capacity 24, so

```text
Delta + Omega = 24 - 22 = 2.                              PROVED
```

Both terms are nonnegative integers. Hence any hypothetical eight-cover has
at most two units of unused dual capacity plus repeated dual coverage. In
particular at least six of its eight responses must have dual load three, and
the weighted overlap on the dual support is at most two. The 20 satellites
collectively rule out every such near-partition. This explains why the linear
dual proves only eight while the integer cover number is nine.

## Pairwise obstruction is far too weak

Join two failures when no realized response covers them together. This is an
intrinsic incompatibility graph, not a tournament: many pairs are compatible,
and forcing orientations would discard the relevant ties.

An exact Bron-Kerbosch census gives

```text
full 100-point graph: 740 edges, clique number 4, 50 maximum cliques;
core 41-point graph:  175 edges, clique number 4,  7 maximum cliques.
```

A frozen proper four-coloring of the core graph, checked edge by edge, combines
with its four-clique to give `chi(core)=4` exactly. Since graph coloring is the
optimal cover of vertices by pairwise-compatible classes, the entire core
pairwise-incompatibility abstraction sees only four, while the response
hypergraph cover number is nine. The obstruction is genuinely higher-order.
For the full 100-point graph only the clique census is asserted here; its
chromatic number is not inferred from clique number four.

## Geometry visible in the core

Every core body avoids labels `1,2,4,5,8`; label 26 occurs in 40 of 41 core
bodies, label 20 in 36, and label 12 in 31. These near-core coordinates may be
useful for classifying the 117 dual-saturated maximal traces. They are not by
themselves symmetries of the active response universe: activity remembers the
labelled endpoint geometry, which the failure-body projection discards.

## Map audit and limits

- **Source:** the complete 100-point, 14,368-signature response hypergraph.
- **Target:** the 41-point trace hypergraph, then its 395-member inclusion
  quotient.
- **Preserved predicate:** existence or nonexistence of a `k`-response cover
  of `U`.
- **Destroyed information:** behavior on the other 59 failures, response-mask
  multiplicity and rank, labelled margin geometry, and all simultaneous
  deletion safety.
- **Needed sidecars:** the independently verified complete response atlas and
  the endpoint-590 arithmetic replay that produced it.
- **Cheapest decisive tests:** the exact no-eight/no-seven searches and the
  direct criticality ledgers frozen here.

The packet establishes no exchange, no physical entry, no terminating
descent, and no proof of LRC(14).

## Next mathematical pulls

1. Classify the 117 load-three traces by their dual footprints, then seek a
   short finite family of satellite inequalities replacing the depth search.
2. Determine whether the `Delta + Omega = 2` regime admits a human-sized case
   split based on the near-universal labels 26 and 20.
3. Search earlier endpoint response hypergraphs for similarly critical cores;
   persistence of a common core type would be evidence for a reusable
   endpoint-lowering lemma rather than an isolated finite certificate.
