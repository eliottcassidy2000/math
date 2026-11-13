# The two requested Collatz graphs

**PROVED elementary reductions / FINITE-EXACT experiments; universal
reachability remains OPEN.** During this session, incoming commit
798965f7a5 supplied an independently developed
[guarded Collatz note](../fill_tree_matrix_20261113/collatz/NOTE.md), including
the same literal-graph reduction and forced-excursion family. We credit
that packet for the shared proof and do not claim these as a second new
discovery. Our additional receipt computes exact minimax heights in a
fixed box and checks an explicit greedy strategy.

## Ordinary edges in both directions

Let C(n)=n/2 for even n and3n+1 for odd n. Allowing each edge to be
traversed in either direction does not weaken root reachability: n is
connected to1 iff its ordinary forward orbit reaches1. Along every edge
the forward orbits coalesce, and coalescence is transitive.

There is a stronger finite-height version. For every cap B>=1, n has a
bidirectional path to1 entirely inside[1,B] iff its ordinary forward path
to1 stays inside[1,B]. The induced graph is a partial functional graph.
Every weak component contains at most one terminal cycle or boundary
sink; the component of1 has the ordinary root cycle when B>=4, and has
the terminal vertex1 for B<4. A vertex connected to1 cannot instead
escape through a different boundary sink. Consequently the minimum
possible maximum visited integer is identical in the two models.

The native inverse predecessors are2n, and(n−1)/3 only for n=4 mod6.
Deleting1→4 gives a rooted inverse tree on the basin of1. It spans all
positive integers exactly if Collatz holds. On odd vertices, accelerated
preimages form infinite fibres (2^a n−1)/3 with parity/congruence guards;
binary branching applies to the native step clock, not that acceleration.
These are standard guarded formulations; see
[Lagarias's overview](https://arxiv.org/html/2111.02635v1) and the incoming note.

## The literal parity-restricted extra moves

Retain C and add even n→2n, and odd n→(n−1)/3 when n>1 and n=1 mod3.
The latter is new, e.g.7→2; it is not the inverse of an ordinary odd
Collatz edge. Added even doubling can be deleted from any successful
path without increasing its peak: between odd vertices it only moves
up and down one dyadic chain, whose odd exit is unchanged.

Thus the exact odd-vertex reachability model has the two choices

    oddpart(3u+1), with native peak3u+1;
    oddpart((u−1)/3), for u>1 and u=1 mod3, with native peak u.

The parity-restricted graph differs from the older repository E-graph,
which adds even FORWARD3n+1 edges. Its Q1/Q2 claims cannot be transferred
here by changing arrow labels. Nor may a full dyadic ray be contracted
as a strongly connected component: its odd endpoint does not have an
allowed doubling edge.

For even K>=2, every successful literal path from2^K−1 must reach at
least2·3^K−2. The forced odd states are
3^j2^(K−j)−1, j=0,...,K−1. The first is0 mod3 and every subsequent state
is2 mod3, so the new inverse move is unavailable throughout. The peak
ratio tends to infinity. This rules out a universal constant-factor
height buffer, not eventual reachability. The numerical obstruction at
1023 has exactly the predicted peak118096; the next forced family member
4095 already needs1062880, above the one-million box.

## Exact experiments

`collatz_graphs.py` runs minimax Dijkstra on reversed edges of each
induced graph on1..1,000,000. It audits every reverse guard against the
forward definition on1..300 and stores named simple path witnesses.
For starting integers1..100,000:

| Model | Reaches1 inside the box | Not connected inside this box |
|---|---:|---:|
| ordinary forward |92,097|7,903|
| every ordinary edge bidirectional |92,097|7,903|
| literal extra moves |98,545|1,455|

The complete minimax arrays for the first two models agree at every
vertex in the box, as proved above. Missing box connections are not
divergent orbits. In the literal graph the first box misses are4095,
6143,8190,9215,10239.

A separate explicit strategy chooses the new inverse move whenever it
is available at an odd vertex, otherwise the ordinary move. It reaches1
for EVERY starting integer1..100,000 in the exact calculation. Its
memoized certificates contain161,649 visited states with strictly
decreasing remaining path length. The largest observed peak is86093440
from65535; the longest path has104 edges, from81407. These paths can
leave the one-million box, explaining the different totals. No claim
of termination for every positive integer follows from this finite test.

## What a spanning-tree analogy can transfer

For a countable digraph, universal root reachability is equivalent to a
natural-number rank with an available decreasing outgoing edge at every
nonroot vertex, and hence to a spanning in-arborescence. Such a rank
exists on the known basin; defining it everywhere is the substantive
proof obligation. A local binary branching rule does not supply it.

For tournament constructions the analogous useful object is a tree of
REALIZABLE endpoint states with a decreasing inverse target. The
selector in LANGUAGE supplies exactly that tree for its sparse language.
It does not establish that the tree's projection covers every allowed
odd count. This is the precise shared obstruction: well-defined local
transitions versus surjectivity onto the intended integer set.
