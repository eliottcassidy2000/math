# A literal bridge: common blocks measure successor loss

**PROVED elementary corollary.** This connects the two graph mechanisms,
not their arithmetic or geometric realizability. No universal1–2–4 law is
claimed; square-sum graphs are not identified with Euclidean unit-distance
graphs.

Let S be a nonempty compatible family of directed Hamilton paths on the
same n>=1 vertices. Let D be its directed arc union, and G the underlying
simple undirected union. Compatibility means no pair has both orientations,
so G and D have the same edge count m. The union contains a spanning path,
hence is connected. Let b be its number of maximal common directed blocks
and J the set of distinct ordered block joins. The common-block theorem gives

    m=n-b+|J|.

Consequently its cycle rank and its minimum edge-deletion loss when arbitrary
reorientation is allowed are

    beta(G)=|J|-b+1,
    epsilon(G)=max(0,|J|-b).

**Proof.** Connectedness gives beta=m-n+1. If G is a tree then m=n-1 and
its pseudoforest deletion defect is zero. Otherwise it has one cyclic
component and no tree component, so the exact m-n+t formula gives epsilon
=m-n. Substitute the common-block count in both cases. This also proves
that |J|>=b-1. The statements include a one-path family (b=1,J empty).

When contracting common blocks, retain the quotient as an undirected
multigraph. Opposite ordered joins between nontrivial blocks can use
different original endpoint pairs and become **parallel edges**, not one
simple quotient edge. Deleting that multiplicity would change cycle rank
and the functional deletion defect. Reversed joins between singleton blocks
are excluded by compatibility, as the common-block criterion requires.

## Two exact application witnesses

| Compatible family | n | b | Number of joins | m | beta | epsilon |
|---|---:|---:|---:|---:|---:|---:|
| An order25 compatible pair |25|3|4|26|2|1|
| The displayed order26 triple |26|5|8|29|4|3|

The first row uses the inherited UAB/UBA pair. The second uses the new
five-block three-path carrier. Their nonzero defects make an exact statement:
no single partial successor system can retain all union edges, and precisely
one or three edge deletions, respectively, suffice if reorientation is allowed.
This is a graph-capacity transfer, not a claim that the original directed
paths have one common successor choice or that these graphs are unit-distance
extremizers.

## Preserving the arrows changes the deletion budget

Let z be the number of vertices of D with outdegree zero. If retained edges
must keep their current directions, the exact minimum deletion budget is

    epsilon_arrow(D)=m-n+z.

Indeed, at each vertex v one may retain at most one of its outgoing arcs.
Keeping one whenever its outdegree is positive retains exactly n-z arcs;
these choices are independent because each arc has exactly one source.
Equivalently, sum_v max(0,outdegree(v)-1)=m-(n-z). A Hamilton path already
has outgoing arcs at every vertex except its terminal, so z<=1 for its union;
z=1 exactly when every member of S has the same terminal.

This distinction survives the concrete examples:

| Family orientation | z | Arbitrary-reorientation loss | Fixed-arrow loss |
|---|---:|---:|---:|
| Displayed n25 pair starting at18 |0|1|1|
| Its reversed pair, both ending at18 |1|1|2|
| Displayed n26 triple starting at18 |0|3|3|
| Its reversed triple, all ending at18 |1|3|4|

Thus reversal preserves compatibility, union size, polynomial weights and
undirected cycle/defect data, but need not preserve the fixed-successor
compression cost. The required sidecar is terminal/sink information. In
particular, saying merely that the n25/n26 budgets are1/3 without specifying
reorientation or the displayed starting18 orientation would discard a real
one-edge distinction.
