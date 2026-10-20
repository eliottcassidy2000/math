# The Map of All Tournaments

**Corrected 2026-10-20.** This explainer distinguishes two different graphs:
the single-arc-flip graph and the blue/black complement-line multigraph.
Earlier numerical counts, the mixed-triangle prohibition, and a suggested
global greedy-optimization guarantee were incorrect or unsupported.

A tournament has one directed edge between each pair of distinct vertices.
For a fixed size n, identify tournaments that differ only by relabelling.
The number of resulting classes is A000568(n). For n=3,4,5,6,7,8 the counts
are respectively 2,4,12,56,456,6880.

## The single-arc-flip graph

Each vertex of G_n is one tournament isomorphism class. Two distinct
classes are adjacent when reversing one arc of a representative can move
between them. Reversals that preserve the class may be retained as loops
when the chosen convention requires them.

This graph is connected: the labelled orientation cube is connected by
single-arc reversals, and its projection onto isomorphism classes preserves
walks. Consequently a search starting at the transitive tournament and
canonicalizing all one-arc reversals enumerates every class. This avoids
Burnside's formula, but still visits every class; it is not an efficient
scalar counting formula at large n.

The Hamiltonian-path count H(T) is a useful vertex statistic. A graph
picture or finite experiment does not prove that greedy H-increasing moves
reach a global optimum. No such guarantee is asserted here.

## Converse merging

Reversing every arc sends T to its converse T^op. Some classes are
self-converse; all other classes occur in pairs. If a_n is the ordinary
class count, s_n the self-converse count and m_n the merged-node count,
then

    m_n=(a_n+s_n)/2,       a_n=2m_n-s_n.

Equivalently, add the number of paired merged nodes to m_n. At n=7,
a_n=456, s_n=88 and m_n=272. At n=8 the corresponding counts are
6880,176 and3528. Self-converse counts are A002785.

Calling SC--SC edges the spine, SC--NS edges the ribs, and NS--NS edges
the sea is a classification of edge types; those names do not establish
connectivity or other geometry. The rib-only graph is bipartite because
its two parts are SC and NS. Every closed walk crosses that cut an even
number of times. A mixed triangle with two ribs IS possible, already at
n=4: take the transitive class, a source above a directed triangle, and
the strongly connected four-vertex class. Each pair is one-arc adjacent.

## The blue/black line multigraph is a different object

Fix a directed Hamiltonian path. Its n-1 arcs are prescribed; the remaining
binom(n-1,2) bits form a tiling cube X_n. Thus the TOTAL number of tilings
is 2^binom(n-1,2). Over an ordinary class T the fibre has size
H(T)/|Aut(T)|, which varies with T.

There are two commuting operations on tilings:

- Reflection of the triangular grid represents tournament converse after
  relabelling. A reflection-fixed tiling is blue.
- Complementing every FREE bit preserves the prescribed path. Its two
  tilings form one complement line, blue when its endpoints are blue and
  black otherwise. This is generally neither one arc reversal nor full
  tournament converse.

Map both ends of every line to their converse-merged tournament nodes.
Keep loops and parallel lines: their multiplicities contain fibre data.
This produces the blue/black line multigraph. Complementing a tiling does
not determine a unique partner of its node; different path presentations
inside one node can lead to different partner nodes.

## What changes with size

The faithful recursive object retains tilings, lines and nodes together.
A tiling is reconstructed from its two endpoint-deleted faces, their common
interior, and one apex bit. Bare face lines additionally lose a coherent
endpoint phase. Deletion gives a weighted correspondence between nodes,
not a unique parent-node function.

These are proved all-size statements in
[THM-796, recursive incidence](../../01-canon/theorems/THM-796-three-sorted-recursive-tiling-line-node-incidence.md).
The fibre dictionary and class search are described in
[THM-1430](../../01-canon/theorems/THM-1430-the-tiling-class-metagraph-dictionary-and-which-tricks-pay.md).
The [2026-10-20 growth study](../../05-knowledge/results/metagraph_growth_20261020/README.md)
adds checked multiplicative constructions and explains the remaining
obstacle to a compressed recurrence for A000568.
