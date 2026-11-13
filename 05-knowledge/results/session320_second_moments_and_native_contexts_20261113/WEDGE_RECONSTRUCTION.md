# First moments lose chirality; wedge co-usage recovers it

An undirected wedge consists of two different labelled pairs {u,v},{v,w}.
Its co-usage bit says whether both pairs occur consecutively somewhere
in the same spanning Hamilton path. Retain also the support of single
edge usage (necessary as a convention for order two). The observer
forgets the order of the path and the directions of all arcs.

## A one-slot reconstruction theorem

For any tournament this observer determines the tournament up to global
converse. This uses two older results: [legal single-vertex relocation
connectivity](../hamilton_path_relocation_connectivity_20261020/THEOREM.txt)
and [reconstruction from directed positive support](../session70_attachment_guard_alignment_20261026/session70-guards/REPORT.md),
the latter using the classical [Busch–Jacobson–Reid theorem](https://doi.org/10.1002/jgt.20173).

Make a graph on the positive edges; join them when they are adjacent in
a positive wedge. Every Hamilton path's edges lie in a connected chain
in this graph. Two Hamilton paths related by one legal relocation share
an edge. For n>=5, the common word after removing the moved vertex has
n−2 pairs; at most two of these are split by the two insertion positions,
leaving at least n−4 common pairs. At n=4, the only way to split both
pairs would be inserting the moved vertex into two adjacent interior
gaps, impossible since their common vertex would have to beat and lose
to the moved vertex. At n=3, the only internal gap cannot be used
simultaneously with an end gap: either end would contradict the
orientation to its adjacent vertex. If the paths are different and
both use end gaps, the common edge survives. Orders one and two are
trivial. The common word in distinct relocations is itself a directed
path: if an edge in it failed, it would have to be split by both
insertion positions, contradicting distinctness.

It follows from relocation connectivity that the wedge graph on all
positive edges is connected. At a positive wedge, exactly one of its
two edges points into the common vertex: their relative directions
are fixed. Choose an orientation of one edge and propagate along the
connected graph. The data come from a tournament, so the constraints
are consistent; the only two possible positive directed supports are
mutual reversals. Apply the inherited directed-support decoder to
recover the dark arcs. This proves injectivity up to converse even
from the **positivity bits**, not the co-usage counts.

## A two-slot version with a fixed optimal query list

If two optional slots are available, the wedge bit is positive exactly
when the local wedge is a directed 3-vertex path. Sufficiency follows by
putting those three vertices in one slot and any Hamilton path on the
remaining induced tournament in the other slot, allowing empty.
Necessity is local because two edges sharing a vertex cannot belong
to distinct vertex-disjoint slots. This version does not need any
relocation or dark-arc theorem.

Let M=binom(n,2). Choose in advance a spanning tree of the line graph
of the complete graph. Query its M−1 wedges (for n>=2). Each answer,
even the negative ones, determines the relative orientation of the two
edges: a negative bit says both point in or both out at their common
vertex. Integrating along the tree recovers all arcs up to one global
flip. This number of yes/no queries is optimal in worst case for
labelled tournaments up to converse, because there are 2^(M−1) such
labelled states. This count is not a query lower bound for unlabelled
isomorphism classes, or for queries returning full integers.

The contrast with the [strong first-moment collision](STRONG_SHADOW.md)
is exact: a tower of unary slot moments is insufficient, but an
appropriately selected binary relation restores the missing signs.
