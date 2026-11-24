# The two Moore-spine tournament conventions have different exact triangle counts

Take a simple d-regular graph G with a specified spanning path, and
label its vertices 1,...,n in that order. All consecutive arcs point
forward. The sparse convention points an off-spine arc backwards
exactly on graph edges. The dense convention points it backwards
exactly on graph nonedges. They are complements on the free tiles,
not converses of each other. All statements on Moore(57,2) are
conditional on the graph and the specified spanning path.

## Triangle counts and girth

For any undirected graph F on this order, reverse the edges of F
from the all-forward tournament. Put
b_i=d_F,left(i)−d_F,right(i). Then the exact cyclic triangle count is

    C3 = sum_(ij in F, i<j)(j−i) − (1/2) sum_i b_i^2.

Indeed the new outdegrees are n−i+b_i; the number of transitive
triangles is the sum of their binomial coefficients choose(outdegree,2).
The b_i sum to zero, and each edge ij contributes i−j to
sum_i (n−i)b_i. Subtract from the all-forward total.

Take F=G for the auxiliary all-edge reversal. The dense convention
is exactly its converse, so has the same cyclic triangle count.
The sparse convention instead takes F=G minus the consecutive
spine edges. This reduces the arrangement sum by n−1; it changes
no internal b_i, because one left and one right edge are removed
at an internal vertex. The two endpoint values change from ±d
to ±(d−1), reducing sum b_i² by 4d−2. Consequently

    C3(dense) − C3(sparse) = n−2d;
    at Moore size n=d²+1 this equals (d−1)^2.

In particular the exact differences for the orders 10, 50, 3250
are 4, 36, 3136. The sparse number of backward chords is
nd/2−(n−1), respectively 6, 126, 89376. Each such chord, together
with the forward spine segment, forms a directed cycle. The
dense convention also has backwards chords (every graph nonedge).
Every cyclic tournament has a directed triangle, by taking a
shortest directed cycle and using a chord if its length exceeds
three. Thus both conventions have directed girth exactly three;
the girth-five property only governs the undirected carrier.

For degree 57, the linear arrangement sum for G is the sum of cut
sizes along the order. By the Moore Laplacian spectrum, cut size
for a set of k vertices lies between 50k(n−k)/n and 65k(n−k)/n.
The arrangement sum thus lies between 50(n²−1)/6 and
65(n²−1)/6. Since |b_i|≤57 for the all-edge reversal, sparse
C3 is at least 82,738,064 and dense C3 is at least 82,741,200.

## A transitivity bound from inversion forests

For an induced transitive subtournament on a chosen set S,
the graph of backwards arcs in the spine order is the inversion
graph between that order and the induced transitive order.

**Lemma.** A triangle- and quadrilateral-free inversion graph
is a forest. Orient each inversion from earlier to later in
the first order. This orientation is transitive (two decreasing
comparisons in the second order imply the third inversion).
Triangle-freeness therefore implies that every nonisolated
vertex is a pure source or a pure sink. In each of these classes
the two vertex orders agree. Order the sinks by either order:
each source has an interval of sink neighbors, starting at the
first sink that is later in the first order and ending at the last
sink that is earlier in the second. Both interval endpoints
are nondecreasing along sources in their common order.
If a shortest cycle existed, it would be induced and of length
at least six. Take its earliest source and its sink neighbors
b<c on the cycle. The other cycle source adjacent to b is later;
its interval starts no later than b (as it contains b), and
ends no earlier than c by monotonicity. Hence it is adjacent to
c, a chord unless the cycle was a quadrilateral. Contradiction.

In dense convention, one can pass to the auxiliary converse,
which has exactly G as its inversion graph. Hence a transitive
set of size s satisfies e_G(S)≤s−1.

In sparse convention the inversion graph H is G[S] minus the
selected consecutive spine edges, so is again a forest. Split S
into r consecutive-spine runs, of which r_1 are singletons.
Along each run both orders increase since all its spine arcs
are retained. An inversion between an earlier run and a later
run must use, respectively, the last and first vertex of the
two runs. Otherwise the next vertex of the earlier run or
the previous vertex of the later run forms a second inversion
and hence, with the within-run spine edge, a triangle in G.
There are at most 2r−r_1 run endpoints available to the forest.
The other graph edges are the s−r spine edges, giving

    e_G(S) ≤ s−r + (2r−r_1−1)
            ≤ s+floor(s/2)−1.

For a Moore(57,2) graph the induced edge lower bound is

    2e_G(S) ≥ 65s²/3250 − 8s,

by decomposing the characteristic vector into its constant and
zero-sum parts and using least adjacency eigenvalue −8.
Comparison gives s≤499 in the dense case and s≤549 in the
sparse case.

## A lower bound on Hamiltonian paths

Greedily remove vertex-disjoint directed triangles until none
remains. The remaining subtournament is transitive, so the
preceding bounds give at least 917 removed triangles in the
dense convention and 901 in the sparse convention.

Hamilton paths can be merged two at a time, preserving both
of their relative orders. At each step choose the head that
beats the other head. The previously chosen vertex reaches
both available heads: upon picking a head, it beats the other
head by choice and beats its own successor in its directed path.
This is a maintained invariant, with the one-head case trivial.
Iterate this deterministic merge on any number of vertex-disjoint
paths. The final path retains each original path order, so the
choice of original constituent paths is recovered from its
restrictions. Each removed directed triangle has three choices
and the remainder has a directed Hamilton path. Hence

    H(sparse) ≥ 3^901,       H(dense) ≥ 3^917.

These are structural lower bounds on these particular encodings,
not a general constraint on arbitrary tournaments.
