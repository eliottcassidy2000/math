# The hard-core conflict graph: induced universality, full-graph rigidity

For a tournament T, let Ω(T) have a vertex for each directed odd cycle
of T, and join two whenever their vertex sets intersect. An independent
set in Ω(T) is a pack of disjoint odd cycles; the independence polynomial
A_T(z) of Ω(T) has A_T(2)=H(T) by the
[Grinberg–Stanley formula](https://arxiv.org/abs/2307.05569).
The hard-core analogy itself is inherited (see
[the earlier conflict-graph note](../../../07-reflections/real-rootedness-of-independence-polynomial.md)).
Two all-order statements delimit the analogy.

## Every graph is a selected conflict graph

Every finite simple graph G occurs as an induced subgraph of Ω(T) for
some strong tournament T of order at most 3|V(G)|+|E(G)|+2, for
nonempty G. To construct it, give every edge of G one shared vertex.
For each v in G, use the d(v) shared vertices corresponding to incident
edges and enough private vertices to make an odd cycle of length at
least three. Choose the least odd length at least max(3,d(v)+1), so
there is a private vertex even if all incident vertices could fit
without one. This length is at most d(v)+3. Any two prescribed cycles
meet in exactly one vertex if their G-vertices are adjacent, and in
none otherwise. Distinct prescribed cycles therefore do not prescribe
orientations on the same unordered pair. Complete all unassigned
orientations arbitrarily, then add two singleton portals making a
directed triangle on the two portals and the whole old tournament as
one block. This makes the host strong and leaves the selected cycles
unchanged. The sum of selected cycle lengths is at most
3|V(G)|+2|E(G)|, but each shared vertex was counted twice, giving the
claimed host bound. The empty graph is induced on the empty selection
in any host, e.g. the one-vertex tournament.

The least host order is exactly 3d for G=K_{1,d}, d≥2, even among
strong hosts. The d leaf cycles must be disjoint and use at least
three vertices each. Conversely, make d disjoint directed triangles.
A central directed cycle can use one vertex from each triangle if d
is odd; if d is even, use one additional vertex in one triangle as a
consecutive pair oriented consistently with its triangle. Complete
other arcs arbitrarily. The union of the specified cycles is already
spanning and strongly connected, and their conflicts form exactly
the prescribed star. For d=1 the minimum is four (two intersecting
triangles); for d=0 it is three (one triangle).

## The full conflict graph cannot be arbitrary

Every nonempty connected component of Ω(T) has a universal vertex
(one adjacent to every other vertex in the component), and hence has
diameter at most two. In fact these components correspond exactly to
the nonsingleton SCCs of T.

Any directed cycle is inside one SCC, so cycles in different SCCs
are disjoint. Conversely, a strong tournament on n≥3 vertices has
cycles of every length 3 through n. An elementary proof starts
with a directed triangle, which exists because a triangle-free
tournament is transitive. Given a cycle not yet spanning: if an
outside vertex has mixed adjacency to the cycle, insert it at a
transition. Otherwise the outside is partitioned into W, beating
the entire cycle, and L, losing to the entire cycle. Both are
nonempty by strongness, and some l→w crosses from L to W (if not,
no arc leaves L). Replace a consecutive segment u→v→t of the
cycle by u→l→w→t to increase length by one. Consequently there
is an odd cycle of length n or n−1. Every other odd cycle has
at least three vertices and therefore intersects it. This proves
universality within each nonsingleton SCC.

Thus the selected conflicts can realize an arbitrary hard-core
model if unselected cycles are explicitly assigned activity zero.
They do not realize arbitrary G as the entire conflict graph:
for example, a four-vertex path is connected with no universal
vertex, so cannot be the whole Ω(T) or one of its components.
Nor does the induced embedding set the uncontrolled cycles'
activities to zero in the unweighted Hamiltonian evaluation.
It gives no scalar inverse for Chan–Pak.
