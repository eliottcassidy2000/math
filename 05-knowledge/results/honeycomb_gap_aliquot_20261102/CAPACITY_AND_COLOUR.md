# Honeycomb capacity, information and unit distance colourings

**PROVED elementary finite information bound and scoped deductions;
CITED inherited geometry; FINITE-EXACT graph controls.**

## Correct capacities and what the faithful network preserves

The [root5 certificate](../collatz_hexagonal_growth_20261102/GEOMETRY.md)
has 704 distinct positive integers reaching5 within21 ordinary Collatz steps.
The honeycomb ball has 1+3*21*22/2=694 vertices. The deficit is **10, not20**.
At radius26 the triangular-lattice comparison is2263 versus2107, a deficit156.
These finite inverse graphs are trees; deleting the1,4,2 cycle does not help.
The new checker independently regenerates the inverse levels and both lattice
balls. No geometric placement is assumed in this counting obstruction.

The [faithful honeycomb theorem](../session140_collatz_hexagonal_geometry_20261102/collatz_hex/PROOF.md)
instead replaces each ordinary Collatz edge by a finite route, with unbounded
route lengths. It preserves every original vertex and transition after
contracting marked routes, including unknown components. It uses
[Kuehn's grid-minor theorem](https://web.mat.bham.ac.uk/D.Kuehn/unvacc3.pdf)
with its topological hypotheses, plus a square-grid minor of the honeycomb.
It is a subgraph model, not an induced or spanning model. Markers and route
clocks are essential: a source within21 ordinary steps need not lie within
21 honeycomb steps.

## A sharp finite sidecar entropy theorem

Let a uniformly distributed source S have N possible values. Suppose its
visible representation f(S) has at most M values, with 1<=M<=N. A sidecar
that reconstructs S must carry at least H(S|f(S)) bits of conditional entropy.
Write N=qM+r, 0<=r<M. Then

    H(S|f(S)) >= [(M-r)q log2(q)+r(q+1)log2(q+1)]/N.       (1)

This is sharp for unrestricted maps f. Its equality fibres have M-r sizes q
and r sizes q+1. Proof: for fibre sizes a_j,
H(S|f(S))=sum_j a_j log2(a_j)/N. Splitting a non-singleton into an unused bin
decreases this value, and the strict convexity of t log t shows that moving
one element from a fibre at least two larger than another decreases it.
Thus the minimum uses all M bins and balanced sizes. These moves also prove
the equality description. Exact reconstruction implies the sidecar bound
by H(S|f(S))=I(S;sidecar|f(S))<=H(sidecar|f(S)).

For N=704,M=694, the minimizing pattern is684 singleton fibres and10
doubletons. Equation (1) becomes

    H(S|f(S)) >= 20/704 = 5/176 bits.                      (2)

Here20 counts the source occurrences in those10 ambiguous pairs. This is
a legitimate secondary appearance of20, with an explicit counting map;
it is NOT the vertex deficit and gives no relation to aliquot20. A fixed-width
sidecar needs at least one bit in the worst case. The entropy bound concerns
uniform sources and does not assert that a nearest-neighbour map attaining
these fibre sizes exists; adjacency may force strictly more information.

The inherited all-size capacity theorem says that no fixed finite number of
sheets repairs exponential inverse growth in a fixed polynomial-growth lattice.
Equation (1) is a finite, distribution-specific refinement, not a new proof
of that all-size obstruction. Negative curvature changes the available
spatial capacity; the [gap construction](GAP_GEOMETRY.md) additionally retains
coarse LCA distance.

## Subdivision is not a colouring or metric equivalence

A hexagon C6 has a proper two-colouring. Mark alternate vertices and contract
the three intervening length-two routes: the recovered graph is C3, which
needs three colours. In the usual unit-edge hexagon those marked vertices
have Euclidean separation sqrt(3), not1. Thus faithful subdivision can
preserve a directed computation while changing both original edge lengths
and chromatic number. Triangle inequality only bounds endpoint distance by
route length; it never declares a contracted route to be a unit segment.

The [proper induced Euclidean realization theorem](../session140_collatz_hexagonal_geometry_20261102/hadwiger_nelson_refinement/PROPER_UNIT_DISTANCE_REFINEMENT.md)
for countable locally finite pseudoforests is a different repair. It preserves
exactly which pairs are at Euclidean distance1, but allows arbitrarily close
nonadjacent vertices and straight-edge crossings. Those concessions distinguish
it from a fixed lattice and from a planar subdivision drawing.

Hadwiger-Nelson concerns all Euclidean unit-distance pairs and their chromatic
number, not routing capacity. For example [de Grey's finite graphs](https://arxiv.org/abs/1804.02385)
already require five colours. A binary tree has chromatic number2; the ordinary
Collatz pseudoforest has chromatic number3 because of its known triangle.
Neither growth nor a honeycomb subdivision produces a larger chromatic bound.

## Small unit-distance counts and arithmetic guards

The primary [Alexeev-Mixon-Parshall table, version2](https://arxiv.org/html/2412.11914v2)
reports u(10)=20,u(14)=33,u(15)=37,u(20)=54; its Theorem1 gives
72<=u(25)<=78 and76<=u(26)<=84. Here u(n) is the maximum NUMBER OF EDGES
among planar unit-distance graphs on n points. These are the stated source's
results, not claims of newly resolved values at25 or26. The carrier33 counts
observer states and is a different invariant from u(14)=33.

For the square-sum graph on1..15, leaf stripping leaves exactly the cycle

    1--3--6--10--15--1.

Successive triangular numbers sum to squares; the last sum is16. The whole
graph has15 edges on15 vertices and is connected. It is therefore unicyclic,
and the inherited pseudoforest realization theorem supplies an induced
Euclidean unit-distance realization. Its chromatic number is exactly3.
This is an actual graph-preserving geometric transfer, not an equality of counts.

Exact three-colouring certificates also exist for the square-sum graphs at
20,24,25,26,33. They all contain this odd cycle, so their chromatic number is3.
This computation does NOT supply their Euclidean unit-distance embeddings.
G25 has32 edges and degree-four vertices, already obstructing any faithful
SUBDIVISION of it inside a degree-three honeycomb. Splitting such a vertex
into a cluster changes the operation and requires its ordered ports and
allowed traversals. A minor model alone does not preserve those path responses.

The first triangle in these interval square graphs appears at30, and the
first four-cycle at46; the latter is1--3--46--35--1. These are complete
finite checks through46. Thus G25 and G33 have no four-cycle. The rhombus
constraints used to initialize the cited unit-distance solver supply no
equations on either graph. This identifies a specific algorithmic boundary,
not a claim that the entire solver or geometric realization must fail.
