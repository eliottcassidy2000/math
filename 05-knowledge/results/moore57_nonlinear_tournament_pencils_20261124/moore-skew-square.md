# Moore edge orientations and their signed complement squares

Let G be a Moore graph of degree d and diameter two, and let T be
an arbitrary tournament on the same named vertices. Write K for
the skew matrix obtained by restricting the tournament signs to the
edges of G (zero elsewhere). Conditional results at degree 57 do
not establish the existence of G.

## Exact signed-complement encoding

The matrix Q=KKᵗ−dI is supported precisely on graph nonedges,
and for distinct nonadjacent u,v it has entry

    Q_uv = K_uw K_vw,

where w is their unique common neighbor. If they are adjacent
there is no common neighbor. Consequently the complementary
edges are partitioned into complete graphs on the d neighbors of
each middle vertex w. Each such complete graph has a rank-one
signing: the sign product on each of its triples is positive.

There is an exact converse. Given a symmetric sign Q_uv on
each complementary edge, first require all star-triangle sign
products to be positive. Then for each w choose incidence signs
a_(u,w) for its neighbors such that
Q_uv=a_(u,w)a_(v,w) on that star. This is possible by fixing one
incidence sign and using its edge signs to determine the others.
Any second such factorization differs by one star sign t_w.
The wanted skewness condition is precisely

    t_u t_w = −a_(u,w)a_(w,u)    for each edge uw of G.

This is soluble exactly when the right sides multiply to +1
on a cycle basis of G, by fixing t at a spanning-tree root,
propagating along the tree, and testing the remaining edges.
Connectedness then gives exactly two sign choices, differing by
simultaneous converse of all original graph edges. Thus the
sign-square map is a 2-to-1 map onto an explicitly affine subset
over F₂ with 2^(e(G)−1) elements. Its local star tests and global
cycle tests are both necessary. Switching Q by signs at vertices
preserves the image since diagonal switching K↦DKD switches
KKᵗ by the same conjugation.

## Every switching is substantially nonbalanced

At each star, its negative signs count the product of incoming
and outgoing degrees in the orientation. If m_w is the minority
of these two degrees, then the number of negative unordered
complement edges is

    F = Σ_w m_w(d−m_w).

For d=57 each summand is at least 29m_w. Partition vertices by
their majority direction. Every original edge inside either part
contributes a minority incidence at one endpoint. By the Laplacian
top eigenvalue 65, any cut of the graph has at most
floor(65n/4) crossing edges, so it leaves at least
ceil(n(57−8)/4)=39813 within-part edges. Hence for every signing
and every switching,

    F ≥ 29·39813 = 1,154,577.

In particular the complement signs cannot be made all positive.
A second spectral explanation of this last fact is that an
all-positive switched KKᵗ would have a simple largest positive
eigenvalue d², whereas the nonzero singular-square eigenvalues
of a real skew matrix occur in pairs.

An optimizer of F can always be sought among total orders:
indeed put b_w=outdeg(w)−indeg(w), so minimizing F maximizes
Σ b_w². Flipping an arc u→v changes that sum by
4(b_v−b_u+2). At an optimum it follows that b_u≥b_v+2
on every arc. Hence every optimum is acyclic and extends to a
total order. Its majority cut is oriented consistently on all
crossing edges. A concrete open question is whether some optimum
order can always have a maximum-cardinality majority cut for
the exceptional Moore graphs; the maximum-cut lower bound by
itself does not prove this.

## Exact six-cycle census and signed bounds

Each triangle in the complement either lies in one of the stars
(positive sign product) or its three pairwise common neighbors
are distinct. In the latter case they are also distinct from the
triangle vertices, and the resulting alternating walk is a simple
six-cycle in G. Conversely each six-cycle gives exactly two
alternating complement triangles. For the cycle sign sigma defined
as the product of the six entries K along a traversal, these two
triangle products are both −sigma. The sign is independent of
the traversal direction since the cycle has even length.

The complement degree is n−d−1, and each complementary edge
has n−2d−1 common neighbors in the complement. Thus

    C6(G) = [n(n−d−1)(n−2d−1)/6 − n·choose(d,3)]/2.

At d=57 this is 2,662,660,000. There are many cycles of each sign,
regardless of the orientation:

* **Positive six-cycle sign.** Fix a vertex r and switch at its
  complement neighbors so that all their edges to r are positive.
  The count of negative edges among those neighbors is exactly
  the number of negative complement triangles incident to r.
  Choose the signs at the d original graph neighbors of r randomly.
  Every remaining eligible edge is negative with probability 1/2.
  Their number is d(n−d−1)−choose(d,2): the original neighbors are
  independent, so their internal complement edges were counted
  twice by summing complement degrees. Thus some choice costs at
  most 90,174 additional negatives at d=57. The universal bound
  on F implies that each r lies in at least 1,064,403 negative
  triangles. Summing over r and dividing by six (three vertices
  per complement triangle, two triangles per original six-cycle)
  gives C6(+)≥576,551,625.
* **Negative six-cycle sign.** Put M=KKᵗ. Since G has girth five,
  a closed walk of length six either reduces by immediate
  backtracks to the empty walk or is a simple six-cycle. Indeed
  a nonempty reduced core would have length at least five, and
  removing backtracks changes length by two; at length six a
  repeat would yield a shorter reduced cycle. The five Dyck
  shapes of tree walks give d³+2d²(d−1)+2d(d−1)² choices per
  root, and each contributes a negative sign to K^6 (three
  opposing arc pairs). Thus

      tr M³ = nd(5d²−6d+2) − 12(C6(+)−C6(−)).

  On the other hand tr M=nd and tr M²=nd(2d−1), using its fixed
  diagonal and nonzero pattern. Cauchy on its nonnegative
  eigenvalues gives tr M³≥nd(2d−1)². So
  C6(+)−C6(−)≤nd(d−1)²/12, which at degree 57 implies
  C6(−)≥1,307,124,000.

## Calibration

The standard-library checker constructs C5, the Petersen graph
as KG(5,2), and Hoffman–Singleton via five pentagons, five
pentagrams, and affine cross matchings. It independently checks
all common-neighbor conditions and the complement-triangle
bijection for the increasing-index orientation. Its counts
(F,C6(+),C6(−)) are (3,0,0), (6,3,7), and (270,1929,3321),
respectively. Exhausting all orientations at degrees two and
three yields exact minima 1 and 6 for F. These are finite
calibrations, not assumed features of the unknown graph.
