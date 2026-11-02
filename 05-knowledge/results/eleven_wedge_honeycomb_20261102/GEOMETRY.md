# A gap repairs the hyperbolic tree's inverse distance control

Status: **PROVED elementary refinement**, using the standard upper-half-plane
hyperbolic distance formula. No claim of literature novelty is made. The
[incoming dyadic construction](../session140_collatz_hexagonal_geometry_20261102/curvature_control/PROOF.md)
already supplies a proper induced constant-edge binary tree but proves it is
not quasi-isometric. Its final accepted audit is recorded in
[checkpoint2](../session140_collatz_hexagonal_geometry_20261102/CHECKPOINT2.md).
The refinement here retains gaps between sibling intervals.

## All levels, with the original addresses retained

Fix an integer b>=2. For a binary word w=(e_1,...,e_k), let

    x_w=sum_(t=1)^k (b-1)e_t b^(-t) +(1/2)b^(-k),
    y_w=b^(-k).

The empty word is at(1/2,1). Use the hyperbolic metric
ds^2=(dx^2+dy^2)/y^2, whose standard distance formula is

    cosh d(z,z')=((x-x')^2+y^2+y'^2)/(2yy').           (G1)

The standard formula is the same cited input as the predecessor, sourced to
[Felikson's university notes, section6.2](https://www.maths.dur.ac.uk/users/anna.felikson/Geometry/Geometry23/notes.pdf).
All construction-specific conclusions below are derived here.

**Theorem.** The parent-child geodesics give a proper noncrossing embedding
of the rooted binary tree. Its distinct vertices have minimum distance

    L_b=arcosh(1+5(b-1)^2/(8b)),                       (G2)

attained exactly on parent-child pairs. For b>=3 it also satisfies

    d_T(w,v) log b - log(2b^2/(b-2)^2)
      <= d_H(z_w,z_v) <= d_T(w,v) log b + log3.        (G3)

Thus it preserves the tree metric after scaling by log b, up to a bounded
additive error. At b=3 the edge equation is

    cosh L_3=1+20/24=11/6,
    d_T log3-log18 <= d_H <= d_T log3+log3.

Here20 and24 are explained by the horizontal-plus-vertical displacement in
(G1), not identified with an aliquot source or a path count.

**Proof of separation and planarity.** A parent(a,h) has children
(a +/- c h,h/b), c=(b-1)/(2b), giving(G2). Distinct same-level centres
are separated horizontally by at least(b-1)h, so their cosh distance is
at least1+(b-1)^2/2, strictly larger than cosh L_b. If levels differ by
at least2, the vertical term gives cosh d>cosh L_b: after multiplying
the difference by8b^2 it equals

    (b-1)^2(4b^2+3b+4)>0

at level difference2, and increases thereafter. At adjacent levels, a child
of a different parent has horizontal displacement at least(b+1)h/(2b),
larger than ch. Hence equality at L_b occurs exactly at tree edges.

The circle centre of the edge to the child on side s is
a-sDh, where D=(3b+5)/(4b)>0, and its squared radius is h^2(D^2+1).
Both endpoints lie on the same side of that centre. The selected geodesic
segment is monotone, stays in its horizontal interval and height band[h/b,h],
and meets the band boundaries only at its endpoints. Outgoing intervals of
different parents are disjoint, as are the open bands at different levels.
The only intersections are genuine shared endpoints. A compact hyperbolic
set has a positive minimum height, so it meets only finitely many levels
and edges. This proves the proper embedding claim.

**Proof of(G3).** Let k,l be word lengths, h their longest-common-prefix
length and t=k+l-2h=d_T(w,v). Both centres lie in the same depth-h interval,
so |x_w-x_v|<=b^(-h). Each of the three numerator terms in(G1), divided
by y_w y_v, is at most b^t. Thus cosh d<=3b^t/2 and d<=log(2cosh d)
gives the upper bound. If neither word is an ancestor of the other, their
distinct child intervals are separated by(b-2)b^(-h-1). Hence
cosh d>=(b-2)^2 b^t/(2b^2), and d>=log(cosh d) gives the lower bound.
For an ancestor pair, d>=|log(y_w/y_v)|=t log b gives it directly. QED.

At b=2 the gap vanishes. The rightmost and leftmost depth-k points on opposite
sides of the root split have tree distance2k but constant cosh distance3/2.
The predecessor's counterexample is therefore the exact limiting boundary
of this proof, not an obstacle removed by assuming a missing inverse bound.

## What the triangle inequality does and does not say

The triangle inequality puts all depth-k vertices within radius kL_b of the
root. It does not bound how many separated vertices that ball can contain.
Euclidean separated packing has quadratic area growth; hyperbolic balls have
area2pi(cosh R-1), which grows exponentially. The metric still obeys the
triangle inequality. Its angular geometry and volume growth differ.
The Cantor gap adds a second property: separated branches cannot secretly
approach each other at bounded hyperbolic distance as their tree distance
grows. This is a control on distance from below, absent from the triangle
inequality's path-length bound from above.

This remains an abstract binary-tree construction. Applying it to a selected
Collatz inverse tree retains its word addresses and ordinary-step clocks;
it does not embed every Collatz edge, identify shared native halving tails,
or show universal forward entry. Those obligations remain distinct.

## Honeycomb capacity, small unit-distance maxima and coloring

The independent controls reproduce704 sources reaching5 by ordinary time21,
versus694 sites in the honeycomb radius21 ball, and2263 versus2107 in the
triangular-lattice comparison at26. The first deficit is20 for this ROOT,
CLOCK and METRIC. Changing any of them changes the comparison. Removing
the1–4–2 cycle has no effect on this already acyclic subtree.

The older June unit-distance work asks a different question: maximize the
number of unit pairs among a FIXED number of planar points. AMP's dated
[Theorem1 and Table1](https://arxiv.org/html/2412.11914v2) give u(10)=20,
u(20)=54 and u(21)=57. These are cited results, not re-proved here.
The present deficit20 counts missing hosting SITES, not unit EDGES.

Incoming [Session141](../session141_metric_aliquot_square_paths_20261102/unit_distance/PROOF.md)
makes the useful exact bridge: a graph admits k partial successor layers
iff every induced W satisfies |E(W)|<=k|W|. Its proof uses endpoint slots
and Hall's theorem. Consequently the first u(n)>kn is exactly the first
unit-distance graph exceeding k functional layers. This is inherited
classical capacity machinery, not a discovery of this session.

A honeycomb subdivision retains reachability only with marked original
vertices and route clocks. It stays bipartite and degree<=3, unlike dense
unit-distance extremizers. Conversely, Cartesian powers of an equilateral
triangle have unbounded edge density and remain three-colorable, as the
incoming faithful-product theorem shows. Density, metric capacity and
chromatic forcing are three different preserved predicates.

## An actual geometric restriction on the stability calculation

The sharp formal eight-vertex bracket is attained at2K4. Four points cannot
be mutually equidistant in the Euclidean plane. In fact enumerating all
6,431 K4-free graphs on8 vertices gives the stronger endpoint bound

    512+33w-32m2+7m3 >=240,  hence B8>=15/32.          (G4)

It follows that **if the positive-variance support graph of A is K4-free**,
the order-eleven theorem retains triangle coefficient
15/32-7/40=47/160 instead of153/640. All other terms remain unchanged.
Separate affinity is applied only inside the fixed K4-free support, so
every endpoint graph remains K4-free. In particular this applies if the
support edges are realized as unit distances among distinct plane points.

This is a precise restricted-family bridge from Euclidean geometry to the
weighted tournament inequality. Arbitrary tournament variance supports are
not unit-distance graphs, and this corollary gives no HN coloring bound.
The equality witness for(G4) is two K4-minus-edge diamonds joined at their
degree-two ports by two independent edges. `geometry_controls.py` verifies
the complete restricted census; K4-freeness is necessary, not sufficient,
for geometric realizability.

This particular equality support IS geometrically realizable. Take the
diamond with vertices(+-1/2,0),(0,+-sqrt(3)/2), its translate by(0,1),
and the two vertical unit edges joining the degree-two ports to their
translates. All eight vertices are distinct and every selected edge is
unit length. Extra unit pairs are allowed: the support need only be a
subgraph of a unit-distance graph. Thus unit realizability of the support
alone cannot improve(G4) beyond15/32. The code checks this over Q(sqrt(3))
and tests the stronger theorem on exact rational weighted matrices. This
does not assert sharpness of the final triangle coefficient47/160.
