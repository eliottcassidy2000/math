# Scale, flags, and a common refinement that retains both parents

2026-10-04 (America/Denver).

**PROVED elementary:** the two mixed face-count formulas, one-level common
refinement, intrinsic nonisomorphism at degree two, and the intersection
construction. **FINITE-EXACT:** the explicitly bounded Fraction-mesh checks
below. **REFUTED:** commuting the two operations, or iterating the one-level
common-refinement statement without a new containment proof. No historical
novelty is claimed for barycentric subdivision or polyhedral intersections.

## Inheritance and precise objects

The closest proved mechanisms are the cumulative lattice coordinates in
[edgewise coordinates](edgewise_cumulative_coordinates_20261004.md), the typed
incidence reconstruction in [medial pair tiles](medial_pair_tiles_20261004.md),
and [noble flag refinement](noble_subdivision_flags_20261004.md). The hostile
is a pair of complexes with the same counts and different incidences. The
corrected near miss is extrapolating a one-step containment to every depth.
The least-used sidecar is a pair of parent-cell labels, with local weights.

Let Delta be the affine triangle with vertices e1,e2,e3. A point has
barycentric coordinates x+y+z=1, all nonnegative. Its vertices are Euclidean
unit vectors; its other points generally are not. Scaling to x+y+z=r puts
the integer lattice on the face. If a spherical realization is wanted,
radial normalization is a separate operation after forming this mesh.

Write E_r for degree-r edgewise subdivision and S for barycentric
subdivision. Composition is read from right to left. For triangles E_r
does not depend on a vertex ordering. Both constructions extend affinely
to any nondegenerate triangle and preserve its boundary as a subcomplex.

The cumulative rule in the companion note is the standard common-sign
condition on C(x,y,z)=(x,x+y,x+y+z). Taking coordinatewise absolute values
instead produces crossing chords already at r=2; that separate graph is
not used in this note. The primary definition is
[Athanasiadis, section 4](https://arxiv.org/pdf/1310.0521).
The classical multiplicative edgewise refinement context is
[Edelsbrunner--Grayson](https://www.graysonfamily.org/dan/Papers/Esd-scan.pdf).

## 1. Two parameters with different composition laws

For a triangle,

    f(E_r Delta)=((r+1)(r+2)/2, 3r(r+1)/2, r^2),
    boundary_edges(E_r Delta)=3r,
    E_s(E_r Delta)=E_(rs) Delta.

The equality is equality of geometric meshes. Triangle edges lie on the
three integer barycentric coordinate-line families. Refining an r-cell by s
adds precisely the same parallel lines at denominator rs. In dimensions
at least three the inherited local order must be retained; arbitrary
reordering does not support this equality (see the companion audit).

For k>=0, let V_k,E_k,F_k,B_k count S^k Delta. Every triangle gives six new
triangles and every boundary edge gives two boundary edges. Euler's disk
identity and 3F=2E-B then give

    F_k=6^k, B_k=3*2^k,
    V_k=(6^k+3*2^k)/2+1,
    E_k=(3*6^k+3*2^k)/2.

Thus both procedures can be iterated. The distinction is multiplicative
lattice scale versus a sequence of incidence flags, not recursion versus
an inability to recurse. The same face count does not make them equivalent:
S^2 Delta has (25,60,36), whereas E_6 Delta has (28,63,36).

## 2. Even the entire f-vector can miss the operation order

**PROVED.** For every positive integer r,

    f(E_r(S Delta)) = f(S(E_r Delta))
                   = (3r^2+3r+1, 9r^2+3r, 6r^2).

Both sides have 6r^2 triangles and 6r boundary edges. The same Euler and
edge-incidence identities give the remaining counts. Nevertheless the
meshes need not be abstractly isomorphic. At r=2 their graph degree
multiplicities are:

| Order | Degree 3 | Degree 4 | Degree 6 | Degree 7 |
|---|---:|---:|---:|---:|
| E_2 after S | 6 | 6 | 7 | 0 |
| S after E_2 | 9 | 3 | 4 | 3 |

These can also be obtained directly. In S(K), an old vertex has degree
deg_K(v)+number of incident faces, an edge center has degree three on the
boundary or four inside, and every triangular face center has degree six.
For E_2(S Delta), each old vertex keeps its degree; each new edge midpoint
has degree four on the boundary and six inside. The table separates the
two abstract graphs without metric assumptions.

The quotient to the f-vector therefore loses order-sensitive incidence.
The cheapest separating invariant here is the vertex-degree multiset;
it is not claimed to classify arbitrary subdivisions.

## 3. A useful one-level common-refinement theorem

**PROVED.** For every r>=1, S(E_r Delta) is a geometric subdivision of
both E_r Delta and S Delta, for any affine nondegenerate triangle Delta.

The first containment is the definition of barycentric subdivision. For
the second, each coordinate transposition acts affinely on the
S3-invariant edgewise mesh. Its fixed line is one of the three medians.
In a barycentric subdivision, a simplex is a chain of faces of distinct
dimensions. A group element stabilizing this chain must fix each of its
vertices, since it preserves those dimensions. The action consequently
has no inversions, and its fixed set is a subcomplex: the minimal simplex
whose relative interior contains a fixed point is stabilized and fixed
pointwise. Here each median is therefore a union of edges and vertices
of S(E_r Delta). The three medians divide Delta into the six triangles of
S Delta. No small triangle can cross a median, so each lies in a single
one of those six chambers. This proves the claim.

This proof needs affine coordinate permutations, not literal Euclidean
reflection symmetry of an arbitrary physical triangle. It also exposes
why a dimension change needs a fresh audit: cumulative edgewise meshes
in higher dimension do not generally retain the full coordinate symmetric
group.

**Hostile 1: the opposite order.** E_2(S Delta) does not refine E_2 Delta.
One of its triangles has vertices

    (0,1/4,3/4), (0,1/2,1/2), (1/6,5/12,5/12).

These are not all contained in any one coarse E_2 triangle. The first
vertex lies above z=1/2 and the last below it.

**Hostile 2: blindly iterating the theorem.** S^2(E_2 Delta) does not
refine S^2 Delta. An exact witness triangle is

    (0,1/4,3/4), (1/18,11/36,23/36), (1/12,5/24,17/24).

The script exhausts the 36 coarse S^2 triangles and checks that none
contains all three points. The failed implication is
"K refines L, therefore S(K) refines S(L)" for arbitrary geometric
refinements. Barycenters of different parent cells need not align.
The strongest survivor is the one-level symmetry theorem above.

## 4. Repair: intersect first and keep two parent labels

For two finite geometric triangulations K,L of the same triangle, form
all their nonempty simplex intersections, with duplicates identified.
This is a polyhedral complex: an intersection of two cells is a face of
each, since supporting hyperplanes for the corresponding parent faces
restrict to supporting hyperplanes after intersection. Conversely each
face of an intersection is obtained by activating parent facet inequalities,
so it is itself an intersection of parent faces. Its cells cover
the triangle. Each full-dimensional cell is a convex polygon, with a
unique pair of full-dimensional parents (sigma,tau).

Barycentrically subdivide this intersection complex. A polygon's vertex
average lies in its relative interior; join it to the midpoint of each
edge and both edge endpoints. Adjacent polygons use the same midpoint
and endpoint on their shared edge. The result is a triangulation refining
both K and L. Degenerate intersections are their shared lower-dimensional
cells, not positive-area triangles.

Retain the two parent labels on each child, and its exact affine local
coordinates. These yield two independently checkable containment maps.
The labels alone do not reconstruct the full ancestral meshes, so the
record also retains their cell tables. This is a faithful shared carrier,
not an assertion that the two operations commute or that the chosen
common refinement is minimal.

For K=E_2 Delta and L=S^2 Delta the exact intersection complex has 54
polygonal cells: 42 triangles and 12 quadrilaterals. Its barycentric
triangulation has (187,534,348), boundary length 24, and a saved parent
pair for every triangle during verification. These counts show the cost
of a uniform repair; a more economical common triangulation is possible
and is not needed for the containment theorem.

**General symmetry repair, all finite dimensions.** Let K be any finite
geometric triangulation of a d-simplex. Overlay all its coordinate-permuted
copies gK, for g in S_(d+1), using the same intersection construction, now
with convex polytopes. The resulting complex O is invariant under every
coordinate permutation and refines every gK. Barycentrically subdividing O
removes inversions exactly as in section 3. Each fixed hyperplane x_i=x_j
is consequently a subcomplex of S(O). The chambers of these hyperplanes
inside the simplex are exactly the simplices of S Delta. Thus S(O) is a
common refinement of every gK and S Delta, and carries the full affine
coordinate-permutation group. This is a finite constructive theorem,
not a claim that the overlay or its subdivision is minimal.

For a degree-two tetrahedron, the three distinct coordinate-order meshes
have the [12-cell center completion](tetrahedral_pair_refinement_20261004.md)
as their overlay. Its barycentric subdivision has 12*24=288 tetrahedra and
also refines S Delta. The smaller 12-cell complex itself does not refine
S Delta: its corner cells are larger than the 1/24-volume flag tetrahedra.
This distinguishes restoring a group action from removing simplex inversions.

## 5. What transfers to route certificates, and what does not

The [rational-chart audit](anchor_coverage_refinement_20261004.md) subtracts
overlapping dyadic source cells exactly. Dyadic cylinders have the simpler
property that two intersecting cells are nested; geometric triangle cells
need not. Both calculations benefit from retaining the two inputs to a
refinement, but their maps and preserved predicates differ:

| Transfer | Source -> target | Preserved | Required retained data | Cheap hostile |
|---|---|---|---|---|
| Scale | ordered mesh -> E_r mesh | points, incidence, parent containment | parent chart and local order in higher dimension | reordered tetrahedron child |
| Flags | map -> order complex | typed incidence and group action | rank; geometric realization for metric claims | degree-4 edge centers |
| Common refinement | two meshes -> intersection subdivision | both containment relations | both parent tables and child coordinates | S^2(E_2) versus S^2 |
| Source-cell subtraction | two cylinder banks -> disjoint remainder | exact source set and density | residues, modulus, original-source descent clock | removing a whole partially overlapped parent |

No general Collatz proof, canonical tournament orientation, or noble
polyhedron follows from these representation transfers. The useful
recommendation is concrete: store addresses and incidence maps at the
level where a desired operation is defined, and separately audit what a
coarser statistic discards.

## Reproduction and boundary of computation

    python -X utf8 -B 04-computation/experiments/mixed_triangle_refinement_20261004.py
    python -X utf8 -B -O 04-computation/experiments/mixed_triangle_refinement_20261004.py

[Script](../../04-computation/experiments/mixed_triangle_refinement_20261004.py)
and [saved output](mixed_triangle_refinement_20261004.out) use Fraction
coordinates, exact signed-area containment, and active checks under
optimized Python. Universe: r=1..12 for mixed counts/common refinement;
r,s=1..5 for pure composition; k=0..4 for barycentric counts; r=1..4 and
k=0..2 for the intersection repair. Every tested mesh satisfies exact
total area, nondegeneracy, disk Euler characteristic, and edge incidence.
The two failed-containment controls remain expected failures. Finite
verification supports, and does not replace, the general proofs above.
