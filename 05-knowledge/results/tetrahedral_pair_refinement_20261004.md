# Three pair-matching charts and a symmetry-preserving tetrahedral refinement

2026-10-04 (America/Denver).

**Status:** PROVED elementary geometry, reconstruction, symmetry, and
minimality statements with their hypotheses below; FINITE-EXACT for all
enumerations and rational-coordinate checks. No new general subdivision
theory or arithmetic dynamics is claimed.

## Inheritance and the precise bridge

The [pair-tile note](medial_pair_tiles_20261004.md), sections 1--2 and 4,
already distinguishes the six unordered pairs of four labels, their
octahedral overlap graph `J(4,2)`, and its 48 automorphisms from the 24
permutations induced by the original labels. Remembering which octahedral
faces are old vertex-stars removes the extra duality. The
[cumulative-coordinate note](edgewise_cumulative_coordinates_20261004.md),
sections 3--4, shows that a degree-two tetrahedral edgewise subdivision
selects one central diagonal and retains only eight coordinate symmetries.

These are the same six objects, through an explicit affine map: an
unordered pair `{i,j}` goes to the midpoint `m_ij=(e_i+e_j)/2`. Pair
overlap is adjacency of the midpoint octahedron, and pair complement is
opposition through its center. This supplies an actual connection between
the pair graph and subdivision, rather than a match of their vertex counts.

The hostile is retaining a selected diagonal while claiming all original
vertex permutations remain symmetries. The corrected near miss is the
distinction between an invariant family of three charts and an invariant
single mesh. The least-used sidecars are the chosen perfect matching,
the old-vertex/old-face color of octahedral facets, and the sign omitted
by one chart. The board consists of pair labels, matching charts, central
cut coordinates, common refinement, volume, and retained parent maps.

## 1. The six pairs give three possible edgewise charts

Work in the standard tetrahedron `Delta=conv(e_1,e_2,e_3,e_4)`, with
nonnegative barycentric coordinates summing to one. Its six midpoints
form an octahedron O. There are three opposite-pair diagonals:

    m_12--m_34,  m_13--m_24,  m_14--m_23.             (1)

They are exactly the three perfect matchings of the original four labels.
All meet at `c=(1/4,1/4,1/4,1/4)`.

Each degree-two edgewise mesh contains four corner tetrahedra

    T_i=conv(e_i, m_ij for j != i).

The remaining octahedron is triangulated into four tetrahedra using one
diagonal from (1). Its other four vertices form a cycle around that
diagonal; adjoining each cycle edge to the two diagonal endpoints gives
the four cells. All eight tetrahedra have relative volume 1/8. Their
f-vector is `(10,25,24,8)`, with 16 triangular boundary faces.

The cumulative chart with coordinate order 1,2,3,4 selects `13|24`.
Permuting the coordinate order gives exactly three embedded meshes, each
arising from eight of the 24 permutations. The action on the three
matchings has Klein-four kernel; the stabilizer of a chosen matching has
order eight and is dihedral. No matching is fixed by the entire S_4.

This retains an S_4-equivariant *family* of charts if the matching label
is kept. A single chosen chart has a smaller symmetry group. The repair
below instead constructs one common invariant mesh.

## 2. Cut coordinates retain the lattice condition

For general coordinates x_1,...,x_4, introduce

    s=x_1+x_2+x_3+x_4,
    u=x_1+x_2-x_3-x_4,
    v=x_1+x_3-x_2-x_4,
    w=x_1+x_4-x_2-x_3.

The transformation matrix is

    H = [1  1  1  1
         1  1 -1 -1
         1 -1  1 -1
         1 -1 -1  1].

It satisfies `H^2=4I` and `det(H)=-16` in the displayed row order. The
inverse is H/4, explicitly

    x_1=(s+u+v+w)/4,  x_2=(s+u-v-w)/4,
    x_3=(s-u+v-w)/4,  x_4=(s-u-v+w)/4.               (2)

This is not a unimodular integer shear. An integer tuple `(s,u,v,w)`
comes from integer x_i if and only if the four coordinates have the same
parity and `s+u+v+w=0 mod4`. Equivalently, all four numerators in (2)
are divisible by four. Necessity follows directly from Hx; sufficiency
follows because the four numerators then agree modulo four. This is an
index-16 sublattice. Nonnegative simplex coordinates additionally require
the four numerators in (2) to be nonnegative.

**Exact Fourier interpretation.** Label the four vertices by the Klein
four group `V_4=(Z/2Z)^2` in this order:

| Vertex | Group label | u character | v character | w character |
|---|---|---:|---:|---:|
| 1 | (0,0) | +1 | +1 | +1 |
| 2 | (0,1) | +1 | -1 | -1 |
| 3 | (1,0) | -1 | +1 | -1 |
| 4 | (1,1) | -1 | -1 | +1 |

Then H is literally the character table: its rows are the functions
`1, (-1)^a, (-1)^b, (-1)^(a+b)`. The rows are closed under pointwise
multiplication, with `u*v=w` and each nontrivial row squaring to 1.
Their orthogonality gives `H H^T=4I`, so H/2 is orthogonal, whereas the
inverse of the unnormalized histogram transform H is H/4. A barycentric
histogram on the four vertex labels is thus stored as total mass s and
the three nontrivial Walsh coefficients u,v,w. These three characters
are exactly the three pair-matching cuts, rather than numerical analogues
of them.

Every permutation of the four labels is an affine automorphism of V_4:
the faithful affine group has `4*|GL(2,2)|=4*6=24` elements. Translations
produce the even sign changes of the three nontrivial characters; linear
automorphisms permute those characters. This explains the Klein-four
kernel of the three-matching action and the signed-permutation subgroup
used below. This Fourier carrier is V_4; it is not an identification with
an arithmetic valuation clock or an unrelated cyclic Fourier group.

At s=1, the six pair midpoints have cut coordinates `(u,v,w)` equal to
the six signed coordinate unit vectors, so

    O = {|u|+|v|+|w| <= 1}.

The four original corners are the four cube-sign triples with product
`uvw=+1`. The eight triangular octahedral facets correspond to all eight
sign triples. Four are old vertex-stars and have positive sign product;
four lie on old tetrahedron faces and have negative sign product.

The common center has cuts `(s,0,0,0)`. At degree two its barycentric
integer-address coordinates would be `(1/2,1/2,1/2,1/2)`: the cuts
`(2,0,0,0)` have equal parity but fail the mod-four condition. At degree
four it is `(1,1,1,1)`, with cuts `(4,0,0,0)`. These denominators follow
from this exact coordinate inverse; they imply no identification with
unrelated quadratic parameters such as -7/4 or -29/16.

## 3. The center gives the unique minimum common refinement

The chart selecting the u-axis diagonal has four central tetrahedra

    O intersect {sign(v)=epsilon_v, sign(w)=epsilon_w},

where the inequalities are closed and each epsilon is +/-1. Equivalently,
each is the convex hull of both u-axis endpoints and one chosen endpoint
on each of the other two axes. Thus the u-chart forgets the sign of u.
The v-chart forgets the sign of v, and the w-chart forgets the sign of w.

Intersecting cells from any two distinct charts fixes all three signs.
Their full-dimensional common cells are precisely the eight tetrahedra

    conv(c, m_u^{epsilon_u}, m_v^{epsilon_v}, m_w^{epsilon_w}),       (3)

where `(m_u^+,m_u^-)=(m_12,m_34)`,
`(m_v^+,m_v^-)=(m_13,m_24)`, and `(m_w^+,m_w^-)=(m_14,m_23)`.
In cut coordinates these are
`conv(0,epsilon_u e_u,epsilon_v e_v,epsilon_w e_w)`. Adding the four
unchanged corner tetrahedra gives a mesh K with

    f(K)=(11,30,32,12).

The boundary is exactly the original 16-face midpoint boundary. The four
corner tetrahedra each have relative volume 1/8, and the eight central
tetrahedra each have relative volume 1/16. Their volumes sum to one.
The symmetry repair therefore sacrifices equal cell volume.

**PROVED minimal common-refinement theorem.** K is the common refinement
of any two distinct matching charts and hence of all three. Every geometric
simplicial common refinement of those charts has at least 12 tetrahedra;
equality forces K itself. It must contain the center as a vertex.

Indeed the eight central octants and four corner cells are the nonempty
full-dimensional intersections of coarse cells. Every fine tetrahedron
must lie in one such intersection, and all 12 interiors must be covered.
At least one fine tetrahedron per intersection is necessary. With exactly
12, each intersection is itself that tetrahedron, proving uniqueness.
The center is also forced directly by the crossing of the two selected
coarse diagonals: a conforming common refinement must subdivide both
edges, whose intersection must be a vertex. K attains this lower bound
with just one added vertex. This is a minimality assertion for these
specified embedded charts, not for every conceivable tetrahedral mesh.

**PROVED ten-vertex obstruction with fixed boundary.** Prescribe the
standard degree-two subdivision on all four outer faces and allow only
the four original corners and six edge midpoints as mesh vertices. Then
the only tetrahedral triangulations are the three matching charts. In
particular there is no S_4-equivariant eight-cell equal-volume mesh under
these hypotheses, or any S_4-equivariant tetrahedral mesh under them.

To prove this, a corner's only permissible neighbors are its three
incident edge midpoints. An edge to another original corner skips an
existing midpoint; an edge to a nonincident midpoint lies in an outer
face and violates its prescribed subdivision. Thus the four corner
tetrahedra are forced. In the central octahedron, a nondegenerate
four-vertex tetrahedron must contain an opposite pair: otherwise it has
at most three vertices, one from each pair. It cannot contain two opposite
pairs because those four vertices are coplanar. Therefore it uses exactly
one of the three diagonals. Two distinct chosen diagonals would cross at
the absent center, so all central tetrahedra must use the same one.
There are exactly four such nondegenerate tetrahedra, and all four are
needed to fill O. This proves the classification and obstruction.

The finite verifier separately enumerates all 16 permitted nondegenerate
candidate tetrahedra: four corner cells and 12 central cells, four per
matching. The three resulting triangulations are directly compared with
the cumulative-coordinate rule and all 24 coordinate permutations.

## 4. Three views recover the missing sign

Label the eight central fine cells by
`epsilon=(epsilon_u,epsilon_v,epsilon_w)`. The three parent maps are

    pi_u(epsilon)=(epsilon_v,epsilon_w),
    pi_v(epsilon)=(epsilon_u,epsilon_w),
    pi_w(epsilon)=(epsilon_u,epsilon_v).              (4)

Each individual map is two-to-one. Any two distinct parent maps determine
all three signs, and the third then follows. The joint three-chart carrier
has eight consistent states inside the 64 formal triples of four-state
chart labels. The repeated sign coordinates must agree. The script
verifies (4) by independent tetrahedral containment calculations.

This is a concrete information-storage connection: the source is the
central common-refinement cell; the target is its three coarse parent
cells; the maps preserve geometric containment; each coarse projection
loses one bit; a second chart restores it. It describes cells, not exact
points. Local barycentric weights are still needed to recover a point
inside a cell, and points on common faces need their minimal-face address
or a retained incident-cell convention.

## 5. Which symmetry was restored, and which was never present

The midpoint overlap graph is the octahedron `J(4,2)`. Its full group of
48 acts by signed permutations of the three cut axes. Its six vertices
alone forget whether an octahedral triangular facet represents an old
vertex-star or an old face. Pair complement is central inversion and
exchanges these two types.

The original S_4 is the subgroup with an even number of cut-coordinate
sign flips, preserving the product-sign color on the eight facets. It
is not the orientation-preserving octahedral subgroup: odd axis
permutations are allowed. Of the 24 original actions on cut coordinates,
12 have determinant +1 and 12 have determinant -1.

The exact group counts in the finite verifier are:

| Carrier | Full graph automorphisms |
|---|---:|
| Six midpoint vertices and octahedral edges | 48 |
| Same carrier with one selected opposite-pair edge | 16 |
| Central octahedron coned to its center | 48 |
| Full ten-vertex edgewise tetrahedral mesh | 8 |
| Full eleven-vertex center repair K | 24 |

For K, original corners have degree three, the center has degree six, and
midpoints have degree seven. These classes are intrinsic. Each midpoint
is adjacent to exactly its two original endpoints, so an automorphism
injects into S_4; all 24 permutations are realized. The full complex has
three tetrahedron orbits, each of size four: corner cells, central cones
on vertex-star facets, and central cones on old-face facets. It is not
tetrahedron-transitive. The bare central star, after forgetting the old
corner attachments, has one orbit of eight tetrahedra under its larger
group.

These are abstract complex/graph and affine-coordinate symmetries. They
are Euclidean symmetries for the displayed regular simplex. An arbitrary
affine image need not preserve Euclidean lengths or angles, although the
incidence, containment, and relative-volume results persist. The midpoint
graph's extra 24 symmetries are not silently transported to the original
tetrahedron.

## 6. A separate 288-cell completion also refines the original flags

The center repair K has 12 tetrahedra. Its barycentric subdivision S(K)
has 24 children per tetrahedron and

    f(S(K))=(85,420,624,288).

It has 96 boundary triangles. The relative-volume histogram is 96 cells
of volume 1/192 and 192 cells of volume 1/384. This is a separate flag
completion; no claim says 288 is a minimum.

**PROVED.** S(K) refines both all three matching charts and the original
tetrahedron's barycentric subdivision S(Delta). The first assertion
follows since K already refines the three charts. For the second, each
coordinate transposition acts simplicially on the invariant mesh K.
After barycentric subdivision, a stabilized simplex has all its vertices
fixed: its vertices encode a strict chain of faces of distinct dimensions,
and these dimensions cannot be permuted. If a point is fixed, its unique
minimal simplex is stabilized, so the entire minimal simplex is fixed.
Thus the fixed plane section `x_i=x_j` is a subcomplex of S(K), lying in
its two-dimensional skeleton. The six such planes separate Delta into
the 24 coordinate-order chambers, precisely its barycentric tetrahedra.
No full-dimensional child can cross a wall, proving the assertion.

The argument needs an affine simplicial action, not an isometric one.
The [mixed-refinement note](mixed_triangle_refinement_20261004.md)
develops the corresponding orbit-overlay repair in general dimension.
Here every one of the 288 cells is checked against its actual K parent,
all three edgewise parents, and an original barycentric chamber, with
independent exact determinant containment tests.

## Verification

Run:

    python -X utf8 -B 04-computation/experiments/tetrahedral_pair_refinement_20261004.py
    python -X utf8 -B -O 04-computation/experiments/tetrahedral_pair_refinement_20261004.py

The [saved output](tetrahedral_pair_refinement_20261004.out) retains exact
counts, volumes, matching labels, all eight joint chart codes, automorphism
counts, lattice conditions, and the 288-cell completion. Its finite universe
is all three matching meshes, every four-subset of the ten original mesh
vertices, all 24 original coordinate permutations, exhaustive graph
automorphisms of the five small carriers in the table, all 2,401 integer
cut tuples in [-3,3]^4, and every barycentric child of K. The lattice-image
test includes the degree-two parity-only hostile and degree-four repair.

All geometric computations use Fraction coordinates and determinant
ratios. Normal and optimized outputs agree; all acceptance checks remain
active under Python -O. Minimality and general group statements rely on
the proofs above, with these finite computations as independent controls.
