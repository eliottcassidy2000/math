# Cumulative coordinates, edgewise refinement, and the missing common sign

2026-10-04 (America/Denver).

**Status:** PROVED elementary adjacency, counting, symmetry, refinement,
and decoding statements below; FINITE-EXACT for the explicitly bounded
controls; CITED classical edgewise-subdivision context. No novelty claim
is made for edgewise subdivision or its coordinate construction.

## Inheritance and the object being tested

The closest repository mechanism is the distinction between faithful
addresses and their quotients in
[the pair-tile and typed-flag note](medial_pair_tiles_20261004.md), especially
sections 1 and 5. Its tetrahedral hostile distinguishes the transported
symmetry group from a larger graph automorphism group. The
[earlier tiling note](tiling_modular_atoms_20261003.md) similarly preserves
loop/path/free-pair positions rather than identifying equal scalar labels.
The historical `04-computation/triangle_matrix_codec.py` retains positions
within diagonal layers; it is background provenance, not a subdivision
theorem. Targeted searches found no earlier repository definition of this
cumulative-coordinate edgewise construction.

Here the corrected near miss is replacing a vector with its componentwise
absolute value before testing membership in a positive binary cube. The
smallest hostile is a scale-two triangle; the missing sidecar is whether
all nonzero cumulative differences have the same sign. The working board
is lattice addresses, signed cut differences, actual simplicial faces,
coordinate symmetry, refinement degree, and local reconstruction weights.
The anchor is the owner's proposed adjacency; the niche is higher-simplex
ordering; the wildcard is a lossless hierarchical point address.

For primary context, Athanasiadis,
[Edgewise subdivisions, local h-polynomials and excedances, arXiv:1310.0521v3](https://arxiv.org/pdf/1310.0521),
section 4, defines edgewise faces using cumulative-coordinate differences
with one common sign, together with a support-face condition. Edelsbrunner
and Grayson,
[Edgewise Subdivision of a Simplex, author conference manuscript](https://www.graysonfamily.org/dan/Papers/p24-edelsbrunner.pdf),
sections 2--3, describe an ordered abacus construction, equal-volume cells,
face compatibility, and multiplicative refinement. These are classical
dependencies/context, not newly named constructions. The explicit tests
and elementary proofs below make our precise conventions reviewable.

## 1. Literal rule, first failure, and repair

For an integer r>=1 let

    V_r={(x,y,z) in Z_{≥0}^3 : x+y+z=r},
    C(x,y,z)=(x,x+y,r).

These are integer coordinates in the dilated triangle. Divide by r to
obtain points in the original triangle with vertices e_1,e_2,e_3. Degree
zero collapses the triangle and is not a subdivision. The zero binary
vector represents equality, so it is excluded when making a simple graph.

The literal componentwise-absolute rule connects distinct u,v whenever
each coordinate of `C(v)-C(u)` lies in {-1,0,1}. Its last coordinate is
always zero; only four of the eight absolute binary triples can occur.
At r=1 the rule happens to work. At r=2 take

    u=(0,2,0), v=(1,0,1), C(v)-C(u)=(1,-1,0).

The literal rule admits uv. Its squared length in the dilated coordinate
plane is 6, whereas an ordinary mesh edge has squared length 2. It crosses
the legitimate edge joining p=(0,1,1) and q=(1,1,0), since
`(u+v)/2=(p+q)/2`. All four vertices form a clique under the literal rule,
but they are coplanar: its abstract 3-simplex is geometrically degenerate.
Thus the literal graph does not give the proposed geometric triangulation.

**Repair:** connect distinct u,v exactly when

    C(v)-C(u) in {0,1}^3 OR C(u)-C(v) in {0,1}^3.       (1)

The choice of sign is shared by the whole vector. Equivalently, the two
cumulative vectors must be comparable coordinatewise and at distance at
most one in every coordinate. In the triangle the three positive nonzero
possibilities are `(1,0,0),(0,1,0),(1,1,0)`. Applying the inverse map

    C^-1(a,b,c)=(a,b-a,c-b)

shows that (1) is exactly `v-u=e_i-e_j` for distinct i,j, with either
orientation allowed. It is the ordinary triangular lattice graph.

For every r>=1, its vertex/edge/triangle counts are

    f=(C(r+2,2), 3r(r+1)/2, r^2).                     (2)

There are C(r+1,2) upward triangles and C(r,2) downward triangles.
For example, their vertices are respectively `b+e_i` with sum(b)=r-1,
and `b+e_i+e_j` with sum(b)=r-2. Each of the three unoriented root-edge
directions occurs C(r+1,2) times. This proves (2), including boundary cases.

The literal rule adds exactly C(r,2) edges, of differences
`+/- (1,-2,1)`. Each is a second diagonal in one of the C(r,2) complete
cumulative-coordinate unit squares. Its clique counts are therefore

    (C(r+2,2), 2r^2+r, 2r^2-r, C(r,2)).               (3)

The last entry counts four-vertex cliques, not geometric tetrahedra; zero
entries are omitted. Every added diagonal creates two extra triangle
cliques and one four-clique, and a clique cannot span more than one unit
square. This proves (3). At r=2, (2) is `(6,9,4)` whereas (3) is
`(6,10,6,1)`. These are clique counts, not faces of a planarized drawing.

## 2. The graph itself remembers triangular addresses

For the repaired triangle graph,

    dist(u,v)=1/2 sum_i |u_i-v_i|.                    (4)

Each root transfer reduces this quantity by at most one; transferring a
unit from any surplus coordinate to any deficit coordinate attains the
bound and stays in V_r. Its three corners are exactly the degree-two
vertices. If c_i=r e_i is the chosen i-th corner, then

    dist(v,c_i)=r-v_i, hence v_i=r-dist(v,c_i).         (5)

Thus the unweighted graph plus ordered corners reconstructs every integer
barycentric address. Without the corner order, precisely the S_3 relabeling
ambiguity remains: every graph automorphism permutes the corners, (5)
then fixes all vertices, and all six coordinate permutations are actual
automorphisms. This proves the full triangle graph automorphism statement.
It does not assert such a statement for every higher-dimensional graph.

The map C is an invertible integer shear, with determinant one. It loses
no information. Its differences can be read as signed cumulative cut
flows between ordered coordinate bins. Taking componentwise absolute
values forgets whether those cut flows agree in direction. Retaining only
an unordered multiset of r coordinate labels instead forgets the original
ordering of those labels, with fibre size `r!/(x!y!z!)`. These are distinct
losses and should not be conflated.

## 3. Higher simplices: exact faces, counts, and coordinate symmetry

For dimension d, replace triples by d+1 nonnegative integer coordinates
summing to r, and use all d+1 cumulative sums. A set of vertices is a face
when every pair satisfies (1) in this larger space. For a full simplex
this is a flag complex, so its graph determines its faces. For an arbitrary
original simplicial complex one must additionally retain the condition
that the union of supports lies in an original face. A graph alone can
lose that condition: a triangle boundary and a filled triangle already
have the same graph at r=1.

The first d cumulative coordinates identify the dilated simplex with

    0<=s_1<=...<=s_d<=r.

Partition integer unit cubes by the order of their fractional coordinates.
The resulting simplices have vertices forming coordinatewise chains,
each successive difference being a distinct standard basis vector. Their
pairwise differences are exactly those in (1). The chamber walls
`s_i=s_(i+1)` and its outer integer walls are unions of faces in this
partition, so restricting to this chamber is a conforming triangulation.
Each top simplex has determinant +/-1; the cumulative map is unimodular.
Consequently there are r^d top simplices, all of equal volume, and
`C(r+d,d)` vertices.

An exact formula for every face count, proved from this description, is

    f_k(d,r)=sum_(j=0,...,k) (-1)^(k-j) C(k,j)
                                      * C(r(j+1)+d,d), 0<=k<=d.       (6)

Indeed every face is unimodular in its affine integer lattice: differences
between successive selected chain vertices have disjoint nonempty
coordinate supports, each containing an integer pivot of size one. A
k-face dilated by an integer t>=1 thus has `C(t-1,k)` lattice points in
its relative interior. The disjoint relative interiors give

    C(rt+d,d)=sum_(k=0,...,d) f_k(d,r) C(t-1,k).

Taking successive finite differences at t=1 proves (6). The tetrahedral
case r=2 has `(10,25,24,8)`.

**Coordinate-symmetry theorem.** For d+1=m>=3 and r>=2, the coordinate
permutations preserving this subdivision are exactly the dihedral group
of the circularly ordered m coordinate positions. The count is 2m; at
r=1 every coordinate permutation works. For m=2 the group has size two.

For sufficiency, a cyclic shift changes the cumulative difference sequence
by subtracting the value at its new cut and cyclically reindexing. A binary
sequence with one common sign consequently remains binary with one common
sign. Reversal negates and reverses such cumulative differences. Thus all
dihedral coordinate permutations preserve (1).

For necessity, mark four distinct coordinate positions i,j,k,l. After
adding the same `(r-2)e_a` to both, the endpoints `e_i+e_k` and `e_j+e_l`
are adjacent precisely when their positive and negative marks alternate
around the coordinate circle. Their cumulative difference then visits
only 0 and one of +/-1; every nonalternating pattern visits either both
signs or magnitude two. This crossing relation recognizes consecutive
pairs on the circle: such a pair crosses no disjoint pair, whereas every
nonconsecutive pair has a crossing pair. Any preserving coordinate
permutation must therefore preserve the cycle graph, whose automorphisms
are dihedral. For m=3 every permutation is already dihedral.

In a triangle D_3=S_3, explaining why the visible cumulative chart does not
break triangle symmetry. In a tetrahedron, r=2 chooses the internal
octahedral diagonal joining `(1,0,1,0)` to `(0,1,0,1)`. Swapping coordinate
positions 2 and 3 changes this to the nonedge joining `(1,1,0,0)` and
`(0,0,1,1)`. Only eight of the 24 coordinate permutations preserve the
chosen subdivision. The literal absolute rule already has only two of
the six coordinate symmetries in a scale-two triangle.

These are coordinate-permutation groups of the specified embedded mesh.
No unproved identification with a higher graph's full abstract
automorphism group or with the Euclidean symmetries of a scalene simplex
is being made. A common global vertex order, restricted to every face,
is a sufficient way to make higher-dimensional charts agree on overlaps.

## 4. Refinement degree multiplies; recursion depth is a separate parameter

Let E_r denote the repaired edgewise subdivision. Order the vertices of
each child simplex by their cumulative-coordinate chain. Then

    E_s(E_r(Delta^d))=E_(rs)(Delta^d)                 (7)

as embedded complexes. Here a local integer vertex `w=(w_0,...,w_d)` with
sum(w)=s in a child with integer vertices `v_0,...,v_d` maps to

    sum_i w_i v_i,

an integer address of total degree rs. In cumulative space a full child's
vertices are `b,b+e_(pi1),...,b+e_(pi1)+...+e_(pid)`. For chain-ordered
local barycentric coordinates, each global coordinate is an integer offset
plus one minus one local cumulative coordinate. Thus this local map takes
the finer coordinate-chain simplices to the same fine cube triangulation,
with a permutation of coordinate directions and an overall reversal.
Their unions fill the original simplex, proving (7). The same restriction
argument applies to boundary faces.

Child ordering matters for d>=3. With r=s=2 in a tetrahedron, swapping the
two middle vertices of every inherited child chain changes 32 of the 64
fine tetrahedra. Both constructions still have the same top-cell count;
they are not the same embedded refinement. The saved output retains a
specific changed cell. In dimension two every vertex permutation preserves
the local edgewise subdivision, so this particular ambiguity disappears.

Repeating a fixed edgewise factor a for k rounds gives degree `r=a^k`
and `a^(kd)` top cells, not degree k. More generally the final degree is
the product of the factors. The flattened mesh forgets their chosen
factorization and the sequence of refinement operations.

Barycentric subdivision S is different: it uses chains of nonempty faces.
For a triangle its counts evolve as

    (V,E,T) -> (V+E+T, 2E+6T, 6T).

At recursion depth k from one triangle they are

    V=(6^k+3*2^k)/2+1,
    E=(3*6^k+3*2^k)/2,
    T=6^k.                                         (8)

At k=1 these are `(7,12,6)`, compared with `(6,9,4)` for E_2. No k>=1
iterated barycentric triangle has the same f-vector as an E_r triangle:
equality of top counts would require k even, k=2h and r=6^h; then the
vertex counts differ by `3(6^h-4^h)/2>0`. At k=2 the comparison is
`(25,60,36)` versus E_6's `(28,63,36)`. Mixed operations and their common
refinements are treated separately; formula (7) concerns only compatible
edgewise operations.

## 5. An exact address that retains the represented point

There is a useful carrier stronger than a mesh-cell label. For any point
with barycentric coordinates lambda and any degree r, put

    q_i=r*(lambda_1+...+lambda_i)=b_i+f_i,
    b_i=floor(q_i), 0<=f_i<1, 1<=i<=d.

Sort the distinct values among `{0,1,f_1,...,f_d}`. For each consecutive
interval `(a,b)`, choose any interior t and form the cumulative integer
vertex `s_i=b_i+1_(f_i>t)`, with last coordinate r. Give its inverse-C
vertex the positive weight b-a. These vertices form a chain in a unit
cube and remain in the ordered chamber. Moreover,

    sum_(intervals) (b-a)*s_i = b_i+f_i=q_i.

Applying C^-1 proves exact reconstruction of lambda after division by r.
Strictly positive weights identify the unique minimal mesh face containing
the point, so boundary points require no arbitrary incident-top-cell choice.
For the hierarchical operation use this face's inherited chain order.
This retains the point even on shared edges and vertices.

The resulting record is

    (degree, ordered mesh-face vertices, positive local weights).

Its source is a simplex point; its target is a geometrically realized
face address with local barycentric data; its decoding is the weighted
sum above. It preserves the exact point and its smallest containing face.
The cell address alone discards the local weights. For example, at r=2
the distinct points `(1/10,1/5,7/10)` and `(1/5,1/10,7/10)` have the same
three-vertex cell and different weights. Integer cumulative coordinates
alone do not create this loss; discarding the weights does.

Two successive records at degrees r and s flatten by weighted affine
composition to the direct degree-rs record. Retaining the two-stage
record additionally remembers the chosen factorization and parent path.
Thus mesh refinement provides a concrete structured-storage operation,
with its reconstruction data stated explicitly. It supplies no general
compression bound for unrelated data and no arithmetic dynamics merely
from resemblance to a triangular table.

## Reproduction and precise finite universe

Run from the repository root:

    python -X utf8 -B 04-computation/experiments/edgewise_cumulative_coordinates_20261004.py
    python -X utf8 -B -O 04-computation/experiments/edgewise_cumulative_coordinates_20261004.py

The [saved output](edgewise_cumulative_coordinates_20261004.out) records:

- Both triangle clique censuses for every r=1,...,12, all unit-transfer
  equivalences, and independent graph distances between every pair.
- Every face count and top-simplex determinant for d=1,...,4, r=1,...,4.
- Full embedded top-cell comparisons for d=1,2,3 and r,s=1,2,3; the
  incompatible tetrahedral child-order hostile at r=s=2.
- Every coordinate permutation for d=1,...,5 at r=1,2,3.
- 4,212 exact point decodes for d=1,2,3, all barycentric numerator
  compositions with denominators 1,...,8, and r=1,...,6; 2,106 two-stage
  comparisons at factor pairs (2,2),(2,3),(3,2), including all boundaries.
- The barycentric recursion counts through depth five and explicit
  literal-edge, ordering, and lost-weight controls.

The graph/clique enumeration, determinant calculation, integer transfer
criterion, independent distance search, permutation checks, and rational
point reconstruction are separate verification paths. All checks remain
active under optimized Python; normal and optimized outputs agree.
