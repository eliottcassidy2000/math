# The seven-color lattice quotient selects the exceptional Campedelli orbit

**Status:** explicit mathematical deduction and exact finite certificate;
[independent audit passed](../audits/lattice/AUDIT.md). Classical torus/Fano ingredients, new-to-route
application to the Session60 hourglass. No plane chromatic lower bound.

## 1. Keep the objects and markings

Let omega^2+omega+1=0 with positive imaginary part, Lambda=Z[omega],
and Q(a,b)=a^2-ab+b^2. The six unit directions are
R={+/- (1,0), +/- (0,1), +/- (1,1)}. Identify the seven *dual labels*
of the hourglass branch lines E,E+/-x,E+/-z,E+/-(x+z)
with {0} union R. This does not identify the primal central branch E=0
with a point in the affine hexagon: that branch is the line at infinity.

The ideal I=(2-omega) has index seven. The map
c_2(a,b)=a+2b mod7 identifies Lambda/I with F7. Its values on the
seven labels are all distinct. Use those values as the branch labels.
The three actual triple crossings are then
{0,1,6}, {0,2,5}, {0,3,4}.

## 2. Retaining triangles creates a torus, not a unit K7 in the plane

Give Lambda its usual triangular tiling, with triangles
{v,v+1,v+1+omega} and {v,v+omega,v+1+omega}.
The translation quotient by I is a triangulated torus. Its vertices are
F7, its graph is K7, and its two types of triangles are

U_t={t,t+1,t+3},  D_t={t,t+2,t+3},  t in F7.

There are seven vertices, twenty-one edges and fourteen triangular faces.
Every edge is incident to two triangles; each vertex link is a six-cycle.
Orientation is inherited from the plane because the deck group consists of
translations. Thus it is a closed connected orientable surface of Euler
characteristic zero, independently confirming that it is a torus.

The source triangular graph is three-colorable, while its seven-vertex
quotient is K7. The quotient is locally bijective on neighbors (a graph
covering), but is not a chromatic equivalence. Arbitrary source colorings
need not be I-periodic. This is a particularly explicit hostile control:
even a covering map with surjective vertices and edges cannot pull a
chromatic lower bound backwards. The quotient torus is not embedded as a
unit-distance K7 in the Euclidean plane.

## 3. Two genuine Fano structures, after adding the quotient and triangles

Each of U={U_t} and D={D_t} is a Steiner triple system on F7.
Indeed the ordered differences of {0,1,3} are exactly F7 minus {0},
each once; the same is true of {0,2,3}. Therefore every unordered pair
lies in exactly one triple of each system. Each is a Fano plane;
U and D are disjoint as triple sets.

An explicit identification of U with nonzero vectors of F2^3 is, in
residue order 0,1,...,6, the list of binary integers
(1,2,4,3,6,7,5). A triple is a line iff its three labels xor to zero.
For D one may use (1,2,4,5,7,3,6).
These are markings, not additive maps F7 -> F2^3. In particular residue0
maps to a nonzero inertia vector, not to binary0.

This gives a precise repair of the old unqualified assertion that the
centered hexagon 'is' the Fano plane. The bare seven-set or its Euclidean
collinearity does not carry this incidence. The lattice quotient, triangular
faces, a choice of triangle type, and the marking are essential sidecars.

## 4. Transport to the hourglass quotient geometry

Neither U nor D contains any of the three actual branch-crossing triples.
For U, the triples through0 are {0,1,3},{0,2,6},{0,4,5}; for D they are
{0,1,5},{0,2,3},{0,4,6}. Hence their Session60 dependent-crossing count is d=0.
The sum of the seven nonzero F2^3 vectors is zero, so each marking defines
a surjective map from the projective sign deck group
F2^7/<(1,...,1)> to F2^3. Its kernel has order eight. By the independently
audited Session60 free-action classification, each kernel acts freely on
the resolved hourglass and gives a numerical Campedelli surface,
p_g=q=0, K^2=2.

This is a substantive transfer: the torus seven-color quotient has supplied
actual free quotient codes of the K^2=16 hourglass. It is not a map from a
plane coloring to the surface, nor does it preserve a coloring obstruction.
The geometric consequence uses the prior quotient theorem, not merely the
cardinality seven.

## 5. Both chiral ideals and both face types give exactly four codes

The conjugate ideal gives c_4(a,b)=a+4b mod7. Pull back each ideal's two
triangle-type Fano systems to the same seven branch labels. These yield
four distinct Fano structures. This is the complete orbit under the
hourglass's realized D12 symmetry group, not all sixteen free codes.

One may check the action with
r(a,b)=(a-b,a) (60-degree rotation), s(a,b)=(a-b,-b) (conjugation).
For c_2, r acts on residues by multiplication3 and interchanges U and D;
r^2 preserves each. Conjugation exchanges I and its conjugate and cannot
act on the single fixed quotient Lambda/I. On the marked seven-set it
still gives a permutation, carrying the first pair of designs to the
second pair. The four-design orbit has stabilizer C3 in D12.

Session60 proves that the d=0 code orbits have sizes4 and12, and that these
are exactly the two complex Campedelli isomorphism classes. Therefore this
construction selects precisely the orbit4 class, with full automorphism
group C2 x A4 of order24. The other twelve free codes cannot arise from
these two chiral index-seven ideal quotients and their triangular face types.
This exclusion is scoped to that construction, not to arbitrary torus maps.

Central inversion r^3 interchanges U and D. It fixes one branch label and
swaps three opposite pairs. It preserves no Fano incidence on seven labels,
as Session60 already showed. There is no conflict: the construction yields
a pair of exchanged structures, not one inversion-invariant Fano plane.

## 6. Reproduction and boundaries

Run `python3 check_torus.py` and `python3 -O check_torus.py`.
The finite universe is all5040 markings of seven labels, all30 Fano
structures, all12 realized D12 permutations, and the exact14 listed faces.
Checks include pair multiplicities, links, boundary ranks, explicit binary
markings, free-code membership, orbit uniqueness, and hostile controls.

Inherited geometry: Session60 `seven-lines/RESULT.md`, sections3-4,
and `seven-lines/check.py`, at commit a144b255fce03bc2afb4414cd9804584670e3c5a.
The classical triangle/Fano literature is recorded in [SOURCE_STATUS.md](../SOURCE_STATUS.md).
No new K3 Picard rank, integer-square existence claim, arbitrary plane
six-coloring, or seven-chromatic unit-distance graph is obtained here.
