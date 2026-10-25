# The projective eight behind the hexagonal seven

**Status: CITED classical construction; PROVED explicit transport and loss
statements; FINITE-EXACT group and incidence controls.** No new plane bound,
surface automorphism or literature-wide novelty is asserted.
[Independent audit: PASS](../audits/PROJECTIVE_AUDIT.md).

## Inheritance and source

[Session61's torus/Fano map](../../session61_hadwiger_six_seven_20261025/torus-fano/TORUS_FANO_CAMPEDELLI.md)
already gives the norm-seven quotient of the Eisenstein lattice and its
two face designs on F7. The selected design is

    U_t={t,t+1,t+3},   t in F7.

Its marking by the nonzero binary vectors is

    lambda(0),...,lambda(6)=(1,2,4,3,6,7,5).

Binary integers here encode three coordinates, least significant first.
Neither this map nor its extension is an additive map between fields.
The [earlier cardinality packet](../../cardinality_bridges_20261018m.md),
section5, already identifies affine eight-point structures, extended
Hamming codes and their fourteen planes. The
[incoming metric packet](../../ports_regular_actions_hn_20261025/geometry/README.md),
section6, already computes the conic direction basis (7,3,5).

The missing connection is classical: Caldero–Germoni,
[*Equilateral Triangles and the Fano Plane*](https://kam.mff.cuni.cz/~spring/2017/papers/equilateral_triangles_in_Fano.pdf),
American Mathematical Monthly123(8),789–801(2016), Lemmas3–4 and
Proposition9, extends the triangles to a Steiner system on P1(F7).
This source was read at the statement and proof level. We rederive the
finite design below and apply it to the marked conic component set.

## P1. From seven triangles to eight affine points

Adjoin infinity and set lambda(infinity)=000. Form fourteen four-sets:
the seven {infinity} union U_t and their complements in P1(F7).
Under lambda these are exactly

    {x in F2^3 : s.x=b},  s nonzero, b in F2.                 (1)

Indeed each U_t has xor zero, so its three nonzero vectors together with
zero form a two-dimensional subspace; all seven such subspaces occur.
The complements are their seven other affine cosets. Every three distinct
binary points determine a unique affine plane: translate the first to
zero, and the other two are distinct nonzero, hence independent. Thus the
fourteen blocks form S(3,4,8), and deriving at infinity recovers U.

Projectively, these blocks are one PSL2(F7) orbit of {infinity,0,1,3}.
Equianharmonic means cross-ratio in {3,5}; all28 such four-sets split into
two disjoint14-block systems. The other orbit derives to the other torus
face design. This is Proposition9 of the cited paper, also independently
checked here over every four-subset and projective transformation.

**Hostile / needed sidecar.** Taking all equianharmonic four-sets gives
two completions of every triple, not a Steiner system. The choice of
orbit (the original triangle type) must survive extension to eight.

## P2. An explicit map to the eight actual conic components

For the hourglass root labels colored by c2(a,b)=a+2b mod7, the three
opposite pairs are (1,6),(2,5),(3,4). Their binary differences are

    d=(lambda(1)+lambda(6),lambda(2)+lambda(5),
       lambda(3)+lambda(4))=(7,3,5).

They form a basis. In the inherited root-cover theorem a simplex deck
element a translates the conic label by (a.d1,a.d2,a.d3).
Choose the dot-product identification of the inertia space with its dual,
and choose a reference conic. These are additional choices, not a canonical
duality. Then the marked correspondence is

    t -> lambda(t) -> e(t)=(lambda(t).7,lambda(t).3,lambda(t).5),

    t:     0  1  2  3  4  5  6  infinity
    e(t):  7  3  5  4  6  1  2     0.                       (2)

The conic C_e is the actual component with equations

    A_-=i(-1)^e1 A_+, B_-=i(-1)^e2 B_+,
    C_-=i(-1)^e3 C_+, C_+^2=A_+^2+B_+^2.

The invertible linear map in(2) takes the fourteen blocks in(1) to all
fourteen affine half-sets of the conic labels. Complementary block pairs
correspond precisely to the seven nonconstant Walsh characters, up to
their sign. This transports incidence on the component SET. It supplies
no map between a projective line over F7 and a complex K3 surface.

For the K3 the intersection matrix is M=-2I+2A(Q3), and its Walsh
eigenvalues are 4(1-|s|). Exactly three directions, s=1,2,4, give the
actual half-cube hyperplane relations. The other four nonconstant
directions have negative eigenvalues. Therefore the Steiner structure
alone does not recover the divisor relations: it treats all seven
parallel classes symmetrically. A selected cube basis restores that data.

For example, the shear (e1,e2,e3)->(e1+e2,e2,e3) preserves all fourteen
planes but sends an intersecting pair C_000,C_010 to the disjoint pair
C_000,C_110. This is a minimal explicit failure of incidence transport
from the affine design to the surface's curve intersection graph.
The [palette lane](../fano_eight/RESULT.md) computes the full lost-basis
fiber and its significance for six paired colors.

## P3. An affine action can preserve every block and have no origin

This is a further exact warning relevant to frame gluing. In the lambda
coordinates, the projective transformations T:z->z+1 and S:z->-1/z act as

    T(x)=L_T x,        columns(L_T)=(2,4,3),
    S(x)=L_S x+1,      columns(L_S)=(1,4,2).                  (3)

Thus S exchanges the last two bits and flips the first. It has no fixed
binary point, just as z^2=-1 has no solution in F7. Both actions preserve
the fourteen affine planes. The matrices T,S generate PSL2(F7): powers
of T give all upper unipotents and conjugation by S gives the lower ones;
elementary row operations generate SL2(F7), modulo its scalar center.
Consequently the whole168-element projective action is affine on this
eight-set. It is transitive, and has no common fixed origin.

Write g(x)=L_g x+a_g. Then

    a_(gh)=a_g+L_g a_h.

Changing the chosen origin to o removes all translations iff
a_g=(I+L_g)o for every g, equivalently iff there is a common fixed point.
That fails here already for S: its translation1 is outside im(I+L_S).
The action is therefore a concrete nontrivial crossed cocycle, not merely
an abstract group-isomorphism slogan. No characteristic-zero or continuous
cohomology assertion is intended. Its168 distinct linear parts, and the
six elements preserving the standard cube graph, are FINITE-EXACT checks.

The abstract exceptional isomorphism PSL2(F7)~=GL3(F2) does not identify
this transitive eight-point action with the usual point-fixing linear
action. Nor does preserving the affine design make all168 elements
automorphisms of the marked conic intersection graph.

## P4. A cheap hostile to confusing two twenty-eights

The prior root-character formula gives p_g=28+3=31 for eight distinct
spanning branch lines with no four concurrent. Its28 one-dimensional
six-root character labels are the six-subsets, equivalently the omitted
pairs of the eight labels. They are not the28 equianharmonic quadrangles
in the projective completion above.

There is no PSL2(F7)-equivariant bijection between these two28-sets.
Indeed z->4z has two fixed points,0 and infinity, and two3-cycles on the
remaining six points. It fixes exactly one pair, {0,infinity}, hence one
six-subset. It fixes four equianharmonic four-sets: each3-cycle together
with either fixed point. Their fixed-point counts differ, obstructing
even an isomorphism of the corresponding rational permutation modules.
This is a formal action on branch labels; no PSL action on a fixed
eight-line arrangement or its surface is assumed. The repaired connection
is the marked affine eight of P2, with its explicit choices and losses.

## Connection contract, validation and stopping boundary

Source: the marked norm-seven triangular quotient and one face class.
Target: the affine block system on the eight conic labels. Map: extend to
projective quadrangles, apply lambda, chosen duality and the conic-action
matrix. Preserved predicate: block membership and complementary Walsh
partitions. Lost: Euclidean distances, chosen cube directions, intersection
multiplicity, curve classes and complex geometry. Required sidecars:
chirality, reference conic, duality, ordered root pairs and original M.

The checker covers all70 four-subsets,336 projective transformations and
168 square-determinant transformations, verifies each of56 triple
completion counts, both chiral designs, the exact map(2), eight Walsh
eigenvectors and five hostile controls. It uses only integer arithmetic
and explicit failure gates. Run from this directory:

    python3 -B check.py
    python3 -B -O check.py

The construction merges old HN and surface themes at a concrete finite
carrier. It does not give a six-color extension theorem, a seven-chromatic
unit graph, or a new action on the complex surface. The next geometric
obligation is to retain actual boundary color relations and metric seams;
the full affine symmetry cannot replace them.
