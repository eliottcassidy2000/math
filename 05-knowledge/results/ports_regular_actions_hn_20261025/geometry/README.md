# What a regular eight-point action can preserve in the Euclidean plane

2026-10-25. **PROVED elementary** statements below; **INHERITED PROVED**
K3 and root-cover inputs are linked explicitly. No claim of historical
novelty, six-chromatic unit-distance construction, or tournament connectivity.

## 1. Inheritance: the free action really is detected on eight conics

The [audited branch-code criterion](../../six_seven_geometry_20261025/fano_codes/RESULT.md),
section6a, gives the precise carrier. Label the seven branch lines by the
seven nonzero vectors of F2^3, with central label c and opposite pairs
(v1,v2),(v3,v4),(v5,v6). Put

    d1=v1+v2, d2=v3+v4, d3=v5+v6.

The simplex subgroup element a translates the eight conic labels by
(a.d1,a.d2,a.d3). Their sum is c. A crossing is bad iff one d_j=c,
equivalently the other two differences coincide. Thus freeness on the
resolved seven-root hourglass is equivalent to d1,d2,d3 forming a basis,
which is equivalent to a regular action on the eight conic components.
Ranks3,2,1 give respectively one orbit of8, two orbits of4, four orbits of2.
This is a group action on the component SET, not a free action on every
point of the K3 and not a Euclidean isometry action.

The component intersection matrix, proved in the
[root-deletion packet, section4](../../six_seven_geometry_20261025/root_deletion/DELETION_AND_CUBE.md),
is

    M=-2I+2A(Q3).

Its spectrum is4,0,0,0,-4,-4,-4,-8. The three singleton Walsh characters
are actual divisor relations: the two coordinate half-cubes are equivalent
hyperplane sections. The other four nonconstant characters are negative
directions, not additional null relations. The next sections construct an
explicit connection from this relation space to plane unit distances.

Closest mechanism: a translation action with its actual direction labels.
Canonical hostile: regular action alone does not retain geometry or colour
extension. Corrected near miss: a doubled intersection edge is not already
a Euclidean unit segment. Least-used sidecar: all square-face equations.

## 2. Every injective plane unit quadrilateral is a parallelogram

Suppose p0,p1,p2,p3 are four distinct points in R2 and every consecutive
edge of their four-cycle has length one. Both p1 and p3 belong to the
intersection of the unit circles centred at the distinct p0,p2. Since the
intersection contains two distinct points, those are its two points,
symmetric about the midpoint of p0,p2. Therefore

    p0+p2=p1+p3.                                      (1)

This does not require a noncrossing graph drawing. It does require distinct
vertex images and the Euclidean plane. A unit regular tetrahedron in R3
provides a four-cycle violating(1); silently transferring that geometry
would be a false dimension change.

## 3. Cubelike unit-distance graphs are disjoint hypercubes

Let V=F2^r and let S be a set of distinct nonzero vectors. Let Cay(V,S)
have edges x--(x+s), for x in V,s in S. An embedding here is injective
on vertices, with every graph edge represented by a unit segment; crossings
are allowed and no geometric symmetry of the drawing is assumed.

**Theorem U.** Cay(V,S) admits such an embedding in R2 if and only if S
is linearly independent. In that case its connected components are
hypercubes Q_|S| and it is bipartite when S is nonempty. In particular a
regular elementary-2 translation model by itself cannot force six or
seven colours in a plane unit-distance graph.

**Necessity.** For distinct a,b in S, the four distinct vertices
x,x+a,x+a+b,x+b form a unit quadrilateral. Equation(1) gives

    Delta_a f(x+b)=Delta_a f(x),
    Delta_a f(x)=f(x+a)-f(x).                          (2)

If a lies in the span of S minus a, repeatedly apply(2) along a word in
those other generators to get Delta_a f(x+a)=Delta_a f(x). But translation
by a twice is the identity, so the left side is -Delta_a f(x). Over R this
forces Delta_a f(x)=0, contradicting its unit length. Hence every a is
outside the span of the other generators and S is independent.

**Sufficiency.** For independent S, identify each component with {0,1}^m.
Choose unit vectors u1,...,um in the plane and map epsilon to sum epsilon_i*u_i.
They can be chosen so that all subset sums are distinct and no undesired
pair is at distance one. Inductively, when adding a new vector u on the
unit circle, collisions require u to equal one of finitely many existing
differences. An unwanted unit pair between the old copy and its translate
requires |u+d|=1 for a fixed nonzero old difference d; its intersection
with |u|=1 has at most two choices. Avoid this finite set. The infinitely
many rational points on the unit circle even allow rational coordinates.
Different components can be translated far apart. This constructs an
induced unit-distance realization as well. QED.

The ambient group may have arbitrarily large rank. It is the dependence
among chosen directions, not the number of abstract group elements, that
prevents the embedding. Conversely, the graph action need not extend to
Euclidean isometries of the chosen drawing.

## 4. The exact loss is visible in Fourier coordinates

For an arbitrary S, call a scalar function f:V->R square-flat when
Delta_a Delta_b f=0 for every distinct a,b in S. Put H=span(S), t=dim H,
and let b be the number of coloops in S: those a outside span(S minus a).
Then the square-flat function space has dimension

    2^(r-t) (1+b).                                    (3)

On one component, its basis is the constant and one character for each
coloop, negative on that generator and positive on all other generators.
**Proof.** A Walsh character chi has Delta_a Delta_b chi=0 unless both
chi(a) and chi(b) are -1. Thus an allowed character is negative on at most
one member of S. It is either constant on H or singles out a coloop. The
converse is immediate; the choices on different H-cosets are independent.
All noncoloop generator directions are collapsed by every square-flat map.
This gives an exact quotient-loss statement, stronger than a count of modes.

For Q3, precisely the three coordinate characters survive after centering:

    f(epsilon)=sum_(j=1)^3 a_j (-1)^epsilon_j,
    a_j in R2.                                       (4)

Every injective unit-distance embedding has this form, by(1) and(3).
Every edge has length one iff |a_j|=1/2 for each j, and injectivity requires
no nonzero {-1,0,1}-combination of the three a_j to vanish. Consequently
the two centered coordinate functions of the drawing lie exactly in ker M
from section1, with these explicit norm and noncollision conditions.

This is the promised typed bridge: the three K3 divisor-relation coefficient
vectors supply the coordinate function space for all plane cube embeddings.
It transports the cube adjacency, while forgetting intersection multiplicity,
curve classes, complex structure, and the other four negative Walsh modes.
It does not embed the K3 into the real plane or identify its intersection
form with a Euclidean Gram matrix.

## 5. Eight vertices: exact colour and geometry boundaries

For completeness, every graph Cay(F2^3,S) has chromatic number1,2,4 or8.
S empty gives1, and S all seven nonzero vectors gives K8. Any missing
direction d partitions the vertices into four independent pairs
{x,x+d}, proving chromatic number at most4. If S contains a,b,a+b,
the subgroup they generate induces K4, giving4. Otherwise S is a Fano
cap, of size at most4 and contained in an affine hyperplane avoiding zero;
that hyperplane defines a bipartition. The small cap assertion is also
elementary: an independent triple can be sent to e1,e2,e3, after which
the only possible fourth element is e1+e2+e3.

Triangle-free does not suffice for a plane unit-distance embedding:
S={e1,e2,e3,e1+e2+e3} gives K4,4, which is bipartite but dependent.
Theorem U excludes it. Thus even the correct abstract chromatic number
does not restore the discarded square-face geometry.

## 6. The incoming norm-seven torus has an explicit continuation to the cube

The concurrent [torus-to-Campedelli proof](../../session61_hadwiger_six_seven_20261025/torus-fano/TORUS_FANO_CAMPEDELLI.md)
arrived at commit47eea7cad5 during this session. The two Eisenstein quotients
`a+2b mod7` and `a+4b mod7`, together with the two triangular face types,
select exactly the four-element orbit among the sixteen free codes.
Here is its explicit continuation through the regular-action criterion.

For the first quotient, the opposite branch pairs have residues(1,6),
(2,5),(3,4). The two face systems have binary markings, in residue order,

    U: (1,2,4,3,6,7,5),   D: (1,2,4,5,7,3,6).

Their central label is1 and their three difference vectors are respectively
(7,3,5) and(4,7,2), each a basis. For the conjugate quotient, in the same
ordered geometric directions the differences become(7,5,3) and(4,2,7).
Thus all four torus-selected codes act regularly on the conic cube. This
can be checked directly without re-enumerating the thirty Fano structures.
The inherited classification is still needed to identify their complex
surface isomorphism class; the basis test by itself does not do that.

The maps preserve the chosen face incidence, branch crossing exclusions
and marked translation directions. They are not homomorphisms from the
additive group F7 to F2^3: residue0 maps to a nonzero branch label. Nor does
the torus K7 give a planar unit K7. Its triangular universal cover is
three-colourable; a colouring need not descend through the lattice quotient.
This loss of global compatibility complements the square-face obstruction
in Theorem U and the retained boundary relations in the palette lane.

## 7. What can move the open problems forward

To use the conic cube for Hadwiger–Nelson, a construction must break the
global translation rule, retain boundary lists or add attachments whose
geometry is proved. The [path-list interface](../hn/README.md) takes this
route and gives a different exact6/7/8 mechanism. For tournament repairs,
the analogous requirement is to retain each directed path and its actual
reversal; a graph-chart action alone does not preserve blackness. For
counting, the [future-response state](../counting/README.md) retains the
queries that remain affordable instead of treating a quotient count as
complete arithmetic information.

Literature status was refreshed after the incoming commit. De Grey's
[five-chromatic construction](https://arxiv.org/abs/1804.02385) and
[Haugland's construction](https://arxiv.org/abs/2608.04542) are finite
five-chromatic controls, not the latest unrestricted bound. The incoming
[source audit](../../session61_hadwiger_six_seven_20261025/SOURCE_STATUS.md)
cites the September2026 six-colour lower bound. We separately read the
primary [Lean root theorem](https://github.com/openai/math/blob/adc7f1241b42e322a6451854ab7e4b4c146bf78a/lean/OAI/Geometry/PlaneColoring/Five.lean)
`no_proper_five_coloring`, but did not build it or audit its transitive
axioms. That result has status **CITED** here. None of the elementary
deductions in this packet depends on adopting the external bound, and
we make no literature-priority claim.

Reproduce `python check.py` and `python -O check.py`. The checker considers
all128 direction sets on F2^3, independently solves their square equations
and chromatic numbers, tests exact rational cube coordinates and all5040
Fano labellings' action ranks. These finite controls supplement the proofs;
they are not an infinite-plane colouring certificate.
