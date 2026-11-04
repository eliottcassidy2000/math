# Interacting unit hubs, complete boundary zeros, and a Moser-carrier barrier

**PROVED:** the two-hub extension criterion, explicit interacting selectors,
triangular-lattice growth barrier, and the reduction/list-colouring arguments
below. **FINITE-EXACT:** complete boundary tables, all256 fixed-bank assemblies,
and the35-triple algebraic certificate supporting the Moser contact lemma.
No new six-chromatic graph or five-colour UNSAT certificate is claimed.

## Inheritance and target

The [preceding ear packet](../../two_port_unit_distance_20261103/geometry/PROOF.md)
reconstructed a whole six-port tensor through actual unit geometry, but every
ear had positive extension weights at palettes5 and6. This investigation
leaves that language by allowing two private vertices to compete for the
same remaining colour. The
[full partition interface](../../boundary_palette_compression_20261011.md),
[six-colour staircase](../../hn_metagraph_information_20261026/six_color/RESULT.md),
and [five-interior paired-selector theorem](../../paired_budget_selectors_matrix_20261028/FIVE_INTERIORS.md)
are inherited. The last result concerns eight paired ports at palette6;
the new fixed-Moser statement concerns seven different ports at palette5.
Neither statement is substituted for the other.

Closest mechanism: complete allowed-colour lists on interacting private
vertices. Hostile: a conditional zero with a three-colourable full graph.
Corrected near miss: positive-ear impossibility is not an impossibility for
all unit-distance gadgets. Least-used sidecar: the exact number of old ports
visible from a new point in a *fixed non-lattice* configuration.
Board: unused colours; private interactions; full boundary support; exact
circle contacts; growth language; a genuine four-chromatic port carrier.

## H1. The exact two-hub criterion

Let x,y be adjacent private vertices, with their other neighbours in an
already properly coloured boundary B. Let A_x,A_y be the colours absent
from their respective boundary neighbourhoods. The exact extension count is

    w=|A_x|*|A_y|-|A_x intersect A_y|.                         (1)

In particular w=0 if and only if one list is empty, or both lists are the
same singleton. If either nonempty list has at least two colours, a colour
different from any chosen colour at the other hub can be selected. This
proves both directions, including the degenerate empty-list cases.

At palette5, when each hub has four boundary neighbours, neither list can
be empty. A zero occurs precisely when both neighbourhoods use the same
four distinct colours: the private edge then demands two different colours
from the one remaining colour. This is a genuine private-interaction zero,
not a one-vertex rainbow5 exhaustion. At palette6 the same degree bound
leaves at least two colours at each hub, so every proper boundary assignment
extends. The private degree threshold agrees with the minimum-degree5
necessity for a six-critical graph, but does not establish such a graph.

## H2. An actual non-lattice unit selector with exactly two forbidden patterns

Take hubs x=(0,0), y=(1,0), and the six ordered boundary ports

    p0=(1/2,sqrt3/2),    p1=(1/2,-sqrt3/2),
    p2=(-1,0),           p3=(-3/5,4/5),
    p4=(2,0),            p5=(8/5,-4/5).                         (2)

All eight points are distinct. Checking all28 pair distances gives exactly

    xy;  xp0,xp1,xp2,xp3;  yp0,yp1,yp4,yp5.                 (3)

In particular the six ports are independent. There are no hidden boundary
edges or cross-unit contacts. Coordinates and exact incidence are saved in
`coordinates.json`; all computations use Q(sqrt3,sqrt11), with zero sqrt11
coefficients here. The point p3 is outside the triangular lattice fixed by
the unit anchors x,y. Thus the new interaction is not dependent on a
triangular-lattice restriction.

At q=5, H1 shows that the only two forbidden admissible equality partitions
on these independent ports are

    (0,1,2,3,2,3),       (0,1,2,3,3,2).                      (4)

Indeed p0,p1 must have different colours; p2,p3 must use the other two
distinct colours in the common four-colour set; and p4,p5 must use those
same two colours in either order. Both assignments in(4) are proper on
the actual boundary, yet force x=y=4. At q=6 each of them has exactly two
extensions, assigning the two unused colours to the two hubs in either order.

The **full** boundary table is retained, not just these selected zeros.
Of203 set partitions,202 are admissible at q=5; precisely200 extend. This
is15385 of the15625 labelled boundary assignments. At q=6 all203 partitions
and all46656 labelled assignments extend. The complete graph has

    Z(q)=q(q-1)(q-2)^2(q-1)^4,
    Z(5)=46080,             Z(6)=300000.                     (5)

The graph is three-colourable: colour x,y differently, colour both common
ports with a third colour, and colour each leaf differently from its hub.
Hence(4) is a true conditional obstruction, while the global five-colour
relation remains nonempty. There is no six-chromatic conclusion.

## H3. All interacting two-contact centres on a fixed lattice boundary

Write alpha=(1+i sqrt3)/2 and represent a+b alpha by(a,b). Use boundary

    B=((0,1),(1,-1),(-1,0),(-1,1),(2,0),(2,-1)).             (6)

The hubs remain(0,0),(1,0). This changes p3,p5 from(2) to the lattice points
(-1,1),(2,-1). Its complete boundary graph has four edges, and both patterns
in(4) remain proper. The two-hub graph has13 unit edges and the same
private-interaction zeros. Its complete q5/q6 boundary supports have64/67
partition types, with14580/122880 full colourings.

There are exactly eight **external real planar points** at unit distance
from at least two ports of(6):

    (-2,1), (-1,2), (0,-1), (0,0), (1,0), (1,1), (2,-2), (3,-1). (7)

This list is not a lattice-only search claim. For each of the15 port pairs,
the checker supplies exactly two candidate circle intersections when the
squared distance is<4, one when it is4, and none when it is>4. Those are the
maximum possible numbers of intersections of the actual two unit circles,
proving completeness in the full plane. The coincident-port case is absent.

Their boundary contact lists, in the order(7), are

    {2,3}, {0,3}, {1,2}, {0,1,2,3},
    {0,1,4,5}, {0,4}, {1,5}, {4,5}.

The entire private graph is a three-edge path on centre indices2,3,4,5,
with the other four centres isolated. Its central edge is the interacting
hub pair. The two path leaves have at most two boundary neighbours and one
private neighbour; at q>=4 they always extend after colouring the hubs.
Each isolated centre also always extends. Consequently adjoining **all**
eight centres preserves exactly the two-hub boundary support at q5 andq6.

Every one of the256 subsets of these centres, with every induced unit edge
retained, is tested. The192 subsets not containing both hubs preserve the
whole boundary support; the64 containing both remove exactly the two
patterns(4) at q5. At q6 all256 have the same boundary support. No assembly
eliminates all states. `bank_supports.json` stores every full support as a
203-bit mask, with the partition order in `boundary_tables.json`.
The full14-point bank has6031260 five-colourings and322775040 six-colourings,
and an explicit three-colouring is included.

## H4. Why iterative rigid growth cannot escape this particular bank's colour cap

Let L=Z[alpha] be the triangular lattice at unit scale. If two distinct points
of L have a common unit neighbour anywhere in R2, that neighbour also lies
in L. Indeed their displacement has squared norm a^2+ab+b^2<=4. Up to lattice
rotation and reflection, its only possible nonzero types are1,1+alpha,2,
of squared norms1,3,4. Their common unit neighbours are respectively
{alpha,1-alpha}, {1,alpha}, and{1}, translated with the first endpoint.

It follows inductively that starting within L and adding one new point at
a time, each with at least two unit contacts to distinct already placed
points, can never leave L. Every resulting induced unit graph has the
three-colouring a+b alpha -> a+2b mod3. Additional unit edges are included;
all six possible unit displacements have nonzero colour difference.

This is a **growth-language barrier**, not a universal barrier for unit
graphs. A first new point with only one old contact permits a free angle,
as p3 in(2) demonstrates. Simultaneous multi-vertex gluing without such an
ordering is also outside the theorem. Escaping L removes this particular
barrier but, as H2 shows, does not itself create six-chromatic forcing.

## H5. Exact contact exhaustion at a genuine four-chromatic carrier

Use the seven Moser-spindle points

    A=(0,0), B=(sqrt3,0), C=(sqrt3/2,1/2), D=(sqrt3/2,-1/2),
    B'=R(B), C'=R(C), D'=R(D),
    R(s,t)=((5s-sqrt11*t)/6,(sqrt11*s+5t)/6).               (8)

All21 pairs are checked exactly. The unit edges are the two diamonds
AC,AD,BC,BD,CD and AC',AD',B'C',B'D',C'D', plus BB'. In any three-colouring,
the first diamond forces A=B and the second forces A=B', contradicting BB'.
The supplied four-colouring certifies chromatic number exactly4. Independent
literal counts are Z(3)=0, Z(4)=384, Z(5)=5040. This is the familiar spindle
used as a verified carrier, not a newly claimed four-chromatic construction.

**Contact lemma (complete exact algebraic certificate).** Every point outside
these seven points is at unit distance from at most two of them.

For any three noncollinear points there is at most one common unit-circle
centre, their circumcentre. Three distinct collinear points have none.
The35 triples in(8) are exhausted in Q(sqrt3,sqrt11). Exactly10 have
circumradius1; all10 centres are existing spindle vertices. Their triples,
using index order(A,B,C,D,B',C',D'), are

    012,013,045,046,156,234,235,236,256,356.

The other25 have no common unit centre. `spindle_contacts.json` retains
every triple, its exact centre when present, and a second certificate:
for d=P1-P0, e=P2-P0 and D=det(d,e), unit circumradius is equivalent to

    |d|^2 |e|^2 |d-e|^2 - 4D^2 = 0,      D!=0.            (9)

The checker compares this multiplication-only test against solving the
two perpendicular-bisector equations, then checks that each surviving
centre is already an old vertex. Coefficients in the linearly independent
basis1,sqrt3,sqrt11,sqrt33 decide zero exactly; no floating-point tolerance
or sampled plane search enters. This finite certificate proves the contact
lemma for all real exterior points because any third contact would belong
to one of the exhausted triples.

## H6. All five-colour assignments extend across six or fewer private points

Keep the fixed spindle points(8) as boundary. Add any finite set I of new
distinct planar points, retaining the **full induced unit graph**. Every
private vertex has at most two boundary neighbours by H5. Thus every proper
five-colouring of the boundary leaves a list of at least three colours at
each private vertex.

If the private graph is2-degenerate, greedily colour it in reverse of an
ordering that removes vertices with degree at most2. At most two already
coloured private neighbours exclude colours at each step. Therefore:

**Extension theorem.** Every proper five-colouring of this fixed spindle
boundary extends whenever its private unit graph is2-degenerate.

The following elementary refinements cover every |I|<=6.

1. A unit graph on at most five vertices is2-degenerate. A minimum-degree3
   subgraph on four vertices is K4, impossible in the plane. On five vertices
   its complement has maximum degree1, hence is a matching. A matching of
   size0 or1 leaves a K4; size2 leaves K1 join C4. But the unit-neighbour
   graph on a circle is a union of subgraphs of regular C6 orbits: adjacency
   changes angle by exactly60 degrees. It has no C4, excluding that wheel.
2. A six-point unit graph has at most nine edges. With maximum degree<=3
   this is immediate. A degree5 vertex has five neighbours on its unit
   circle, among which there are at most four unit edges. A degree4 vertex
   has at most three unit edges among its four circle neighbours, while
   the remaining vertex has at most two contacts to those neighbours by
   the common-circle-intersection bound. Both cases also give at most9.
3. Consequently a six-vertex private3-core is cubic. Its complement is a
   simple2-regular graph on six vertices, either C6 or two disjoint C3s.
   The original is respectively a triangular prism or K3,3. The latter
   cannot be a unit graph: two vertices would have three common unit
   neighbours. Thus only the prism remains.
4. The triangular prism is3-choosable by an elementary argument. Restrict
   all lists to exactly three colours. The first triangle has at least
   3*2*1=6 proper list colourings. Remove each chosen colour from the list
   at its matched vertex of the second triangle. These three lists still
   have size at least2; they fail to admit three distinct representatives
   only when all three are the same2-set S. Such failure occurs for at most
   one first-triangle colouring: its three distinct colours make the
   original opposite lists pairwise distinct, S is their common pairwise
   intersection, and each chosen colour is the unique element outside S.
   At least five of the initial choices therefore extend. This proof is
   independent of the size of the common palette.

It follows that **every proper five-colouring of the fixed seven spindle
ports extends over any six or fewer private planar points**. A construction
on this carrier that blocks even one such boundary assignment therefore
requires at least seven private vertices and a nonempty private3-core.
The statement retains arbitrary private cross-contacts; they are exactly
the private graph tested in the argument. It is not restricted to a disjoint
union of two-hub fragments, and does not bound arbitrary six-chromatic unit
graphs that lack this prescribed boundary geometry.

The surviving prism is itself geometrically real: translate an equilateral
unit triangle by(3/5,4/5). All15 pairs of the resulting six points give
exactly the nine prism edges. As a separate finite control, the checker
normalizes one three-element list from a five-colour palette to{0,1,2}
and exhausts the100000 choices for the other five lists. All extend. Palette
permutation covers all10^6 unnormalized choices, and larger lists only help.
This corroborates item4; the elementary proof is the unbounded statement.

## Verification and the remaining forcing obligation

Run `python3 verify.py` and `python3 -O verify.py` in this directory. Ordinary
verification compares frozen coordinates, all203 boundary tables, all256
complete bank supports, and the35-triple certificate without overwriting them.
`--write` is explicit artifact-production mode. Checks are active under -O.

The boundary engine is list-colouring dynamic programming on the actual
private forest. Literal enumeration independently checks both two-hub
tables; a separate hub/leaf factorization checks the full bank. Independent
literal graph counting checks the spindle; the inverse and Gram methods
cross-check its geometric contact exhaustion. All100000 prism list witnesses
are verified edge by edge and included in a reproducible rolling hash.

Map/loss ledger: boundary colourings -> allowed lists -> private list
colouring is exact and retains their common palette. A conditional zero
loses all other boundary states if reported alone; the complete tables
restore them. Contact data -> list-size bounds is a one-way relaxation,
so a successful extension proves a no-go, while a failed relaxed test would
still need an actual placement. A private3-core is necessary on the spindle
carrier, not sufficient for forcing; the realized prism explicitly shows
that distinction.

The next genuine forcing candidate must preserve its exact geometry and
whole q5 boundary relation while escaping these proved closures: for the
fixed spindle carrier it needs at least seven private vertices, private
minimum-degree structure surviving peeling, and lists that actually conflict.
Neither the two real zeros(4), a larger lattice bank, nor a three-core by
itself supplies a six-chromatic graph. That final construction remains OPEN.
