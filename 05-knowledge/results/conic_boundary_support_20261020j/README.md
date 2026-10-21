# Ellipse packets preserve six-colour boundary support with two extra exclusions

2026-10-20. **PROVED**, using the **CITED** degree-choosability theorem:
the unit graph of any finite set on one nondegenerate real ellipse is
4-choosable. More strongly, every proper old six-colouring extends when
each new vertex sees at most two distinct old colours, with at most two
exceptions that may each see three. The number of new vertices is
unrestricted. **FINITE-EXACT:** the sharp four-neighbour example and the
bounded abstract-graph audit below. **OPEN:** six versus seven colours for
the plane; this packet supplies an exclusion test, not an equality tester.

## 1. Inheritance and the connection being tested

The closest proved mechanism is the boundary-list reduction and Gallai-tree
argument of [THM-4661, odd-degree carriers](../../../01-canon/theorems/THM-4661-odd-degree-carriers-a-failed-six-colour-exterior-has-at-least-twenty-vertices-and-equality-testers-are-seven-colour-certificates.md).
The present result has no number-field or triangle-free premise; its
replacement is the internal degree bound from a single ellipse. The
corrected near miss is transferring that theorem's carrier-specific
twenty/twenty-two-vertex bound to unrelated placements.

The hostile example is the proper twelve-port colouring of the nineteen
points in [exact palette forcing](../kernels_games_detours_20261020f/forcing.md):
all seven interior lists become the same three-set, and the spindle blocks
extension. Its uncoloured graph is still four-chromatic. The less-used
sidecar is the pairwise common-neighbour limit: two distinct plane points
have at most two common unit neighbours. It upgrades the basic ellipse
extension lemma to tolerate two extra exclusions.

The live concepts are actual boundary lists, finite critical interiors,
intersection degree, pairwise contact limits, and synchronized area clocks.
The source is a finite geometric exterior with its named old neighbours;
the target is a list-colouring instance. The map deletes the old graph
after forming each exact available list. It preserves extendibility of
that specific old colouring, and loses the old coordinates and colours
unless their contact-and-list sidecar is retained. Merely recording the
number of vertices or a cyclic tournament loses the required lists.

## 2. The geometric degree bound is sharp

Let E be a nondegenerate real ellipse, including a circle, and let B be a
finite set of distinct points of E. For p in B, its unit neighbours are
the intersections of E with the circle of radius one centred at p.

If E is not a circle, two conics with no common component have at most
four intersections. Here this can also be seen without a general
intersection theorem: use the rational parametrization

    (x,y)=(a(1-t^2)/(1+t^2), 2bt/(1+t^2))

after an isometry, substitute into |(x,y)-p|^2=1, and homogenize. The
result is a nonzero homogeneous quartic in the parameter pair, so it
has at most four distinct real projective roots. If it were identically
zero, E would be contained in the unit circle about p and hence would
equal it, impossible since p belongs to E and not to that circle.
The projective parameter includes the point missing from the affine
chart. If E is a circle, its intersection with that unit circle has at
most two points; the circles cannot coincide for the same reason.
Thus the unit graph on B has maximum degree at most four.

Four is attained on x^2+16y^2=1. The vertex p=(0,1/4) has the four distinct
unit neighbours

    (+/-2sqrt(6)/5, 1/20), (+/-2sqrt(2)/3, -1/12).

The checker uses exact arithmetic in Q(sqrt(2),sqrt(3)), confirms every
ellipse equation, and checks all ten pair distances. This induced graph
is the four-leaf star. The example establishes sharpness of the internal
degree bound, not sharpness of the list-colouring theorem below.

## 3. Two three-lists among arbitrary four-lists

**Lemma.** Let G be a finite simple graph with maximum degree at most four,
no K4, and at most two common neighbours for every pair of distinct
vertices. Give each vertex a list of at least four colours, except for
at most two vertices whose lists have at least three. Then G is
list-colourable.

**Proof.** Suppose a vertex-minimal non-list-colourable induced subgraph
H exists. Every v has degree_H(v)>=|L(v)|. All ordinary vertices therefore
have degree and list size exactly four. Call these vertices light, and
the at most two exceptions special. If no light vertex remains, H has
at most two vertices with lists of size at least three and is colourable.

Take a connected component C of the light-induced subgraph. Minimality
colours H-C. Deleting the colours used by neighbours outside C leaves
each vertex a list of size at least its degree in C. Since C cannot be
coloured from these residual lists, the classical degree-choosability
theorem says that C is a Gallai tree: every block is complete or an odd
cycle. Since G has no K4, its complete blocks have size at most three.

A singleton C cannot have degree four using only the two special
vertices. If C has one block, an edge has insufficient degree even after
both specials are added. In a triangle or odd cycle every vertex has
internal degree two and must be adjacent to both specials to reach four.
This gives at least three common neighbours, impossible.

If C has multiple blocks, it has at least two leaf blocks. A non-cut
vertex in an edge leaf again has insufficient degree. Every triangle
or odd-cycle leaf has at least two non-cut vertices of internal degree
two; each must see both specials. Two leaf blocks therefore give at
least four distinct common neighbours, again impossible. This exhausts
C and proves the lemma. QED.

The cited input is the Borodin--Erdos--Rubin--Taylor degree-choosability
theorem, already inherited by THM-4661. Its precise statement is also
[Theorem A in Cranston and Rabern, Beyond Degree Choosability](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v24i3p29/pdf/).
No assertion is made that two is the largest admissible number of
three-lists. The common-neighbour argument gives an explicit sufficient
certificate, not a classification of all failed lists.

## 4. Whole boundary support and recursive packet deletion

Let A be any set of old plane points, disjoint from finite B on E, and
let c be a proper six-colouring of the actual unit graph on A. Put

    L(b)={0,...,5} minus {c(a): a in A, |a-b|=1}.

The unit graph on B has maximum degree at most four, no K4 (four distinct
plane points cannot be mutually equidistant), and common-neighbour
counts at most two. If at most two vertices of B see three distinct old
colours and every other vertex sees at most two, the lemma applies.
The old colouring extends. The stronger structural condition "at most
two old neighbours per vertex, with at most two vertices allowed three"
therefore preserves **every** old six-colouring simultaneously.

Consequently a single-ellipse exterior that deletes any old boundary
colouring must have either a vertex seeing at least four distinct old
colours, or at least three vertices each seeing at least three. These
are necessary conditions only. Neighbour counts cannot be substituted
for distinct-colour counts in a purported failure witness.

The rule composes. Delete finite disjoint packets B_k,...,B_1 in reverse
order, recording a carrier ellipse for each packet and its contacts to
the remaining vertices. If every packet satisfies the structural
two-contact/two-exception guard, any colouring of the remaining base
extends by adding B_1,...,B_k. The carrier ellipses may differ. A guard
depending on actual colours must instead be checked for every colouring
encountered; one successful colouring does not establish full support.
This is a block deletion certificate. No example here establishes that
it is strictly stronger than degree-five vertex deletion on ellipse
placements.

Incoming [circle-port palettes](../circle_ports_polygon_palettes_20261021/README.md)
gives a finer special case: circle unit graphs are paths and cycles;
at radius one only paths and even6-cycles occur. The resulting2-list
extension rule is stronger on that carrier. The present contribution is
the ellipse case with two allowed extra exclusions, using its degree4
and common-neighbour bounds. Both certificates retain the actual old
colouring and contacts rather than a palette cardinality alone.

## 5. What the Kepler transfer does and does not preserve

The [regular-polygon area-clock result](../runner_area_clock_20261020j/README.md)
preserves the common runner clock and, for a matching regular polygon,
the exact critical loneliness predicate. A Kepler orbit on an ellipse
instead has constant swept area about a focus, generally nonlinear
perimeter and angular speed, and no reason to preserve unit distances.
An affine map from a circle to an ellipse preserves incidence and area
ratios but sends the unit circle of differences to an ellipse. It is
not a unit-distance graph isomorphism unless its relevant directions
satisfy a separate length condition.

The useful transfer is therefore a constraint: an elegant single-conic
assembly with sparse old contacts is guaranteed harmless to six-colour
boundary support. A coupled construction must exceed the stated contact
gate, retain exact unit edges across its carriers, and then test the
full separating-colouring support. Changing the area clock cannot
create a missing geometric constraint.

## 6. Reproduction and finite scope

    python 05-knowledge/results/conic_boundary_support_20261020j/check.py
    python -O 05-knowledge/results/conic_boundary_support_20261020j/check.py --verify

The checker verifies the exact five-point geometry independently of the
quartic argument. It also enumerates the NetworkX graph atlas (every
unlabelled simple graph through seven vertices) and every choice of at
most two special vertices, testing the necessary degree and Gallai-tree
conditions for a minimal failed instance. No candidate survives all
three graph hypotheses. This finite audit is not the proof for arbitrary
graphs or arbitrary lists; that proof is Section 3. A seven-vertex
near-miss survives when the common-neighbour guard is removed, recording
why that coordinate matters to this proof. The near-miss is not claimed
to be non-list-colourable.
