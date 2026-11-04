# Seven private vertices still cannot block a fixed Moser boundary

**Result: PROVED with the complete FINITE-EXACT seven-graph classification
below.** Every proper five-colouring of the fixed seven Moser-spindle points
extends after adding **any seven or fewer distinct private planar points**,
retaining every unit-distance edge. The new case is seven private vertices.
A blocking extension on this particular carrier now needs at least eight.
This is not a universal order bound for six-chromatic unit graphs, and no
new six-chromatic graph or global plane-colouring claim is made.

The [independent mathematical and executable review](../recovery/GEOMETRY_REVIEW.md)
accepts the complete proof chain and the stated finite scopes without blockers.

## Inheritance, quantified scope, and missing coordinate

The [previous packet, H5–H6](../../two_port_beyond_nine_20261103/geometry/PROOF.md)
proved that every point outside the fixed spindle has at most two unit
contacts to it, and that six or fewer private vertices preserve every
five-colour boundary state. The earlier positive-ear theorem and two-hub
conditional zeros are not being promoted to an unrestricted forcing result.

Use the same seven points, in the order(A,B,C,D,B',C',D'):

    A=(0,0), B=(sqrt3,0), C=(sqrt3/2,1/2), D=(sqrt3/2,-1/2),
    B'=R(B), C'=R(C), D'=R(D),
    R(s,t)=((5s-sqrt11*t)/6,(sqrt11*s+5t)/6).                 (1)

Let I be any set of at most seven additional points, pairwise distinct and
disjoint from(1). The graph throughout is the full induced unit-distance
graph on the union. No old colouring is selected for convenience: the
claim quantifies over **every** proper assignment of five named colours
to(1), including every realizable equality partition.

Closest proved mechanism: exterior-contact bounds become lower bounds on
private colour lists. Hostile: a private copy of the Moser graph itself is
not three-colourable and can defeat arbitrary abstract three-lists. Corrected
near miss: a permissible abstract bad list need not be realizable by contacts
to the old boundary. Least-used sidecar: **reciprocal** contact capacity after
the private core is proved congruent to the old carrier.
Board: private3-core; full boundary partitions; nonuniform lists; rigid
realization; reciprocal incidences; conditional versus global failure.

## S1. The finite core classification is complete

Each private vertex has at most two old neighbours by the inherited contact
lemma, so a fixed proper boundary colouring leaves it a list of at least
three colours. If the private graph has a vertex of degree at most2, delete
it; the remaining at most six-vertex unit graph is three-choosable by the
previous packet's elementary proof. Restore the deleted vertex, which sees
at most two privately used colours. Thus a seven-vertex counterexample
must have minimum private degree at least3.

The complete census `census.py` visits all2^21 labelled simple graphs on
vertices0,...,6. It keeps minimum degree at least3 and imposes only these
necessary geometric conditions:

1. Every pair has at most two common neighbours, by intersection of unit
   circles with distinct centres.
2. Each vertex's induced neighbourhood is a union of paths and six-cycles.
   Indeed its neighbours lie on its unit circle. Adjacency on that circle
   changes the angle by exactly60 degrees, so components lie in regular
   six-point rotation orbits. The pair condition already bounds internal
   degrees by2; cycles of lengths3,4,5 must still be excluded.

These are necessary filters, not an asserted geometric characterization.
There are236926 minimum-degree3 labelled graphs;8190 pass both filters.
All8190 are classified, with full labelled-orbit coverage, into six types:

| Mask | Labelled orbit size | Degree multiset | Description used below |
|---:|---:|---|---|
|120383|420|3^6,6|wheel with six-cycle rim|
|120444|2520|3^6,4|deleting0 leaves a theta|
|120568|1260|3^6,4|deleting1 leaves a theta|
|122616|840|3^4,4^3|triangle with three reflected caps and their common hub|
|378106|630|3^6,4|Moser spindle|
|380278|2520|3^4,4^3|spanning Moser plus one additional edge|

The mask uses lexicographic unordered pairs(0,1),(0,2),...,(5,6), with
bit0 first. `census_receipt.json` supplies every edge set, labelled-orbit
cardinality, and a hash of all8190 surviving labelled masks. The code checks
that the union of those complete permutation orbits equals the survivor
set exactly. No vertex-transitivity assumption, representative heuristic,
unproved edge bound, or truncated graph family enters the census.

### The two remaining false geometric types

In122616, vertices1,2,3 form a unit equilateral triangle. Vertices6,5,4 are
the second common unit neighbours of its sides12,13,23, respectively;
they are distinct from the old triangle vertices. Circle intersections
force those three reflected cap positions. They form an equilateral
triangle of side2. Vertex0 is required to be at unit distance from all
three caps, impossible: two radius-one circles with centres distance2
are tangent, and their three pairwise midpoints differ. The checker verifies
the forced squared distances4 exactly.

Deleting edge(2,3) from380278 gives a spanning copy of378106. The rigidity
lemma S3 below shows every injective unit realization of the Moser graph
is congruent to(1), whose complete induced unit graph has exactly eleven
edges. Hence the additional edge cannot also have unit length. This rules
out380278 without overlooking another realization of its spanning graph.

Consequently every actual seven-vertex private3-core is either one of the
three types treated in S2 or a Moser spindle. It is unnecessary to assert
that the two abstract theta types have geometric realizations: proving them
three-choosable safely handles that possible overinclusion.

## S2. The three non-Moser survivors are three-choosable

Lists may belong to any common palette and may be larger than the minimum
sizes; restrict them to the sizes used below.

**Wheel120383.** Colour its hub arbitrarily. Each rim vertex retains a list
of at least two colours. An even cycle is two-choosable: if all two-lists
are equal, alternate their two colours; otherwise choose an adjacent pair
with different lists, colour the first from its list but outside the second's,
and greedily follow the long path around to that second vertex. Its first
neighbour's colour is absent from its list, so the final step also succeeds.

**Theta120444.** Colour its unique degree4 vertex0 arbitrarily. The remaining
graph has degree3 ends1,2 joined by three paths

    1-2,    1-6-2,    1-5-3-4-2.

The ends retain three colours, because they were not adjacent to0. Every
internal vertex retains at least two. Restrict those internal lists to two.
Choose a colour for end1 outside the list at6. Removing end1 leaves the path
6-2-4-3-5. Vertex6 retains two colours, the opposite endpoint5 has at least
one, and every internal vertex has at least two. Greedily colour this path
from5 towards6. All edges to the removed, coloured vertices have already
been accounted for in the lists.

**Theta120568.** Colour its degree4 vertex1 arbitrarily. The remaining
degree3 ends0,4 are joined by paths

    0-4,    0-5-3-4,    0-6-2-4.

Again the ends retain three colours and the internal vertices at least two.
Choose end0's colour outside the two-list at5, then colour the remaining
path6-2-4-3-5 towards5. This is the same direct argument, with every vertex
and list-size inequality explicit. Thus none of these three abstract cores
can be a counterexample.

These are unbounded list-colouring arguments. The independent common-five-
palette list census in the parent packet is a corroboration, not the reason
they hold for arbitrary lists.

## S3. Every unit realization of the private Moser core is rigid

In a diamond K4 minus an edge, let C,D be the adjacent degree3 vertices.
The other two vertices A,B are distinct common unit neighbours of C,D.
The two intersecting unit circles supply exactly those two points, so A,B
are the opposite equilateral triangle tips and |A-B|=sqrt3. The diamond is
therefore rigid up to isometry and exchanging its two wings.

A Moser graph consists of two such diamonds with a common tip A and an
edge joining the other tips B,B'. Both vectors B-A and B'-A have length
sqrt3, and |B-B'|=1. Their inner product is5/2 and their angle has cosine5/6.
Hence their two possible relative angles differ by reflection. Each diamond's
other two vertices are already determined by its tips. It follows that
every injective unit realization of the seven-vertex Moser graph has the
same point set as(1), up to a planar isometry and harmless wing relabellings.
The exact full-pair test of(1) has no additional unit edges.

In particular the inherited contact lemma applies **reciprocally**: if seven
private points form this core, every old boundary point, being distinct
from all seven private points, has at most two unit neighbours in that
private core. This conclusion would not follow from an abstract graph
isomorphism without the rigidity proof.

## S4. The only bad Moser lists are identical three-colour lists

**Lemma.** A Moser graph with lists of size at least3 is uncolourable from
those lists if and only if all seven lists are exactly the same three-set S.

First suppose all lists have size3, using labels from(1). The graph obtained
by deleting its degree4 vertex A is two triangles BCD and B'C'D' joined by
the edge BB'. Its tips B,B' have degree3 and its four wings have degree2.

We use an elementary degree-list fact: in a connected graph, lists of
size at least each vertex's degree are colourable if one vertex has strictly
more colours than its degree. Root a spanning tree at that vertex and colour
in reverse tree order. Every nonroot has an uncoloured parent and therefore
at most degree-minus-one coloured neighbours; the root has its surplus.

If some c in L_A is absent from any wing list, colour A with c. In the
remaining connected graph the wing lists have size at least2, the tip lists
size3, and that particular wing still has size3. The degree-list fact gives
an extension. Therefore a bad list assignment must satisfy

    L_C=L_D=L_C'=L_D'=L_A=S.                                 (2)

For a choice c in S at A, each pair of adjacent wings uses S\{c}. The
remaining possible colours for B and B' are respectively

    L_B\(S\{c}),       L_B'\(S\{c}).                        (3)

Both sets are nonempty. If L_B differs from S, choose c so that the first
set has at least two elements: two outside-S colours suffice for any c;
with exactly one outside-S colour, choose c from L_B intersect S. Then the
adjacent tips can take different colours, regardless of the second nonempty
set. The same argument holds with B and B' exchanged. Thus failure forces
both tip lists also equal S. Conversely identical S lists fail because the
Moser graph is not three-colourable.

For lists larger than3, choose any three-subsets. If the original assignment
were bad, the chosen three-subsets would have to be uniform by the preceding
argument. A list larger than3 permits changing just its three-subset while
holding the others fixed, making them nonuniform and therefore colourable.
This contradiction proves the lemma in the full stated scope.

## S5. Uniform bad lists cannot come from the old spindle boundary

Fix any proper five-colouring of the old spindle, and suppose a private
Moser core fails to extend it. By S4 every private list is exactly the same
three-set S. Each private vertex has at most two old contacts by the old
carrier's contact lemma. Therefore it must have exactly two such contacts,
with their colours precisely the two excluded colours in[5]\S. Across all
seven private vertices this requires **fourteen** old-private unit edges.

The old spindle has independence number2. Indeed its six noncentral vertices
split into the triangles BCD and B'C'D'; an independent triple would therefore
have to include A. Its only two possible companions are B and B', which are
adjacent. Consequently each boundary colour class has at most two vertices,
and the two excluded colours occupy at most four old points.

By S3's reciprocal contact bound, those at most four old points can supply
at most **eight** unit edges to the private Moser core. This contradicts
the fourteen required edges. Unused excluded colours only reduce the supply.
Thus the residual abstract bad-list system is not geometrically realizable.

Combining S1–S5 proves the announced theorem. Every old proper five-colour
boundary assignment extends for every injective placement of at most seven
private points, with all old-private and private-private unit edges included.
The result strengthens the earlier six-private bound; it does not merely
report another successful finite placement.

## S6. Full boundary spectrum and a connected exact14-point control

The old spindle's complete proper five-colour boundary spectrum consists of
42 equality partitions:16 with four blocks and26 with five. Each represents
120 labelled assignments, giving5040 colourings. All42 are retained. For
each partition and each of the ten excluded colour pairs, `capacity_table.json`
records the actual old vertices available, their reciprocal capacity at most8,
and the required demand14. These420 cases corroborate S5's uniform proof;
no particularly favourable boundary assignment is substituted for the others.

An actual private Moser3-core is also tested, rather than only discussed
abstractly. Put e=B-B', let J(s,t)=(-t,s), and translate a second copy of(1) by

    t=(e - sqrt3*J(e))/2.                                   (4)

Since |e|=1, this is a unit translation. The14 points are all distinct.
Their complete induced unit graph has30 edges:eleven within each spindle,
the seven corresponding cross-edges, and the extra cross-edge from old B
to private B'. Every one of the91 point pairs is checked exactly in
Q(sqrt3,sqrt11). The private B' sees two differently coloured old neighbours
in every proper boundary colouring, so this control actually reaches a
three-element private list; the other private lists have size4.

Literal private colouring enumeration gives a positive extension count for
every one of the42 proper five-colour boundary partitions. Exact weights
and witness colourings are saved in `boundary_tables.json`. The whole graph
has4360320 five-colourings. It also retains all16 four-colour boundary
partitions, has16608 four-colourings, and is exactly four-chromatic because
it contains an old Moser spindle. This is a connected, interacting positive
control, not a six-chromatic witness.

## Reproduction, status and next boundary

Run `python3 census.py` and `python3 verify.py` here, then repeat with `-O`.
The latter freezes the expected census before rerunning it, checks inherited
source hashes, replays both circumcentre tests on all35 old triples, verifies
the canonical Moser mapping and two metric exclusions, and checks all full
boundary and coordinate artifacts without overwriting them. `--write` is
explicit artifact-production mode only. All validations remain active under-O.

The all-seven extension statement uses the complete finite graph census as
one load-bearing input. The list-colouring arguments, rigidity, and reciprocal
contradiction are proofs with no palette cutoff. No list sample is used to
infer universal choosability, and no local graph filter is treated as a
sufficient realization certificate. The parent packet's independent list
search supplies additional controls and proof review separately.

Map/loss ledger: actual geometry -> private3-core -> complete graph census
-> list-obstruction type -> forbidden old colour classes -> reciprocal unit
contact capacity. The graph quotient loses rigidity and old-private incidence;
S3 and S5 restore exactly those coordinates. The Moser uniform list is a
positive abstract hostile and a negative geometric one. Treating an arbitrary
three-list as physically realizable would fail precisely at the14-versus8
incidence inequality.

The next possible obstruction on this carrier needs at least eight private
vertices. Its private3-core, full five-colour boundary relation, and actual
cross-contacts must all be retained. No claim is made that the eight-vertex
case is sufficient, that these lists describe every larger gadget, or that
this carrier must appear in a global six-chromatic construction.
