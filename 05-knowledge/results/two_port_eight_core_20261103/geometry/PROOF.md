# Eight private points preserve the full Moser five-colour boundary

**Status: PROVED with the complete FINITE-EXACT eight-graph classification
in [cores/order8.json](../cores/order8.json).** Every proper five-colouring
of the fixed seven Moser-spindle points extends after adding any eight or
fewer distinct private planar points, retaining every unit-distance edge.
A blocking extension on this fixed carrier therefore needs at least nine
private points. This is a carrier-specific extension theorem, not a new
six-chromatic unit-distance graph or a bound for unrestricted constructions.

There is also an exact small-order failure boundary: **an injectively
realizable planar unit-distance graph on at most eight vertices is
three-choosable if and only if it contains no Moser-spindle subgraph**.
Every such graph containing a Moser spindle has choice number exactly four.
Both statements allow arbitrary lists from an arbitrary common palette in
the three-choosability step; there is no restriction to a sampled palette.

The [independent mathematical and executable review](../recovery/GEOMETRY_REVIEW.md)
accepts these claims within the stated finite and geometric scopes.

## Inheritance and scope

Use the seven old points B, in the order (A,B,C,D,B',C',D'):

    A=(0,0), B=(sqrt3,0), C=(sqrt3/2,1/2), D=(sqrt3/2,-1/2),
    B'=R(B), C'=R(C), D'=R(D),
    R(s,t)=((5s-sqrt11*t)/6,(sqrt11*s+5t)/6).                 (1)

The overloading of B for the boundary set and its second vertex is only
notational; indexed coordinates and all artifacts use the displayed order.
Let I be any set of at most eight distinct planar points disjoint from (1).
The fixed-boundary theorem concerns the full induced unit-distance graph
on B union I. Every old proper colouring is retained, including all 42
equality partitions that occur with five colours, not merely pairwise
boundary marginals or a selected assignment.

Inherited proved mechanisms, with theorem IDs and unambiguous sources:

- **H5**, [two_port_beyond_nine_20261103/geometry/PROOF.md](../../two_port_beyond_nine_20261103/geometry/PROOF.md):
  an exterior point has at most two unit contacts to (1).
- **H6**, the same file: every planar unit graph on at most six vertices is
  three-choosable; the sole possible six-vertex minimum-degree-three core
  is the triangular prism, handled by an elementary list argument.
- **S1–S2**, [two_port_seven_core_20261103/geometry/PROOF.md](../../two_port_seven_core_20261103/geometry/PROOF.md):
  the complete seven-vertex core classification and three-choosability of
  its three non-Moser geometrically surviving types.
- **S3–S5**, the same file: Moser rigidity, the uniform-three-list failure
  criterion, and the reciprocal-contact proof that every five-colouring
  of (1) extends across at most seven private points.

The closest mechanism is exterior contact capacity converted into private
list size. The hostile example is the Moser graph with the same three
colours at every vertex: abstract minimum list size three alone is not
sufficient. The corrected near miss is treating an abstract bad list or a
local metric filter as a realizable unit-distance obstruction. The useful
sidecars here are a graph-polynomial coefficient and a rigid subgraph's
exterior contact capacity. No historical priority claim is made.

Concept board: full boundary relations; low-degree removal; arbitrary
three-lists; exact local metric filters; rigid Moser subgraphs; circle fans.
The map from coordinates to an abstract core preserves its unit edges and
list constraints but loses global geometric compatibility. Sections E3–E4
restore exactly the information needed for the eight residual core types.

## E1. Complete finite core universe

A point outside (1) has at most two old neighbours, hence every proper
five-colouring of the boundary leaves each private point a list of at
least three colours. If an eight-vertex private graph has a vertex of
private degree at most two, remove it. S1–S5 extend the old colouring
over the remaining at most seven private points, and at least one colour
remains when the removed point is restored. Thus an eight-point blocking
extension would need minimum private degree at least three.

The producer [cores/census.cpp](../cores/census.cpp), through
[cores/run.py](../cores/run.py), examines **all 2^28 labelled simple graphs
on eight vertices**. It retains minimum degree at least three and imposes
only these necessary conditions for injective unit realizability:

1. Two distinct vertices have at most two common neighbours, since two
   distinct unit circles have at most two intersection points.
2. Each vertex's neighbourhood induces a union of paths and six-cycles.
   All those neighbours lie on one unit circle. A unit chord changes the
   angle by exactly plus or minus 60 degrees, so components lie within
   regular six-point rotation orbits.

These conditions are necessary, not a claimed realization criterion.
There are 53,470,067 labelled graphs of minimum degree at least three;
724,080 pass both metric conditions. Complete permutation-orbit reduction
gives **38 classes**. The producer checks disjoint coverage of the entire
survivor set by full 8! permutation orbits, with orbit-stabilizer equality.
The JSON records representatives, adjacency rows, degrees, orbit sizes,
automorphism counts, and a checksum of all labelled survivors. Pair bits
are (0,1),(0,2),...,(6,7), in lexicographic order, bit zero first.

Thirty classes have the polynomial certificate in E2. The eight remaining
masks are

    5675891, 5683827, 5691763, 5692274,
    6215909, 6216436, 14587370, 48675698.

E3–E4 show that none of these eight admits an injective planar unit
realization. Consequently every actual eight-vertex unit graph of minimum
degree at least three is three-choosable. It is unnecessary to decide
whether each of the other thirty abstract types has a geometric realization.

The local [verify.py](verify.py) checks every row and independently
recomputes the thirty coefficients by summing signed endpoint choices.
It consumes the pinned complete census; it does not silently substitute
a smaller graph family or claim to rerun all 2^28 graphs itself.

## E2. A nonzero bounded-exponent coefficient certifies arbitrary lists

For a graph G on vertices 0,...,n-1 define

    P_G(x) = product over edges i<j of (x_i-x_j).

Suppose a monomial x_0^d_0 ... x_(n-1)^d_(n-1) has nonzero coefficient,
with 0 <= d_i <= 2 and sum d_i = |E(G)|. Then G is three-choosable.

Proof. Given lists L_i with at least three colours, choose subsets A_i of
size d_i+1. Encode the finitely many colours involved by distinct real
numbers. The product of the univariate Lagrange leading-coefficient
functionals gives

    [x^d] P_G = sum over a in product A_i of
        P_G(a) / product_i product_(b in A_i, b != a_i) (a_i-b).    (2)

For completeness, the univariate functional on powers x^e is zero for
e<d_i and one for e=d_i, directly by interpolation of degree at most d_i.
For a monomial of total degree at most sum d_i, either its exponent vector
equals d or some exponent is smaller than d_i. The product functional
therefore extracts precisely the coefficient in (2), even if another
individual exponent exceeds its d_i. Every denominator is nonzero.
A nonzero coefficient makes at least one P_G(a) nonzero, and that choice
is a proper list colouring. This proves the claim without a palette bound.

The producer expands P_G modulo the ideals x_i^3=0. Since exponents only
increase during multiplication, this truncation does not alter any retained
coefficient. [polynomial_certificates.json](polynomial_certificates.json)
records a nonzero coefficient and its exponent vector for each of the
thirty classes. The independent local computation instead visits all
2^|E| choices of one endpoint of every edge and sums their signs when the
resulting exponents match the certificate. It agrees with all thirty
producer coefficients. A missing coefficient is not treated as a proof
of non-three-choosability; the eight zero-support cases require geometry.

## E3. A Moser subgraph excludes a minimum-degree-three eighth point

Every injective unit realization of the Moser graph is congruent, as a
point set, to (1). Here is the inherited S3 mechanism. In a diamond K4
minus one edge, its two nonadjacent tips are the two distinct common unit
neighbours of its adjacent wings. Thus their separation is sqrt3, and
the diamond is rigid. Two such diamonds sharing tip A have their other
tips B,B' at distance sqrt3 from A. The closing edge BB' has length one,
so the angle between the two long vectors has cosine 5/6. Its two signs
differ by reflection, and each diamond's wings are already determined.

H5 then applies in any isometric copy: an exterior point has at most two
unit contacts to the seven Moser points. The exact H5 computation is
repeated here over Q(sqrt3,sqrt11). It tests all 35 triples of the old
points using both exact circumcentre and squared-circumradius formulas.
Exactly ten triples have unit circumradius, and each centre is already
one of the seven old points. No exterior point has three contacts.
Three distinct collinear points cannot lie on a unit circle, so there is
no omitted collinear triple case.

Six of the eight residual graphs contain a Moser subgraph on seven of
their vertices. The following embeddings map the standard coordinate
order (A,B,C,D,B',C',D') into each residual graph:

| Mask | Moser embedding | Exterior vertex | Its neighbours |
|---:|---|---:|---|
|5691763|0,3,1,6,4,2,5|7|0,1,2|
|5692274|4,0,2,5,6,1,3|7|0,1,2|
|6215909|2,0,1,7,3,4,5|6|0,1,3|
|6216436|4,0,3,5,7,1,2|6|0,1,3|
|14587370|0,1,2,7,3,4,6|5|1,2,3,4|
|48675698|1,2,3,7,5,4,6|0|2,5,6,7|

All eleven Moser edges and the displayed exterior contacts are checked
in [metric_exclusions.json](metric_exclusions.json). The subgraph need
not be induced: its specified unit edges alone force the rigid geometry.
The exterior point has at least three prescribed contacts, contradicting
the bound of two. Therefore these six types are not unit realizable.

## E4. A six-vertex circle fan forces a collision

**General metric obstruction.** Let h,a,b,c,d,z be six distinct vertices.
Require the nine edges

    ha,hb,hc,hd;  ab,bc,cd;  za,zd.                         (3)

There is no injective planar unit-distance realization of (3).

Proof. Normalize h to the origin. The four neighbours a,b,c,d lie on its
unit circle. Each of the three path steps ab,bc,cd changes the angle by
plus or minus 60 degrees. Changing sign between successive steps would
repeat the vertex from two steps earlier, which injectivity forbids.
The three signs are therefore the same. The four neighbours are four
consecutive regular-hexagon vertices, and a,d are antipodal. Their distance
is two. The two unit circles centred at a,d are tangent at their midpoint
h, so a point z at unit distance from both must equal h. This contradicts
the stipulated distinctness. Additional edges cannot repair the collision.

Equivalently, one can normalize the four neighbours exactly as

    (1,0), (1/2,sqrt3/2), (-1/2,sqrt3/2), (-1,0).

The checker verifies these distances in the exact field. This is a
six-vertex, nine-edge forbidden subgraph; no claim of edge-minimality is
needed or made. If injectivity is dropped, z=h gives a unit drawing,
showing precisely which hypothesis supplies the obstruction.

Both residual masks 5675891 and 5683827 contain (3) with

    (h,a,b,c,d,z) = (0,1,7,2,5,4).

Vertex 6 and the remaining eighth vertex play no role. All nine required
edges are verified directly. This excludes the last two residual types.

## E5. The extension theorem and the exact small-order boundary

If the private graph has at most seven vertices, inherited S1–S5 apply.
If it has eight vertices and a vertex of degree at most two, remove and
restore that vertex as in E1. Otherwise it has minimum degree at least
three. The complete E1 census, E2 certificates, and E3–E4 exclusions show
that it is three-choosable. Its lists from any fixed old five-colouring
have size at least three, so that colouring extends. This proves the
announced theorem for **every** actual placement and **every** old proper
five-colouring. No boundary state is discarded by a quotient or selected
to make the proof work.

There is a stronger structural corollary:

> An injectively realizable planar unit-distance graph G on at most eight
> vertices is three-choosable if and only if G contains no Moser subgraph.
> If it contains a Moser subgraph, its choice number is exactly four.

For the forward implication, the Moser graph has chromatic number four,
so a graph containing it cannot be three-choosable. Conversely, take a
graph with no such subgraph and induct on its order. A vertex of degree
at most two can be deleted and greedily restored for arbitrary lists of
size three. If minimum degree is at least three, H6 treats orders at most
six, S1–S2 treat order seven (the sole non-three-choosable realizable core
there is Moser), and E1–E4 treat order eight. The induction is valid because
deleting a vertex neither creates a Moser subgraph nor destroys unit
realizability. This proves the equivalence.

For the last assertion, the Moser graph is three-degenerate and hence
four-choosable by greedy deletion and restoration. Any eighth point has
at most two contacts to its rigid seven-point Moser subgraph, by E3.
There are no extra unit edges among those seven points. Removing the
eighth point thus leaves a three-degenerate graph, so the whole graph is
four-choosable. Its Moser subgraph gives the matching lower bound four.
This corollary also allows unit-edge subgraphs with some actual unit
contacts omitted; all its metric tests used necessary edge requirements.
The fixed-boundary theorem above, and the actual control below, retain
all contacts.

## E6. An actual connected private three-core, with the full boundary table

The positive control uses the old points (1) and the private points

    p_m = (0,1) + sum_(i=0,1,2 with bit_i(m)=1) v_i,
    v_0=(1,0), v_1=(0,1), v_2=(3/5,4/5),  0 <= m < 8.       (4)

All fifteen points are distinct. Exact testing of all 105 unordered pairs
gives 25 unit pairs and 80 nonunit pairs. Among the private points, an edge
occurs exactly when the two bit masks differ in one bit, so the private
graph is the cubic eight-vertex cube, a genuine minimum-degree-three core.
The **only** cross-edges are A-p_0 and C-p_0. Since AC is an old edge,
every proper boundary assignment gives p_0 a list of exactly three colours
at q=5. This is an interacting, connected unit realization, not a remote
disconnected abstract list example.

For each q=4,5,6, the checker enumerates all q^8 private assignments,
counts every proper cube colouring by the colour at p_0, and combines
those exact counts with every proper old boundary equality partition.
[boundary_tables.json](boundary_tables.json) retains each partition's
weight, named-colour multiplicity, and explicit private colouring witness.
All cross-edges and all private edges are used when checking each witness.

| Palette | Old boundary partitions | Cube colourings | Extensions per old assignment | Full graph colourings |
|---:|---:|---:|---:|---:|
|4|16|2,652|1,326|509,184|
|5|42|29,660|17,796|89,691,840|
|6|52|198,030|132,020|4,182,393,600|

The constant extension weight within each row follows here from colour
symmetry at p_0 and its two differently coloured old neighbours. It is a
property of control (4), not an assertion that arbitrary attachments have
constant boundary weights. The full old Moser graph has zero three-colour
assignments, 384 four-colour assignments and 5,040 five-colour assignments,
checked literally. The fifteen-point control is therefore exactly
four-chromatic. Its positive five-colour relation is a proof control, not
a candidate six-chromatic obstruction.

## Replay, status, and next obligation

From the repository root:

```text
python 05-knowledge/results/two_port_eight_core_20261103/cores/run.py --order 8
python 05-knowledge/results/two_port_eight_core_20261103/geometry/verify.py
python -O 05-knowledge/results/two_port_eight_core_20261103/geometry/verify.py
```

The first command reruns the complete labelled universe; the local checker
checks its pinned result and the certificates/geometry described above.
The exact producer command syntax is documented by `cores/run.py`.
Ordinary local runs compare frozen artifacts and receipts; only an explicit
`--write` updates them. All required checks use exceptions, not removable
Python assertions. [SOURCE_PINS.json](SOURCE_PINS.json) pins inherited
proofs, the field implementation, the old coordinates, and the current
complete census and its producer sources. Sources are checked both before
and after computation. Normal and optimized receipts agree.

**PROVED:** the fan obstruction, the interpolation implication, the rigid
Moser exterior exclusion, the small-order three-choosability equivalence,
and the fixed-boundary eight-private-point extension theorem.
**FINITE-EXACT:** the complete 2^28 graph census and its 38 orbit classes;
the thirty coefficient certificates; all 105 distances and full boundary
tables for (4). **OPEN here:** arbitrary nine-or-more-private-point
attachments, and any genuine six-chromatic planar unit-distance forcing.
Nine private points is a necessary next size for blocking this carrier,
not a sufficient construction or a claim that such a blocker exists.
The [nine-point frontier reduction](FRONTIER_NINE.md) gives the exact
two-hub boundary condition when such a private core contains a Moser
subgraph; its realization and blocking outcome remain open.

The failed implication was “abstract core plus locally permissible metric
data implies realizable obstruction.” The repaired mechanism keeps both
the list certificate and the exact geometric sidecar. For an actual
six-colour obstruction, one must still give distinct planar coordinates,
all unit contacts, and a proof that **every** five-colour assignment fails.
