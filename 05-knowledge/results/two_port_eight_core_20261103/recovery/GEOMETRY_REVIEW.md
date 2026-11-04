# Independent audit of the eight-private-point theorem

**ACCEPTED within the stated scope.** Reviewed E1–E6 in
[../geometry/PROOF.md](../geometry/PROOF.md), the producer and deletion-parent
core enumerators, the geometry checker, all eight exclusion witnesses, and
the actual fifteen-point control. No proof gap was found. This review uses
the previously audited H5–H6 and S1–S5 dependencies named in that proof;
their present source hashes are checked by the producer's source pins.

## Complete universe and arbitrary lists

The producer's Gray traversal visits each labelled eight-vertex graph once.
Its maintained adjacency and degree arrays toggle exactly the changed edge.
Every candidate is tested for minimum degree three, pair codegree at most
two, and neighborhood components that are paths or six-cycles. These are
necessary conditions, including for graphs with some unit contacts omitted.
Every unit neighbor lies on its center's unit circle and has at most two
unit neighbors there, at angular offsets plus and minus 60 degrees.

The independent [../cores/independent.cpp](../cores/independent.cpp) uses a
different complete traversal. Deleting the designated vertex 7 from any
eight-vertex minimum-degree-three graph leaves a seven-vertex graph of
minimum degree at least two. Every degree-two parent vertex must belong
to the new attachment cut; the new vertex needs at least three neighbors.
Those are exactly the degree requirements used in the child enumeration.
The parent filters are hereditary, so applying them early loses nothing.
Every child has a unique deletion parent and cut.

The independent filter uses codegree at most two and bipartite neighborhoods.
These are equivalent to the producer's local filter at orders seven and
eight: codegree bounds the degree inside a neighborhood by two; bipartiteness
leaves even cycles; a four-cycle would give its opposite vertices their two
cycle neighbors plus the center as three common neighbors; there are too
few vertices for an eight-cycle. Thus only six-cycles remain. The independent
traversal exhausts the entire 724,080-element union of the 38 disjoint full
permutation orbits, not a sample. Its frozen normal and optimized receipts
agree and are source-pinned in this review's receipt.

[LIST_CERTIFICATE.md](LIST_CERTIFICATE.md) gives a separate proof of E2 by
finite interpolation on arbitrary finite selected lists. The independently
implemented [check.py](check.py) expands all 548,864 complete endpoint choices
over the 38 representatives and checks the complete capped support, not only
one selected coefficient. It confirms all thirty signed certificates and
exactly eight zero-support residuals. A zero coefficient bank is not treated
as a list obstruction or a geometric realization.

## The eight metric exclusions

The six Moser-containing witnesses specify all eleven Moser edges, seven
distinct embedded labels, and a distinct eighth vertex with at least three
neighbors in those seven. Mandatory Moser edges already force the rigid
two-diamond realization; no induced-subgraph hypothesis is required. The
rigidity argument fixes the long-tip distance at square root of three and
the angle between the two long vectors by cosine 5/6. The two signs differ
by reflection. H5 applies after the resulting isometry.

The independent [geometry_review.py](geometry_review.py) checks all 35 old
triples using the squared-circumradius determinant identity in an independently
implemented rational four-dimensional field. Exactly ten triples lie on a
unit circle, each with its unique center among the seven old points. The
collinear case cannot hide a circle through three distinct points. Therefore
no exterior point has three contacts, excluding all six witnesses.

For the other two masks, the nine required edges are the four hub spokes,
a path through the four neighbors, and a separate cap joined to its endpoints.
Three successive steps on the hub's unit circle have angle plus or minus
60 degrees. A sign change revisits the preceding vertex, so injectivity
forces all three signs equal. The endpoints are antipodal. Their two unit
circles meet only at the hub, forcing the distinct cap to coincide with it.
All nine edges and all six distinct labels are checked directly. This is a
mandatory-edge obstruction; additional edges cannot evade it. Setting the
cap equal to the hub gives the explicit hostile when injectivity is dropped.

## Quantifiers in E5

For the fixed boundary, every private list has at least three allowed colors
by H5. If an eighth private vertex has private degree at most two, first
apply the inherited seven-private-point theorem and then restore that vertex.
This restoration forbids at most two colors from its private neighbors.
Otherwise the complete minimum-degree-three classification applies. The
thirty list certificates and eight metric exclusions exhaust it. Therefore
every proper old five-coloring extends for every actual placement of at most
eight private points. The lower bound of nine is a necessary size for this
fixed carrier, not a sufficient construction or an unrestricted chromatic
bound.

The stronger small-order equivalence is valid. On at most eight vertices,
Moser-free unit graphs are three-choosable by vertex deletion and the audited
core classifications at orders at most six, seven, and eight. Deletion
preserves both unit realizability and Moser-freeness. Conversely a Moser
subgraph has chromatic number four, so its containing graph cannot be
three-choosable. With at most one additional vertex, the rigid Moser subgraph
has no added internal unit edges and the outside vertex has degree at most
two. Deleting it leaves the three-degenerate Moser graph. Hence the containing
graph is four-choosable and has choice number exactly four. The argument
also applies when some unit contacts are omitted.

## Actual full-boundary control and independent checks

The independent field implementation checks all 105 unordered pairs of the
fifteen listed points, their distinctness, all 25 unit edges, and all 80
nonunit pairs. The private graph is exactly the cube; the only cross contacts
are the same private vertex to adjacent old vertices A and C. Hence every
proper old coloring forbids exactly two different colors at that private
vertex.

The producer counts all private assignments. The independent route uses
edge-subset inclusion–exclusion over all 4,096 subsets of the cube's twelve
edges: the chromatic polynomial at q is the sum of (-1)^|F| q^components(F).
Color permutation symmetry then gives exactly (q-2)/q of these colorings
after excluding the two boundary colors. This independently reproduces
2,652, 29,660, 198,030 cube colorings and 1,326, 17,796, 132,020 extensions
per old assignment at q=4,5,6. Every stored boundary partition and explicit
private witness is checked against all edges, and named multiplicities are
recomputed. In particular all 42 five-color equality states remain present.

Reproduction commands:

```powershell
python 05-knowledge/results/two_port_eight_core_20261103/recovery/check.py
python 05-knowledge/results/two_port_eight_core_20261103/geometry/verify.py
python 05-knowledge/results/two_port_eight_core_20261103/recovery/geometry_review.py
python -O 05-knowledge/results/two_port_eight_core_20261103/recovery/geometry_review.py
```

The review checker compares the separately replayed producer receipt with
its frozen receipt, checks all dependency/artifact pins, and emits its own
receipt. Normal and optimized independent receipts agree. It does not claim
to independently rerun the producer's full 2^28 traversal; the distinct
deletion-parent traversal and its pinned receipts establish that cross-check.
