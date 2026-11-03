# Two-anchor reflection phases and an infinite unit-distance recurrence

**PROVED elementary family classification; FINITE-EXACT independent controls.**
Local labels GP1--3. The objects are supplied metric configurations, not
abstract unit-edge graphs assumed to be rigid. No improved bound for the
chromatic number of the plane is asserted.

## GP1. The exact information lost by independent reflection quotients

Fix named anchors A=(0,0),B=(a,0), a>0. For each i=1,...,k a supplied
three-point metric piece has its off-axis point at

`P_i=(x_i,epsilon_i y_i), y_i>0, epsilon_i in {+1,-1}`.

Assume the x_i are pairwise distinct, so every choice gives distinct points.
Each piece separately is unchanged by its anchor-fixing reflection. In a
common plane, however,

`|P_i-P_j|^2=(x_i-x_j)^2+y_i^2+y_j^2-2 epsilon_i epsilon_j y_i y_j`. (G1)

Because y_i y_j is nonzero, a cross-distance determines the relative sign.
A spanning tree of k-1 such probes determines all signs up to simultaneous
reflection, and k-1 binary probes are necessary to distinguish2^(k-1)
possibilities. This is exact metric gluing, not a numerical analogy.

For unit adjacency alone, form the **eligibility graph** E: ij is an edge
when exactly one of the two distances in (G1) is1. Both cannot be1. Mark
that edge by the required product for unit distance. Edges outside E never
give unit distance and carry no sign information. On each component of E,
the complete observed unit graph determines every relative sign. Conversely,
an arbitrary flip of an entire component changes none of those observations.
If E has c components, every attainable unit graph has exactly

`2^(c-1)` metric phase classes modulo simultaneous reflection. (G2)

The image is characterized by the usual product-one condition on every
cycle of relative signs. Equivalently, observed unit/nonunit bits are an
affine cut-space coset over F2. The kernel is the space of componentwise
constant signs. Adding a new piece can connect two components and detect
an ambiguity that was formerly invisible. The retained sidecar is relative
reflection phase; individual anchored isometry types do not supply it.

For larger multi-vertex pieces one must also check coincidence and all
cross-distances. The unrestricted fibre formula is not asserted without
those guards.

## GP2. A four-point chromatic hostile

Take A=(0,0),B=(1,0),P=(0,1),Q=(sqrt(3)/2,+/-1/2). The separate anchored
pieces are the same for both signs. Their full induced unit graphs differ:

| Relative sign | Squared PQ distance | Full unit graph | Chromatic polynomial |
|---|---:|---|---|
| + |1|triangle APQ with leaf B at A|q(q-1)^2(q-2)|
| - |3|star K1,3|q(q-1)^3|

All other squared distances are checked exactly. Thus forgetting the relative
reflection can even change the required number of colors from3 to2. This
does not realize a high-chromatic Hadwiger--Nelson witness; it identifies
the missing geometric coordinate before a coloring argument is attempted.

## GP3. Complete regular-circle family and generating function

For every integer d>=1 put A=(0,0),B=(1,0) and

`P_j=(cos(j*pi/(3d)), epsilon_j sin(j*pi/(3d))), 1<=j<=3d-1`.

The x_j are distinct, so GP1 applies. There are3d+1 points including the
anchors. Two points on the unit circle are a unit distance apart exactly
when their angular difference is +/-pi/3. Consequently eligibility edges
are precisely

* j-i=d, with required relative sign +;
* i+j=d or5d, with required relative sign -.

The eligibility graph has the complete component classification

`K2 disjoint-union [P3 if d is even] disjoint-union floor((d-1)/2) C6`. (G3)

**Proof.** Partition indices by their residues modulo d. The zero residue
has only d,2d and gives K2. If d is even, residue d/2 has three indices,
linked in a path. Every other pair r,d-r has two three-vertex paths; the
two sum conditions link their ends to form C6. Distinct such residue pairs
do not meet. Each C6 has exactly two negatively marked eligibility edges.
This accounts for every vertex and every possible unit chord.

There are c=floor(d/2)+1 components. Every distinct marked unit graph has
exactly2^floor(d/2) phase orbits, and there are2^(3d-1-c) such graphs.
Every point is a unit distance from A, and BP_d is always a unit edge:
the number of fixed edges is3d+1. The K2 permits either eligibility edge
status; P3 permits four patterns independently. On a C6, the attainable
sets of present eligibility edges are exactly the even subsets: the six
relative products multiply to1 and the required-sign product is1. Thus

`E6(z)=1+15z^2+15z^4+z^6`

counts its distinct edge patterns. If U_d(z) counts the full marked unit
graphs by their total edge count, then

`U_d(z)=z^(3d+1)(1+z)^(1+2[d even]) E6(z)^floor((d-1)/2)`,   (G4)

and in particular

`U_(d+2)(z)=z^6 E6(z) U_d(z)`.                              (G5)

The initial values are U_1=z^4(1+z), U_2=z^7(1+z)^3. Counting metric phase
orbits instead multiplies U_d by2^floor(d/2), and introduces a factor2 in
(G5). These are marked graph counts, not unlabelled isomorphism counts.

All graphs in this family have chromatic number exactly3. Their outer-circle
graphs are induced subgraphs of disjoint six-cycles, hence bipartite; use a
third color for A. The fixed triangle A,B,P_d gives the matching lower bound.
In particular, growing geometric phase ambiguity alone does not force an
increasing chromatic number.

[geometric_phases.py](geometric_phases.py) independently enumerates every
phase orbit through d6 (65536 at d6), verifies uniform fibres and every
coefficient, and checks unit chords by exact cyclotomic arithmetic through
d12. The all-d proof is the residue decomposition above, not this finite
test. The coefficients15 come from binomial(6,2)=binomial(6,4), rather than
from identifying these graphs with the fifteen symplectic Lagrangians.

## Recovered geometric mechanisms

The [older Hadwiger synthesis](../hadwiger_six_seven_eight_20261025/SYNTHESIS.md)
and [metric spindle audit](../hadwiger_six_seven_eight_20261025/metric_spindle/RESULT.md)
already separate full metric conditions from combinatorial shadows.
The [eight-frame packet](../hadwiger_six_seven_eight_20261025/fano_eight/RESULT.md)
retains affine transport and pair-permutation monodromy discarded by scalar
parity. GP1 is the concrete planar C2 instance, with actual cross-distance
probes and GP3 supplying an infinite classified family. It does not
transport the old palette proper-coloring predicate without its guards.
