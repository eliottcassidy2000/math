# Even counts do not force a symmetry of the actual line diagram

**FINITE-EXACT complete classification at orders3,4,5; proof from the exact
class-and-line census.** The conclusions concern the literal merged
complement-line multigraph, retaining parallel lines and doubled loop entries.
They do not refer to the auxiliary arc-flip graph or to an unweighted support.

Let X_n be the ordinary tournament classes, R their converse involution,
and G_n the merged multigraph on X_n/R. A pairing J is **compatible with the
diagram** when it commutes with R and its induced permutation of X_n/R
preserves the full dart adjacency matrix A. A colored-compatible pairing
also preserves both blue and black matrices. This is a precise requirement;
even cardinality of X_n by itself supplies only an abstract matching.

The previous [pairing theorem](../paired_recursions_20261103/PAIRINGS.md)
allows an H-preserving J commuting with R, with minimum fixed-point count
sum_h(c(n,h) mod2). At orders3,4,5 those minima are2,2,4. Adding actual
diagram compatibility changes the exact minima to **2,2,8**. The same
minima hold with the colored diagram guard.

## A general algebraic enumerator for compatible pairings

**PROVED IN NOTES; standard permutation-group lifting, applied with the
literal diagram guard.** This lemma works for any finite set X with an
involution R, not only tournament classes. Write Y=X/R, partitioned into
singleton fibres Y_s and double fibres Y_c. Let G be the group of weighted
graph vertex automorphisms of Y preserving these two sectors. If H or colors
are required, also restrict G to preserve those labels or matrices.

After choosing one name in each double fibre, the full R-equivariant
compatible permutation group is noncanonically

```text
C2^(|Y_c|) semidirect G,
```

where G permutes the double-fibre coordinates. For an involution g in G,
let f_s(g) and f_c(g) count its fixed singleton and double fibres, and let
p_c(g) count its two-cycles of double fibres. The polynomial enumerating
all compatible involutions J by their number of fixed ordinary elements is

```text
I_G(x) = sum_(g in G, g^2=1)
             2^p_c(g) x^f_s(g) (1+x^2)^f_c(g).          (L)
```

The identity is included as an involution. In particular the minimum fixed
count is `min_(g^2=1) f_s(g)`. If G is trivial, the complete polynomial is
`x^|Y_s| (1+x^2)^|Y_c|` and the unique minimizer is R.

*Proof.* A lift sends each R-orbit to the orbit prescribed by g. A double
fibre has two equivariant bijections to another double fibre; a singleton
has one. These choices give the displayed semidirect product. For a lift
to square to the identity, g must square to the identity and the two
transport bits around each exchanged double-fibre pair must agree. Thus
each such pair supplies two choices and no fixed ordinary element. On a
fixed double fibre the identity fixes two elements and the flip fixes none,
giving `1+x^2`. Each fixed singleton supplies one fixed element; exchanged
singletons supply none. Choices are independent on disjoint g-orbits,
proving (L) and its minimum. QED.

This formula classifies the lift once the actual graph symmetry group is
known. It is not an algorithm recovering that group or its wiring from
scalar class totals. Choosing arbitrary independent fibre flips gives an
ordinary-class action; only the global converse lift is asserted below to
have the specified physical reflection map. The [independent lift review](metagraph/PAIRING_LIFT_REVIEW.md)
includes synthetic exchanged-double-fibre controls and the finite census.

## The order-five diagram is vertex-rigid

The ten merged nodes have the following data. Codes use bit1 for i->j at
the corresponding lexicographically ordered pair i<j; each code is the
minimum over relabelings and converse. Degree counts all literal darts.

| code | H | automorphism order | dart degree | self-converse |
|---:|---:|---:|---:|---|
|0|1|1|1|yes|
|2|3|3|2|no|
|4|5|1|10|no|
|5|3|3|1|yes|
|8|9|1|9|yes|
|10|9|1|9|yes|
|12|11|1|11|yes|
|40|15|3|5|yes|
|41|13|1|13|yes|
|76|15|5|3|yes|

Only the degree1 bucket {0,5} and degree9 bucket {8,10} could be
permuted by a graph automorphism. Node0 has just one neighbor,10, whose
degree is9. Node5 has just one neighbor,41, whose degree is13. Therefore
0 and5 cannot be exchanged, and both are fixed. The now-fixed node0
distinguishes10 from8, so both degree9 nodes are fixed too. Every other
degree occurs once. **The weighted merged graph has trivial vertex
automorphism group**, even without imposing H preservation.

This statement allows permutations of parallel edges or darts fixing all
vertices; it is not a claim that every automorphism of the incidence object
is trivial. The exact neighbor rows are in the independent
[small census](audit/small.normal.json).

## All compatible pairings, including the self-converse classes

There are eight self-converse classes at order5. Vertex rigidity forces
each one to be fixed by any R-equivariant permutation compatible with the
diagram. Each of the other two merged nodes represents a converse pair:
ordinary codes {2,11}, of H3, and {4,9}, of H5. Within either pair one
may independently exchange the ordinary classes. These four choices give
the complete compatible permutation group C2 x C2; every element is an
involution and preserves H. The possible fixed-point counts are8,10,12,
with multiplicities1,2,1. The unique eight-fixed minimizer is R itself.

Converse has an actual presentation lift: reverse the displayed Hamilton
order and reverse every arc. On fixed-spine tilings this is reflection,
which commutes with chord complement and preserves both line colors.
Thus the minimum is attained with a physical dart map, rather than only
with an abstract permutation of class names.

Without the diagram guard there are16 H-preserving R-commuting involutions
at order5, with fixed-count histogram {4:1,6:4,8:6,10:4,12:1}. The unique
four-fixed one exchanges the two H9 self-converse classes and the two H15
self-converse classes, as well as both chiral pairs. Its H9 exchange
destroys the transitive node's incidence; its H15 exchange destroys the
degrees5 and3. These are the two missing constraints.

At order4 the merged degrees1,2,5 are distinct, so its vertex group is
already trivial. Its two self-converse classes must stay fixed; converse
exchanges the one remaining chiral pair. At order3 the graph is a single
edge, whose endpoints are the H1 and H3 self-converse classes. It has a
free diagram symmetry, but that symmetry changes H. On the nondegenerate
n>=3 complement-line tower, a free diagram-compatible pairing first fails
at order4; the first additional
obstruction to the **optimal H-preserving** pairing occurs at order5.

## Reproduction and the next interface

[check_small.py](audit/check_small.py) independently enumerates all labelled
orientations through6 by adjacent-transposition orbits, then enumerates all
physical Hamilton paths and builds the actual dart matrix. At3,4,5 it tests
every H-preserving R-commuting involution and every degree-compatible merged
vertex permutation. Its checks import no producer code. The metagraph
agent independently checked the same finite classification using its
different frozen generator before this note was written.

The producer generator uses a different, score-sorted canonical gauge.
Its merged labels32,68,97,258 correspond here to5,41,12,10, and its ordinary
partners256,64 correspond to11,9. The original draft combined those labels
with the unrestricted-minimum caption; the final figure audit caught and
repaired that mismatch. [The independent review](metagraph/PAIRING_REVIEW.md)
records the full transport and correction lineage. No incidence or group
claim changed.

This gives a specific obstruction to interpreting the evenness of A000568
as a free symmetry of the actual line diagram. The strongest survivor is
the abstract pairing already proved, or a diagram action with fixed
vertices and explicitly transported darts. Higher-order vertex symmetry
groups and their optimal compatible fixed-point counts are not classified
here. A future recurrence must specify whether it retains class cardinality,
H fibres, weighted wiring, or full marked incidence.
