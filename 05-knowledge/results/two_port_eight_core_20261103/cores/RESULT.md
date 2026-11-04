# Complete eight-vertex core universe and independent coverage

**FINITE-EXACT COMPLETE.** The universe is all2^28=268435456 labelled simple
graphs on vertices0,...,7, with pair bits in lexicographic order. There are
53470067 graphs of minimum degree at least3. Necessary planar unit-distance
filters leave724080 labelled graphs, in exactly38 isomorphism classes.
This is not a claim that all38 are geometrically realizable.

The filters retain full adjacency. Any two distinct points have at most
two common unit neighbours. Every vertex neighbourhood lies on its unit
circle; its induced edges move by60 degrees, so its graph is a union of
paths and six-cycles. Allowing paths longer than can actually fit is harmless
overinclusion; no such path is used as a realization certificate.

`census.cpp` visits all labelled graphs by binary reflected Gray code,
updating the unique changed edge and its two degrees. It checks every
eligible adjacency. After sorting the survivor masks, it expands all8!
permutations of each least remaining representative, checks orbit-stabilizer,
and deletes its entire orbit from the survivor set. Overlap, a missing image
or an unclassified survivor is fatal. Complete representative rows and orbit
sizes, not just a hash, are saved in `order8.json`.

## A different exhaustive enumeration

`independent.cpp` instead visits every one of the2^21 graphs on the named
vertices0,...,6 in ordinary integer order and attaches vertex7. Deleting a
vertex from an eight-vertex minimum-degree3 graph leaves minimum degree at
least2; every resulting degree2 vertex must be adjacent to7 when restored.
The metric filters are hereditary, so rejecting a failing seven-vertex
parent cannot discard a later admissible graph. The new vertex needs at
least three neighbours. These are the only inherited pruning conditions.

The alternate metric implementation checks codegree at most2 and
**bipartiteness of every neighbourhood**. Codegree bounds each degree
within a neighbourhood by2. A four-cycle in a neighbourhood would make
its opposite vertices share three neighbours including the centre, so
is excluded already. With at most seven neighbours, the only remaining
possible cycle in a bipartite degree2 graph is C6. Thus this alternate
test is equivalent to the producer's filter at the stated orders.

The audit visits1052793 minimum-degree2 parents;217290 pass the hereditary
filters. It tests5392215 admissible attachments and finds exactly724080
survivors. Independently expanded permutation orbits of the38 representatives
are disjoint and equal that complete survivor set. This checks the universe
by a different enumeration order, pruning argument and local graph test.

## List certificates and exact geometric exclusions

Thirty representatives have a nonzero coefficient of their graph polynomial
product_(uv in E,u<v)(x_u-x_v) with every exponent at most2. Finite Lagrange
interpolation proves this sufficient for arbitrary three-lists. The producer
uses multiplication modulo(x_0^3,...,x_7^3). The independent
[coefficient check](../recovery/LIST_CERTIFICATE.md) expands every choice of
one endpoint from each edge, checking all capped coefficients and support
sizes:548864 choices across the38 graphs.

The eight remaining types have no such certificate, which alone is not a
list obstruction. Six contain a rigid Moser spindle with an exterior vertex
having at least three prescribed contacts; two contain a six-vertex fan
forcing two distinct vertices to coincide. [The geometry proof](../geometry/PROOF.md)
therefore excludes all eight exactly. No approximate distances or unproved
realization test is used. The30 certificates may safely include additional
nonrealizable graphs.

## Reproduction

```text
python3 run.py --order 7
python3 run.py --order 8
python3 independent.py
python3 ../recovery/check.py
```

Repeat with `python3 -O`; all validation uses explicit failures. Both complete
eight-vertex enumerations and the full coefficient audit have matching
normal/optimized receipts. The order7 producer control reproduces the
previous six representatives and all8190 labelled survivors. Compiler
optimization is-O3 in both Python modes. Sources and result inputs are
hashed before and after each run. The proof of completeness and of the
list-colouring implication remain necessary alongside these computations.
