# Independent checkpoint-two audit

**ACCEPTED:** exact frozen producer identities in `REVIEWED_PRODUCERS.json`.
All construction-specific proofs were read in full. This checker is written
separately, imports no producer code, and has no sibling-directory dependency.
The separate checkpoint-one audit remains frozen; its proof graph is not
silently extended to these files.

## Curvature, capacity, and endpoint geometry

The graph-edge upper bound gives only d_X<=L d_G. Uniform separation (or
separated occupied sites plus finite multiplicity) is essential for packing.
Disjoint OPEN radius-delta/2 balls and containment in radius Lr+delta/2 prove
the stated constant-curvature area inequality. Taking exponential rates proves
kappa L>=log b, or log2/5 for the explicitly certified five-native-step supply.
This is necessary, not sufficient for a Collatz realization or convergence.

The dyadic outer radius formula, the strict additive log(5/4) upper bound,
and its asymptotic equality are correct. Radial efficiency does not imply
pairwise quasi-isometry: cross-root near-neighbor pairs have unbounded tree
distance and constant ambient distance. Properness does not fix that loss.

For 0<q<1/2,0<W<=1 put g=W(1-2q),s=log(1/q). Across the LCA split,
|dx|>=g q^h and heights q^(h+a),q^(h+b); the asinh distance formula gives
d>=s(a+b)+2log g. Ancestor pairs instead have d>=s(a+b), which is stronger.
The upper bound follows from2cosh d<=(W^2+2)q^(-(a+b)). Since g<1,
C=-2log g is positive, so the lower estimate covers every case.

The edge-circle center sign beta^2+q^2-1<0, beta=W(1-q)/2, proves monotone
height and horizontal endpoint confinement. Disjoint outgoing horizontal
intervals at one level and disjoint open height bands across levels establish
noncrossing without relying on a false general geodesic-height assertion.
Positive minimum height of a compact set gives properness. The same-level
nearest-center formula and vertical separation prove the positive delta bound.
Small-W siblings can be closer than edges: equal edges do NOT preserve the
old induced-minimum-distance property. The producer explicitly retains this
boundary. Taking q up to1/2 and W down to0 gives edge length tending tolog2;
each fixed model remains separated, but separation is not uniform across
the family. The prescribed-length inversion is valid for every L>log2;
attainment at equality remains open in the packet.

Independent rational all-pair checks run through depth6 for nine models,
and exact dyadic extreme cosh radii through depth12. The native inverse set
is regenerated through radius26 and every vertex's first-hit clock is forward
replayed. Separate lattice BFS confirms counts704 versus694 (difference10),
2263 versus2107 (difference156), and first violations in the checked prefixes.
These finite tests do not prove the all-level statements; the proofs do.

## Square paths, common blocks, and the cubic response

The general compatible-subset expansion and n25 polynomial are inherited.
The common-directed-arc intersection is a directed linear forest. Its maximal
blocks are traversed whole in every bank member. Ordered block joins map
injectively to distinct endpoint arcs, giving n-b+|J| union size. A reverse
connector requires each endpoint simultaneously first and last in its own
fixed block, hence both blocks singleton. Shared internal arcs cannot be
reversed by any member. This proves both directions of the singleton test;
nontrivial opposite block joins use different original vertex pairs.

The source-forced blow-up lemma is valid: indegree-zero first vertex forces
the global start, so no block interior can start a Hamilton order. Every
internal predecessor is therefore mandatory. The directed-cycle hostile
shows why the source condition cannot be removed. The leaf18 condition
excludes Hamilton cycles in the stated n18..30 range, giving the d>=2 pair
bound there, not for arbitrary guards.

Independent unpruned all-start DFS enumerates EVERY directed path for every
n1..30. Independently enumerating compatibility cliques with their arc unions
recovers every coefficient vector, including n25 matching structure and
n26's two triples. At n26 the factor is2^291 and guard coefficients are
17179857952,10240,928,64. Pair unions are27:8,29:2,30:2; both triple unions
have29 arcs. Orders27..30 and the nonmonotone maxima5,5,3,3 also match.
The producer's second active-orientation algorithm is additional independent
computational structure, not silently treated as part of this checker.

Abstract normal-form controls exhaust one/two/three-path banks through n4
and one/two-path banks at n5. Source-forced blow-up is checked for all16
three-vertex loopless digraphs whose specified vertex0 has no incoming arc,
using block lengths1,2,3. Source-free directed cycles give hostile counts3,4.
All counts concern GUARDED Hamilton paths, not the ordinary tournament count.

## Replay and limits

Run `python3 -B verify.py`; it checks file/source hashes and replays normal
and optimized computation from an unrelated temporary cwd. Proof snapshots
and exact source pins are included. No absolute scratch path, network fetch,
Git database, external package, or producer checker is needed. Explicit
exceptions survive `-O`. No all-size square-sum Hamiltonicity, complete Collatz
embedding, global aliquot classification, or Hadwiger-Nelson improvement is
claimed. Primary area/geodesic notes were independently inspected; the
producer's separate metric citation is not misreported as an audit fetch.
