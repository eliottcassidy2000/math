# Independent audit: the regular-eleven correction

24 October 2026. Status: PROVED structural arguments; FINITE-EXACT independent
computation. The frozen production packet passes. This audit does not publish.
Ordinary oriented classes are retained; converse classes are never merged.

## Findings

The independently authenticated almost-regular ten-vertex atlas has 13,333
classes, with complete automorphism-order histogram

    1: 13,282; 3: 38; 5: 8; 9: 5.

Its labelled orbit mass is 48,251,508,480. An independent residual-degree
coefficient calculation gives exactly the same entire labelled universe.
Distinctness, exact full group actions, and this equality prove completeness
without treating 13,333 as an expected-value input.

Every one of the 146,663 deleted-root destinations agrees with production.
The diagonal trace is 13,453, giving 1,223 regular-eleven classes and regular
error 133,210. Independently enumerating the FULL automorphism group of one
regular representative per component gives

    1: 1,205; 3: 6; 5: 7; 9: 1; 11: 3; 55: 1.

Those groups agree with every parent stabilizer and diagonal orbit size.
The error histogram is 0:4, 22:7, 66:5, 88:2, 110:1205.
There are 229 self-converse parents and 229 self-converse regular classes.

Combining only with the separately audited nonregular correction gives

    E_s,11 = 133210 + 9166 = 142376;
    E_11   = 133210 + 19480 = 152690.

This does NOT evaluate U11 or a11. The identity a11=(U11-152690)/11 remains a
conditional route once U11 is supplied. There is no full order-ten or eleven
census, no target-eleven canonicalizer, no imported target count, and no claim
that an unmarked metagraph or scalar count alone carries this information.

## Independent construction and computation

`independent.cpp` is not a copy of the production cyclic normalizer. It starts
from one cell, refines by INCOMING counts into ordered cells, and enumerates
ALL within-cell permutations. It minimizes an upper-triangle bit encoding,
not the producer's outgoing-row tuple and not intrinsic rotations. This is a
complete normal form because isomorphisms preserve the intrinsic ordered
partition, and every remaining label permutation is searched.

Full automorphism lists use another algorithm: depth-first partial bijections,
choosing the unmapped vertex with fewest compatible images. Exact adjacency,
canonical color, and common-outneighbor counts are necessary filters only.
Every surviving full bijection is returned, so no mere subgroup is mistaken
for the full group. Every parent's supplied permutation set equals this exact
independent set. For the regular children this routine finds isomorphisms to
self; it does not compute or compare any target-eleven canonical form.

The generator scans each source parent whose degrees permit an extension,
then tests EVERY one of its 2^9 attachment masks by the actual resulting degree
sequence. It does not reuse production's must/optional mask construction or
subset-orbit representatives. It obtains 10,252 possible source parents and
135,418 eligible masks. Every resulting independent key lies in the audited
atlas, and every output key is generated. This rechecks the entire selected
universe, not just a chosen sample.

The order-eleven run makes 308,747 parent normalizations and examines
27,741,098 unrestricted within-cell labelings. It performs 181,334 nodes of
full-automorphism backtracking and 30,636,390 explicit checks. Every reroot
index, multiplicity, component, converse index, parent action, and pointed
source orbit count is checked, not only the aggregate answer.

The default-compiler, O2, and O3/NDEBUG reports match byte-for-byte. The same
code gives the lower controls:

    regular order7: parents5, trace21, classes3, error14, labelled mass2640;
    regular order9: parents85, trace135, classes15, error630, mass3230080.

`degree_dp.py` is a separate maximum-degree deletion recursion using literal
INDEXED incoming-neighbor subsets. Production instead deletes a minimum-degree
vertex and compresses equal-degree subsets by binomial coefficients. No Landau
filter, unlabeled count, or group calculation enters this independent recursion.
It gets 191,474,240 for one fixed assignment of five degree4 and five degree5
labels, then multiplies by binom(10,5). The all-degree5 vector on eleven labels
independently has the same count. Every labelled degree-vector count at orders
three and five is also checked against literal exhaustive edge enumeration.
Normal and optimized Python reports agree exactly.

## Structural proof audit

Deleting a specified vertex from an almost-regular ten-tournament gives degrees
in {3,4,5}. Old degree3 vertices must beat the new root, old degree5 vertices
must lose to it, and degree4 vertices are free; the incoming root subset has
size4 or5. These conditions are necessary and sufficient, with five vertices
of each resulting degree forced by the total edge count. Actual source-group
orbits give pointed types exactly, not unpointed types.

Stable refinement only splits the initial two five-vertex degree cells. Each
stable cell induces a regular tournament and hence has odd size1,3,or5. The
unique regular three- and five-vertex tournaments have intrinsic directed
successor cycles and rotation-only automorphisms. For order5, if the root's
outneighbors are a->b and inneighbors c->d, degree2 forces the cyclic order
root,a,b,c,d. The successor test in the code agrees with this proof.

An isomorphism preserves the intrinsically ordered stable cells and their
successors, hence appears among the product of cyclic rotations. Each initial
five-cell can split as5,3+1+1,or1+1+1+1+1, so at most25 rotations suffice.
The FULL crossing matrix must be retained; a gain kernel or crossing-degree
shadow alone is not an isomorphism certificate. Production correctly compares
all rows, and its minimizers give the full group by transporter composition.

An almost-regular even parent has exactly one regular extension. Parent
isomorphisms extend over its distinguished root because adjacency is forced
by degree. Thus almost-regular parent types are exactly regular pointed types.
Rerooting recovers their vertex-orbit components. Each row contains the orbit
sizes and each component has diagonal trace11. Therefore R11=tr(D)/11 and
E_regular11=11(M10-R11). Orbit-stabilizer gives the exact child group order
from each parent, independently checked here against all target automorphisms.

## Analytic exclusion of error44

The new orbit-partition lemma also passes an independent proof audit. Actual
Aut-orbits are odd, internally regular, and mutually biregular; crossings of
coprime orbit sizes must be uniform. Of the twelve odd partitions of eleven:

- 9+1+1, 7+3+1, and 7+1^4 cannot give a fixed vertex outdegree5.
- 5+3+3 forces the five-orbit to beat exactly one triple; that triple then has
  internal degree1 plus five outgoing edges, already exceeding degree5.
- 5+3+1^3 forces each fixed vertex's degree in its induced triangle to be0 or2,
  contradicting total triangle degree3.
- 3+3+3+1+1 gives a fixed degree equal to a multiple of3 plus0 or1, never5.
- 5+1^6 forces each fixed vertex's internal degree to be0 or5. Total degree15
  would require three degree5 vertices, although at most one can beat all five.

Only 11, 5+5+1, 3+3+1^5, 3+1^8, and 1^11 survive. Thus the numbers of orbits are
1,3,7,9,11, analytically excluding five orbits and regular error44. The audit
does not mistake the observed empty histogram bin for this general proof.

## Converse and scope

The converse involution on full parent matrices is independently normalized;
all indices and induced components match production without merging them.
For a self-converse odd tournament, automorphisms have odd order, while the
extended group of automorphisms and anti-automorphisms has twice that order.
An involution in the extended group is an anti-automorphism. It fixes exactly
one vertex: parity forces a fixed vertex, and two fixed vertices would reverse
their mutual arc. Each invariant odd Aut-orbit contains such a fixed vertex,
so exactly one orbit is invariant. These are precisely the self-converse
pointed types. This proves the producer's general odd-order self-converse
pointing equality, and the regular-extension bijection gives the parent equality.
No novelty claim is made.

## Hostiles, lineage, and portability

Ten independent controls run normally and optimized: pristine and arbitrarily
relabelled positive cases, plus rejection of a closed proper subgroup, a
nonliterally identical isomorphic duplicate, a one-sided arc toggle, a missing
parent, a changed reroot multiplicity, a wrong converse, a wrong component,
and a missing state row. All pass.

An early producer SOURCE_PINS metadata field named an obsolete nonregular
manifest. The archived source bytes were already correct. This audit requested
repair; the frozen producer now pins the actual final nonregular manifest and
checks its theorem and numeric receipt entries. No mathematical or numerical
correction was required. `SOURCE_PINS.json` identifies the exact final freeze.

`replay.py` copies ONLY manifest-enumerated files into a fresh temporary directory,
then reruns all three independent compilation modes, both Python modes, semantic
hostiles, production's complete two-mode replay, both frozen validators, and
byte-corruption manifest controls from an unrelated working directory. The
replay receipt is deliberately outside its own hash manifest to avoid recursion.
