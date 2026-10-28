# From rooted metagraphs to an exact, pruned class generator

Session87, 2026-10-28. Source pin:
`b8672af00bcc6d440a9a00a6915325faf39fcdff`.

**PROVED:** canonical-root selector lemma and its attachment specialization.
**FINITE-EXACT:** the accompanying standard-library implementation and bounded
checks reported below. **INHERITED:** rooted-metagraph completeness and canonical
prime/composite decomposition. No claim of method novelty, polynomial complexity,
wall-clock record, or a new scalar formula for A000568.

## 1. What is actually improved

The scalar Burnside formula counts tournament classes without constructing them.
Session79's root-switch operator constructs them, retaining all pointed states
and all their rerooting rows. The present implementation keeps that exact source
information but avoids computing most of those rows:

1. Canonicalize each parent once, retaining **all attaining transports**.
2. Enumerate binary attachment masks modulo the actual parent automorphism action.
3. Discard roots outside an isomorphism-invariant maximum-signature locus.
4. On the surviving locus, compute aligned parent-deletion keys only until a
   smaller key rejects the candidate.
5. Cache repeated literal deletion matrices together with all their canonical
   transports. Never cache just their degree sequence or H-value.

For child order seven, the degree-only selector reduces the full-row comparison
budget from 21,280 aligned-key requests to 600. Refining degree by the degrees of
out- and in-neighbors reduces it to **90**, while returning the same 456 ordinary
classes. These are measured operation counts, **not a claimed 236-fold runtime
speedup**: parent preparation, mask orbits, signature construction, cache storage,
and the unequal costs of individual canonicalization searches also matter.

The generic canonical construction-path method is established literature, not
new here: Brendan D. McKay, *Isomorph-free exhaustive generation*, Journal of
Algorithms 26 (1998), 306–324,
[author PDF, including errata](https://users.cecs.anu.edu.au/~bdm/papers/orderly.pdf).
The contribution here is a tested specialization and pruning of this repository's
aligned root-switch pipeline.

## 2. Inheritance and source comparison

All repository paths below are relative to the pinned root.

* `05-knowledge/results/session67_metagraph_enumeration_information_20261026/session67-root/COUNT_FREE_ENUMERATION.md`
  and `enumerate.cpp`: exhaustive fixed-Hamilton-path cube, relocation unions,
  then reflection/complement line bookkeeping. Exact and count-free, but its
  `2^((N-1)(N-2)/2)` presentation universe is different from the current parent
  attachment universe. No timing comparison with that C++ code is claimed.
* `05-knowledge/results/session79_tournament_enumeration_20261027/ROOT_SWITCH.md`
  and `check_root_switch.py`: `(parent,mask)/Aut(parent)` is precisely the pointed
  class bank. Each child contributes `Q_U=1_r(m_1,...,m_r)`, with `Q²=NQ` and
  `rank_F2 Q=a_N`. The existing code constructs the state bank by canonizing the aligned deletion pair of
  every raw extension and computes all N rerooted keys per pointed state.
* `05-knowledge/results/session79_tournament_enumeration_20261027/PRIME_CYCLE_INDEX.md`
  and `05-knowledge/results/session86_relations_attachment_carriers_20261028/A000568_COMPOSITION.md`:
  recovered aggregate prime indices recompose ordinary counts, but do not
  reconstruct prime representatives or individual transporter actions.
* `05-knowledge/results/session86_relations_attachment_carriers_20261028/81_27_TRANSPORTERS.md`
  and `CLOSED_ATTACHMENT_CARRIERS.md`: even identical literal nonempty pointwise
  stabilizers and all named subset-H values can miss between-port transport.

Closest mechanism: aligned rooted deletion. Hostile: the 81/27 decorated C3
pair. Corrected near miss: treating a minimum unrooted deletion card, mask
cardinality, or separate parent/mask minima as a pointed key. Restored sidecar:
all transports attaining the parent's canonical form.

Concept board: complete binary relation; pointed versus ordinary class;
invariant root selection; canonical child colors; action on attachment masks;
exact memoization versus observer collisions. The anchor is actual generation;
the niche is safe pruning; the wildcard is transfer to other relational classes.

## 3. Canonical-root selector lemma

Let `U` be a finite relational object on N vertices. Suppose deletion and the
attachment alphabet encode each pointed object `(U,v)` without loss. Let
`k(U,v)` be a complete, totally ordered invariant of pointed isomorphism.
Let `sigma(U,v)` be any totally ordered vertex invariant, with

    D(U) = {v : sigma(U,v) is maximal}.

Assume a construction stage supplies **one representative of every pointed
class**, as parent-attachment orbits do for tournaments. Keep its representative
`(U,r)` iff

    r is in D(U), and k(U,r) = min_{v in D(U)} k(U,v).          (1)

**Lemma.** Exactly one pointed class survives over each ordinary child class.

**Proof.** The nonempty set D(U) is invariant under isomorphism. The minimum of
its pointed keys is therefore intrinsic. Complete pointed keys are equal
exactly when the corresponding vertices are related by an automorphism of U.
Consequently the minimizing vertices form one vertex orbit, hence exactly one
of the supplied pointed classes. Distinct ordinary children cannot share a
pointed class. QED.

A candidate can be rejected immediately when any smaller competitor is found.
No full rerooting row, its matrix rank, or child canonical form is needed.
This is a selector of a root-switch support block, not recovery of that block
from the scalar rank. The proof does not use odd automorphism order, so it is
more general than the special binary-rank identity.

### Exact tournament key

Represent the new root's out-neighborhood in its deletion parent T by a bitmask
S. Minimize `(code(T), transported(S))` lexicographically under parent labelings.
The implementation first orders vertices by outdegree, then exhausts permutations
within equal-degree cells. Such orderings suffice because isomorphisms preserve
degrees. It retains **every** ordering attaining the least parent code and
minimizes the transported mask among that set.

For an already canonical parent, those attaining orderings are exactly its
actual automorphism permutations. They act on masks, so visiting one mask from
each orbit supplies exactly one pointed class. Independently minimizing the
mask under arbitrary permutations would destroy the alignment and is invalid.

Every production canonical search in the N-vertex generation stage has at most
N−1 vertices. A runtime guard enforces that bound. The representatives are raw
N-vertex matrices; when they become parents at the next stage, canonicalizing
them is legitimate parent preparation, not an omitted free operation.

## 4. Two inexpensive selectors, and a hostile to overcompression

The first selector is `sigma(U,v)=outdegree(v)`. For an attachment mask S on
parent T, it needs only

    |S| >= outdegree_T(w) + 1_{w not in S}    for every w in T.

The second selector is the lexicographically ordered tuple

    sigma(U,v) = ( d(v),
                   sorted(d(w) : v -> w),
                   sorted(d(w) : w -> v) ).                 (2)

Both are isomorphism-invariant, not complete pointed invariants. For example,
through order five the refined selector happens to leave one pointed orbit per
class, but at order six it leaves **59** candidates for **56** classes. The
exact key therefore rejects three additional candidates. At order seven it
leaves 481 for 456. The receipt includes a concrete adjacency code and two
vertices with identical maximizing signatures but different rooted keys.
The first refined witness is order-six raw upper-triangle code 8340, roots
5 and 3: both have signature `(3,(2,2,3),(2,3))`, but their rooted keys are
`(68,22)` and `(40,22)`. Thus they are not automorphic. This is why replacing
the exact tie breaker with the cheap signature is unsound.
Regular tournaments show a structural worst-case boundary: every vertex has the
same signature (2), so neither cheap filter separates their roots. No universal
asymptotic improvement follows from the successful small-order counts.

## 5. Finite-exact measurements

Universe: ordinary tournament classes, not converse-merged classes. All masks
are allowed; no Hamiltonian, prime, regularity, or arithmetic filter is imposed.
Known A000568 values are consulted only after construction.

| Child order N | Raw parent masks | Pointed mask orbits | Degree survivors | Degree competitor keys | Refined survivors | Refined competitor keys | Classes |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 1 | 1 | 1 | 0 | 1 | 0 | 1 |
| 2 | 2 | 2 | 1 | 0 | 1 | 0 | 1 |
| 3 | 4 | 4 | 2 | 2 | 2 | 2 | 2 |
| 4 | 16 | 12 | 5 | 4 | 4 | 2 | 4 |
| 5 | 64 | 48 | 15 | 12 | 12 | 6 | 12 |
| 6 | 384 | 296 | 82 | 70 | 59 | 18 | 56 |
| 7 | 3,584 | 3,040 | 717 | 600 | 481 | 90 | 456 |
| 8 | 58,368 | 54,256 | 10,934 | 9,117 | 7,164 | 693 | 6,880 |

The full Q row budget is N times the pointed column: 48, 240, 1,776, 21,280,
and 434,048 at orders four through eight, before the old raw state-bank key
requests. These comparisons do not equate a key request with a constant-time
operation. For the degree selector at seven, exact literal-matrix caching turns
600 competitor requests into 298 new canonical searches; parent preparation is
separately recorded. At order eight the refined selector makes 693 competitor
requests, of which 574 require a new canonical search; parent preparation takes
912 requests and 630 searches. The cumulative construction still visits every
parent-mask orbit, so this is not sub-output-size enumeration.

A controlled local timing test uses the same exact canonicalizer and
literal-matrix cache in both algorithms. Three alternating runs include parent
preparation but exclude the class-set audit and earlier levels. The full-Q
reference constructs the raw pointed bank and every Q row; it is **not** the
uncached historical executable. Median stage times were:

| Child order | Full Q, same cache | Refined selector | Ratio of medians |
|---:|---:|---:|---:|
| 7 | 0.2506 s | 0.04516 s | 5.55 |
| 8 | 8.933 s | 0.9316 s | 9.59 |

Both output class sets were independently compared after timing. These are
local measurements, not hardware-independent factors or asymptotic bounds.
Reproduce with `python3 -B benchmark_generation.py`; individual runs and the
comparison scope are in `benchmark.json`.

The tests independently: (a) compute every full Q row and verify constant
rank-one support blocks; (b) canonicalize full children only in an isolated audit
engine and compare exact class sets; (c) partition every labelled tournament
through six into full S_n orbits, covering all `2^15=32,768` labelled order-six
tournaments. Explicit exceptions survive Python `-O`. The production-only
mode disables all three audit paths to check that none secretly supplies state.

### Four vertices and six edges

For N=4 the independent labelled universe has `2^6=64` tournaments and four
ordinary classes. The rooted construction instead has 16 raw parent masks,
12 pointed classes, five maximum-degree candidates, and four accepted classes.
The neighbor-degree refinement leaves four candidates already, but still checks
ties inside their vertex orbits. The fixed-Hamilton-path metagraph uses three
free edges and eight presentations, yet another carrier. Equal small numbers
or equal bit dimensions do not identify these different quotients.

## 6. The 81/27 hostile now changes the actual attachment workload

Let A0 be `1 -> C3`, A1 its literal converse on the same four labels, and put

    T = C3[A0,A0,A0],       U = C3[A0,A0,A1].

The inherited proof identifies their full groups as `C3^3 semidirect C3` and
`C3^3`. They have the same 67 **labelled** cyclic triples, identical H on all
4,096 named induced subsets, and full-set `H=462915`. Their literal pointwise
stabilizers agree on every nonempty named tuple. But their binary attachment
mask orbit counts are now explicitly checked to be **176 and 512**.

There is a transparent calculation. A0 has two choices on its singleton and
four subset orbits on its triangle, hence q=8 subset colors. For T the outer
C3 rotates the three child slots, giving `(q^3+2q)/3=176`. In U, the different
ordinary child color kills that outer rotation, giving `q^3=512`. For a single
marked port the corresponding counts are 2 and 6. The transporter from named
singleton 0 to singleton 4 has size 27 in T and zero in U.

Thus a memo table keyed by the full subset-H observer, cyclic-triple hypergraph,
abstract child groups, or local point stabilizers would already merge unequal
attachment workloads. The actual child colors and their common-frame transport
are not optional. The checker verifies the listed permutations preserve arcs,
direct mask orbits agree with their Burnside sum, all subset-H entries agree,
and literal single-point stabilizers agree. The completeness of the listed
groups uses the inherited Session86 degree/incidence proof; a separate Session87
audit checks full groups by unrestricted adjacency backtracking.

An all-order variation uses `A_m=T_m -> C3`, m>=1. It has q=`2^(m+2)` subset
orbits, and its converse has the same action shape but a different ordinary
child type. The chain vertices are individually fixed and the triangle rotates, so its
automorphism action is C3. Its SCC order distinguishes it from its converse.
Canonical prime-root decomposition therefore gives attachment counts
`(q^3+2q)/3` versus `q^3` for the corresponding two blowups. Reversing the last
module preserves every named induced H by reversing maximal path runs in that
module. This explains the family mechanism; the finite receipt checks m=1 only.

### Direct four-vertex scaffold bridge (same-triangle variation)

Keep the same singleton root and the same **named directed triangle**, but change
`A=1 -> C3` to `Aprime=C3 -> 1`. A scaffold remembering only that root and
triangle forgets precisely this source/sink bit. This is not the literal
converse convention used in the finite receipt: let t swap two triangle labels
in the changed leaf. Then `Uprime=t(U)` is the desired same-triangle variation.
The transposition normalizes the child's C3, so the displayed literal full
groups remain G and H; attachment-orbit counts remain 176 and 512.

For any named subset S, its triangle part and its image under t have the same
cardinality. C3 is transitive on the k-subsets of its three vertices for every
k=0,1,2,3. Consequently tS lies in the Aut(T)-orbit of S, and

    H(Uprime restricted to S) = H(U restricted to tS)
                             = H(T restricted to tS)
                             = H(T restricted to S).

The cyclic-triple sets also remain identical. Thus the four-vertex scaffold's
lost source/sink bit really changes ancestor slot transport and the number of
inequivalent binary attachments; this is a structural connection, not merely a
coincidence among dimensions or small cardinalities. This paragraph is an
all-subset proof for a relabeled variation. A separate Session87 audit also
checks this variant by flipping only spokes 8–9, 8–10, 8–11, verifying the full
group by adjacency backtracking and every named subset H; see the sibling
`audit/scaffold-variant.json`. The generator receipt itself retains the Session86
literal-converse convention.

## 7. Prime/composite grammar and safe transfers

Canonical SCC/strong-module decomposition can construct composites from actual
prime representatives: retain the prime quotient relation, the ordinary child
type assigned at every slot, and the quotient action. Enumerate assignments
modulo the actual slot action. At a transitive root, retain the intrinsic SCC
order, not an unordered multiset. Canonical decomposition avoids arbitrary
parenthesization duplicates. This can be combined with the attachment generator,
but it does not manufacture the missing prime representative bank.

Aggregate prime cycle indices suffice for ordinary composition **counts**.
Session86's recovered index through twenty was obtained from Burnside, so its
recomposition is not an independent enumeration discovery. Likewise a scalar
rank of Q counts blocks but cannot regenerate the block representatives needed
for the next stage. The present speedup removes redundant comparisons while
retaining those representatives and their exact actions.

The selector lemma transfers to finite graphs, digraphs, hypergraphs, and
colored relational structures when their deletion/attachment encoding and
isomorphisms are correctly specified. Use the correct attachment alphabet and
action; the root invariant need only be isomorphism-invariant. It does not
transfer an unlabeled quotient to problems whose integers or port roles are
fixed. Square-sum path states must retain actual integer endpoint names and
connectivity; graceful-labeling states must retain named label/difference
resources and legality. Abstract tree shape symmetries cannot erase these
coordinates unless they preserve all permitted future guards.

**Reproduction:** run `check_generation.py --max-order 8` and then with
`--selector neighbor_degrees`; repeat each command with `python3 -B -O` and
compare JSON. `--production-only` skips the full-Q, child-canonicalization and
labelled-universe audit engines. Default maximum order is seven; hard limits
are eight for generation and six for the full labelled-universe check.
