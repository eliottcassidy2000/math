# Independent review of the order-ten census and joint response controls

Status: **FINITE-EXACT — accepted within the inherited complete bank**.
Reviewed by the independent rank_audit lane, 2026-11-03. No blocking defect
found in the mathematical reduction, canonicalization, integer bounds, or
frozen receipts. Universal two-port reconstruction remains OPEN.

## What the complete computation establishes

The completed [receipt](../census/full10.json) processes all9,733,056 accepted
order-ten classes in the pinned inherited bank. It finds no equality of the
H,S,E,C sidecar between different classes, and no sidecar automorphism which
fails to preserve the tournament. Therefore that sidecar is injective on
**labelled** order-ten tournaments. Together with the earlier audited census
through nine, any universal two-port counterexample has at least eleven old
vertices. This is a finite lower bound on a counterexample, not an all-order
proof and not an order-eleven exhaustive computation.

The implication from observations is the exact Session27 tomography in
[proof.md, sections7-8](../../session27_repunits_lattices_20261022/session27-response/proof.md):
actual at-most-two-port queries recover H,S,E and the **ordered** C matrix.
The census does not replace ordered C by its symmetric part. Nor does it
claim HSEC contains every coefficient of the full two-port response.

The bank's completeness and one-representative-per-isomorphism-class property
are inherited from
[consecutive_flags_order10_20261028](../../consecutive_flags_order10_20261028/verification10.json).
They are not re-established merely by matching A000568, the automorphism
histogram, or the labelled mass. The new wrapper verifies both compressed
and uncompressed bank hashes; decoding checks strict increasing flag
addresses, valid parent indices, exactly9,733,056 entries, and end of input.
An explicit limited run is marked incomplete and cannot certify order ten.

## Why the two-stage shortcut is exact

The first-stage key is the complete sorted multiset of the ten pairs(S_u,E_u).
It is stored as integers, with no hashing or possible hash collision. Any
HSEC isomorphism preserves this key. Consequently distinct endpoint-key
groups cannot contain an HSEC collision.

A group containing exactly one tournament class can bypass C only when
its ten endpoint pairs are also distinct. Class uniqueness excludes an
isomorphism to any other class; distinct pairs force every sidecar
automorphism to fix every vertex. Thus both required conclusions are
already certified in this case. If either condition fails, the program
constructs the full C matrix and canonicalizes it. In particular a unique
class with repeated endpoint pairs is not incorrectly treated as rigid.

The endpoint counts use the exact subset recurrence for Hamilton paths
ending at a vertex. Starting counts are ending counts in the converse.
For each nonempty proper subset A, the C contribution is

    C_uv += E_u(T[A]) S_v(T[V\\A]).

Disjoint support and nonempty slots are built into the subset partition;
each ordered two-path cover occurs exactly once.

## Canonicalization and the class-to-labelled implication

Initial vertex colors are the exact ordered endpoint pairs. Refinement
uses the old color together with the sorted multiset of triples
(color(v),C_uv,C_vu). This is invariant under every sidecar isomorphism.
Including the old color prevents cell merging. Sorted integer keys assign
colors consistently between isomorphic inputs; once no cell splits remain,
the numbering is stable. Every permutation within each remaining color
cell is visited, and no permutation preserving the sidecar can move a
vertex between distinct cells. This is exhaustive canonicalization after
an invariant restriction, not an assumption that refinement alone decides
isomorphism.

At every leaf the key contains all endpoint pairs and all ordered C entries
in that vertex order. All leaves realizing the minimum key are retained
through their count and minimum/maximum tournament adjacency code. Relative
to any fixed minimizing labelling, these leaves are exactly the sidecar
automorphisms. Therefore mincode=maxcode is equivalent to every sidecar
automorphism also preserving tournament adjacency. The reverse inclusion
Aut(T) subset Aut(HSEC(T)) holds because all these counts are intrinsic.

Distinct classes must have distinct canonical sidecar keys. Within a class,
the preceding automorphism equality ensures that any labelling ambiguity
preserves T itself. These two properties together establish labelled
injectivity. Checking only distinct class keys would be insufficient.

The receipt records:

- 9,726,272 endpoint-multiset groups;493 contain multiple classes; largest49.
- 227,764 classes with a repeated endpoint pair.
- 227,772 full C canonicalizations and3,272,365 canonical leaves; maximum362,880.
- Zero sidecar class collisions and zero extra sidecar automorphisms.
- Sidecar automorphism histogram identical to the inherited tournament
  histogram, and labelled mass35,184,372,088,832=2^45.
- Zero classes with exact S=E. This last fact is finite evidence about
  even-order converse separation, not a universal theorem.

## Integer and implementation boundaries

For n=10, every endpoint DP value is at most10!. Every C entry and every
nonnegative partial sum is at most the number of permutation/cut pairs,
9*10!=32,659,200. Thus uint32_t arithmetic, including each product in the
subset convolution, is safe. Tournament codes require45 bits; accepted
addresses are checked below2^27. The uint64_t counters, automorphism
counts, canonical-leaf totals, and mass are safe: even the crude total
leaf bound9,733,056*10! is below2^64. These bounds are specific to the
implemented order and must not be silently extrapolated.

The comparison keys are full vectors. The prefix optimization in the
independent brute canonicalizer skips an ordering only when its complete
vertex prefix is already lexicographically greater than the best prefix.
Its edge suffix cannot improve that ordering, so the rejection is sound.
This control still enumerates all10! vertex permutations. The refined and
unrestricted canonical representatives need not have identical integer
keys because their ordering conventions differ; their automorphism count
and the presence/absence of extra tournament codes are the relevant common
invariants, and these are the comparisons made.

## Full response and the three-path information

The formula in [algebra/PROOF.md, A4](../algebra/PROOF.md) was reviewed by
splitting every Hamilton word at its two ports. With x->y, the x-before-y
middle run can be empty; the y-before-x middle run must be nonempty. The
three successive old runs account bijectively for every term. The factors

    (1+e a)(1+(1-a)^t G c)(1+(1-c)^t s)
    +(1+e c)((1-c)^t G a)(1+(1-a)^t s)

retain each old support and all joint endpoints before extracting[V].
Singleton middle runs are included. The full formula has Boolean degree
at most two in each mask. Hence values on every pair of masks of weights
at most two determine the complete response by ordinary subset Möbius
inversion. Equality of that response need not individually recover every
entry of the unsymmetrized three-path tensor D; the proof makes no such
claim.

In [response/check.cpp](../response/check.cpp), the literal path enumerator
walks every permutation and every one- or two-cut position. Its M,C,U,V,D
indices agree with the path endpoints specified in the proof (code V is
proof W). The independent squarefree implementation computes path matrices
on all subsets and multiplies entries only on disjoint nonempty supports.
The response evaluator includes both port orders, both exterior singleton
cases, and all three-path terms. Its values are compared to a separate
Hamilton subset DP on the actual enlarged tournament.

The frozen normal and optimized tensor receipts agree. Their eight
controls are C5 and its converse, four prescribed bank-ten classes, the
transitive ten-vertex tournament, and a universal source over the cyclic
three-by-three wreath tournament. The last two give positive symmetry
controls of orders1 and81. Each ten-vertex control also compares census
S,E,C against literal permutation/cut counts and refined canonicalization
against all10! permutations. The four sampled bank indices are0,100000,
4866528,9733055; these are prescribed controls, not another exhaustive
ten-vertex response census.

Each control evaluates the complete low-weight battery plus nine
high-weight Möbius reconstructions:265 queries per order-five control and
3145 per order-ten control,19,400 actual queries in total. All M,C,U,V,D
entries agree between the literal and squarefree routes. The two C5
orientations have distinct actual query batteries, so the hostile is
detected. The response integers fit uint64_t: at most12! Hamilton words,
and at most choose(9,2)*10! ordered three-path covers. The signed
interpolation intermediates also fit int64_t even under the crude bound
3136*16*12!.

## Frozen artifacts reviewed

The current source hashes match their saved receipts. The two response
receipts are identical; the algebra normal/optimized receipts are also
identical and pin the current checker. The complete long census was not
rerun by this audit: review used its finished receipt, its entire source,
the separate literal controls, and the inherited bank certificate.

SHA256 of reviewed result files:

```text
36cc2f622c5d0fc6e3af14c7f9ba80735c20a1449c8e21430dc936a20614270f  census/full10.json
c37f9d0cdce26d009f28348ccaec52424ba42269850e2f10bfe075f009d7ef3d  response/tensors.json
a4c27c3b78188c20e2fef2c7f8fc3b50631ab97cf6b8babdf022f5e5b090c1e6  algebra/normal.json
```

The strongest supported conclusion is finite labelled faithfulness through
ten, together with the separately proved all-order odd-module obstruction.
Neither statement proves universal two-port reconstruction.
