# Independent review: complete odd self-converse family and faithful sidecars

**Verdict: ACCEPT within the declared self-converse family and orders3..11.**
No mathematical or implementation blocker was found in the complete-source
codec, individualization/refinement canonicalization, or minimum-leaf
automorphism test. This review does not promote the finite result to all
eleven-vertex tournaments or to arbitrary order.

Reviewed versions:

- `census.cpp`: `37312242df9c03eff19d3c9587e8d7a36942bec109e9d6424da5003944e9109a`
- `run.py`: `48ed6869acf5e05c235035c7f2aa37ee2ea2bbdde7070682b20612863ea93fd5`
- `full11.json`: `d9e66644e0a56776ac636602b943f857362e732baf08f41045b6e8be2456896f`

The author of this review independently read the source and the inherited
[Session89 codec proof](../../session89_consecutive_tournament_pairs_20261028/encoding/ODD_SELF_CONVERSE_LIFT.md),
checked the complete small source bank by a separate brute permutation
implementation, compiled and replayed the new engine through order7, and
reconciled the full order11 receipt. The complete order11 producer run is
the parent's recorded run, not a second full replay by this reviewer.
Separate selected-order11 endpoint/automorphism controls are owned by the
parent packet and should be reported separately when completed.

## 1. Why the source enumerates every odd self-converse class exactly once

The inherited proof supplies more than a matching class count. A tournament's
ordinary automorphism group has odd order, since an involution exchanging
two vertices would reverse their mutual arc. For a self-converse tournament,
its signed symmetry group has twice that odd order. Sylow2 gives an
anti-involution; on odd order it fixes exactly one vertex and pairs all
others. The two endpoints of each pair have a unique tail/head order.

After naming those tails p_i and heads q_i, the tournament is completely
specified by:

- the tournament T on the p_i;
- the symmetric cross-pair bits R_ij;
- the root-to-p_i bits d_i.

The code's `lift` implements the inherited five arc rules exactly and
checks the anti-involution on every ordered pair. There is no independent
pair-flip gauge after tail/head normalization. Relabelling the h pairs is
precisely the simultaneous S_h action on(T,R,d).

All anti-involutions of a fixed tournament are conjugate by ordinary
automorphisms. Sylow conjugacy initially gives a signed conjugator; if it
is negative, multiplying by one of the involutions makes it positive
without changing the conjugation. Thus forgetting the supplied involution
does not merge distinct input isomorphism classes. This is the load-bearing
reason the source's pair classes form a complete, nonredundant *unmarked*
self-converse bank; it is not inferred merely from the returned totals.

The producer uses one complete representative bank for T, then enumerates
every(R,d) bit string modulo Aut(T). Its `transform` applies the actual
vertex permutation to the symmetric R and vertex bits d. A key is accepted
only when none of its complete Aut(T) orbit is smaller. For an accepted key,
the stabilizer is counted over the entire group, not a truncated scan.

The source's centralizer order must not be confused with the full output
automorphism group. It is exactly C_Aut(W)(tau); multiplying its source orbit
size by h! gives the correct labelled source mass. The code checks

    sum_(accepted sources) h!/stabilizer = 2^(h^2).

This agrees at h=1,...,5 with2,16,512,65536,33554432. The quotient
|Aut(W)|/|centralizer| counts the conjugate anti-involutions only after the
full automorphism group has separately been established.

## 2. Independent verification of the inherited core bank

The bank read by `run.py` has SHA256
`feffa1e43ffc168af2e4f1fa6a7ac881438535da7876ef8720f4e85bd614668c`.
For each order1..5 the reviewer independently applied every vertex
permutation to each listed tournament, checked that these complete orbits
were pairwise disjoint, and checked that their union contained **every**
labelled tournament of that order. Automorphism orders were independently
recovered as n!/orbit_size and matched every saved `aut` value.

| Order | Listed classes | Labelled tournaments covered |
|---:|---:|---:|
|1|1|1|
|2|1|2|
|3|2|8|
|4|4|64|
|5|12|1024|

Thus the small core bank's completeness is verified independently here,
rather than trusting only its cardinality checks inside `run.py`. The runner
does report the bank hash and check its class counts; those checks alone
would not prove completeness for an arbitrary replacement bank.

## 3. The endpoint data and arithmetic are correctly typed

`ends(n,a)` is the usual subset recurrence for directed Hamilton paths ending
at a given vertex, with singleton value1. Applying it to the converse graph
gives the corresponding starting counts. The matrix entry C_uv sums the
product of left-path endpoint count and right-path start count over all
nonempty proper vertex splits with u on the left and v on the right.
This is the intended ordered two-path endpoint count; no edge u->v is
silently required across the cut. Diagonal entries are zero because the
two parts are disjoint. H is redundant here, being the sum of the endpoint
counts, so its omission from the stored canonical key loses no information.

The32-bit unsigned count storage is safe at the stated maximum n=11.
Endpoints are at most10!, and C_uv is also at most10!: concatenating the
two words gives a permutation with u immediately before v, uniquely
recovering the split. Each product inside the recurrence is bounded by
(|U|-1)!*(11-|U|-1)!, at most9!. Even the crude total H<=11! fits32bits.
These bounds do not automatically authorize increasing n in the fixed-size
implementation. The eleven-vertex tournament code uses55bits and fits64bits.

## 4. Canonicalization really is invariant under relabelling

Initial colours are the ordered endpoint pairs(S_u,E_u). Refinement appends
the old colour and the sorted multiset

    (colour(v), C_uv, C_vu),       v!=u.

The old colour is retained, so refinement only splits classes; colour names
are assigned by sorted complete keys. This is a deterministic, relabelling-
equivariant refinement. Repeated refinement terminates; a mere renumbering
after a split can require one extra pass but cannot cause a merge or cycle.

The first nonsingleton colour class is selected canonically. Every vertex
of that class is individualized in a separate branch; no branch is discarded
before a complete leaf. The doubling convention keeps the singled-out
vertex before the rest of its previous class and preserves all other class
orders. Consequently the search tree is equivariant under every sidecar
isomorphism. Different branch choices cannot produce duplicate final vertex
orders: at their first differing individualization the two vertices retain
opposite relative class order forever.

At each leaf the key includes **all** ordered(S,E) values and the full ordered
matrix C. Taking its minimum over this canonical search tree gives a valid
canonical sidecar key. It need not be the minimum over all n! arbitrary
vertex orders; the latter stronger property is not needed or established.
Equality of two full leaf keys itself exhibits a sidecar isomorphism, and
tree equivariance ensures isomorphic inputs attain the same minimum.

## 5. All minimum leaves, not one representative, establish automorphism faithfulness

Every sidecar automorphism acts on the search tree and permutes its minimum
leaves. Fix one minimum leaf ordering p. Any other ordering q with the same
full key defines a unique sidecar automorphism p_i->q_i. Conversely every
sidecar automorphism maps p to another minimum leaf. The action is free
because fixing a full vertex order fixes every vertex. Thus the number of
minimum leaves is exactly |Aut(sidecar)|.

The code counts **every** minimum-key leaf. If a later smaller key appears,
it resets the count and both tournament-code extrema; equal-key leaves
increase the count and update both extrema. Therefore `result.aut` is the
full sidecar automorphism order, and `lo==hi` holds exactly when every one
of those automorphisms preserves the tournament's arcs. Tournament
automorphisms always preserve the sidecar, so equality establishes
Aut(sidecar)=Aut(tournament), not merely inclusion of one selected subgroup.

Comparing the canonical full keys across the complete source bank separately
tests class injectivity. A repeated key is reported even if the associated
tournament codes agree, so an accidental duplicate source would also be
detected rather than silently accepted. This distinction between class
collisions and extra within-class symmetries is necessary and is correctly
retained in the receipt.

## 6. Replays, receipt reconciliation, and the precise consequence

The reviewer compiled `run.py --max-half 3` and obtained exactly the frozen
receipt's orders3,5,7:2,8,88 classes, with no class collisions or extra
sidecar symmetries. Both source hashes and the bank hash matched
`full11.json`. The complete order11 receipt was also independently checked
for all histogram totals, odd automorphism orders, labelled source masses,
class/sidecar equality, and empty collision arrays.

The parent producer reports at order11:

    279968 source classes = 279968 distinct sidecar classes;
    0 class collisions;
    0 extra sidecar symmetries;
    290128 canonical leaves, maximum81 for one input.

Together with the complete codec and the audited canonicalization, this
supports injectivity on **all self-converse eleven-vertex tournaments**
and automorphism faithfulness for each member of that family. It does not
exclude a non-self-converse eleven-vertex tournament sharing a sidecar with
one of these members, and it does not establish the conjecture for all
eleven-vertex tournaments or all orders. The self-converse promise remains
part of the finite reconstruction theorem.

The full census is a reproducible exact computation, not a Lean proof or a
hosted-CI claim. No further producer change is required by this review.
