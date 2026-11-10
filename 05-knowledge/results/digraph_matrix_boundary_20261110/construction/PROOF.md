# A compact binary construction for strong oriented path totals

**PROVED.** Every integer h>=3 is the Hamiltonian-path total of an explicit
strong oriented graph on `bitlength(h-2)+2 = ceil(log2(h-1))+2` vertices.
This is a size refinement and a direct constructor, not a new claim of
oriented-spectrum completeness. Session235 already proves that completeness
and supplies strong fixed-count examples of unbounded order. No external
priority or minimum-order claim is made.

## C1 The binary deletion family

Fix k>=1. Use a transitive chain `0 -> 1 -> ... -> k-1`, adding every
forward chain arc, and two additional named vertices a,b. Add every
`i -> a`, add `a -> b`, and add `b -> i` exactly for i in a chosen subset S
of the chain. There are no other arcs. Call this graph G(k,S).

It has no loops or reciprocal pair. Its only missing pairs are `{b,i}`
for i not in S. Thus it has exactly `k-|S|` missing pairs. It is strong
if and only if `0 in S`: that condition supplies the spanning directed
cycle `b,0,1,...,k-1,a,b`, while otherwise vertex0 has indegree zero.

The exact count, also when the graph is not strong, is

    H(G(k,S)) = 1 + sum_(i in S) 2^(k-1-i) + [0 in S].            (C1)

To prove this, separate paths according to whether they use `a -> b`.
Those that do compress this ordered pair into one vertex w. The result is
the transitive tournament on `0,...,k-1,w` with sink back-arcs `w -> i`
for i in S. Its count is `1+sum_(i in S)2^(k-1-i)` by
[THM-5071 (4a), weighted-cycle-formula-for-digraph-path-totals](../../../../01-canon/theorems/THM-5071-weighted-cycle-formula-for-digraph-path-totals-and-one-reciprocal-pair-exclusions.md).
The compression and expansion are mutually inverse for these paths.

For a path avoiding `a -> b`, a must be last because it has no other
out-neighbor, and b must be first because it has no other in-neighbor.
The remaining chain has only its increasing Hamiltonian path. Therefore
there is exactly one such path, `b,0,...,k-1,a`, if `0 in S`, and none
otherwise. This proves the additional term and (C1).

Equivalently, begin with the tournament G(k,{0,...,k-1}), of count `2^k+1`.
Delete any subset of the arcs `b -> i`, keeping `b -> 0`. Each deletion
subtracts exactly `2^(k-1-i)`. There are no mixed deletion terms: no
Hamiltonian path can use two outgoing arcs of b. The retained arc to0
is the strongness guard, not a cosmetic choice.

## C2 Direct inversion and its exact range

For h>=3, set `k=bitlength(h-2)` and put i in S when bit `k-1-i` of h-2
is1. The highest bit is1, so `0 in S`. Formula (C1) gives H=h and the
graph is strong. Its number of missing pairs is exactly

    bitlength(h-2) - popcount(h-2).

At the fixed order k+2 the strong subfamily contains exactly one graph
for each integer in `[2^(k-1)+2, 2^k+1]`. Different counts imply different
isomorphism classes. These consecutive intervals cover every h>=3.
They are the range of this family, not the entire fixed-order spectrum.

The graph is a tournament precisely when every bit is1, giving the
inherited values `h=2^k+1`. For the missing tournament targets the required
holes are explicit: h7 gives one missing pair on five vertices; h21 gives
two on seven vertices. The inherited Session235 census supplies a smaller
six-vertex, one-missing-pair H21 example, so this constructor is not optimal.

Complete the nonnegative spectrum with an edgeless pair for0, a singleton
for1, and the inherited [Session235 audited four-vertex rows](../../session235_arc_pair_models_20261109/audit/AUDIT.md)
`(2,4,9,2)` for2.
The positive examples are strong. The construction takes O(log h) vertices
and O((log h)^2) explicitly written adjacency bits.

## C3 Why this does not produce missing tournament values

The deletion operation leaves the tournament category. Filling a deleted
pair with its original arc restores a power-of-two contribution, while
filling it in the opposite direction introduces different paths. Neither
operation preserves the specified total. Similarly, writing each missing
pair as raw11 gives a no00 matrix with the same count only under the
explicit cancellation decoder `C(A)_ij=A_ij(1-A_ji)`. Its ordinary count
is a different observable.

The mechanism is related to the earlier adaptive-root work through named
endpoint response: b is the unique controlled departure port, and a is
the mandatory entry into that port. Restricting each path to one controlled
departure removes mixed terms. In the tournament induction, the extra
orientations needed to fill all holes destroy this linear binary rule.

## Controls and provenance

`check.py` tests all126 words with1<=k<=6, all62 words with1<=k<=5 also
by direct permutation counting, the full interval consequence, strongness,
hole counts and139 selected targets. Nonstrong words, wrong target types
and the raw/cancelled distinction are hostile controls. The independent
implementation in `../history/compact_audit.py` checks the construction
without importing its producer. Saved normal and optimized outputs must
match byte for byte. These finite tests validate implementation; the
preceding bijection proves the all-order quantifier.

The inherited ingredients are THM-5071 (4a) above and
[Session235 PB5-PB6](../../session235_arc_pair_models_20261109/theory/PROOFS.md).
The present contribution is the one-sink split, its extra-path correction,
the logarithmic strong construction and its exact binary hole budget.
