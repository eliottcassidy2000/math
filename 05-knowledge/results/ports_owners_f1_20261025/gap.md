# A second exact band of zero-heavy weighted-core classes

2026-10-25. **PROVED, independently reviewed:** the class formulas below
for every budget k>=3 and 0<=delta<=2k-1. The extension beyond delta=k is
new in this packet; the earlier layers are inherited. **FINITE-EXACT:**
the same formulas for k=2 through delta=3, and bounded independent owner controls.
No unrestricted A000568 recurrence is claimed.

## Inheritance and the next actual boundary

The first unit gap was already closed in
[the high-weight theorem](../charge_palettes_boundary_20261020g/zero_heavy.md).
The first missing gap, at delta=k, was closed in
[the next-layer packet](../zero_heavy_next_20261020h/README.md), also promoted
as [THM-4736, delta-equals-k zero-heavy layer](../../../01-canon/theorems/THM-4736-the-delta-equals-k-zero-heavy-layer-exact-counts-missing-gaps-and-the-adjacent-matched-pair-collapse.md).
Its correction matters: a shared backward-edge endpoint costs k+1, not k.
The canonical hostile is the three cloned records behind 45 versus 42
at budget two. The less-used sidecar is the first-closed-cut decomposition
of ordered matchings, which can be differentiated when precisely one
component changes type. The complete port-word owner from
[the fixed-anchor theorem](../zero_heavy_ports_20261020g/zero_heavy.md)
remains the reference when small weighted records cease to be intrinsically
decodable.

The board is **mass deficit / endpoint overlap / missing gaps / unit cells /
ordered-component decomposition / canonical ownership**. The source is a
pair-free weighted tournament Q[w], 1<=w_i<=k, with weighted feedback <=k.
The target record keeps an ordered defect graph and every gap weight.
It preserves the entire weighted tournament, literal clone predicates,
feedback budget and strongness. Labels alone are quotiented. External defect
partners and the intrinsic order of all weight-above-one vertices are the
essential sidecars. No scalar mass or degree distribution replaces them.

## 1. Statements

Put m_k=(2k-1)!!, M_k=2k^2+3k and S_k=2k^2+k. As before, I_k counts
perfect matchings on 2k ordered endpoints crossing every proper initial
cut. Thus I_1=1 and

    I_k=m_k-sum_(j=1)^(k-1) I_j m_(k-j).

Define

    J_k=(k-1)I_k+sum_(j=1)^(k-1) I_j I_(k-j).                 (1)

This counts ordered graphs consisting of one three-vertex path and k-2
disjoint edges, crossing every proper cut. The edges in this graph are
backward tournament arcs, not Hamiltonian-path edges.

Use binom(N,r)=0 if r<0 or N<r. For k>=3 and 0<=delta<=2k-1, the number
Z_(k,delta) of all such weighted-core classes of mass M_k-delta is

    m_k [ binom(2k+delta,delta)
          -(2k+1) binom(k+delta-1,delta-k-1)
          -binom(k+delta-1,delta-k)
          +(k-1) binom(k+delta-2,delta-k-1) ].              (2)

The number Zs_(k,delta) of nontrivial strong classes at mass S_k-delta is

    I_k [ binom(2k+delta-2,delta)
          -(2k-1) binom(k+delta-3,delta-k-1) ]
      -(2k-3)I_(k-1) binom(k+delta-3,delta-k)
      +J_k binom(k+delta-4,delta-k-1).                     (3)

Every admitted core has exact weighted feedback k. At the first genuinely
shared-endpoint layer delta=k+1, formula (2) simplifies to

    Z_(k,k+1)=m_k [binom(3k+1,k+1)-3k-2].                 (4)

This is proved by the band theorem for k>=3; k=2 is a separately exhaustive
finite case. It gives 81, 2985 and 133665 at k=2,3,4. The corresponding
strong counts are 15,650,33616.

## 2. Completeness from the mass identity

For any order of weighted backward cost <=k, let r be its number of
backward arcs, t their distinct endpoints, W the endpoint weight sum,
g the number of gap vertices, and U their deficit from weight k. The
inherited identity, with q=k-r, is

    delta=(2k+1)q+(k+r-W)+k(2r-t)+k(t+1-g)+U.            (5)

All terms are nonnegative; use t-1-g instead of t+1-g for strong cores.
Since delta<=2k-1<2k+1, q=0. There are k positive backward products
whose sum is at most k, so all products are one. Every endpoint has
weight one, W=t, and the order has exact cost k. Write d=2k-t and h for
the number of unoccupied allowed endpoint gaps. Then exactly

    delta=(k+1)d+k h+U.                                  (6)

Consequently d<=1. If d=0 the defect graph is a perfect matching;
there is at most one missing gap. If d=1, all gaps are occupied and the
degree sequence has one vertex of degree two and all others of degree
one. The graph is precisely a three-vertex path plus k-2 separate edges.
These exhaust the records; no star of larger degree, second shared
endpoint, triangle, or omitted second gap is silently discarded.

The conclusion holds for every good order of the class. In particular
every vertex of weight >1 avoids backward arcs in every such order.
Its induced tournament is transitive, with one intrinsic order. The
smaller-cost possibility would require r<k and contradict (5), proving
that the weighted feedback minimum is exactly k.

## 3. Pair freeness and exact ownership

With all endpoint gaps occupied, a defect graph with no isolated endpoint
gives a pair-free tournament. Endpoint pairs are distinguished by intervening
gap vertices. Gap pairs have an endpoint between them. An adjacent gap and
endpoint are distinguished by a defect neighbour of the endpoint. For a
nonadjacent gap and endpoint, an intervening gap distinguishes them unless
there is only the one intervening endpoint, in which case its defect
incidence distinguishes them. This proves the claim directly, including
the path-of-length-two defect component.

For a matching with one gap missing, the inherited exact criterion remains:
clones occur iff the missing gap lies between endpoints paired to each
other. Other positive gap weights do not affect this literal adjacency
test. Such a record is rejected, not contracted and counted again.

To prove distinctness of all surviving records, fix the intrinsic ordered
vertices of weight >1. Incidences to them place the remaining unit vertices
in intrinsic cells. Consider the three possibilities separately.

* In a missing-gap matching record with delta<=2k-2, U=delta-k<=k-2,
  so every occupied gap has weight >=2. The only nonsingleton unit cell
  consists of the two endpoints surrounding the omission. Its directed
  arc fixes its order. They are not paired to each other, since that
  was rejected. The additional delta=2k-1 case is proved below.
* In a shared-endpoint record, U=delta-k-1<=k-2, so every gap has weight
  >=2. All endpoint cells are singletons. The complete defect graph,
  including its degree-two vertex and both partners, is fixed.
* In a full-gap matching record, at most two gaps have weight one,
  since three would cost 3(k-1)>2k-1 for k>=3. With at most one unit gap, use
  the previously audited local three-vertex argument. With two unit
  gaps the same two-gap argument as the old k=2 boundary applies at
  arbitrary k: separated interior cells are triples, while a unit gap
  at an exterior boundary gives a two-vertex cell with its unique directed
  order. Each interior triple is transitive
  unless its two endpoints are paired internally, when it is a triangle
  module. A transitive triple has a unique order. Rotation of a triangle
  changes a chosen vertex role but leaves exactly the same numerical
  matching and gap record.

For completeness, two adjacent unit gaps at an exterior end give a
four-unit cell. It is transitive, or has a triangle module and one extreme
vertex; both determine the record. Two adjacent interior unit gaps give
a five-unit cell with ordered roles a,g,b,h,c. If no two endpoints are
internally paired, it is transitive. If a,b or b,c are paired, it is a
triangle module joined in order to two singletons. If a,c are paired,
it is the transitive five-tournament with its longest arc reversed.
That tournament is rigid: its scores are (3,3,2,1,1), and the directed
arc inside each repeated-score pair fixes the pair pointwise. Its displayed
roles therefore agree between two copies of the same local shape. The
middle endpoint also retains its external partner. The three local types
have different strong-component decompositions, so cannot be interchanged.
These observations include strong end cells. No automorphism division is
inserted for a triangle rotation.

### The first missing-gap / unit-gap contact

At delta=2k-1 a missing-gap record can have one unit gap: then all its
other gaps have weight k. When the omission and unit gap are separated,
their unit cells are decoded as above. A transitive cell has a unique
order; a triangle cell has its two endpoints paired internally, so its
rotation preserves the numerical record. Exterior two-vertex cells
have their directed order. No partner in a mutable triangle is needed
to identify another cell, because its matching edge stays inside it.

If the missing and unit gaps touch an exterior boundary, the only larger
cell has three vertices. A triangle fixes the one allowed pattern: the
omission is exterior and the unit gap is between the matched endpoints.
Otherwise it is transitive. Its endpoint roles are then identified by
their distinct external matching partners, while the gap has no such
external defect incidence. Those partners are fixed singleton cells.

At an interior contact, the only new cell has four unit vertices. Its
two candidate role words are (a,b,g,c) and (a,g,b,c), with a,b,c endpoint
roles and g the unit gap. All external endpoints lie in fixed singleton
cells. Each local endpoint not paired internally therefore has a retained,
fixed external partner. The possible internal matching edges are exhaustive:

1. No internal edge: the cell is transitive. Its order and the external
   partners distinguish all three endpoint roles from the gap.
2. A short internal pair across the occupied unit gap: the cell is
   `1 -> C3` or `C3 -> 1`. The extreme vertex fixes which side contains
   the triangle; its rotations preserve the same numerical record.
   A pair across the omitted gap was already rejected by the clone test.
3. The long pair a,c: both candidate role words have the same strong
   four-vertex tournament, obtained by reversing the longest arc of a
   transitive four-tournament. It is rigid: scores (2,2,1,1), with each
   repeated-score pair fixed by its internal directed arc. Rigidity
   fixes the four positions, but the internal table alone cannot tell
   which middle position is the endpoint b. The external matching
   partner of b distinguishes it from g. This is the required port
   guard; omitting that partner would incorrectly merge the two records.

   Explicitly the common internal bit rows are (6,12,8,1). Either
   position1 or position2 can be the endpoint role; its matching partner
   outside the cell supplies the differing reversed arc. Retaining only
   this four-vertex induced tournament erases that role distinction.

The three cases have different strong-component decompositions, so their
records cannot cross. These cases also cover strong end contacts. Hence
the first missing-gap / unit-gap contact still has a unique numerical owner.

Finally, the full-matching, missing-gap and shared-endpoint constructions
have distinct total numbers of vertices (4k+1,4k,4k-1 in the unrestricted
case; 4k-1,4k-2,4k-3 in the strong case). Thus their records cannot collide
across types. This proves the weighted-isomorphism quotient exactly.

## 4. Counting the two defect types

For a matching let g=2k+1, or 2k-1 in the strong case. Gap deficits are
integers in [0,k] summing to delta; a deficit k records the one omitted
gap. Since delta<2k, exclusion of an overlarge deficit involves only one
coordinate. Their number is

    binom(g+delta-1,delta)-g binom(g+delta-k-2,delta-k-1).

For an omitted interior gap whose endpoints are matched, distributing
the remaining deficit delta-k over the other g-1 gaps gives
binom(g+delta-k-2,delta-k). In the unrestricted case there are m_k
choices of a matching with a marked adjacent pair. In the strong case
there are (2k-3)I_(k-1), by deleting that marked adjacent pair and
retaining the insertion position. These are precisely the rejection
terms in (2)--(3).

There are (k-1)m_k unrestricted ordered path-plus-matching graphs:
choose the degree-two vertex, its two neighbours, and a matching on
the remaining vertices. Its 2k unrestricted gaps (2k-2 in the strong
case) have total deficit delta-k-1. The cap is inactive in this band.
This gives the final terms of (2)--(3).

To count the strong defect graphs, use ordinary formal series
M(z)=sum m_k z^k, I(z)=sum_(k>=1) I_k z^k and
D(z)=sum_(k>=2) (k-1)m_k z^k. Every matching is a sequence of its
indecomposable initial-cut components, hence M=1/(1-I). A defect graph
with one degree-two vertex has exactly one component of the new type;
the matching components before and after it are independent. Therefore

    D=M^2 J,
    D=z M'-M+1,
    J=z I'-I+I^2.

Taking coefficients proves (1). Strongness of the constructed tournament
is exactly crossing every proper endpoint cut: every consecutive arc of
each admitted displayed word is forward, so its strong components occur
contiguously; a noncrossed cut is exactly a one-way cut. The decomposition
uses actual endpoint positions, not merely component sizes.

## 5. Hostiles, tests and the next boundary

The factor (k-1) alone cannot turn I_k into J_k: at k=2 it gives two,
but all three positions of the degree-two vertex give indecomposable
graphs. The convolution I^2 supplies the missing positional information.
The earlier 45/42 missing-gap hostile survives unchanged as a clone
rejection mechanism, now with a retained deficit distribution elsewhere.

The [checker](gap_check.py) generates literal tournaments, tests every
homogeneous pair, and canonicalizes under every permutation inside intrinsic
unit cells. Its second quotient at k=2 uses only full degree and weight
classes, and exact subset DP independently checks feedback. Reversed and
cyclic label permutations test the old mutable-sort-key failure. No
isomorphism or strongness filter excludes an uncounted source. Exhaustive
record universes and output counts are in [the receipt](gap_check.out).
The older complete compiler supplies separately inherited k=2 coefficient
controls; that receipt is not represented as newly independently rerun.
The separate [owner audit](gap_owner_audit.py) imports the disclosed record
generator but never the unit-cell canonicalizer. It uses full weight and
outdegree cells, exhausts every permutation therein, and checks clone pairs
by literal third-vertex comparisons. It checks every structural record at
k=3, delta=4 and5, and every missing-gap/unit-gap record at k=4, delta=7.
The new full-layer totals at k=3, delta=5 are6060/1046 (all/strong).
The k=4 contact subfamily has6720/2808 classes. See its
[receipt](gap_owner_audit.out) for the raw and rejected counts.

At delta=2k-1 the external-partner lemma above repairs the first contact.
At delta=2k a second missing gap can enter; at 2k+1 fewer
backward edges can enter. These are distinct transitions. No equality or
failure beyond the proved band is inferred from the formulas.

Reproduce from the repository root:

    python 05-knowledge/results/ports_owners_f1_20261025/gap_check.py
    python -O 05-knowledge/results/ports_owners_f1_20261025/gap_check.py
    python 05-knowledge/results/ports_owners_f1_20261025/gap_owner_audit.py
    python -O 05-knowledge/results/ports_owners_f1_20261025/gap_owner_audit.py

The finite k=2, delta=3 owner enumeration covers its exceptional possible
three-unit-gap cells. It establishes that isolated finite extension, not a
uniform proof for arbitrarily many adjacent unit gaps.

The root independently audited the mass classification, the all-order
unit-cell cases including the strong four-vertex external-port guard,
the separation by total vertex count, the capped-composition and clone
rejection terms, and J=zI'-I+I^2. That proof audit accepted the stated band.
Both scripts passed normally and under `-O`, with byte-identical receipts:
7,739,819 active checks in gap_check.py and46,725 in gap_owner_audit.py.
