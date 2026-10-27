# Exact low-order detours and the first weighted endpoint

**2026-10-25. PROVED scoped / FINITE-EXACT controls.** The new structural
result counts the first branch with fewer than `k` backward edges at total
mass deficit `2k+1`. A literal pointed-triangle operation connects this
branch to the preceding owner-collision layer. A separate exhaustive audit
repairs an omitted endpoint filter in THM-4734. Unrestricted A000568 and
all-order black-core connectivity remain **OPEN**.

## Inheritance and the current boundary

The closest proved mechanism is the exact endpoint/gap deficit identity in
[the mass theorem](../charge_palettes_boundary_20261020g/zero_heavy.md).
The first unit-gap layer `delta=k-1` and
[the missing-gap layer `delta=k`](../zero_heavy_next_20261020h/README.md)
are already closed. The subsequent
[two-unit corridor](../zero_heavy_extension_20261021a/two_unit_cells.md)
extends to `2k-1`, and its
[first-owner correction](../zero_heavy_extension_20261021a/first_owner.md)
closes `2k` for `k>=4`; the
[independent audit](../zero_heavy_extension_20261021a/audit/AUDIT.md)
checks its exact port quotient. Those results are inherited, not new here.

The canonical hostile is the pointed `C3` with one real external matching
port: one matching/two-empty-gap record and two `P3`/one-unit-gap records
describe the same class. The corrected near miss is dividing by a guessed
automorphism or treating those records as different classes. The less-used
sidecar is the marked span-two edge whose count was previously given by an
indecomposable-block convolution. Here it has a simple bijective formula.

The concept board is **ordered matching / pointed C3 / full numerical gap
word / weighted endpoint / actual anti-reversal / continuation quotient**.
The first four now connect by an explicit operation; the detour operation
still needs its full directed paths. Neither port collapse nor matching
enumeration is a black `K` move. A scalar equality of counts does not supply
an unrestricted recurrence or global reachability.

**Concurrent recovery:** MISTAKE-646 had already recorded the order-six
exception, and the October 25 [ports/owners checkpoint](../ports_owners_f1_20261025/README.md)
repaired THM-4734's surviving wording while this packet was being audited.
The small-order result below is independent corroboration plus the explicit
order-five, order-seven endpoint-four, and parity controls. The incoming
fixed-rank owner language uses feedback vertex number; this packet's mass
lane uses weighted feedback arc cost. No parameter identification is made.

## 1. The first weighted-endpoint branch

Let `Q[w]` be a pair-free weighted tournament: a homogeneous pair means two
vertices seen identically by every third vertex. Assume integer weights
`1<=w_i<=k` and weighted feedback

    f(Q[w]) = min_order sum_(backward ij) w_i w_j <= k.

Put `M_k=2k^2+3k`, `S_k=2k^2+k`. In the general case require mass
`M_k-(2k+1)`; in the nontrivial strong case require `S_k-(2k+1)`.
Consider precisely those classes admitting a good order with fewer than
`k` distinct backward edges. This is an existence property of the weighted
class; the theorem below shows that every good order then has the same
number of backward edges.

**Theorem.** For every `k>=3`, each such class has exactly this record:

- a perfect matching on `2k-2` ordered endpoints;
- exactly one distinguished endpoint of weight two, all others weight one;
- one weight-`k` vertex in every endpoint gap, including the two exterior
  gaps in the general case and excluding them in the strong case.

Conversely every such record is pair free and has exact weighted feedback
`k`. For the construction omitting exterior gaps, it is strong exactly when
the matching crosses every proper initial cut. General records have an
exterior source and sink and are not strong. Distinct records give distinct weighted-isomorphism classes. Thus the
two branch counts are

    B_k   = 2(k-1)(2k-3)!!,
    B^s_k = 2(k-1)I_(k-1),

where `I_j` counts cut-indecomposable perfect matchings on `2j` ordered
points. These count this branch only, not the complete `delta=2k+1` layer.

### Proof of completeness and exact feedback

For a good order let `r` be its backward-edge count, `t` the number of their
endpoints, `W` the total endpoint weight, `g` the number of nonendpoint
vertices, `q=k-r`, and `U=sum_gaps(k-w_i)`. Pair freeness bounds `g<=t+1`.
The inherited identity is

    M_k-v = (2k+1)q + (k+r-W) + k(2r-t) + k(t+1-g) + U.       (1)

Every term is nonnegative. For strong cores replace `M_k` by `S_k` and
`t+1-g` by `t-1-g`: the first and last vertices of a good order must be
endpoints, or they would be a source or sink.

At deficit `2k+1`, `q>=1` forces `q=1` and all remaining summands zero.
Thus `r=k-1`, `t=2k-2`, the backward graph is a perfect matching,
`W=2k-1`, and all allowed gaps are present at weight `k`. Positivity and
integrality force exactly one endpoint of weight two. Its one incident
matching edge contributes two, and the other `k-2` contribute one each,
giving displayed weighted cost `k`.

The class has `4k-3` vertices, or `4k-5` in the strong case. It cannot have
another good order with `r=k`. Such an order has all endpoints unit, since
its `k` positive integer edge products sum to at most `k`. Write
`d=2k-t` and `h=t+1-g` (strong: `h=t-1-g`). Then (1) reduces to

    2k+1 = (k+1)d + kh + U.

For `k>=2` this allows only `d=0,h<=2` or `d=1,h<=1`. The corresponding
vertex counts are at least `4k-2`, or `4k-4` strongly, a contradiction.
Also `r<=k-2` is impossible directly in (1). Therefore every good order
is of the stated kind and has weighted cost exactly `k`. This also
prevents isomorphism collisions with any `r=k` family at this mass.

### Pair freeness, ownership, and strongness

The ordinary full-gap matching argument applies: endpoints separated by
an occupied gap are distinguished by it; two gaps are distinguished by an
endpoint between them; a gap and an endpoint are distinguished by an
intervening vertex or by that endpoint's matching partner. Endpoint
weights do not alter these unweighted adjacency checks.

For `k>=3` the weight-`k` vertices are exactly the gap vertices. They induce
a transitive tournament, so their individual order is intrinsic under any
weighted isomorphism. Incidences to that fixed order locate every endpoint
in its own singleton interval, including the extremes in the strong case.
Thus the matching and the distinguished weight-two endpoint are recovered.
There is no symmetry division or unrecorded clone contraction.

All consecutive arcs of the displayed order are forward. Its strong
components are consequently contiguous along that Hamilton path. A proper
cut is one-way exactly when the corresponding matching prefix has no
crossing edge. Thus the strong condition is exactly cut irreducibility.
The counts follow by choosing a matching and then one of its `2k-2`
endpoints.

## 2. A matching bijection and a literal pointed-triangle operation

Let `H_k` count perfect matchings on `2k` ordered points with one marked
edge `(i,i+2)`, and `H^s_k` its cut-indecomposable subfamily. Delete the two
ends of the marked edge and distinguish its surviving middle point. This
is a bijection to a matching on `2k-2` ordered points with one distinguished
point. The inverse inserts two points immediately before and after the
distinguished point and pairs the new points.

Cut irreducibility is preserved in both directions. A cut of the reduced
order lifts to a cut of the expanded order with the inserted triple kept
on the same side as its middle point. The inserted edge crosses no such
cut, so a crossing in the expanded matching comes from the reduced one.
Conversely the old cuts retain their crossing and each of the two new
cuts inside the triple is crossed by the inserted edge. Therefore, for
every `k>=2`,

    H_k = 2(k-1)(2k-3)!!,
    H^s_k = 2(k-1)I_(k-1).                                (2)

The second identity replaces the earlier convolution description with a
direct structural formula. For `k=2,...,7` it gives
`2,8,60,592,7060,97944`.

In the actual weighted-core model the marked span-two edge specifies the
previous layer's owner fibre: its two intermediate gaps are absent, and
the resulting three unit vertices form a `C3`. The middle endpoint is
distinguished by its external matching partner. Delete the two unmarked
vertices of this pointed `C3` and change the distinguished port's weight
from one to two. Retain every outside vertex, arc and gap weight.

This is precisely the matching bijection above, now on weighted objects.
It removes two vertices and decreases total mass by one. The internal
backward edge of cost one disappears; the port's external backward edge
changes cost from one to two. The displayed cost, and hence exact feedback,
remains `k`. The gap vertices immediately before and after the triangle
remain separate. Strongness is preserved by the cut argument (both general constructions
retain their exterior source and sink).

For `k>=3` this gives a bijection between that specified pointed-owner
subfamily at deficit `2k` and the new branch at deficit `2k+1`. Its inverse
expands the distinguished weight-two endpoint into the uniquely marked
triangle record. **This is port collapse, not clone contraction:** the
three vertices of a directed triangle are not a homogeneous chain.
No `K` move or black-component predicate is preserved by assertion here.

### Sharp hostile at `k=2`

At `k=2`, the endpoint of weight two is indistinguishable in type from the
weight-two gaps. Both records collapse to one weighted class. Strongly
this is `C3` with weights `(2,2,1)`; generally it is a weight-two source,
that weighted `C3`, and a weight-two sink in order. The two weight
sequences `(2,2,1)` and `(1,2,2)` agree under a triangle rotation.
Thus (2) still counts two pointed matching records, but the new weighted
branch has **one** class, both generally and strongly. Its source has two
distinct pointed-owner classes, so port collapse has a two-to-one fibre.
The condition `k>=3` in the class bijection is necessary.

## 3. Complete endpoint-three audit through order six

Use `P=(0,...,n-1)` and the convention `K_Q(T)=T^op xor E(Q)`, complementing
every nonpath pair. A word is blue exactly when its positional reversal
is an anti-automorphism. The computation independently enumerates every
labelled tournament for `1<=n<=6` and retains exactly the directed, closed,
blue words with endpoint scores `(n-4,3)`. No mixedness, strongness or
isomorphism filter is applied.

| Order | Full labelled universe | Retained presentations | Three-black detours | Target has no black path |
|---:|---:|---:|---:|---:|
|1|1|0|0|0|
|2|2|0|0|0|
|3|8|0|0|0|
|4|64|0|0|0|
|5|1,024|2|0|2|
|6|32,768|12|11|1|

For each of the eleven successes, the receipt supplies three actual words
and a failed anti-equation at every step. A separate directed-arc-set
implementation replays directedness, blackness and the literal final
labelled tournament. The remaining three cases are impossible at any
positive black-path length because their targets have no black Hamilton
word at all. A word is directed at both ends of its own `K` move, and
blue/black status is unchanged there: path reversal preserves the path-edge
set, so complementing all other pairs preserves every anti-equation.
Thus absence of a black word also excludes incoming black moves. This is
stronger than failure of the three-step search.

The two order-five sources and their targets are

    (2,4,9,19,7)   -> (30,28,24,16,0) = TT5,
    (2,12,9,17,7)  -> (30,20,24,18,0) = 1 -> C3 -> 1.

Their targets have respectively one and three Hamilton paths, all blue.
At order six the unique excluded source is the inherited hostile

    (6,4,8,51,39,7) -> (58,60,57,16,32,8) = C3 -> C3,

whose target has nine Hamilton paths, all blue. Thus THM-4734's unfiltered
sentence asserting a detour for every `n>=6` drops a necessary boundary
condition. The strongest inherited repair is: every `n>=7`; at `n=6`,
exactly the eleven presentations outside this pure-blue-target exception.
The source packet's stated strong-mixed filter is also a valid repair.
The order-ten computation is unaffected by this correction.

A separate actual order-seven source,

    rows=(42,124,57,112,33,64,21), I=J={2,4},

has the required blue closed path and endpoint scores `(3,3)`. Therefore
the assertion that overlap two `I=J` can occur only at even order is also
false. A reflected two-element pair is possible in either parity. Only
overlap **one** forces an odd-order fixed middle position.

### The analogous endpoint-four hypothesis also needs its filter

An explicit order-seven source

    (6,4,8,19,103,79,15) -> (122,124,121,112,32,64,16) = C3 -> 1 -> C3

has the standard directed closed blue word and endpoint scores `(2,4)`.
Its target has nine Hamilton paths, all blue, so no positive-length black
detour reaches it. This **REFUTES** the unfiltered `n>=7` claim in
HYP-9348(a); a version excluding pure-blue targets remains **OPEN**. This
is one explicit hostile, not an exhaustive order-seven endpoint-four census.
The self-contained `endpoint_four_seven_control.py` enumerates all 5,040
permutations of this target and checks the source and target conditions.

The common obstruction is the palindromic order-join `C3 -> TT_m -> C3`
for `m>=0`. Each Hamilton path uses the three components in order; the
middle chain has one path and each triangle has three, so all nine paths
are blue under positional reversal. Applying `K_P` to the standard path
gives a closed blue source with scores `(2,n-3)`, `n=m+6`. Thus the two
counterexamples are adjacent members of one exact family, and every
endpoint-`k` law starting unfiltered at `n=k+3`, `k>=3`, meets this hostile.
This is an obstruction family, not a statement about the non-pure-blue core.

## 4. Controls, receipts, and the remaining interface

The branch checker uses no producer imports. It canonicalizes all records
for `k=2,...,6`, generally and strongly, by weight/directed colour
refinement followed by all unresolved permutations. It does not use the
gap decoder being tested. All literal clones, mass, strongness and port
collapses are checked. Independent subset dynamic programming verifies
exact weighted feedback in 178 controls; reverse relabellings test
covariance. All marked span-two matchings through `k=7` are generated and
checked under both inverse maps and both cut predicates.

| k | General branch classes | Strong branch classes |
|---:|---:|---:|
|2|1|1|
|3|12|8|
|4|90|60|
|5|840|592|
|6|9,450|7,060|

The positive controls include all explicit matched constructions; the
hostiles are the `k=2` two-to-one class map, the three pure-blue detour
targets, and odd-order overlap two. Every check uses explicit exceptions
and survives optimized Python. Reproduction from the repository root:

    python3 -B 04-computation/shape_ports_20261025/low_order_detours.py
    python3 -B -O 04-computation/shape_ports_20261025/low_order_detours.py
    python3 -B 04-computation/shape_ports_20261025/weighted_endpoint_boundary.py
    python3 -B -O 04-computation/shape_ports_20261025/weighted_endpoint_boundary.py

Saved receipts:
[low order](tournaments_low_order_detours.json),
[weighted boundary](tournaments_weighted_endpoint_boundary.json).
Normal and optimized receipts are byte-identical, with SHA256 respectively
`518c453b03148d7eac9ccbf5e6b428b4e2f267f010ca2083520d724bb484e750` and
`eba9c1ccfd5a44a432309454c9696d23c9d1f13a197485d167f535678570f965`.

The new operation connects two precise boundary languages while keeping the
external port and the exact weights. It does not resolve the `r=k`
families at deficit `2k+1`, later unit-cell owner collisions, all-order
black-core connectivity, or a closed unrestricted scalar recurrence for
A000568. Those remain distinct live obligations.
