# The high-weight boundary inside the zero-heavy stratum

**PROVED:** sharp total-weight bounds and exact upper layers.
This differs from heavy-coordinate deficiency and the incoming feedback-
vertex parameter. It does not count the entire unrestricted boundary.

Use the canonical pair-free weighted core Q[w] of
[the sparse theorem](../defect_boundary_20261020/sparse_theory.md).
Fix k>=1, weighted feedback f(Q[w])<=k, and 1<=w_i<=k. There are no heavy
coordinates at cap k+1. Write v=sum_i w_i.

## Sharp bounds

Choose an order of weighted backward cost c<=k. Let r be its number of
backward edges, t their endpoint count, W their total endpoint weight,
and g its nonendpoint vertex count. Pair freeness gives g<=t+1. Also
```
W<=sum_(ij backward)(w_i+w_j)<=sum_(ij backward)(w_i w_j+1)
 <=k+r.
```
The inequalities count each endpoint at least once and use
(w_i-1)(w_j-1)>=0. Since r<=k and t<=2r,
```
v<=k+r+k(2r+1)<=2k^2+3k=:M_k.                            (1)
```
In a nontrivial strong core, first and last ordered vertices must be
backward-edge endpoints, since otherwise they would be a source or sink.
Now g<=t-1, giving
```
v<=k+r+k(2r-1)<=2k^2+k=:S_k.                             (2)
```

Both bounds are sharp. Reverse a perfect matching on2k ordered unit
endpoints and place a weight-k vertex in every gap. The core is pair free:
intervening vertices or a matching partner distinguish each possible pair,
including the extreme endpoint/gap pairs. It has displayed cost k and
weight M_k. For strongness omit the two exterior gaps and use a matching
crossing every proper endpoint cut. All consecutive arcs remain forward;
every cut has a backward crossing, so the tournament is strong and has
weight S_k. At k=1 these are 1->C3->1 and C3.

Writing T_j=j(j+1)/2, these are S_k=T_(2k) and M_k=T_(2k+1)-1.
The triangular numbers here arise from the actual endpoint/gap budget,
not a numerical identification: 2k unit endpoints plus 2k-1 or 2k+1
gaps of weight k. In particular k=3 gives the bounds21 and27.

## Exact class counts near maximal weight

For k>=2 and 0<=delta<=k-1, the numbers of zero-heavy weighted core classes
of weights M_k-delta and, respectively, S_k-delta in the nontrivial strong
case are
```
Z_(k,delta)   =(2k-1)!! binom(2k+delta,delta),              (3)
Z^s_(k,delta)=I_k binom(2k+delta-2,delta),                  (4)
I_1=1,
I_k=(2k-1)!!-sum_(j=1)^(k-1) I_j(2k-2j-1)!!.             (5)
```
Here I_k counts perfect matchings crossing every proper initial cut;
the first values are1,2,10,74,706. Splitting at the first closed initial
segment proves (5).

**Completeness.** Set q=k-r and U=sum_(gaps)(k-w_i). The exact mass deficit is
```
M_k-v=(2k+1)q+(k+r-W)+k(2r-t)+k(t+1-g)+U.                (6)
```
Every term is nonnegative. Deficit delta<k forces r=k, t=2k, g=2k+1.
The defect graph is a matching, and all endpoint weights are1 because its
k positive edge products sum to at most k. The gap deficits sum to delta.
For delta<=k-2 every gap weight is at least k-delta>=2. For strong cores replace t+1-g
by t-1-g in (6), giving all2k-1 interior gaps and no exterior ones.

**Exact quotient for delta<=k-2.** This argument holds for every good order of the
class. Thus every weight-above-one vertex avoids all backward edges in
every such order. These vertices have one intrinsic transitive order and
are fixed individually by any weighted isomorphism. Each unit endpoint
lies alone in a prescribed gap between them, including the two extreme
positions in the strong case. Consequently its position is fixed too.
The matching and complete weight word are recovered from the class; no
generic symmetry division, deleted port, or hidden clone merger occurs.

Conversely every matching and weak composition of delta among the gaps
gives one distinct admissible core. Stars-and-bars gives (3). Strongness
is exactly crossing every proper cut: SCCs must be contiguous along the
displayed Hamiltonian path, and a one-way cut is precisely an unmatched
prefix cut. This proves (4).

**The unit-gap boundary delta=k-1.** The same deficit identity still
forces a matching and every gap occupied. If all gaps have weight>=2,
the preceding argument applies. Otherwise exactly one gap has weight1
and every other gap has weight k. The vertices of weight>1 still have
their intrinsic transitive order. Incidences to these vertices fix each
remaining unit interval. All but one such interval are singletons; the
exceptional interval has two vertices at an exterior gap, or three at
an interior gap (also at a strong end when the first/last interior gap
has weight1). Its location identifies the lowered gap in the record.

At an exterior gap the endpoint's matching partner is outside its unit
interval and fixed, so it distinguishes the endpoint from the gap. In
a three-vertex interval, either the two endpoints are matched to distinct
fixed vertices outside it, which distinguish all three roles, or they
are matched to each other. In the latter case the interval is a directed
C3, with identical external incidences. Its three rotations change the
choice of gap vertex but give the very same matching and weight record.
The outside-singleton positions are fixed before this local comparison;
no circular choice of partner is used. Thus distinct records still cannot
be isomorphic, though a C3 automorphism can occur. Every construction is
still pair-free. This extends (3)--(4) through delta=k-1 without dividing
by three or deleting the unit gap.

**A marked-gap bijection.** At delta=1 and k>=2, (3) becomes
Z_(k,1)=(2k-1)!!(2k+1)=(2k+1)!!. Insert a new endpoint in the marked gap
of a matching on2k ordered endpoints, append one endpoint at the end,
and match these two new endpoints. Deleting the last endpoint and its
partner is the inverse. At k=3 this is105=15*7. The map preserves the
ordered matching plus its gap mark. Reinterpreting it as a new tournament
changes the feedback budget and does not preserve K moves or graph minors.
The strong subfamily has its additional proper-cut condition.

Extremality fixes this much boundary information. Arbitrary zero-heavy
vertices do not inherit that conclusion. At delta=k missing gaps can
enter, but shared endpoints require at least delta=k+1, as proved in the
[next-layer classification](../zero_heavy_next_20261020h/README.md).
This corrects the earlier sentence placing both phenomena at delta=k.
The complete
k2 compiler gives42 classes at weight12 (delta2), whereas extrapolating
(3) would give45. The formula is not unrestricted.

## Exact controls and the growing boundary

The complete inherited compiler at k1,2 verifies all zero-heavy mass
coefficients and both maxima. Separate matching generation checks the
stated layers through k5, delta2, testing all actual homogeneous pairs
and all strongness predicates. Sixty subset-DP tests independently check
weighted optimum k on the declared smaller examples. Larger construction
counts supplement the proof and do not independently prove completeness.

The separate unit_gap_boundary.py canonicalizes all weight-preserving
isomorphisms using the intrinsic large-vertex order and every permutation
inside the remaining unit cells. At delta=k-1 it finds15/6,420/150,
and17325/6216 all/strong classes for k=2,3,4. It checks all clone pairs
and strongness, and independently verifies the marked-gap bijection for
k=2,3,4. These are quotient checks, not just presentation counts.

One7-cycle gives a useful budget hostile: all eight invariant tournaments
are regular and pair free, with the same prescribed cycle type. Six have
minimum feedback6 and two have feedback7. The receipt retains exact rows.
Thus unrestricted fixed-permutation counts cannot be used as a feedback
filter. The [signed clone transform and boundary operator system](../defect_boundary_20261020/boundary.md)
already supply stronger unrestricted cycle-index identities than a fresh
derivation here would. They still need a genuine budget observable.

The incoming one-feedback-vertex formula closes a different infinite slice;
neither result contains the other. An unrestricted A000568 recurrence must
handle the growing remaining interface, not extrapolate these top layers.

