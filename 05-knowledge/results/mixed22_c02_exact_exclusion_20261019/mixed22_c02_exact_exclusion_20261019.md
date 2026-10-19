# Exact exclusion of seven-thin mixed-(2,2) roots with at most two common neighbors

## Theorem and scope

Let P be a finite poset in which every vertex is incomparable with at most
seven others. Suppose P has no balanced pair, where balanced means that the
probability of either order in a uniformly random linear extension lies in
[1/3,2/3]. Choose an incomparable pair x,y, oriented x before y in the strict
majority order, of maximum index span in that order.

There is no such root with two lower slaves, two upper slaves, and c common
incomparable neighbors for c=0,1,2.

This is a finite morphology exclusion, not a proof for all seven-thin posets.
The proof uses Peczarski's published six-thin theorem for short-span orders.
The positive generator certificates and their independent replay are included
in this packet. No historical novelty is asserted.

## 1. Majority order and core structure

Under the no-balanced-pair hypothesis every pair has a unique orientation of
probability greater than 2/3, with comparable pairs oriented as in P. These
orientations form a total order: a directed triangle would make the sum of
its three event probabilities exceed 2, while every linear order satisfies
at most two of those events. Every cyclic tournament contains a triangle.

Put Q=I(x) union I(y), including both roots since they are incomparable.
The mixed-(2,2) core has vertices

`x,y,b1,b2,z1,...,zc,a1,a2`.

Each bi is incomparable with x and below y; each aj is above x and
incomparable with y; each zi is incomparable with both roots. The maximum-span
property forces every common neighbor between x and y in majority order.
Indeed, a common neighbor outside that interval would make its span to the
opposite root greater than the selected root span.

No core comparability runs from an upper slave to a common or lower vertex,
or from a common vertex to a lower slave, because that would make the relevant
vertex comparable with a root to which it is declared incomparable. Within
each group, label vertices in a linear extension. Thus every possible induced
core is represented by a naturally labeled poset on the c+4 nonroot vertices,
ordered by groups lower/common/upper, together with bi<y and x<aj.

The counts of these catalogues for c=0,1,2 are 40,357,4824. The independent
verifier obtains them by enumerating all subsets of forward nonroot relations
and deduplicating their transitive closures. Search used a separate ideal-
extension catalogue. Equality with the stored core set is checked for each c.
Every core satisfies the all-vertex degree bound seven.

## 2. Exterior decomposition and boundary generators

Every vertex outside Q is comparable with both roots. It must be below both
or above both; either mixed direction would make the roots comparable. Write
these exterior parts as A and B. In particular A<B.

Fix a core extension w. All A vertices must be inserted before the first root,
and all B vertices after the second root. The number of ambient extensions
restricting to w is a product L(prefix)*R(suffix), where the prefix consists
of core vertices before the first root and the suffix of those after the
second root. Constraints between an exterior vertex and a core vertex outside
these boundary segments are either automatically satisfied or impossible in
the poset; the remaining left and right choices are independent.

The left counts are nondecreasing as one appends available core vertices to
the prefix: append the new core vertex after a completed left interleaving.
Set U(empty)=L(empty) and U(qv)=L(qv)-L(q). Then U is nonnegative. Define
nonnegative right increments V dually. Telescoping expresses the ambient
extension measure as a nonnegative combination of generators G(u,v), where
G(u,v) is the uniform count of core extensions beginning with ordered prefix
u and ending with ordered suffix v. This is an outer relaxation: arbitrary
nonnegative generator coefficients need not be realizable by exteriors.
A valid inequality on every generator is nevertheless valid on every actual
ambient extension measure.

If a core vertex has incomparability degree seven already within Q, it is
comparable with every exterior vertex. Necessarily A lies below it and B
above it; a reversed relation would make that core vertex comparable with
both roots. Hence every left increment at or after a prefix containing it
vanishes, and every right increment at or before a suffix containing it
vanishes. Such generators are omitted. The included saturation lemma proves
this for all later increments, not only for the increment at that vertex.

## 3. Exact analytic screens

For incomparable u,v define

`delta(u,v)=|D(v)\D(u)| + |U(u)\U(v)|`.

The included directed-defect injection proof shows that delta<=1 implies
`p(u<v)<=2/3`, and delta=0 implies `p(u<v)<=1/2`. Thus bounds <=1 in both
directions give a balanced pair; a single such bound forbids that direction
in the hypothetical strict majority order.

Use the core roots and saturated vertices as forced exterior-direction seeds.
Let F_A be their upward closure in Q and F_B their downward closure. Set
alpha(u)=0 on F_A and otherwise `7-degree_Q(u)`; set beta(v)=0 on F_B and
otherwise `7-degree_Q(v)`. The ambient bound is

`delta_P(u,v) <= delta_Q(u,v) + alpha(u) + beta(v)`.

Indeed all A are below each vertex of F_A and all B above each vertex of F_B;
only a missing lower exterior predecessor at u or upper exterior successor at
v can enlarge these directed differences. Each is bounded by the appropriate
remaining incomparability slack. This also proves alpha=0 on upper slaves and
beta=0 on lower slaves. The included side-specific note gives the full proof.

The verifier checks the resulting balanced-pair witnesses, majority-direction
exclusions, and any forced directed cycle. Comparabilities, the root direction,
and common-neighbor directions x<zi<y are also included in those cycle checks.

## 4. Published-theorem short-span screen

No exterior vertex can lie between the roots in the majority order, so their
projected core index span is their full span. If the selected maximum span is
at most three, every vertex is incomparable with at most three preceding and
three following vertices. P is then six-thin.

Marcin Peczarski, *The Gold Partition Conjecture for 6-Thin Posets*, Order 25
(2008), 91–103, DOI 10.1007/s11083-008-9081-9, proves Gold partition, and hence
balance, when each element is incomparable with at most six others. The
publisher's primary abstract explicitly states this scope:
https://link.springer.com/article/10.1007/s11083-008-9081-9
(checked 19 October 2026).

Consequently projected maximum-span orders of span at most three are excluded
using that published theorem. We neither reprove nor claim novelty for it.

## 5. Finite positive certificates and exhaustive replay

For every remaining projected majority order, the stored certificate consists
of positive integer weights w_i and forward incomparable comparisons E_i.
For every retained generator, independently count its total E and its comparison
numerators N_i. The verifier checks, using unbounded Python integers,

`3*sum_i(w_i*N_i) <= 2*E*sum_i(w_i)`.

It follows that the weighted sum of the corresponding actual probabilities
is at most `(2/3)*sum_i(w_i)`. They therefore cannot all exceed 2/3, contradicting
the hypothesized majority order. Equality in a generator bound is sufficient
because every majority probability is strictly greater than 2/3.

The standard-library verifier regenerates every core extension, every retained
generator, every admissible maximum-span projected order, and all certificate
numerators. It verifies that every order is analytically excluded, removed by
the published short-span theorem, or covered by a forward positive certificate.
The record order-assignment digest is recomputed exactly. NumPy and SciPy were
used only for candidate discovery; neither is needed for the proof replay.

The final counts are:

| c | cores | analytic balanced cores | directed-screen orders | generator-certified orders | short-span orders |
|---|------:|------------------------:|-----------------------:|---------------------------:|------------------:|
| 0 | 40 | 0 | 806 | 42 | 67 |
| 1 | 357 | 0 | 18,464 | 862 | 734 |
| 2 | 4,824 | 115 | 703,635 | 18,500 | 1,257 |
| total | 5,221 | 115 | 722,905 | 19,404 | 2,058 |

These are disjoint counts within the recorded screen pipeline, not counts of
all orders that could satisfy each separate screen. The first 29 c=0 records
were reused from the earlier exact search without the short-span screen; those
records already certify all their admissible orders directly.

There are 3,465 distinct stored positive certificates: 2,762 discovered by a
single-row screen, 581 by a two-row screen, and 122 by rationally exactifying
LP suggestions. Their support has size at most seven. Some integer weights
are large (maximum 4,460,564,226,458,160), so exact unbounded arithmetic matters.
The discovery labels do not confer validity: every stored inequality is
recomputed independently.

## 6. Method limitation preserved

The accompanying `mixed22_core_realization_and_span_dependency_20261019.md`
and standard-library verifier exhibit an actual 16-point, five-thin poset
whose six-point mixed-(2,2) core has all eight forward comparisons biased by
more than 2/3. It nevertheless has eight balanced exterior pairs. Thus core-only
boundary-cone failure is not automatically a conjecture counterexample, and
product-structured or genuinely realizable exterior weights need not repair
such a core screen. Its candidate root order has span three, so Section 4
resolves its relevance to the present theorem.

## Reproduction

From this directory, using Python 3 only:

```sh
python -O verify_mixed22_c02_exact_exclusion_20261019.py
python -O verify_mixed22_core_realization_20261019.py
```

The first command writes `independent_replay_receipt_20261019.json`; the second
writes `mixed22_core_realization_receipt_20261019.json`. Explicit checks remain
active under `-O`. The frozen certificate gzip is hash-checked before replay.
The manifest records every proof, data, and verifier dependency in this packet.
