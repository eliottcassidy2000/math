# Equal path-cover data do not determine the one-birth spectrum

**FINITE-EXACT sharp order-five witness, with elementary PROVED mean law.
This refutes a universal transition kernel on H or unmarked path-cover
profiles. It does not exclude a recurrence for the one fixed sequence F_n.**

The inherited corrected [THM-4803, cyclic-triple realization fibres](../../../01-canon/theorems/THM-4803-intransitivity-is-the-relation-modulo-its-modular-tree-and-a000568-from-the-prime-cycle-index-census.md)
preserves the full ordered path-cover profile under changing the order of
strong components. The [Session146 response-kernel result](../residue_290_transport_20261102/INCOMING_RECONCILIATION.md)
concerns a different, named attachment observer: its order-ten kernel is
2520-dimensional and propagates upward. Neither a scalar H nor an unmarked
path-cover list is that observer. We use no obsolete all-order minimality claim.

## G1. A complete, minimal-order hostile

Let A be a directed triangle and use the dominance join, all arcs pointing
from the left block to the right. Set

    U = 1 -> 1 -> A,             V = 1 -> A -> 1.

These are ordinary five-vertex classes. Both have H=3, automorphism order3,
and ordered path-cover profile

    (P1,P2,P3,P4,P5) = (3,42,162,240,120).

An ordered k-cover is an ordered tuple of k nonempty vertex-disjoint directed
paths covering all old vertices. The two tournaments can be placed on the
same named triangle and two singleton vertices; their cyclic-triple
hypergraphs are then identical. Their dominance order is the lost coordinate.

Append a new named vertex w and independently choose each of its five
incident orientations. Each of the 32 cuts has weight1. The exact histograms
of H after this birth are:

| child H | number of cuts for U | number of cuts for V |
|---:|---:|---:|
|3|4|4|
|5|6|6|
|9|4|6|
|13|3|6|
|15|4|2|
|17|3|0|
|23|3|0|
|25|3|0|
|27|1|2|
|29|0|6|
|33|1|0|

Thus the attainable child supports already differ, as do the full integer
histograms. No statistic consisting only of (n,H,|Aut|,P1,...,Pn) determines
this one-birth response for every tournament. The census of every ordinary
class through order four shows that H-equal parents there have equal
histograms: the only nontrivial collision is the order-four converse pair
1->A and A->1, whose cuts correspond by reversing every arc. Hence five is
the sharp first order for the H-only histogram obstruction.

This is a statement about **all labelled cuts on one fixed source**, not a
uniform law on child isomorphism classes. Repeated child classes retain
their source multiplicities. Quotienting cuts by parent automorphisms or
forgetting the new root changes the measure and requires separate weights.

## G2. Why the mean survives and the distribution does not

For every nonempty tournament T, under independent fair orientations of a
new vertex,

    E[H(T+w)] = H(T) + P2(T)/4.                         (G2)

Every child Hamiltonian order either has w at one end, or has one nonempty
old path on each side of w. Each of the two endpoint choices for each old
Hamiltonian order succeeds with probability1/2. Each ordered two-cover
succeeds with probability1/4, since its two required arcs have distinct old
endpoints. These cases partition the child paths, proving (G2) by linearity.
For the empty parent the child has H=1; that boundary is excluded from (G2).

For U,V both means are 27/2. The unnormalized first moments are432,432;
the second moments are8064,8640, respectively. Thus their variances under
the same 32-cut measure are279/4 and351/4, a difference of18. The mean only
needs one path-cover count; a product of path indicators needs their joint
named endpoint incidence in the **same** cut. Replacing that shared cut
by independent fresh cuts discards the covariance.

This is the exact connection to [toothpick blocked ports](toothpick.md):
both quotients can preserve a total and lose the common-incidence predicate
that governs the next operation. The source objects are different, and no
toothpick-to-tournament bijection is claimed. The sidecar is actual endpoint
incidence, with source weights and vertex names retained.

## Reproduction and limits

Run `python3 -I -B joint_audit.py` in this directory; compare with `-O`.
The source imports no repository computation. It exhausts all labelled
tournaments at orders0..6, canonicalizes by every degree-compatible
permutation, checks orbit mass n!/Aut, and computes H with integer subset DP.
Literal permutation enumeration checks H through order five and all64
six-vertex children of the hostile. The JSON includes the exact arc lists,
bit convention, ordinary F_n coefficients, converse-fixed coefficients,
and positive/hostile controls. No strongness or symmetry filter is imposed.

The support-image problem is still to realize desired H using legal source
states and operations. An exact response formula, a mean with the desired
value, or an arithmetically admissible inverse state does not supply that
realization. This witness does not contradict the explicit scalar A000568
recurrence recovered in the preceding session, and does not rule out
nonlinear or history-dependent recurrences for the actual global sequence.
