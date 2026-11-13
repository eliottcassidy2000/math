# Same Hamilton count, different dark arcs, for any nontrivial reservoir

The two [dyadic switch](../session318_reservoir_descent_20261113/DYADIC_SWITCH.md)
charts are blow-ups of vertex t=3 in the following seven-vertex quotients
(outgoing masks, labels 0,...,6):
```
A: (54,64,122,115,98,66,1)
B: (54,32,122,115,98,64,3).
```
Replace t by any nonempty tournament U with its outside incidences uniform.
Label 2, called x, is the only vertex that can enter U; U beats every
other exterior vertex. Both quotients are strong, hence both blow-ups
are strong. The [reservoir note](RESERVOIR_POLYNOMIAL.md) gives the full
path-cover separation. Here we compute the whole support and the weight
of any arc by a small recursive interface.

## The marked layout recipe

Compress each maximal run of U in a Hamilton path to a token t. There
are at most two runs, because a return must use x. In the two-run
case, the first run starts the whole path (otherwise the first entry
would also have used x). A one-run layout is a spanning path in the
seven-vertex quotient. A two-run layout is a word with each of the six
exterior labels once and t twice, all consecutive transitions required
to be quotient arcs. Neither can have consecutive equal tokens. Expand
the token(s) by one Hamilton path of U, or by an ordered pair of nonempty
disjoint directed paths covering U. These choices and compression are
inverse. In both A and B there are 27 one-run and 9 two-run layouts.
Thus H=27H(U)+9F_2(U). In particular H=9(2^m+1) for U=TT_m.

To compute an exterior arc's weight, count layouts containing that arc
and weight them by H(U) or F_2(U) according to the run count. To compute
an arc between U and the exterior, impose the marked endpoint on the
appropriate run, and count only those path covers. An internal arc
has weight 27 times its count in one-run paths plus 9 times its count
in two-run covers (in either run). This is an exact recursive recipe,
rather than a numerical fit at several orders.

## Exact dark arc inventory

For every |U|>=2, A has precisely **one** dark arc, x→6, and B has
precisely **three**: 0→5, x→5, x→6. To see this, the finite
seven-label layout replay checks that these are the only exterior
arcs absent from both layout types. Every possible exterior contact
with U occurs in a two-run layout, for each chart. To mark any specific
contact vertex of U, assign that vertex alone to the relevant run and
take any Hamilton path through the nonempty remainder. For each
internal arc of U with |U|>=3, take it as a two-vertex run and take any
path of the remainder as the other run. At |U|=2 the internal arc
already occurs in one-run layouts. This proves the complete inventory,
including when U itself has dark arcs. The unweighted support of
the blow-up actually **restores every internal arc** of such a U.

For U=TT_m, m>=2, the entire internal weight of i→j, i<j, is identical
in both charts:
```
27 [j=i+1] + 18 ( 2^(m-j+i-1) - [j=i+1] ),
```
with indices 1,...,m. In a two-run cover put i,j in the same run
(two choices), put every strictly intermediate site in the other,
and distribute all outside sites freely. Only for j=i+1 must the
uniform assignment be excluded. The observer can therefore miss the
whole change if it only sees the growing module. It takes just the
*number of zero weights* in the exterior to distinguish the charts,
even without directions, and this distinction persists at every size.
