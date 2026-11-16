# The even-contact shadow: faithful cut totals, unfaithful continuation

A tournament itself is complete, so planarity of the original graph
would exclude every order at least five. This note instead puts
planarity on the literal support of the rooted Hamilton contribution.
No structural deduction is made merely from a picture of the quotient.

Delete a named root r to obtain D. Let Q be the root's outcontacts in
D, and L its incontacts. Let W_ij count directed even-order paths in
D from i to j, weighted by the Hamilton count of their omitted
subtournament (the empty one has count 1). The **even-contact shadow**
is the weighted bipartite support graph on Q and L with edge iff
W_ij>0, oriented for reading from Q to L. By the odd-cycle packing
identity at the root,
```
H(T)=g+2 sum_(i in Q,j in L) W_ij,  g=H(D).
```
A weighted quotient obtained by merging vertices within each side,
summing weights over each pair of classes, preserves this scalar
exactly. It need not preserve individual endpoints, contact counts
inside Q, E/O, or a subsequent rooted extension. The terminology
here is the projection defined above, not planarity of any other
metagraph or of the full path graph.

## Ferrers classification and small scalar capacity

When D is transitive, with vertices in forward order, W_ij>0 iff
i<j (the edge i→j itself is an even-order path, and a reverse path
is impossible). The outcontact neighborhoods are therefore nested:
this is a Ferrers or chain graph.

Such a shadow is planar iff no K3,3 occurs, equivalently iff at
most two rows have size at least 3. Necessity follows by choosing
three common neighbors of the third largest row, if it exists.
Conversely, order the columns so every row neighborhood is an
initial segment. The at most two large rows make K2,q on their
common columns, with possible additional degree-one columns for
the larger row. Any small row touches only the first two columns.
Place these columns consecutively in the plane K2,q embedding
(choose the corresponding pair of two-edge paths next to each
other). All degree-two small rows are internally disjoint paths
between that same pair of columns, and fit in their common face;
add degree-one leaves and isolated vertices freely. If fewer
than two large rows, the same embedding after adding a dummy
large row covers the graph as a subgraph.

A strong one-apex tournament with transitive deletion is encoded
by an n-bit word starting 1 and ending 0, with the ones the
outcontacts. Planarity allows at most two ones before the third
last zero. If there are at least three zeros, put that zero at
zero-based location t. Before it is the initial one and at most
one more one, t choices. After it is exactly the penultimate
zero, at any of n−t−2 positions, then the last zero at the end.
These give
```
sum_(t=1..n−3) t(n−t−2)=binom(n−1,3)
```
words. At most two zeros contributes n−1 more words. Thus the
exact total is n−1+binom(n−1,3). For n≥3 the first and last
contacts alone contribute H≥1+2^(n−2): every even subset of the
n−2 intervening vertices gives a different even path between
them, each with at least unit weight. Hence the number of
absolute counts at most X from strong planar-shadow one-apex
words is only O((log X)^4). These counts are not all equal
to their words, so this is a capacity upper bound.

More generally, the same Ferrers graph is induced on the singleton
SCCs of an arbitrary deletion, in their forward SCC order. They
have positive W only in that order, and every forward pair has
its direct path. The unrestricted binary words of length s with
planar Ferrers support number O((s+1)^5): strip leading zeros and
trailing ones (isolated in this shadow), then use the cubic
enumeration for the occupied segment and sum over lengths and
over placements. Zero-inversion words add at most s+1.

## A quantitative deletion-budget barrier without readiness assumptions

For any strong T, the first deletion SCC meets Q, and the last
meets L. If they are distinct and there are j≥1 singleton SCCs
strictly between them, choose any even subset for a directed path
from a first-component outcontact to a last-component incontact.
Then H≥1+2^j, directly by the shadow formula. There are at most
two additional singleton SCCs at the endpoints, so the total
singleton order is O(log(X+2)) for H≤X. The one-component case
already has at most one singleton.

Suppose g≤G and H≤X; replace G by min(G,X) by induced monotonicity.
Let n_i,g_i be the orders and Hamilton counts of nonsingleton
deletion SCCs. Their product is g, each g_i≥3, and Busch gives
n_i−1≤3 ln(g_i)/ln 5. Thus the number of such SCCs, and total
order of their vertices, are both O(ln(G+2)), and
```
sum_i binom(n_i,2) ≤ [9/(2 ln² 5)] (ln G)² + O(ln(G+2)).
```
All between-component arcs are prescribed. Enumerate internal arcs
over their binary slots; ordered component size lists and root
contacts on nonsingleton vertices add only exp(O(ln(G+2))). With
s=O(ln(X+2)) singletons and k=O(ln(G+2)) nonsingleton SCCs, their
interleaving has at most (s+1)^(k+1) possibilities. The planar
singleton contact words have O((s+1)^5) choices, independent of all
the enumerated internal arcs. Summing over s gives at most
```
exp(κ (ln G)^2 + O(ln(G+2) ln ln(X+3))),
κ = 9 ln 2/(2 ln² 5)
```
presentations, hence at most that many distinct counts. Therefore
if G≤exp(c sqrt(ln X)) with 0≤c<1/sqrt(κ)≈0.911, these strong
planar-shadow readers cover a density-zero subset of odd targets.
There is **no** bounded-defect assumption here. This does not
constrain arbitrary nonplanar-shadow tournaments, or even their
planar quotients, which need not have the same information as
the support.

## Same exact weighted planar quotient, opposite continuation behavior

For q≥3 take the cycle of blocks
```
r → C3 → TT_q → r
```
and its converse, equivalently r→TT_q→C3→r (reverse the directed
triangle labels). Both are strong and their shadow supports are
the nonplanar K3,q, with the out/in shores exchanged. Contract
each shore to a single vertex; both weighted planar quotients
are the same single edge of weight (9·2^q−6)/2. They have the
same H=9·2^q−3 and g=3, but their normalized profiles are
```
                      E                 O                    Δ=H−2E+g
triangle outcontacts  9·2^(q−1)         E+3                  0
chain outcontacts     6(2^q−1)          E+3                  12−3·2^q
```
A direct proof uses the rooted block transfer: append an
incontact block U to a root r with all nonroot vertices beating
U, and U beating r. Its readout is h_U I+q_U Z and its E'=q_U Z/2,
where h_U=H(U), q_U=P_U(2), I=h−2O+2g, Z=E+O−g on the old root.
If the old root beats an entire block B then h=g=h_B,
O=q_B/2 (alternating two-slot partition), Z=q_B−h_B and
I=3h_B−q_B. Now C3 has (h_B,q_B)=(3,12); TT_q has (1,2^q).
Substitution gives the displayed profiles in both orders.
At q=3 these are two readings of H=69 with E=36 or 42 and
defect 0 or −12. One further near-sink reads H+2E, obtaining
141 or 153, so the failure is observable by a literal next move.
This is a whole nonplanar ready family with an arbitrarily
large deleted transitive side, not an isolated counterexample.

The same computation classifies the fully colored block
r→B→TT_q→r, q≥1: the normalized defect is 4H(B)−P_B(2), so it
is ready iff B is TT2 or C3. Indeed every Hamilton path can be
cut in |B|+1 different positions to give distinct ordered
two-slot covers; concatenation recovers the original path.
Hence P_B(2)=4H(B) requires |B|≤3, and direct checks of the
orders 1,2,3 leave exactly TT2 and C3. Only the triangle
choice crosses the nonplanarity threshold for q≥3.

### A planar baseline too

The elementary one-apex words 1^p0^q and 1^q0^p equally give
H=(2^p−1)(2^q−1)+2, and g=1. Their weighted quotient is again
one edge but the first normalized E=(2^p−1)2^(q−1) and defect
4−2^p; swap p,q to get the other reading. Taking p=2 and q≥3
makes both FULL supports planar, yet only one reader is ready.
Taking p,q≥3 makes both supports nonplanar with the same
weighted planar quotient. The missing information is the
root's actual divided boundary profile, not a planarity flag.
