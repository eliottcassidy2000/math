# Trees, square-sum ports and the distinction between two enumeration problems

Status: PROVED elementary bridge identities and CITED inherited trijections;
exact finite checks and same-agent review complete. No independent audit.

A000055 counts unlabelled unrooted trees; A000568 counts tournament
isomorphism classes. They are not equal. The useful correspondence is on
marked, oriented objects, before the automorphism quotient.

The audited [THM-4835](../../../01-canon/theorems/THM-4835-euler-characteristic-criterion-for-marked-excess-minus-one-carriers-in-tournaments.md)
already gives an Aut(R)-equivariant bijection for every tree R on n vertices:

    Tour(V(R)) <-> Or(R) x Tour(E(R)).                  (T1)

The free-pair map sends nonadjacent vertices to the two end edges of their
unique connecting path, retaining the orientation of the pair. Tree edges
carry Or(R). This is the existing end-edge trijection, not a new theorem.
It is equivariant for Aut(R), not automatically for all of S_n; forgetting
the tree marking has nonconstant fibres. Its general excess-minus-one
criterion requires fixed-subgraph Euler characteristic1 for every subgroup.

At order 4, A000055 has the path and star. The third3-edge graph is a
triangle plus an isolate, already outside the tree species. Its complement
is the star; P4 is self-complementary. THM-4836's pair-or-triple carrier law
already fails from order 6 at C3=>C3. No all-order trijection of unmarked
class sets follows from the small picture.

## T2. The exact smaller carrier for a single bridge attachment

For a nonempty finite simple graph G and a named port p, let e_G(p) count
ordered Hamilton paths starting at p. Let H_und(G) count all ordered
Hamilton vertex words. These undirected counts are separate from raw
directed tournament H. Singletons have e=H_und=1 by convention.

Take disjoint nonempty G,H with ports p,q and add exactly the bridge pq.
Then

    H_und(G --pq-- H)=2 e_G(p)e_H(q).                  (T2)

Every Hamilton path must cross the only bridge exactly once; its two
pieces are Hamilton paths in the components with the indicated endpoints.
The two global directions give factor2. Conversely every such pair glues.

Thus a single integer e_G(p), or its nonzero bit for existence only, is
an exact carrier for THIS attachment test. It can be much smaller than
G, but unbounded counts are not a fixed number of bits. An unmarked
H_und is insufficient: an endpoint and an interior port of P4 both belong
to the same graph, but have e=1 and 0. An attached leaf distinguishes them.

If G itself is a tree of size>=2, e_G(p) is 1 precisely when G is a path
and p is an endpoint, and otherwise 0. This collapses many A000055 types
for the bridge-existence query; it does not count those types. Enumeration
must keep fibre multiplicities and their automorphism actions.

## T3. Two ports force a joint response

Let e_G(p,q) count ordered Hamilton paths from p to q, for distinct ports.
Adding a private leaf at each port gives exactly 2e_G(p,q) ordered Hamilton
paths. Endpoint marginals do not determine this joint response.

Explicit labelled hostile on vertices 0,1,2,3:

    G: cycle 0-1-2-3-0;
    H: cycle 0-1-3-2-0.

Both have e(v)=2 for every v and H_und=8. But e_G(0,3)=1 and
e_H(0,3)=0. Leaves at ports 0 and 3 therefore give2 paths versus 0.
This is a literal source/port-preserving attachment separation, not just a
collision of graph invariants. It explains why independent boundary
marginals cannot generally replace joint component-pairing states.

## Square-sum application and current boundary

For Q_n on1..n, adjacency means distinct vertices sum to a square. The
known Q15 path is

    8,1,15,10,6,3,13,12,4,5,11,14,2,7,9.

Its endpoint responses occur only at 8 and 9. Vertex16 attaches as a leaf
at 9;17 then attaches as a leaf at 8. Vertex18 attaches as a leaf at 7,
creating three degree-one vertices 16,17,18 and obstructing Hamiltonicity.
This recovers a marked-port mechanism behind the small exceptional run;
it is inherited in THM-4505/4510 and is not a new square-sum existence proof.

The stronger existing all-n>=25 path construction and n>=32 cycle
construction remain inherited. Neither the bridge lemma nor the tournament
completion supplies their proof. The next useful extension is a prescribed
port or two-port response across the native Q_n->Q_(n+1) update, which may
add several edges. Formula(T2) applies only when the new attachment really
is a single bridge.
