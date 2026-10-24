# Actual blue wiring on the order-ten positive-error slice

24 October 2026. **FINITE-EXACT**, with an exact carrier argument: the induced
blue multigraph on the24 positive-error self-converse order-ten nodes is
recovered completely. This is not the full metagraph, and boundary target
classes are not identified or counted as distinct nodes.

## Frozen input and output

The input is the independently audited converse-error packet with manifest
SHA256 `30b62a8fdf93065b19f071050d9fcfe5ad08c5f79cfcfbf30afcdae53babde5a`.
It classifies all160 positive-error ordinary order-ten classes, all92 merged
nodes, and all24 self-converse members. Those24 are rigid, strong two-C5-cell
classes, each with error40; four have crossing degree1, twenty degree2. Every
one is mixed, not an all-blue point of symmetry. Their4200 blue presentations
are independently verified in the base packet.

**Exact result.** The induced blue multigraph on these24 nodes has:

- 90 blue edges: four loops and86 nonloop edges;
- 180 internal blue ports and4020 boundary blue edges;
- one connected component on the20 crossing-degree-two nodes, and four
  isolated crossing-degree-one nodes;
- 70 occupied unordered node pairs including loops, with multiplicity
  histogram `1:52, 2:16, 3:2`.

All4020 boundary blue edges end at zero-error self-converse nodes. This counts
boundary edges, not distinct boundary nodes. No statement about their mutual
identifications or further adjacency is made.

Node IDs are the actual representative indices in the frozen converse receipt,
not arbitrary new atlas labels. The loops are at nodes9,47,74,75, once each.
The two multiplicity-three pairs are (28,74) and(31,65). The complete adjacency
list, every internal edge's two presentation words and source paths, and
per-node boundary counts are in `normal.json`.

## 1. The exact presentation carrier

Fix the base directed Hamiltonian path0->1->...->9. The36 non-path pairs are
(i,j) with0<=i<j<=9 and j>=i+2. For an actual Hamiltonian path
P=(v0,...,v9) of a representative T, define its word w(P) by the bit
T(vi,vj) on each of these pairs. Together with the fixed base-path arcs this
word determines the normalized labelled tournament exactly.

If two paths of T yield the same word, mapping the entries of one path to
the corresponding entries of the other is an automorphism of T. Because all
24 source classes are rigid, the paths are equal. Words from two different
ordinary source classes cannot coincide either. Therefore the4200 blue paths
give4200 **distinct actual presentation words**, not a multiset or a chosen
coordinate representative from each fibre.

Let rho(i)=9-i. A word is blue exactly when its bit at(i,j) equals its bit at
(rho(j),rho(i)), by the proved grid-symmetry/anti-automorphism equivalence in
THM-644(a), inherited through the frozen base packet. All enumerated words are
checked literally against this relation, independently of their construction
from the source anti-involution.

## 2. The line operation and membership test

The native line involution kappa fixes the base-path arcs and reverses every
non-path arc. On the36-bit presentation carrier it is exactly

    kappa(w) = (2^36-1) XOR w.

It is fixed-point-free, squares to the identity, and commutes with reflection.
Thus it takes blue words to blue words. A blue line is the unordered pair
{w,kappa(w)}. The checker verifies the encoded complement against a separate
literal flip of the non-path arcs in the normalized tournament rows.

Let W be the4200-word set and owner:W->{the24 source nodes} its exact map.
Then the complete induced positive-error blue multigraph is obtained by:

1. For each w in W, calculate kappa(w).
2. If kappa(w) is in W, record its actual owner and count the line once by
   requiring w<kappa(w). Equal owners produce a loop.
3. Otherwise record one boundary edge at owner(w).

There is no target canonicalization and no tournament-class guess in this
operation. It is exact because W contains every blue presentation of every
positive-error self-converse node, not merely selected examples. Blue nodes
are self-converse; conversely every self-converse node is blue-incident by the
proved anti-Redei theorem. Therefore a blue partner outside W cannot have
positive refinement error: otherwise it would be one of the24 known source
classes and its actual presentation would occur in W. Its error is zero.

Ordinary versus converse-merged classes introduce no additional identification
here: every blue presentation belongs to a self-converse class, and reflection
fixes every blue word individually. The native complement-paired blue lines
are exactly the counted lines on the merged node set.

## 3. Conservation and finite certificate

The receipt verifies the port accounting

    4200 = 2*90 + 4020.

Each induced loop contributes two internal ports at its node. Each internal
nonloop edge contributes one at each endpoint. Every boundary edge has exactly
one positive-error endpoint, so it is counted once. The four degree-one nodes
contribute492 boundary edges and no internal port; the twenty degree-two nodes
contribute3528 boundary edges and all180 internal ports.

`check.py` recomputes every blue path from the actual unique anti-involution,
checks the previously certified per-node path counts, injectivity, literal
reflection, actual normalized rows, complement-versus-row-flip equality,
internal partner rows, the involution, and every per-node conservation law.
It stores both original source paths on every internal edge, giving direct
labelled witnesses rather than only an adjacency matrix. Normal and optimized
Python produce byte-identical outputs, and a real-copy relocation reproduces
the result.

The base packet remains immutable. This extension reconstructs one actual
induced blue wiring slice; it neither predicts black edges nor identifies the
zero-error boundary nodes, and it is not a full order-ten metagraph recovery.

## 4. Inherited endpoint transport explains the four isolated nodes

The exact actual-path transport identity is inherited, not newly claimed:
Session40 `marked-ports/PATH_TRANSPORT_AND_FOREST_RULE.md`, Section1, archived
here as `sources/PATH_TRANSPORT_AND_FOREST_RULE.md`. For a tournament T and an
actual directed Hamiltonian path P with distinct start s and end t (n>=2),
let K_P(T) retain its path arcs and reverse the other arcs. Then

    d_K(v) = n-1-d_T(v) + 1_(v=s) - 1_(v=t).

One direct derivation takes the complete converse, restores every path arc,
and observes that each interior vertex has one restored outgoing and one
restored incoming arc, whereas only one remains at each endpoint. The actual
predecessor/successor incidences, not independent endpoint histograms, are the
retained data in the stronger Session40 subset-count identity.

This immediately explains the isolation in this new slice. Every positive
self-converse order-ten target has degree multiset either3^5,6^5 or4^5,5^5.
For a degree-one source, all eight interior target degrees remain3 or6; hence
any positive target would have to be degree-one. But its start has target
degree10-d_T(s), which is7 or4, neither3 nor6. No such target exists. Therefore
the four degree-one positive nodes are isolated in the induced positive-error
blue graph, without using the lookup count.

For a degree-two source, the interior target degrees are4 or5, so any positive
target must also be degree-two. The start condition forces d_T(s)=5 and the
end condition forces d_T(t)=4. This is necessary, not sufficient: the actual
lookup determines which eligible ports remain positive-error. Every one of
the180 internal ports is explicitly checked against these conditions.

Two further scalar consequences of the same endpoint identity, for n>=2,
are

    sum_v d_K(v)^2 - sum_v d_T(v)^2 = 2(d_T(t)-d_T(s)+1),
    c3(K_P(T)) - c3(T) = d_T(s)-d_T(t)-1,

where c3 counts directed cyclic triples. The first follows by squaring the
complement degrees and changing only the two endpoints. For the second use
c3(T)=binom(n,3)-sum_v binom(d_T(v),2). The checker verifies the degree identity,
squared-score identity, and the triangle identity by actual triple enumeration
on all4200 source/partner pairs. In particular internal edges in this slice
preserve c3, but no sufficiency claim is made from that scalar condition.

For any regular odd-order source with n>=3 and degree q=(n-1)/2, the target
has endpoint degrees q+1 and q-1 and interior degree q. Thus it is not regular
and c3 decreases by one. The n=1 coincidence of the two endpoints is excluded
from these two-endpoint scalar formulas. These are inherited-transport
consequences and new scoped applications, not a general wiring recovery law.

The endpoint condition's exact insufficiency is visible without a new class
search. Of the3708 blue ports on the degree-two source nodes,1671 have source
endpoint degrees4,5 and are rejected by the score obstruction. The other2037
have endpoint degrees5,4 and preserve both the4^5,5^5 degree multiset and c3.
Only180 of these are internal ports;1857 still leave the positive-error slice.
Thus even these preserved scalar data do not replace the literal complement
membership test. The four degree-one nodes have162 blue paths from degree3
to6 and330 from6 to3; all492 leave the slice.

## 5. Exact topology of this induced multigraph only

The induced graph has24 vertices,90 edges (including its four loops), and
five connected components. Its graph-theoretic cycle rank is therefore

    90 - 24 + 5 = 71.

The checker removes each **actual edge** individually and recomputes connected
components, so parallel edges are handled correctly and loops are never
mistaken for bridges. It finds no bridge. It also removes each vertex and its
incident edges and recomputes the number of components, finding no articulation
vertex. These are properties of this induced slice, not assertions about
bridges or articulation in the full graph with its4020 outgoing blue edges.

With loops contributing two to degree, the induced degree histogram is

    0:4, 2:1, 4:1, 6:2, 7:1, 8:1, 9:5, 10:1, 11:4, 12:4.

These finite topology values are derived from the explicit certified adjacency
and have no asserted connection to unrelated occurrences of the same numbers.
