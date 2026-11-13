# Hamiltonian arc-usage shadow: port products, separators, and collision

Provenance: THM-050 already concerns counting specified consecutive positions; the basic reversal bijection and ordered-SCC product for H are inherited. Additional claims here concern the whole weighted shadow and its lost orientation.

For an n-vertex tournament, let H be the number of directed Hamiltonian paths and W(u,v) the number whose consecutive pair is the unordered pair {u,v}. Let s(v), e(v) be number starting and ending at v. For n>=2, the exact identities are
```
sum_(u<v) W(u,v) = (n-1)H,
sum_(u!=v) W(u,v) = 2H-s(v)-e(v).
```
Each path contributes n-1 edges; it contributes degree 2 except at its two endpoints. Converse leaves W unchanged pointwise. A self-converse witness gives an automorphism of W, but that automorphism need not recover the direction. The n=1 case has H=1 and no weighted edges.

## Full ordered-SCC product

Write T=C_1→...→C_k for its unique strong-component sequence, H_i=H(C_i). Every Hamiltonian path is precisely one concatenation of independently chosen paths through the components in order. Thus, with h=product H_i:

* In one block C_i: W_T(u,v)=W_i(u,v) product_(j≠i) H_j.
* Adjacent blocks C_i,C_(i+1), u in the earlier and v in the later: W_T(u,v)=e_i(u)s_(i+1)(v) product_(j≠i,i+1) H_j.
* Nonadjacent blocks: W_T(u,v)=0.

Each s_i, e_i is positive (singleton trivial; otherwise cut a directed Hamiltonian cycle at the desired vertex). Adjacent seam is therefore a complete bipartite positive support with weight matrix of rank exactly one. After division by h its weights are the outer product of the component endpoint probabilities. Intra-component normalized weights are unchanged under ordinal padding. In particular if component membership is supplied, W/h determines each interior component's end distribution on its outgoing seam by row sums and start distribution on its incoming seam by column sums. It cannot generally see the unused first-component start or last-component end separately from their local sums; these two can also be recovered if their component decomposition and local W are both supplied, by the degree identity. Without supplied blocks the rank-one structure alone does not identify which partitions are SCCs.

If the direction on each positive edge is retained, deleting all zero-use arcs preserves exactly the SCC partition: each nontrivial strong component has a directed Hamiltonian cycle whose every arc lies in some Hamiltonian path of that component, hence in a global one. Intercomponent arcs all point forward. In the undirected support graph, articulation vertices are exactly the singleton SCCs with an earlier and later SCC. Proof: removing an internal singleton breaks the two adjacent seams with no nonadjacent edges to replace them; removing a vertex from a nontrivial SCC preserves its connectedness (a spanning undirected cycle minus the vertex is connected), and its adjacent complete joins remain present. An endpoint singleton has no two sides to separate. If there are no such singleton separators, the support can still be missing an edge inside a strong block: on five vertices, assign the lexicographic pair bits 10 (1 means increasing direction; bit positions start at zero). This is strong with H=9; only {1,4} has zero weight. Example therefore refutes a clique characterization of strong components by positive support.

One can compute the full weights without enumerating paths: for S subset V, let F(S,u) count Hamiltonian paths of T[S] ending at u and B(S,v) those starting at v. The usual add-one DP gives F(S,u)=sum_(w in S-u, w→u) F(S-u,w), with singleton 1, and B dually. For the arc u→v, split each path at the arc:
```
W(u,v)=sum_(S contains u, omits v) F(S,u) B(V\S,v).
```
The factorization is exact because every path through the arc has a unique prefix vertex set. This also handles induced subgraphs or enlargements (reuse old DP cells). A zero means no partition with both positive factors, not merely that the arc does not belong to a directed Hamilton cycle.

## Explicit loss of converse-class information

Let G=C7⋊C3 with pairs (a,b) in residues (mod 7, mod 3) and multiplication (a,b)(c,d)=(a+2^b c,b+d). Let
```
S={(0,1),(1,0),(1,1),(1,2),(2,0),(2,1),(2,2),(3,0),(4,1),(4,2)}.
```
Define U by x→xs for s∈S. The set S takes exactly one element from each inverse pair of nonidentity elements, so U is a tournament. Left translations act vertex-transitively. It is strongly connected because (1,0),(0,1) are in S and their positive powers generate G. It is not self-converse: any anti-isomorphism can be made to fix the identity by a left translation. In the induced out-neighborhood at the identity, the unique vertex of internal outdegree 7 beats vertices of internal degree multiset [3,4,4,4,5,5,5]. In the converse, the corresponding unique vertex beats multiset [3,4,4,5,5,5,5]. These elementary ten-vertex fingerprints differ. See cayley_certificate.py for literal verification.

Every vertex of U has s=e=H(U)/21, by transitivity and the facts that each path has exactly one start and one end. Hence changing U to U^op does not change either boundary profile pointwise; it also does not change its internal weighted shadow. Now T=U→U and T'=U^op→U, on the same two labelled vertex blocks, have identical W by the SCC product. They are not isomorphic (compare first SCC), nor is T' isomorphic to T^op=U^op→U^op (compare second SCC). This is a certified collision at order 42; no minimality claim. More generally, k blocks of U with arbitrary independent converse choices yield 2^k distinct ordered isomorphism classes with identical W. For all k, reversal complements and reverses the bit word. For odd k it fixes none and the merged count is 2^(k-1); for even k it fixes 2^(k/2), so the count is (2^k+2^(k/2))/2. This gives exponentially many distinct merged classes in a single weighted-shadow fibre along multiples of 21.

### A smaller and stronger collision with a labelled endpoint probe

Append one terminal singleton to U and U^op instead of a second copy of U.
Their shadow is again identical pointwise, since the final seam uses the same
flat end vector. They are not isomorphic because their first SCCs are not,
and they are not converses up to isomorphism because converse puts the
singleton first. This gives a collision at **order 22**, with no minimality
claim. Both individual endpoint count vectors agree too: only the final
singleton can end a global path, and each of the 21 first-component
vertices has the same start count. Thus even the complete labelled W
together with the separate labelled start and end vectors can lose more
than global converse. With k 21-vertex blocks and one final singleton all
2^k orientations are distinct even up to global converse at order 21k+1,
with the same W and endpoint vectors. This is an exact lower bound on
information loss, not a search extrapolation.

## Finite controls and scope

replay_shadow.py exhaustively chooses one representative of each relabelling+converse orbit through order 7, then independently canonicalizes W by all vertex permutations. There are no extra collisions up to that order. Image counts for n=2..7 are 1,2,3,10,34,272. A zero-use strong arc first occurs at order 5. The order-22 and order-42 collisions are proofs by product and a local invariant, not a numerical enumeration of its Hamiltonian paths. No construction here claims a Chan–Pak proof.

## A sharper reconstruction test with endpoint polarization

Keeping endpoint labels separate (but still not the arc directions) gives a necessary signed divergence constraint. For any vertex v, outgoing used weight is H-e(v) and incoming used weight is H-s(v). Thus their difference is s(v)-e(v). If two orientations of the same labelled W realize the same endpoint imbalance, orient their changed positive edges as in the first tournament. At every vertex the changed incoming and outgoing *weights* balance, since reversing a changed edge alters divergence by twice its weight. Consequently the changed positive edges form a weighted circulation; each edge must lie in a directed cycle in that changed subgraph (flow decomposition), and expanding every positive integer weight into unit parallel copies decomposes the difference into directed cycles. In particular if W has no nonempty equal-sum pair of disjoint incident positive-edge subsets at any vertex, no edge direction can change at that vertex. If the vertices at which such a collision is possible induce a forest in positive support, the entire positive orientation is determined by W and endpoint imbalance. This is a sufficient reconstruction criterion on the used support only; it does not decide orientations of zero edges or claim that every balanced circulation switch is physically realized by the same shadow.
