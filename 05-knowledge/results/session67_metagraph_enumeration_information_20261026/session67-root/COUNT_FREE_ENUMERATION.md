# Count-free local enumeration of tournament classes and the metagraph

Status: PROVED as a corollary of the inherited all-order relocation-connectivity theorem; FINITE-EXACT implementation through order eight. No novelty claim for connectivity, no polynomial complexity claim, and no black-core connectivity conclusion.

## Exact source and target

For n>=3 let X_n be the 2^((n-1)(n-2)/2) tournaments on positions 0,...,n-1 with every consecutive arc pointing forward. A free bit x_(i,j), j>=i+2, is one iff i->j. These are marked Hamiltonian-path presentations modulo simultaneous relabelling, not tournament classes. The inherited theorem in `inputs/THEOREM.txt` proves that the components of legal single-vertex relocation on X_n are exactly ordinary tournament isomorphism classes. Its proof is all-order; the new finite computation does not replace it.

A rightward relocation takes position a to immediately after b, where b>=a+2. It is legal exactly when:

- x_(a,b)=0;
- if a>0, x_(a-1,a+1)=1;
- if b<n-1, x_(a,b+1)=1.

The new word is (0,...,a-1,a+1,...,b,a,b+1,...,n-1). Reindex the entire binary relation in this word. Both endpoints of each resulting edge belong to the same ordinary component. Leftward moves are inverses, so enumerating all rightward moves and taking undirected unions includes every edge needed by the theorem. Adjacent swaps are never legal. Deletion alone is not presumed to preserve a Hamiltonian path.

## Count-free recognition theorem

Suppose a candidate partition of X_n has nonempty blocks. It is exactly the ordinary class partition iff:

1. Every block is connected by legal relocation edges within the block.
2. Every legal relocation edge has both endpoints in the same block.

For the converse-merged partition, include reflection sigma in both conditions, where sigma reverses positions and takes the converse, so sigma merely permutes free bits.

Proof. Connectivity puts each proposed block inside one actual move component. Edge closure puts that entire actual component inside the block. Thus they coincide. Conversely actual components satisfy both conditions. The inherited connectivity theorem identifies these components with the intended isomorphism classes. This argument requires all cube points, not a selected family subset. Checking equality of the block labels at both endpoints of every rightward edge is sufficient; it is an undirected edge-closure condition, not a one-sided reachability claim.

The previous local-forest criterion in `inputs/LOCAL_FOREST_ADDENDUM.txt` used a supplied correct number of nodes K_n to prevent splitting true fibres. Here complete edge closure replaces that number. It costs checking the move boundary rather than knowing A000568 in advance. A sparse spanning forest alone permits arbitrary splitting and is insufficient; closure alone permits merging distinct true classes and is insufficient. No circular count input remains.

## Full metagraph reconstruction

After relocation unions, reflection unions give converse-merged nodes. Pair each free-bit word x with its complement tau(x). Each pair is one literal metagraph line, including loops and repeated lines. Its endpoints are the merged components; the line is blue iff sigma(x)=x. Complement and reflection commute, so both half-edges agree on colour. Count each complement pair once. A loop contributes two to node degree. Consequently each merged node's degree, with multiplicities and loops, equals its number of cube presentations.

Ordinary node size equals tc(T)=H(T)/|Aut(T)|, not H(T). If A is the number of ordinary nodes and M the merged number, the number S of self-converse classes is S=2M-A. These identities are consequences/checks, not inputs to generation. Do not confuse auxiliary relocation edges with literal complement lines or with the black-core subgraph.

## Algorithm and information ledger

The enumerator visits each x, evaluates the exact three-bit guards, computes the full relabelled relation using precomputed bit lookup maps, and unions its legal neighbours. It never calls graph isomorphism, Burnside, a precomputed class atlas, an automorphism routine, or a supplied A000568 value. The whole binary relation is retained. Guard bits alone would only recognize moves, not determine their destinations.

Time is O(n^4 2^((n-1)(n-2)/2)) up to union-find overhead with a direct quadratic-bit word encoding, improvable by bit lookup tables as implemented; space is exponential in the cube dimension. The source executable is deliberately bounded to 1<=n<=8 to prevent accidental unbounded memory allocation. The mathematical algorithm works at each finite order, but no practical enumeration at arbitrarily large n is claimed. For n=0 the empty class gives a_0=1 separately. For n=1,2 the implementation records one class and no lines; the fixed-point-free complement-line convention applies only for n>=3.

The generated node tables enumerate ordinary classes by a smallest fixed-path-cube representative, give exact tc and converse-merged representative, and flag self-converseness. The default saves all line ledgers through n=7. Supplying a fourth argument 8 exports the complete order-eight aggregated coloured endpoint ledger as well. Multiplicity retains every literal line even when its endpoints coincide with another line.

## Inheritance and boundaries

Source pin: b4234f728c05 (full pin in provenance). Closest proved mechanism: the October 20 relocation-connectivity packet. Corrected near miss: a connected forest plus an incorrect node count need not be complete. Canonical controls: C3 has three labelled Hamilton words but one normalized ordinary tiling; transitive tournaments can have nonzero real additive curl; ordinary and converse-merged class counts differ at order four. Least-used sidecar restored here: all legal boundary edges, rather than a guessed global node count.

Live board: complete binary relation; lawful path re-marking; canonical weighted cores; class/merged distinction; information lost at quotients. The count-free local certificate and the weighted-core recurrence are independent exact organizations of the same ordinary class universe, not a claim that all coarse metagraph statistics determine it.

## Binary relation algebra, with the correct rank

Let V be the vector space with basis e_x for all x in X_n. For any field K, let W be spanned by e_x-e_y for every legal relocation edge. Then dim(V/W)=A000568(n). In particular over F2 the relations are e_x+e_y. Proof: in a connected component, a spanning tree identifies all its basis vectors, and the coordinate-sum functional proves that exactly one independent class remains. Distinct components have disjoint supports. Thus rank of the edge-incidence matrix is |X_n|-A000568(n), in every characteristic. Reflection adds the converse-merging relations.

Dually, the algebra of functions f:X_n->F2 satisfying f(x)=f(y) on every relocation edge is the product of one copy of F2 per tournament class. Its primitive idempotents are precisely class-indicator functions. The preceding V/W is only a vector-space quotient, not a quotient of the pointwise algebra: W need not be an ideal. For an edge xy, multiplying e_x+e_y by e_x gives e_x, which need not belong to W. This is a lossless binary formulation of the class relation, not a small description: its domain has exponentially many presentations. Do not substitute the mod-two graph Laplacian, whose rank need not equal the incidence rank. Nor is this the cut/cycle decomposition of the original n-vertex tournament; the vertices here are full presentations.

The smallest nontrivial Laplacian hostile occurs at n=4. In the bit order used by the executable, one relocation component has masks {0,1,2,4,5} and simple edges {01,04,12,15,24,45}: it is K_(2,3) with shores {1,4} and {0,2,5}. Its incidence rank over F2 is 4, but its Laplacian rank is 2. The full cube has four ordinary classes, while the mod-two Laplacian has nullity six. The indicator of {0,2} is in the latter kernel but is not constant on the component. This is why the exact binary quotient uses incidence relations, not the Laplacian. The independent audit supplies direct construction and row-reduction checks.
