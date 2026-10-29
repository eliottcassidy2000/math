# Avoiding words have unit-triangular native two-path observers

Status: **PROVED / INDEPENDENTLY AUDITED, all-order argument.** Independent audits: [audit 1](../audit/ALL_ORDER_SPECTRAL_PROOF_AUDIT.md) and [audit 2](../audit2/REPORT.md). Relative to the actual abcde grammar at d3132eaf, not a classification claim for other tournament families. No publication or canon promotion performed here.

## 1. A vertex-sized certificate

Let C be a directed Hamilton cycle of a digraph T; write s(u) for successor on C. Define D_C on the same vertices by an arc s(u) -> v for each arc u -> v of T not belonging to C. It has no loops. A vertex v is a source of D_C exactly when its only incoming T arc is its cycle predecessor's arc.

**Lemma.** If D_C is acyclic, the integer ordered-two-path spanning-cover observer M2 is permutation-equivalent to a unit upper-triangular matrix, and its bipartite support has exactly one perfect matching.

Choose integer weights f strictly increasing along every D_C arc. For any spanning two-path cover with ordered starts (a,b) and ordered ends (c,d), sum over its used arcs u -> v:

    sum [f(v)-f(s(u))] = f(s(c))+f(s(d))-f(a)-f(b).

Each summand is zero for a cycle arc and strictly positive for a noncycle arc. Reindex columns by the row pair (s(d),s(c)); equivalently row (a,b) matches column (pred(b),pred(a)). Every off-diagonal nonzero entry strictly increases pair weight f(a)+f(b): equality forces all used arcs to be cycle arcs, and the only two-path cover using cycle arcs with these starts is obtained by deleting the two predecessor arcs. Its ordered ends are (pred(b),pred(a)), and it contributes exactly one. Order row pairs by pair weight, breaking ties arbitrarily, and use the same order for the reindexed columns. The matrix is unit upper triangular. Such a support has exactly one perfect matching by successive forced choices. This proof allows singleton paths and requires no determinant inference.

## 2. Exact insertion rule for D_C

Insert x after parent p in C, let q=s(p) be its old successor, and let x beat every old vertex except p. In D_C':

* x has no incoming arcs, hence is a source;
* for old vertices other than q, all outgoing arcs are unchanged;
* x has an arc to q and to every former out-neighbor of q in D_C;
* q's outgoing arcs are replaced by q -> v for every old v outside {p,q}.

Indeed the old row p is now sourced at x rather than q, and its former cycle arc p->q becomes the noncycle arc x->q in D_C'. The new row x is sourced at q, giving the last bullet. These exhaust all T arcs.

**Preservation lemma.** Suppose D_C is acyclic. Then D_C' is acyclic if either (i) q is a source of D_C, or (ii) the only incoming arc of q is p->q and p is a source of D_C.

After removing the new source x, any new cycle must contain q, since all other arcs form a subgraph of D_C. In case (i) q still has no incoming arcs. In case (ii), q's only possible predecessor is p, and p still has no incoming arcs: its old incoming set was empty and the replacement row q omits p. Thus no cycle exists.

## 3. Grammar proof

The seed T has cycle 0->1->2->3->0 and chords 0->2,1->3. Its D_C has just arcs 1->2 and 2->3, and is acyclic.

For operations a,c,d, the old successor q of the chosen parent is an indegree-one maximum: for a,c this is the selected first/unique maximum; for d it is the second maximum. Hence q is a source of D_C and case (i) applies.

For operation b, the current type is A, so the previous letter exists and was a or c. In that preceding insertion, the newest vertex p was inserted immediately before an old indegree-one maximum q. Before that insertion q's only incoming neighbor was its predecessor r. Afterwards its incoming neighbors are exactly r and p, with pred(p)=r and s(p)=q. The newest p has indegree one. Therefore in the present D_C, q has the sole incoming neighbor s(r)=p, and p is a source. Case (ii) applies.

For operation e in a word avoiding be and ee, either the word is the one-letter word e or its preceding letter was d (the only letters producing type B are b,d,e). At the seed, p=1, q=2, pred(p)=0, and q's incoming neighbors are precisely {0,1}; p has indegree one, so case (ii) applies. If the previous letter was d, that preceding insertion inserted the current newest/second maximum p immediately before the old second maximum q. Exactly the same incoming-neighbor argument as for b proves case (ii).

Every prefix of an avoiding word avoids be and ee. Induction using these exhaustive cases proves D_C acyclic for every avoiding grammar word, at every order. Section 1 proves the claimed unit-triangular observer, unimodularity, and unique perfect matching. No reverse implication is used or claimed in this argument.

## 4. Stronger and weaker boundaries

The vertex-DAG lemma applies to any directed graph with a marked Hamilton cycle, not just tournaments. It also implies that the adjacency matrix itself, after cycle-column matching, is unit upper triangular: C is the unique cycle cover. The two-path observer conclusion is stronger than determinant absolute value one.

Finite checks through order 12 also find D_C cyclic for every forbidden word. Section 6 proves that reverse structural implication using the strong-bridge-cut theorem; it is not needed for the unit-triangularity proof. Cyclicity alone does not imply singularity of M2 in arbitrary digraphs; the grammar’s singular direction uses the separately audited bridge argument.


## 5. General ordered k-path observers

**PROVED LEMMA (conditional only on D_C acyclic).** For every 1 <= k <= n, let M_k have rows indexed by ordered k-tuples of distinct starts and columns indexed by ordered k-tuples of distinct ends. Its entry counts ordered spanning covers by k nonempty vertex-disjoint directed paths, including singleton paths. If D_C is acyclic, M_k is permutation-equivalent to a unit upper-triangular integer matrix. Consequently all avoiding grammar words have this property simultaneously for all k, by the induction above.

For an ordered start tuple a=(a_1,...,a_k), cut C immediately before every vertex of the start set S. The k resulting directed paths are the unique all-cycle-arc cover with this start tuple. Order them by their starts a_i, and let Phi(a)=(t_1,...,t_k) be their ordered ends. Equivalently, t_i is the predecessor on C of the next member of S encountered strictly after a_i around C (for k=1, the next member is a_i after one full turn).

The map Phi is bijective. Given an ordered end tuple t=(t_1,...,t_k), let E be its underlying end set. The path ending at t_i begins at the successor on C of the previous member of E in cyclic E-order; for k=1 that previous member is t_i itself. This specifies the ordered start tuple uniquely and reverses the cut construction. This also covers k=n, when all paths are singletons and Phi is the identity.

For any cover counted by M_k[a,t], exactly the same cancellation as in Section 1 gives

    sum_used_arcs [f(v)-f(s(u))] = sum_i f(s(t_i)) - sum_i f(a_i).

All summands are nonnegative, and equality holds exactly when every used arc belongs to C. Thus equality holds exactly for the single canonical cover t=Phi(a); its entry is one. Moreover the sets {s(t_i)} and {a_i} agree for a canonical cover, so its column weight equals its row weight. Reindex columns by Phi, sort row tuples by sum_i f(a_i), and use the same order on reindexed columns. Every noncanonical nonzero entry strictly increases that weight. The resulting matrix is unit upper triangular. This argument proves both uniqueness of the support perfect matching and unimodularity, without large determinants or a rational-sector decomposition.

The all-k statement is a graph-general lemma, not an assertion that all unique-Hamilton-cycle tournaments satisfy its DAG hypothesis.


## 6. Equivalent adjacency permanent certificate

**PROVED elementary corollary.** For any loopless digraph T with a marked directed Hamilton cycle C and 0–1 adjacency matrix A,

    D_C is acyclic  iff  per(A)=1.

Here the permanent counts spanning directed cycle covers, with their directed edges (not cyclic starting-point labels). Normalize the Hamilton-cycle matching to the diagonal by setting B[i,j]=A[pred_C(i),j]. Row permutation preserves the permanent. B has diagonal entries one, and its off-diagonal support is exactly D_C. If D_C contains a directed cycle, select the entries of that cycle and every other diagonal entry to obtain a second permutation term of value one in per(B). Conversely, any nonidentity nonzero permutation term decomposes into permutation cycles; at least one has length at least two and gives a directed cycle in D_C. Thus the diagonal is the only nonzero permutation term exactly when D_C is acyclic. In particular, this equivalence concerns the permanent, not the determinant: determinant one alone would not establish it.

Within the actual grammar, the following three predicates are therefore equivalent:

1. the generating word avoids be and ee;
2. D_C is acyclic;
3. the Hamilton cycle C is the unique directed cycle cover, equivalently per(A)=1.

The implication 1⇒2 is Section 3; 2⇔3 is the preceding elementary argument. For the converse use the [strong-bridge-cut theorem](../singularity/FORBIDDEN_WORDS_ARE_STRONG_BRIDGE_CUTS.md): a forbidden word has a partition into two strongly connected induced tournaments R,Q, each of size at least three. Each block has a Hamilton cycle; their disjoint union is a cycle cover distinct from C. Hence per(A)≥2.

For completeness, the needed strong-tournament fact has a short elementary proof. A strong tournament of size at least three has a directed cycle. Choose a longest one Z. If an outside vertex has both incoming and outgoing arcs to Z, a consecutive pair on Z permits inserting it, contradicting maximality. Thus every outside vertex is either dominated by all of Z (class X) or dominates all of Z (class Y). If there are outside vertices, strong connectivity forces both classes to be nonempty and forces a directed path from X to Y wholly outside Z: take a path from Z back to Z that exits to X and first returns from Y. Replacing any arc of Z by its initial endpoint, that outside X-to-Y path, and its final endpoint gives a longer cycle. This contradiction proves Z spans the block.

No permanent recurrence is asserted or needed. This corollary also does not assert that cyclic D_C by itself forces M2 singularity in arbitrary digraphs; the grammar's singular direction uses its separate bridge argument.
