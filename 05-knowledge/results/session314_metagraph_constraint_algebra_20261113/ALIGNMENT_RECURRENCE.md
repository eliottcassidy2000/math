# A score-pruned algebraic evaluator for the actual merged line graph

The objects here are complement lines in the fixed-Hamiltonian-path tiling cube, not the single-arc-reversal graph. A node is an isomorphism class modulo converse. The inherited [three-sorted incidence](../../../01-canon/theorems/THM-796-three-sorted-recursive-tiling-line-node-incidence.md) gives the complete two-face pullback on labelled presentations, with an overlap face and one apex bit. [THM-4729](../../../01-canon/theorems/THM-4729-closed-and-open-tilings-hamiltonian-cycle-bipartition-of-the-metagraph-lines-and-the-forced-arrangement.md) proves that reversing path arcs instead of all other arcs gives the same merged graph. [THM-4731](../../../01-canon/theorems/THM-4731-tiling-line-metagraph-self-lines-by-twisted-burnside-sigma-twins-depth-and-the-tribonacci-partner.md) computes self-line witnesses using twisted affine systems. Here the path-reversal convention yields a determinant and an endpoint-specific recurrence for **every** adjacency entry, not merely the diagonal or its marginal. The procedure accepts already constructed class representatives; it is not a small scalar recurrence for unlabelled nodes.

## Disagreement criterion

Fix labelled representatives T and U at order n ≥ 3. For an alignment f:V(T)→V(U) (a bijection), let D_f be the directed graph of pairs on which T and f*U disagree, using T's arc directions. Then f is an actual half-edge witness if and only if D_f is a spanning directed path. Indeed, in that case reversing exactly those arcs gives f*U. Conversely every physical Hamiltonian path P of T whose reversal lies in U admits exactly |Aut U| such f, and the path is then recoverable from f.

For color c, let W_c(T,U) count these alignments, taking blue if the map sending the ith vertex of the recovered path to the (n+1−i)th is an anti-automorphism of T, and black otherwise. This is the fixed-path reflection criterion, and path reversal preserves the color. Write R(u) for the one or two ordinary class representatives constituting the merged node u. Its oriented color adjacency entries are
$$
 A_c(u,v)=\sum_{T\in R(u),\,U\in R(v)}\frac{W_c(T,U)}{|\mathrm{Aut}(T)|\,|\mathrm{Aut}(U)|}.
$$
They count all actual outgoing half-edges. For distinct nodes this is the number of lines, and on the diagonal it is twice the loop count. Equivalently, use any one T from R(u) and multiply by |R(u)| by simultaneous converse symmetry. The product action of Aut(T) and Aut(U) on valid alignments is free: an automorphism of T in a stabilizer must preserve the unique directed disagreement path and hence fixes every vertex, after which the target automorphism is also the identity. So the formula can retain actual orbit representatives rather than just their numbers. Path reversal and reversing the recovered vertex order give the partner half-edge, using inverse alignment and the indicated label transport. This does not merge parallel edges.

## An exact endpoint recurrence and determinant filter

Preselect source a and end b of the discrepancy path. Its required indegrees are 1−[i=a], and its outdegrees are 1−[i=b]. The target must satisfy
$$
 s_U(f(i))=s_T(i)-[i=a]+[i=b].
$$
Reject the ordered endpoint pair if those shifted scores have the wrong multiset. Otherwise recurse over a partial injection f. When assigning a new vertex, add every disagreement arc to the earlier assigned vertices. Reject if any in/out degree exceeds its requirement, or if a directed cycle closes. At completion require every deficit to vanish. The remaining graph then has n−1 arcs, no directed cycles, at most one incoming and outgoing edge per vertex, and is necessarily a spanning directed path. A next vertex can be chosen by smallest available score bucket, and degrees and cycles can be maintained using disjoint directed path fragments. Even without this second pruning, the alignment count to test for a given endpoint pair is at most the product of factorials of target score bucket sizes. This is effective all-order enumeration, not a complexity bound uniform over score profiles.

There is a complementary algebraic formulation. For f in the full symmetric group form the degree monomial
$$
 \prod_i x_i^{\operatorname{out}_{D_f}(i)}y_i^{\operatorname{in}_{D_f}(i)}
$$
weighted by det(I−D_f), where D_f now also denotes the adjacency matrix. Retain only the squarefree monomials of total degree 2(n−1) (or the exact degree pattern for a,b). For these monomials the graph is a disjoint union of directed paths and directed cycles. If it has a cycle then det(I−D_f)=0; otherwise the adjacency matrix is nilpotent and the determinant is 1. Thus evaluating the retained part at ones gives the uncolored W(T,U), with the blue coefficient obtained by the reflection check on its surviving paths. The determinant is an exact cycle guard, not a spectral approximation.

**Hostile.** Take T=TT2→C3 and U obtained by reversing the internal TT2 arc and all three triangle arcs, keeping the ordinal cross arcs. The identity alignment satisfies all the score and degree tests and gives n−1 discrepancies, but the discrepancy is a two-vertex directed path plus a directed triangle. The determinant vanishes. In fact U is isomorphic to T and there is no loop at this nonstrong node. This illustrates why the recurrence must retain path connectivity, even if every score condition has been imposed.

## A factored bivariate score projection

A faster generating function exists after summing over class pairs with prescribed score words. This refines the ordinary score factorization in [THM-554](../../../01-canon/theorems/THM-554-tiling-score-partition-function.md) to the **joint** ordered scores and both colors. It is a coarse projection, not a replacement for the alignment recurrence.

Fix P=(1,...,n), write z_i=x_i y_i, and put E={{i,j}:j>i+1} for the free pairs. In the path reversal convention, fixed path arcs contribute
$$
 M=\prod_{i<n}x_i\prod_{i>1}y_i.
$$
Each free pair has the same orientation in both endpoints, giving factor z_i+z_j. Thus the polynomial of all oriented half-edges, weighting the two aligned score words by x and y respectively, is M times the product of these factors. For the blue polynomial put r(i)=n+1−i. Each reflection-fixed free pair contributes z_i+z_j, while a two-element orbit with representative {i,j} contributes z_i z_{r(j)}+z_j z_{r(i)}. This is because a choice i→j forces r(j)→r(i). Multiply over free-pair orbits and by M. The black polynomial is the difference of the two. Multiplying one new endpoint strip (all) or adjoining the new reflection orbits in a two-end growth (blue) updates the products. Collecting monomials by sorted score words, then identifying score words with sorted complements, gives the exact score-block adjacency sums; no class canonicalization is needed for those sums. Distinct classes within a score block can have different entries, and no claim of reconstructing the finer graph from these sums is made.

## Independent replay

The companion [Python replay](alignment_replay.py) builds the colored matrix by both fixed-path tiling enumeration and the score-pruned disagreement recurrence. They agree entrywise, including loop doubling. It also multiplies both score polynomials independently and checks every resulting coefficient. At n=3,4,5 the numbers of merged nodes are 2,3,10; half-edges 2,8,64; blue half-edges 2,4,16; and distinct merged node pairs in support 1,3,20. All calculations are exact. The algorithms and determinant proof hold at all orders.
