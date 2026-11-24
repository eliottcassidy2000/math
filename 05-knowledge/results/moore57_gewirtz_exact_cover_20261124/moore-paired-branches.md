# A two-branch exact-cover test, and the edge-rooted Moore grid

This is conditional on a proposed sharply transitive transport row on the 56 Gewirtz labels. No unrestricted sharp row or Moore graph is asserted. The companion replay exercises the same construction on the Hoffman–Singleton graph.

## The grid obtained by changing the root across an edge

Normalize all matching transports \(f_{0i}\) from a pivot fiber to the identity, using independent relabelings. For any fixed nonpivot branch \(i\), the 56 maps
\[
 \{\mathrm{id}\}\cup\{f_{ij}:j\ne0,i\}
\]
are sharply transitive, by the original matching condition for branches 0 and \(i\). Index the remaining leaf vertices as \((i,x)\), with \(i\) one of the 56 nonpivot branches and \(x\) one of the 56 labels. Delete the old root, all branch vertices, and all pivot leaves. The residual graph has 3136 cells and degree 55. Each edge changes both coordinates. Across any two row fibers there is a perfect matching by definition; across any two column fibers there is also a perfect matching by sharpness, with inverse provided by reciprocity. The interchange of coordinates is geometric: it is the same decomposition rerooted at the other endpoint \(v_0\) of \(r v_0\), using \(r\) as its distinguished first neighbor.

More precisely, for two distinct cells, the number of common neighbors **in the grid** is one precisely when the cells have different rows, different columns, and are not adjacent. Cells in a common row or column share a unique deleted vertex, and mismatched rows/columns have no common neighbor outside the grid. Consequently its adjacency \(C\) satisfies the exact matrix identity
\[
 C^2=55I-C+(J_{56}-I)\otimes(J_{56}-I).
\]
Conversely a 56-by-56 grid with perfect row and column matchings and that equation restores the Moore graph by adding the 114 deleted vertices, with their incidence recovered uniquely from the row and column coordinates. This is a compact exact acceptance test for a proposed completed array; testing sharp rows alone is far weaker.

The coordinate-marginal zero-sum spaces both have eigenvalue \(-1\), since each row or column contacts all others exactly once. The mixed space (dimension \(55^2\)) has eigenvalues 7 and \(-8\), as follows immediately from the quadratic identity. The all-ones eigenvalue is 55. These determine the spectrum without any orthogonality assumption on a proposed Gewirtz labelling.

## Literal extension of one sharp row to a pair of branches

Let \(S\subseteq K\) be a sharp row, containing the identity, in a candidate permutation group \(K\) (e.g. the Gewirtz automorphisms). Choose one nonidentity \(g\in S\) as \(f_{ab}\) for two chosen nonpivot anchors \(a,b\). The remaining 54 columns have known \(p_j=f_{aj}\) varying over \(S\setminus\{1,g\}\). Their unknown second-row values are \(q_j=f_{bj}\). The two further sharp tables are precisely
\[
 T=\{1,g^{-1}\}\cup\{q_j\},\qquad
 R=\{1,g\}\cup\{q_j^{-1}p_j\}.
\]
The table \(T\) is sharpness from pivot 0 to anchor \(b\); \(R\) is sharpness from \(a\) to \(b\), as the two-step routes through a third branch are \(f_{jb}f_{aj}=q_j^{-1}p_j\). Thus these are exactly the Moore constraints between pairs of known anchor rows, including all triangles or quadrangles whose matching entries are already specified.

Encode each option \((j,q)\) by one slot and its two permutation supports in the two ordered \(56\times56\) position sets for \(T\) and \(R\). Drop options intersecting the four fixed permutation supports. Then select one option per slot with disjoint support. Because the count is fixed, disjointness covers exactly the remaining cells in both tables. This is finite exact cover, not the assumption that an arbitrary second sharp row can be independently pasted alongside the first. On the other hand it leaves all equations involving pairs of later unknown rows open.

Store the candidate group and candidate sharp row separately, since there need not be any regular subgroup of the group. The [standalone replay](moore-paired-branches.py) includes the exact-cover generator and exercises it on the classic 50-vertex Hoffman–Singleton graph. It verifies the grid common-neighbor identity there and recovers the actual two-anchor choices from eight feasible individual options at each of the four unknown slots. The Gewirtz branch extension is a precise next finite test *if* an unrestricted sharp 56-set is found.
