# Faithful tournament representations of the four diameter-two Moore degrees

Let \(G\) be a diameter-two Moore graph of degree \(d\in\{2,3,7,57\}\), so \(m=d^2+1\) and \(A^2=(d-1)I-A+J\). The last case is conditional on existence; the other cases are the pentagon, Petersen graph and Hoffman–Singleton graph. The [asymptotic flag construction](https://arxiv.org/abs/2608.03965) is not an exact degree-57 input.

## A faithful, optimally short edit-metric chart

Map each graph vertex \(v\) to the binary row \(X_v=\mathbf1_{N[v]}\) of \(X=A+I\). Its weight is \(d+1\), and two distinct rows intersect in exactly two coordinates if their vertices are adjacent, or one otherwise. Therefore their Hamming distance is exactly
\[
 d_H(X_v,X_w)=
 \begin{cases}2d-2&v\sim w,\\2d&v\not\sim w.\end{cases}
\]
Choose \(m\) free arcs in a tournament frame with prescribed Hamilton spine, and put the code on those arcs, fixing unused ones. The minimum-distance graph of these literal rooted tournaments is precisely \(G\). A graph automorphism simultaneously permutes bit coordinates and codewords, giving a faithful action by physical edit isometries of this family. Conversely any isometry permuting these codewords induces a graph automorphism; this describes the induced action exactly, disregarding ambient extra symmetries fixing this family pointwise.

In fact **\(m\) bits is the minimum** for these prescribed distances. Subtract one fixed anchor row, and let \(D\) be the resulting matrix of the other \(m-1\) row differences. The pair distances determine its Gram matrix by polarization. From \(XX^\mathsf{T}=dI+A+J\), its determinant is
\[
 \det(DD^\mathsf{T})=m\det\big((dI+A)|_{\mathbf1^\perp}\big)>0.
\]
Indeed the anchored contrast basis \(e_v-e_0\) of the centered space has Gram \(I+J\), determinant \(m\), and the displayed centered operator is positive definite. Any binary realization in \(m-1\) bits can first XOR by its anchor (an isometry) to send it to zero. Its anchored difference matrix then has square size, integral entries, and would have a square determinant of its Gram. But for \(d=3,7,57\), the two centered eigenvalues of \(dI+A\) are \((4,1),(9,4),(64,49)\), respectively, all squares, while \(m=10,50,3250\) is not a square. For \(d=2\), the two conjugate eigenvalues both have multiplicity two and product 1, so the Gram determinant is 5. This excludes \(m-1\), and rank excludes all smaller lengths. Consequently the smallest fixed-spine chart orders for this exact prescription are \(5,6,12,83\), using \(\binom{n-1}{2}\) free positions.

One can preserve the exact Hamilton count and also make the unmarked tournament classes distinct. Use two strong height-nine cells with rows \((16,1,3,7,14)\) and \((20,1,2,7,14)\). They share the spine \(0,4,3,2,1\), are nonisomorphic, and differ by one nonspine arc. Replace each bit by its cell and use the ordinal join across the \(m\) cells. The ordered SCC list reads the codeword, the height is exactly \(9^m\) at every graph vertex, and the exact edit distances above are unchanged. To make each encoded state blue, append the opposite cells in reverse order, with their spines reversed: reflection exchanges the two halves. Distinct view classes are still recovered from the first half of the ordered SCC list, height is \(9^{2m}\), and the distances are doubled. These are actions on labelled edit families, not automorphisms of a single tournament. Nor are Moore edges one-arc moves in a single fixed-spine cube: that cube is bipartite, whereas every \(G\) here has a pentagon (extend an edge at both ends and use the unique common neighbor of the two new endpoints).

## One strong, blue tournament with an exact observer Gram

Take source vertices \(u_v\) and observer vertices \(w_x\), one of each for every Moore label. Orient \(u_v\to w_x\) precisely for \(x\in N[v]\), and orient the reverse otherwise. Use any total order in the source layer and its reverse in the observer layer. The source order followed by the observer order is a Hamilton spine; the joining arc uses a positive diagonal entry. Reflection swaps \(u_v,w_v\) and reverses every arc, because \(A+I\) is symmetric. Thus this is a blue presentation of a self-converse tournament. It is strong: if \(v,z\) are nonadjacent in \(G\), there is a directed path \(u_v\to w_v\to u_z\), and the complement of \(G\) is connected (minimum degree \(m-d-1\) is at least \(m/2\) for \(d\ge3\), and for the pentagon the complement is a pentagon). All observers then connect through the sources as well.

The cross-layer signed arc matrix is \(M=2(A+I)-J\), and exactly
\[
 MM^\mathsf{T}=4A+4dI+(m-4d)J.
\]
Hence the graph is recovered by a two-valued off-diagonal observer Gram **inside one literal tournament** on \(2m=10,20,100,6500\) vertices. It uses outside witnesses, not the full same-vertex skew Seidel Gram. The observer width \(m\) is optimal for this exact literal two-distance prescription by the determinant argument.

In all four cases the layers can even be decoded from the unmarked tournament after choosing the order appropriately. Index sources \(i=0,\ldots,m-1\), observers by matching labels, and recall that the latter run in reverse order. Their scores are
\[
 s(u_i)=m+d-i,\qquad s(w_j)=m-d-1+j.
\]
Equality across layers occurs precisely when \(i+j=2d+1\). For \(d\ge3\) let \(F=\{2d+2,\ldots,m-1\}\), the set of degree-unique source and observer anchors. Every ambiguous \(u_i\) beats all \(u_z\) for \(z\in F\), and every ambiguous \(w_j\) loses to all \(w_z\) for \(z\in F\). The observer \(w_j\) fails the first test exactly when \(N[j]\cap F\ne\varnothing\); the source \(u_i\) fails the second exactly when \(N[i]\cap F\ne\varnothing\). It is therefore enough that each ambiguous pair has at least one such index.

- For \(d\ge7\), arrange that \(F\) contains some entire closed neighborhood, which dominates \(G\) by diameter two. It fits since \(m-2d-2\ge d+1\).
- For Petersen (\(d=3\)), choose the two anchors adjacent, so their closed neighborhoods cover six vertices. Put the four uncovered labels at indices \(0,1,2,3\), and the four other covered labels at \(4,5,6,7\). Each pair \(i+j=7\) then has at least one covered index.
- For the pentagon, the degree-unique vertices are \(u_0\) (score 7) and \(w_0\) (score 2). Place the two neighbors of label zero at indices 1,2. The ambiguous pairs satisfy \(i+j=5\) and cannot both lie in \(N[0]\). The tests \(u_0\to v\) and \(v\to w_0\) distinguish each pair: an ambiguous source always passes the first, an ambiguous observer always passes the second, and the opposite test fails exactly for labels not in \(N[0]\).

Thus the source and observer layers are canonical from scores and arcs, the transitive order of each layer is read intrinsically, and reverse ranks match the labels. The full unmarked tournament is rigid and faithfully determines the original graph. It deliberately breaks the native vertex symmetry, while the edit-family construction above preserves it.

## Native symmetry on the original vertices

A permutation group \(K\le S_X\) is contained in the automorphism group of some tournament on \(X\) iff it has odd order. A tournament automorphism cannot be an involution (it would exchange endpoints of some edge), while an odd group acts on each unordered pair orbit without any stabilizer interchanging the two ends, so one may choose and transport an orientation on each orbit. There is an equally precise criterion if global converse is allowed: a character \(\epsilon:K\to\{\pm1\}\) must extend the sign of the action of each unordered-pair stabilizer on its two points. Necessity is immediate; sufficiency follows by choosing the orientation at one pair in each orbit and transporting it with sign \(\epsilon\). In particular any involution fixing two vertices and moving some other vertex makes this impossible.

For \(C_5\), orient \(i\to i+1,i+2\bmod5\): its dihedral rotations preserve and reflections reverse the orientation, and a color on the \(\pm1\) pairs records the original cycle. For Petersen, view vertices as two-subsets with disjointness adjacency. A ground transposition induces an involution violating the criterion. For Hoffman–Singleton the same obstruction occurs explicitly in [Robertson's five-pentagon/five-pentagram coordinates](https://citeseerx.ist.psu.edu/document?doi=83191f1e7be2deb6026e7dd0d7c9a6cdee450591&repid=rep1&type=pdf): put a cross edge \(u_{i,j}\sim w_{k,l}\) iff \(l=ik+j\) in \(\mathbb F_5\), and take \(u_{i,j}\mapsto u_{i,-j}\), \(w_{k,l}\mapsto w_{-k,-l}\). This fixes every \(u_{i,0}\) and is a nonidentity involution. At degree57 the criterion itself is conditional: if the full graph automorphism group is odd, it admits such an orientation, and coloring tournament pairs by Moore adjacency makes the preserved group exactly the graph automorphism group. No assertion is made here about whether a degree57 graph exists.
