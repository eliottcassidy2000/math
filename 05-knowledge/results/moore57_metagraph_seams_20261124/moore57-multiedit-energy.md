# Moore defect under several simultaneous edge edits

For a simple graph A on n=d²+1 labeled vertices define
F(A)=A²+A−(d−1)I−J and energy ℰ(A)=||F(A)||_F². We specialize
distance statements below to d=57. All statements taking an exact
Moore graph as a reference are conditional, not existence claims.

## Exact energy and interaction

For each distinct pair let c_ij=(A²)_ij be its common-neighbor count,
a_ij its adjacency, and d_i its degree. Then F_ij=c_ij+a_ij−1
off the diagonal and F_ii=d_i−d on it. Summing the squares gives

    ℰ(A) = 8·C₄(A) + 12·C₃(A) + 2d(nd−2e(A)).

Indeed ∑_(i≠j)(c_ij+a_ij)=∑_i d_i²; this cancels its diagonal
counterpart in the expansion. The term ∑ c_ij(c_ij−1) counts eight
orientations of each ordinary four-cycle, and 2∑ a_ij c_ij counts
twelve orientations of each triangle. Thus at fixed Moore edge count,
even through irregular intermediate graphs, the defect energy is
exactly this weighted short-cycle penalty.

For an arbitrary symmetric edit E,

    F(A+E) = F(A) + (A+I/2)E+E(A+I/2)+E².

Writing E as a sum of separate edits, the only cross terms in the
matrix defect are E_i E_j+E_j E_i. Independent edit supports therefore
do not interact at the matrix-defect level when no vertex is in
both supports (the scalar squared norm can still have a cross term).

## Real flat directions and curvature

At an exact d=57 Moore graph, on the orthogonal complement of 1
the adjacency eigenvalues are 7 and −8, of multiplicities 1729 and
1520. A symmetric degree-preserving E annihilates 1. In the
orthogonal eigenspace splitting, write its blocks X,Y,Yᵗ,Z.
The derivative of F maps it to diagonal blocks 15X,−15Z
and zero off-diagonal blocks. Its squared Frobenius norm is
225(||X||_F²+||Z||_F²). Thus the linear flat space among real
symmetric degree-preserving perturbations consists exactly of X=Z=0,
of dimension 1729·1520. If zero diagonal is imposed additionally,
the flat subspace loses at most 3249 dimensions, as its trace
already vanishes. Along a straight flat line,

    F(A+tE)=t² E²,    ℰ(A+tE)=t⁴||E²||_F².

This curvature is positive for every nonzero symmetric E.
Also no nonzero degree-preserving actual 0–1 edge edit can itself be
a linear flat direction: the derivative's diagonal entry at i
is minus twice the number of old edges deleted at i.

The exact same-degree real symmetric zero-defect matrices, before
imposing zero diagonal or Boolean entries, correspond to orthogonal
projectors of various ranks, with adjacency on 1⊥ equal to 15P−8I.
For two actual Moore graphs the ranks agree. In one fixed rank
component the cross blocks above are its usual Grassmannian tangent
directions. For two actual exact endpoints the difference E/15 is
a difference of projectors; its nonzero eigenvalues occur in opposite
pairs by the principal-angle decomposition. In particular all its odd
traces vanish, ||E||_op≤15, and if m is the number of unordered edge
flips then the sum of squared sines of the principal angles is m/225.

## Quantitative locality around an exact reference

Let E=B−A be degree-preserving, with S its set of touched vertices,
k=|S|, and m its number of unordered edge flips. The untouched/
touched block of F(B) is A_(S^c,S) E_S exactly. Therefore

    ℰ(B) ≥ 2||A_(S^c,S) E_S||_F²
         ≥ 2(49−ρ(A_S)^2)||E||_F²
         ≥ 4m(50−k).

Here every column x of E is supported on S, and
||Ax||²≥49||x||² by the full adjacency spectrum. The induced
subgraph A_S has no triangles or four-cycles, so every row sum
of A_S² counts a set of distinct vertices within distance two
(other than the root column's repeated returns) and is at most
k−1: indeed the sum of degrees of the neighbors of a vertex
equals its number of neighbors plus its number of distinct
distance-two vertices. Hence ρ(A_S)²≤k−1. Since x has zero sum,
the exact identity ||Ax||²=56||x||²−xᵗAx gives the stronger
alternative coefficient 2(56−ρ−ρ²) before ||E||².

If B is also exact and E≠0 then k≥51. Indeed k≤49 is ruled out.
At k=50 equality in both norm bounds is necessary, and each
nonzero E-column must be a 7-eigenvector of A and its restriction
a 7-eigenvector of A_S. A component of A_S with spectral radius
7 must have order at least 50 by the same girth bound, so A_S
is connected. Its 7-eigenvector is unique up to scale and strictly
one-signed by Perron–Frobenius, impossible for a nonzero zero-sum
column.

## At least 167 edge changes separate distinct exact labeled graphs

Write t_i for the number of deletions at a touched vertex; the number
of additions there is equal and m=Σ_i t_i. If a changed column x
has t positive and t negative entries, its image Ax has support at
least t(115−3t). Indeed the union of 2t original neighborhoods has
size at least 114t−C(2t,2) since any pair intersects in at most one
vertex; cancellations occur at at most t² distinct common neighbors
between a positive and a negative vertex. Since Ax vanishes outside
S when B is exact, k is at least this bound for the chosen column.
If the minimum t_i is 2 or 3 this forces m≥436 or 954 respectively.
If the minimum is at least 4, k≥51 gives m≥204.

It remains to consider t_i=1. Write E e_i=e_u−e_v, where i–u
is added and i–v deleted. Set T={i}∪N_B(i) and
Y=supp(A(e_u−e_v))\T. The support before excluding T is 112
when u,v are nonadjacent in A, or 114 when they are adjacent.
In either case at most two support vertices lie in T. In both
cases i itself lies there. The old neighbors of i in T cannot
be adjacent to v (no triangle); at most one old neighbor is
adjacent to u (unique common neighbor of i,u). If u,v are
adjacent, this one old neighbor is v, excluded from T, but u
itself is counted; otherwise u itself is not counted.
Thus |Y|≥110 or 112.

The exact identity A E=−E(B+I) implies every y in Y has a
changed edge to T. There are at least |Y| such edges, each
counting only once toward the total incidence at Y. Every
vertex of Y needs at least two changed incidences by degree
balance, so at least ceil(3|Y|/2) distinct changed edges are
incident to Y. The changed edge i–u is internal to T, not
counted; in the nonadjacent case i–v is also uncounted since
v does not lie in the support. Thus m≥167 (nonadjacent case)
or m≥169 (adjacent case). These bounds exclude only small edit
trades; the number 167 is not claimed optimal.
