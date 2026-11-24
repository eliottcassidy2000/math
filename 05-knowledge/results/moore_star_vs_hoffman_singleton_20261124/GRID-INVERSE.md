# Two-block attachment grid: constructive local factors and a unique exterior inverse

This continues [the rigid closed-decagon branch](CLOSED-DECAGON.md). We have a prescribed fixed Hoffman–Singleton copy \(F\), a prescribed residual Hoffman–Singleton component \(F'\), and 2500 attachment sites \(C=F\times F'\). The two copies are anticomplete, and the order-five action \(t\) fixes \(F\) and acts on \(F'\) with ten independent orbits. The remaining \(650\) sites \(Z'\) have not been constructed.

## 1. Both coordinate bijections can be satisfied together

Take the simple complement of \(F\), and the voltage multigraph of the simple complement of \(F'\) modulo \(t\). Each is 42-regular; the latter has ten vertices, with two loops of shifts \(\{\pm1\}\), \(\{\pm2\}\) at each vertex, counting twice for degree. Orient both multigraphs by closed Euler walks. Their arcs give 21-regular bipartite multigraphs on source and target copies of the vertices. Successive Hall matchings decompose each into 21 directed factors. In the left factor these give permutations \(\sigma_i\); in the right factor they lift to \(t\)-equivariant permutations \(\psi_i\), retaining the phase on each matched arc.

Pair the factors in any order and put an undirected edge \((f,h)\sim(\sigma_i(f),\psi_i(h))\). The edges and their inverses project bijectively to the 42 complement-neighbors of \(f\) and to the 42 complement-neighbors of \(h\) at every site. In particular there are no duplicate edges, and \(t\) preserves the grid graph. The replay constructs all 52,500 edges and checks both bijections at all 2500 vertices. Its first paired factors give 12,875 triangles and 384,215 squares. This is a baseline, **not** a girth-five seed or a lower bound on cycles. Pairing the factors is just one subclass of all possible matchings.

### At least 550 of the row matchings must use non-automorphisms

For each edge of the complement of \(F\) there is a reciprocal matching of the two corresponding 50-site rows. Its bijection must commute with \(t\) on \(F'\) and send each site to a distinct nonneighbor there (it need not preserve any adjacency of \(F'\)). The centralizer of \(t\) in \(\operatorname{Aut}(F')\) consists of five even rotations of the doubled ring and five uniform phase shifts. In fact a commutant preserves each doubled-edge color. The color-preserving reflections interchange parity and fail to preserve the distance-two versus distance-four simple edges. Rotations preserving color are the five even rotations, and the zero centers force the change of phase origin to be constant on the ring. These 25 maps all are automorphisms, but only the 20 with nonzero shift send every site to a distinct nonneighbor. At any one row all 42 row-matching maps are distinct. Therefore at least 22 use non-automorphisms, counting to at least \(50\cdot22/2=550\) of 1050 reciprocal matchings.

A library of equivariant involutions is also impossible: a nonzero translation inside a five-orbit cannot be an involution, but every site needs its four other orbit sites as second-coordinate destinations. At each fixed orbit of the second coordinate, these destinations form a 4-regular subgraph on the 50 first-coordinate sites made of the \(+1\) and \(+2\) directed permutations and their reverse edges.

## 2. Exact-cover handoff once the grid has girth five

Let \(A\) be the adjacency of a 42-regular grid satisfying both coordinate bijections and having no triangles or squares. Let \(N\) be the two-coordinate incidence with the 100 fixed/block symbols. Define the remaining common-neighbor demand on grid sites by
\[
D=56I+J-NN^\top-A^2-A.
\]
Its diagonal is 13 and its off-diagonal is 0 or 1. To see that off-diagonal demand cannot be negative, two grid sites sharing a coordinate cannot also be adjacent or share a grid neighbor, by the coordinate bijections. Distinct sites with no shared coordinate can be adjacent or have at most one grid common neighbor, not both, by girth. The rows of \(D\) sum to \(56+2500-100-42^2-42=650\).

Consequently the symmetric off-diagonal part is a **637-regular missing-pair graph**. Any additional vertex \(z\in Z'\) must have a grid neighborhood \(B_z\) containing one site in every first-coordinate row and every second-coordinate column: the ambient Moore equation gives one common neighbor with each symbol in each copy. Its 50 sites form a clique in the missing-pair graph. The 650 such cliques must decompose the edge set, not merely cover it. Conversely any edge decomposition into 50-cliques yields a transversal at each clique, since a missing pair never repeats a coordinate. Every grid site lies in \(637/49=13\) cliques, and their grid incidence matrix \(T\) satisfies \(TT^\top=D\). No exterior adjacency is assumed at this step.

## 3. There is no separate search over the 650-site adjacency

**Pinned Laplacian inverse.** For any such clique decomposition there is at most one exterior adjacency compatible with the Moore conditions. Compute it by small rational linear solves and test it.

Let \(S\) be the graph on the 650 cliques joining disjoint ones, and let \(R\) denote the desired exterior adjacency. Distinct cliques meet in at most one site. A clique intersects exactly \(50(13-1)=600\) others, so \(S\) has degree 49 and
\[
T^\top T=49I-A_S+J=L_S+J.
\]
The three still-relevant conditions on \(R\) are
\[
R\text{ simple, symmetric, 7-regular, supported on }S,\qquad
TR=J-T-AT,\qquad R^2+R=6I+A_S.
\]
The middle equation is the grid/exterior block of the Moore identity. The last is the exterior/exterior block, after inserting \(T^\top T\). Together with the exact cover, the two original copies and the coordinate bijections, these equations are also sufficient: the copy/copy and copy/grid blocks follow from those bijections and the HS identity; the copy/exterior blocks follow from transversality.

For column \(j\), write \(Y_j=\mathbf1-T_j-AT_j\). If a solution exists its column is supported on the \(S\)-component of \(j\), its own diagonal is zero, and
\[
L_S R_j=T^\top Y_j-7\mathbf1.
\]
Solve this on that component with the coordinate \(j\) pinned to zero: remove its row and column, solve over \(\mathbb Q\), then reinsert zero and extend by zero to other components. The pinned matrix is positive definite, by the sum-of-edge-squares identity on a connected component. This determines every column. Reject unless the resulting matrix is binary, simple, symmetric, 7-regular, supported on \(S\), and satisfies the two displayed block equations directly. Thus no guessing of radius-two graph roots or exterior edges is required **after** selecting the clique decomposition. If the grid and decomposition are \(t\)-invariant, uniqueness makes any accepted output \(t\)-invariant too.

In the all-closed orthogonal-array case, cliques in different residual components intersect exactly once; then \(S\) is the disjoint union of thirteen \(K_{50}\)'s and the solves each have dimension 49. The same rule is tested on the known Hoffman–Singleton graph using two anticomplete pentagons for coordinates and the other three as the exterior. In that replay, all 15 columns of the three leftover pentagons are recovered from fixed grid edges and the missing-pair clique decomposition. Nothing here constructs the girth-five 2500-site grid or the degree-57 graph; it splits the remaining work into a girth constraint and a clique exact cover with a deterministic last check.
