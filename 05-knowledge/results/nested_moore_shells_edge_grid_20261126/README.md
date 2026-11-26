# Two kinds of nested Moore shells, and the degree-57 boundary equations

2026-11-26. This is a construction analysis, not a degree-57 graph or nonexistence proof. The degree-57 edge grid is forced for any hypothetical graph; the 50-vertex seed below is conditional on an induced Hoffman–Singleton copy. [Godsil’s survey of problems](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v2i1f1/pdf/) discusses the Petersen-in-Hoffman–Singleton inclusion and the lack of a corresponding forced inclusion in degree 57. The existing [marked chart census](../moore57_quadratic_tensors_20261125/nested_core_check.py) examines root-containing copies in a particular partial chart, not arbitrary copies.

## 1. General shell law

Suppose a diameter-two degree-$k$ Moore graph contains an induced diameter-two degree-$d$ Moore graph $H$, of order $m=d^2+1$. Every vertex outside $H$ has at most one neighbor in $H$: two adjacent seed neighbors would make a triangle, and two nonadjacent seed neighbors already have their unique common neighbor in the seed. Put $F_u$ for the outside neighbors of $u\in H$, $F=\bigsqcup F_u$, and let $Z$ be the rest. The following is an exact consequence of the one-edge-or-one-two-step property.

- Every fiber $F_u$ has size $L=k-d$. Between $F_u,F_v$ there is no edge for adjacent $u,v$; for nonadjacent $u,v$ it is a perfect matching. Indeed an outside neighbor of $u$ needs exactly one route to $v$ and cannot use a different seed vertex. There are no edges within a fiber.
- Every $z\in Z$ picks one neighbor from every $F_u$; its internal degree is $k-m$.
- Every $x\in F$ has exactly $d(d-1)$ neighbors in $F$ and $q=k-1-d(d-1)$ neighbors in $Z$. Consequently $|Z|=Lq$. Distinct rows $z\in Z$ agree in at most one chosen fiber coordinate.

This exposes two numeric regimes, $k_-=d^2-d+1$ and $k_+=d^2+d+1$. At $k_-$ there is **no** final shell; at $k_+$ the final shell has degree $d$, just like the seed, and has $2dm$ vertices, while $F$ has $m^2$ vertices.

| Seed and outer graph | $H$ | $F$ | $Z$ |
|:--|--:|--:|--:|
| Pentagon in Petersen | 5 | $5\times1$ | 0 |
| Petersen in Hoffman–Singleton | 10 | $10\times4$ | 0 |
| Pentagon in Hoffman–Singleton | 5 | $5\times5$ | 20 of degree 2 |
| Hoffman–Singleton in a hypothetical degree-57 graph | 50 | $50\times50$ | 700 of degree 7 |

The literal control checks all 12 pentagons of Petersen and all 1,260 pentagons of Hoffman–Singleton. In **every** Hoffman–Singleton case the 20-vertex final shell is four pentagons; the middle 25 vertices are five more pentagons. The five seed-side pentagons give five coordinates on those 25 middle vertices and every pair of coordinates gives a bijection to the $5\times5$ grid. Each middle pentagon covers the complement of the initial pentagon bijectively, so joining it to the initial one gives a Petersen copy. This recovers the nested picture from the routes, not just from the drawing.

## 2. What a Hoffman–Singleton seed would force at degree 57

Let $X$ be the adjacency matrix on the 2,500 middle vertices, $Y$ the adjacency matrix on the 700 last-shell vertices, and $B$ the 2,500-by-700 matrix of cross-incidences. The middle degree is 42, and $B$ has row sum 14 and column sum 50. Let $P$ be the seed-by-middle fiber indicator. The full Moore identity $A^2+A=56I+J$ gives these equations:

    Y²+Y+BᵀB = 56I+J
    XB+B(Y+I) = J
    X²+X+PᵀP+BBᵀ = 56I+J.

The top shell identities are already imposed by the perfect matchings and one-per-fiber condition. Thus these three symmetric/cross equations, with the incidence and degree conditions above, are also a sufficient exact completion test.

**Local exact cover.** For a middle vertex $x$, its 14 chosen last-shell vertices have pairwise disjoint closed $Y$-neighborhoods of size 8. The cross equation requires the 42 middle neighbors of $x$ to partition the **other 588** last-shell vertices into their 14-element incidence rows. Those 42 neighbors must also pick exactly one vertex in each eligible seed fiber. This gives a concrete solver for a proposed $Y$ and $B$: reject a proposed neighbor by comparing physical incidence rows, then branch on the least covered of the 588 elements or on a required fiber. The analogous test was replayed on every pentagon shell of Hoffman–Singleton. Solving each vertex independently would not suffice; the resulting choices must agree as a symmetric $X$ and pass its quadratic identity.

**An exact spectral split.** If $Y$ has $c$ components, $1\le c\le14$, because each is 7-regular of girth at least five and has at least 50 vertices. On centered $Y$ eigenvectors the squared singular value of $B$ is $56-\theta-\theta^2$. The only possible zero occurs at $\theta=7$ since $\theta\ge-7$. The constant maps nontrivially. Hence

    rank(B)=701−c ≥ 687.

The 49 centered fiber directions of $X$ have eigenvalues $-3$ with multiplicity 28 and $2$ with multiplicity 21, via the complement of the seed. The $700-c$ nonprincipal component directions of $Y$ transfer, through $B$, to eigenvalues $-1-\theta$ of $X$, orthogonal to the fiber directions because each last-shell vertex meets every fiber once. On the remaining space of dimension $1750+c$, the last equation imposes $X^2+X=56I$. Its trace is $700-8c$: the constant and fiber directions cancel in trace, and the transferred trace is $-700+8c$. Thus the remaining eigenvalues are **$7^{(980)},(-8)^{(770+c)}$**. In particular the middle graph is connected: 42 occurs only along the constant.

**Maximum replication is a narrower hypothesis.** If all 14 components of $Y$ attain size 50, they are Hoffman–Singleton copies. A middle vertex cannot attach twice to one of those components (diameter two) and has 14 last-shell neighbors, so attaches once to each. The 15 copies, including the original seed, therefore give an **orthogonal array OA(2500,15,50,2)** on the middle vertices: two vertices in different copies have exactly one common middle neighbor. Even under that saturation hypothesis, the middle graph is connected, not a disjoint collection of small graphs as the five pentagons were in Hoffman–Singleton. The array alone does not prescribe $X$, and neither the seed nor saturation is assumed as a necessary fact about a degree-57 graph.

## 3. A forced edge grid, without a Hoffman–Singleton seed

Fix an edge $uv$. Removing its two closed neighborhoods leaves $(k-1)^2$ vertices. Label each residual vertex by its unique neighbor in $N(u)\setminus\{v\}$ and its unique neighbor in $N(v)\setminus\{u\}$. These labels fill the entire grid, because every pair of opposite anchor labels has one common neighbor. Write $R,C$ for the matrices of equality of row and column, including the diagonal, and $X_e$ for the residual adjacency. Then

    X_e² + X_e = (k−1) I + J − R − C.

In particular each vertex has one neighbor in every different row and in every different column. Two distinct vertices in the same row or column have no route of length at most two inside the grid; vertices different in both have exactly one. Conversely a symmetric zero-diagonal $0/1$ solution to this single identity, on the specified grid, reconstructs the Moore graph by adding the two heads and their anchor vertices. The diagonal gives the degree, the route property gives the distinct neighbor columns and rows, and the added paths restore exactly one route for every missing pair. The general reduction is consistent with [Jurišić–Vidali’s Sylvester graph construction](https://www.sciencedirect.com/science/article/pii/S0195669818300313). The small replays are two disjoint edges at each of the 15 Petersen edges, and a connected 36-vertex grid at each of the 175 Hoffman–Singleton edges.

The two centered row and column subspaces each contribute eigenvalue $-1$; on their joint orthogonal complement after discarding the constant, the quadratic is $X_e^2+X_e=(k-1)I$. At degree 57 the exact spectrum is $55^1,(-1)^{110},7^{1617},(-8)^{1408}$. For the analogous Hoffman–Singleton subgraph it is $5^1,(-1)^{10},2^{16},(-3)^9$.

**A local obstruction to simple nesting by translations.** Fix an arbitrary nonconstant tag map $\phi$ from the $m=k-1\ge3$ columns to an abelian group of exponent dividing $m$ (for example residue labels modulo a divisor of 56). Call a matching tag-translational if $\phi(P(x))-\phi(x)$ is constant in $x$. Every full row has **at least three incident matchings that are not tag-translational**, even if all other rows are unrestricted.

Fix rows $a,b$ and a vertex $x$ in row $a$ with column $i$ whose neighbor in row $b$ has column $j$. For each remaining row $c$, let $x_c=P_{ac}(i)$, and let $z_c=P_{bc}^{-1}(x_c)$ be the endpoint in row $b$ after two steps. The sets of $x_c$ and $z_c$ both enumerate all columns except $i,j$. Thus $\sum_{c\ne a,b}(\phi(z_c)-\phi(x_c))=0$ for every $i$.

- If row $b$ has only tag translations, for every $a$ the sum of shifts other than $b\to a$ vanishes. Every outgoing tag shift is then the same, impossible because some input column has at least two tag values among the other $m-1$ columns.
- Exactly one exceptional matching is impossible because at each input column $i$ the sum of all outgoing tag offsets is $\sum_x\phi(x)-m\phi(i)$, constant by the group exponent assumption. The lone remaining offset must be constant too.
- With exactly two exceptional matchings, choose $a$ to be one exceptional neighbor of $b$. The other is the sole nonconstant inverse tag offset in the zero-sum identity. It is evaluated on all columns by the bijection $P_{ac}$, so must be constant after all.

This excludes a single regular group row at degree 57: if its maps with the identity formed a group, the complete row property makes it regular on 56 columns. Label by that group and tag via its abelianization, which is nontrivial since the group of order 56 is solvable. All its maps would be tag translations. This is a local constraint beyond the [whole-chart derangement-group obstruction](https://www.mdpi.com/2227-7390/14/23/4318). In the Hoffman–Singleton controls all matchings have cycle type $2^3$, so only one can be a cyclic translation in any common column relabeling, as the replay also checks.

**Existing chart.** In the actual marked fixture, use the edge from the root to completed branch 56. The surviving layer really is the 56-by-56 grid (row = ordinary branch, column = its sheet in the identity gauge). It currently has 6,897 residual edges. Rows 0 and 28 are complete, yet no column-pair matching is complete. All 1,540 column pairs already have 2–6 of their 56 entries, while 1,378 row pairs have zero entries. Thus the two coordinate views are not equivalent search frontiers: a column reconstruction can propagate globally from existing row data. Neither complete row star is a group under composition in these columns, so no common relabel turns it into a translation group. Row 0 has 18 mod-two tag defects and row 28 has 16; in the mod-seven quotient they have 44 and 45 respectively. The replay lists the other divisor projections.

## Replays

- `python nested_shells.py`: exact small control shells, orthogonal projections and all edge grids.
- `python edge_chart.py`: run from the [pinned chart directory](../moore57_quadratic_tensors_20261125/) or with that directory on `PYTHONPATH` and as the working directory. Reuses the marked contacts file.
