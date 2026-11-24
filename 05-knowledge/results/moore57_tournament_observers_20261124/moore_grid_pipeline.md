# Root-and-anchor permutation words for an exact Moore completion

Put $d=k-1$. Fix a root and one distinguished neighbor (the anchor). Label the anchor's $d$ second-layer sites by columns. Each other branch $i$ has a perfect matching to the anchor branch, so label its sites by the corresponding columns. There are $d$ such branches (rows) and $d^2$ grid sites. The two boundary vertices are the root and anchor; they are adjacent. The root is joined to $d$ row heads, the anchor to $d$ column heads; every grid site is joined to its row and column head. No other edge outside the grid can occur by girth.

Every pair of distinct rows is joined by a matching: for a cell outside row $j$, its unique length-two path to row head $j$ has its intermediate site in row $j$. The same argument applies to columns. Distinct cells in a common row or column already share their head, hence no internal edge or grid common neighbor is possible. Cells in different rows and columns have an edge or one common grid neighbor, since they have no common boundary neighbor. Thus the grid adjacency $B$ satisfies
$$B^2+B=(d-1)I+(J_d-I_d)\otimes(J_d-I_d).$$
Conversely symmetric zero-diagonal binary $B$ satisfying this equation reconstructs a Moore graph via the fixed boundary attachment. Diagonal identifies grid degree $d-1$, the off-diagonal same-row/column zero excludes edges and common neighbors there, and thus each cell has at most one neighbor in each other row/column and exactly one by degree. This resolves all grid/boundary distance-two tests.

## The simultaneous word check

To each grid cell $u=(i,x)$ assign a word $f_u$ on the row alphabet: $f_u(i)=x$, and $f_u(j)$ is the column of its neighbor in row $j$ when $j\ne i$. Each word is a permutation (neighbors use every other column). In a complete deck, require reciprocal arrows:
$$f_{(i,x)}(j)=y \Longrightarrow f_{(j,y)}(i)=x.$$
The full agreement constraints for distinct words are:
- common row or column: zero agreements;
- sites joined by a reciprocal arrow: exactly two agreements, at their two home rows;
- otherwise: exactly one agreement, at neither home row.

These checks are also sufficient. So every two queried rows, and all the other projected coordinates, fit one residual constraint, not independently chosen sharp sets.

## Parallel-clique residual acceptance

Specify $s$ rows and their internal adjacency $C$. Between each two specified rows require a matching, and no repeated neighbor column incident to a site and no same-column edge. Put
$$M=(d-1)I+(J_s-I_s)\otimes(J_d-I_d)-C^2-C.$$
In a completion, $M=HH^\top$ where $H$ is the incidence matrix from known to unknown grid sites. Thus $M$ is positive semidefinite, has diagonal $d-s$, and has 0/1 entries off the diagonal. Its off-diagonal edges must decompose into $d-s$ parallel classes of $d$ cliques of order $s$ (each clique is one new cell's set of known neighbors). A class partitions the known sites. Each clique meets every old row exactly once, in distinct columns, and receives a column absent from its sites. Each class uses every column once. At an old site, the new column tags of its containing cliques must exactly exhaust the old closed-neighborhood column complement. These are usable multiway exact-cover constraints even without solving the remaining new–new block.

## Last-row inverse

For $s=d-1$, the factorization forces the off-diagonal graph of $M$ to be $d$ disjoint cliques of order $d-1$. Include the diagonal 1 in each clique's incidence. Each clique must meet every old row in distinct columns and receives its unique missing column; these missing columns must differ across cliques and each old site's assigned tag must be its unique old closed-neighborhood missing column. These conditions are both necessary and sufficient. Indeed, a clique has no internal old edge or common old neighbor, as indicated by its residual edges, so $(C+I)H$ is binary. If an old site has the same column as the tag, it cannot neighbor a clique member because that member is missing the tag; it cannot be a member either. Each column of $(C+I)H$ has $(d-1)^2$ ones by the known regularity, exactly the number of old sites not in its new cell's column. This is the full cross-block equation; the new–new block is $H^\top H=(d-1)I$.

For two missing rows, once a factorization and tags are chosen, the new–new matching is uniquely determined on its off-diagonal row block: there $B_{\rm new}^2=0$ and $B_{\rm new}=(J_d-I_d)-(H^\top H)_{12}$. Checking the matching, columns, and the cross block then gives an exact completion check.

The Python replay reconstructs the last row across all 30 rooted/anchored Petersen choices and all 350 Hoffman–Singleton choices. It builds the latter from five pentagons and five pentagrams linked by $y=x+ij$ modulo five, and checks the Moore and grid equations directly.

### Two missing rows: line-graph inverse and 2-SAT

For $s=d-2$ and $d\ge5$, the residual graph with adjacency $M-2I$ must be the line graph of a simple bipartite graph $L$. Its $2d$ vertices are new grid cells, its two parts their rows; each old grid site gives the edge joining its two new neighbors. Simplicity follows from $M_{uv}\le1$, and every new site has $d-2$ old neighbors. The grid equation says that $L$ is exactly $K_{d,d}$ with the equal-column matching and the actual new–new matching deleted: other new-row pairs need an old common neighbor, and no within-new-row-length-two path is possible across the two rows. This bipartite graph is connected: any two vertices in one part share at least $d-4>0$ neighbors. Its stars are the maximal cliques of the line graph, since it is triangle-free and degree at least three. Explicitly each residual edge's two endpoints plus their common neighbors give exactly one star. Thus their star decomposition and missing row partition are unique up to swapping the rows.

After checking this reconstruction, give each star two possible labels, the columns it does not use among its $d-2$ old sites (which must lie in distinct rows and columns). Each old cell needs its two adjacent stars' labels to exhaust its missing columns. Require all labels distinct within each side of $L$, and no equal label across any edge of $L$. For a new star $v$, exactly one of its two nonneighbors $a,b$ on the opposite side now has the same label, say $a$; its new matching partner is then $b$. For each possible common label $\ell$ of $v,a$, the precise old/new column equation is
$$ (C+I)H_v+H_b=\mathbf1_{\text{old cell column}\ne\ell}.$$
If it fails, forbid the pair of assignments $\operatorname{tag}(v)=\operatorname{tag}(a)=\ell$. Check the symmetric choice with $b$ in the role of $a$. All constraints in this paragraph exclude unary or pairwise Boolean choices, giving 2-SAT. The reconstructed new–new block already satisfies its Gram equation because precisely the equal-column and matching pairs are missing from $L$. This yields an exact polynomial-time existence test when only two rows are missing.
