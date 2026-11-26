# The degree-57 grid as coupled Latin matrices

2026-11-26. Continues the [nested Moore shells](../nested_moore_shells_edge_grid_20261126/README.md). All statements about a complete degree-57 graph are conditional. The computations use the actual three frozen stars and fourteen planted neighborhoods, not a presumed full graph.

## 1. Each full grid row gives a Latin square

Write $m=56$, and $P_{ac}$ for the matching from the columns of row $a$ to those of $c$, taking $P_{aa}=I$. Put

    T_a(c,j) = P_ca(j).

This is Latin: for fixed $c$ the matching is a permutation, and for fixed $j$ the complete neighborhood of $(a,i)$ meets every column other than $i$ exactly once. Define $N_a$ as the one-hot incidence matrix on all grid cells versus the symbols of $T_a$, with cells as rows. The concatenation is $N=[N_0|\cdots|N_{55}]=I+X$. In particular

    N_aᵀ N_a = mI,         N_aᵀ N_b = J−I+P_ab   (a≠b),

with permutation matrices indexed by source and destination. The second identity counts joint labels: each distinct-column pair occurs once by direct adjacency or a route through a different row, and the direct matched pairs are counted twice, once on each anchor row. This is a coupled near-orthogonal family, not ordinary mutually orthogonal Latin squares. It recasts the [edge-deletion grid of Jurišić–Vidali](https://www.sciencedirect.com/science/article/pii/S0195669818300313) in a form usable in the marked fixture.

On the actual full rows 0 and 28, the Gram profile is **56 absent, 3,024 single, 56 doubled**. They also determine a different complete Latin square: put the pair transport through intermediate row $c$ in its row; use the direct matching and identity in the two exceptional rows. This has **1,950 intercalates**, counted by relative-row 2-cycles. Of its 56 row permutations, 52 commute with the physical shift by 4 and 54 commute with shift by 8 (and by 28). So a period-four pattern is visible even after compiling the two physical stars into one Latin square. This is not evidence for an exact global period.

## 2. Four simultaneous matching coordinates, plus a partial fifth

For an incomplete vertex $x$, every prospective neighbor $y$ consumes a missing row slot, a missing column slot, and one missing label from each of $T_0,T_{28}$. Each of these four sets has exactly $55-\deg_X(x)$ open slots. A valid neighborhood must be a **joint** exact cover. The known complete vertices in row 1 yield a fifth, partially labeled axis: each remaining label in its completed subset must be used once, and the unlabeled choices fill a single wildcard of known capacity. This treats existing two-hop routes as boundary conditions rather than as rough residue information.

On the actual **marker-only** start (6,847 residual edges, 126 complete vertices), there are 3,989,339 undirected candidates after the simple row/column slot test and **3,607,365** after short-route exclusions. The two Latin tags do not remove more at this individual-edge stage: any duplicate tag would already close a short route. I tested perfect-matching support on every unfinished matching on both axes, all six pair projections at all 3,010 incomplete vertices, and all four projections with the partial fifth axis at the 2,968 applicable vertices. **No candidate was deleted.** The alternating-cycle support check is exact; these are only pairwise projections of the joint problem.

The extra completed vertices do give independent simultaneous equations. At $(13,0)$ and $(33,0)$ the original 212 slot equations have rank 209; twelve missing routes into the completed part of row 1 raise the rank to **221**, with a consistent right-hand side at all three tested primes (including 2). This is exact rational rank too: the three dependencies between axis totals give the matching upper bound. The two lists contain 2,405 and 2,403 actual contacts. At the marked $(1,1)$ the four axes instead give the familiar 1,741 options and rank 205 out of 208 conditions, with no fifth-axis routes because the planted complete vertices are in its own row. A bounded separate joint search (four axes, the extra literal routes, and pairwise compatibility of proposed contacts) reached partial choices of 42, 41 and 43 new neighbors at $(1,1),(13,0),(33,0)$ in 80,000 nodes per site, respectively. It found no full star and is not an exclusion. Reinstating all 50 reference contacts already fails the local row/column Hall test at $(1,1)$; this agrees with the previous frozen two-hole bound, not a contradiction for the marker-only fixture.

## 3. The partial graph fixes exactly 359 spectral directions

Let $W$ be the orthogonal complement of row and column indicator vectors, and let $U$ be spanned by those indicators and by $e_x,Xe_x$ for the 126 complete grid vertices. The latter vectors are fully known, and $X^2e_x$ is fixed by the edge-grid quadratic. Thus any completion has the same restriction to the invariant subspace $U$.

To separate its positive and negative parts, write the integer-scaled centered vector

    w_x=m²e_x−m(r_x+c_x)+1,
    Xw_x=m²Xe_x+m(r_x+c_x)−(m+1)·1.

The two lists $(X+8I)w_x$ and $(7I-X)w_x$ respectively lie in the 7 and −8 eigenspaces. Each has two exact dependencies, obtained by summing over either complete row. The rank is **124 in each list**, certified modulo 1009 and 1013 with the matching rational upper bound. Including the 111 row/column directions gives **359 fixed directions**. In the remaining 2,777-dimensional common orthogonal complement, any completion is $15\Pi-8I$, with $\Pi$ an orthogonal projection of rank 1,493 (and negative multiplicity 1,284). The hard condition is that adding the fixed part gives actual zero-one adjacency entries and agrees with all remaining known edges.

## 4. Four modular defects per row, not just three

For a balanced nonconstant abelian column tag, call a row matching exceptional when its tag offset is not constant. The [three-defect argument](../nested_moore_shells_edge_grid_20261126/README.md) already excludes fewer than three. Suppose exactly three are exceptional. At each input from that row, the three exceptional forward offsets form the same multiset, since the full row visits all other columns, the tag frequencies are balanced, and the other offsets are constant.

For each of the three exceptional neighbors, omit it as the starting row in the edge-grid two-step sum identity. The other two inverse offset functions, evaluated on the corresponding matching images, sum to the **same** constant $K$ in all three identities (the sum depends only on the nonexceptional translations). Because each evaluation ranges over all columns, their frequency distributions are pairwise reflections by the same involution $v\mapsto K-v$. Three such relations force all three distributions equal. Therefore for any value in the fixed multiset with multiplicity $r$, each function has exactly $mr/3$ instances of it. If $3\nmid m$, then $r$ must be divisible by three; a multiset of size three can then only be constant, contradiction. At degree 57, $m=56$, so there are at least **four** exceptions in every row for every balanced nontrivial residue tag. In the $m=6$ Hoffman–Singleton grid a balanced binary cut instead has exactly three exceptional matchings: its 9 crossing pairs decompose across the five involutory row matchings, each crossing one or three times.

## 5. Further nested-shell consequences

The shell equation $B^TB=(k-1)I+J-Y^2-Y$ has an exact use for a conditional Hoffman–Singleton seed. If a last-shell component has size $s$, normalized component indicator contrasts are centered 7-eigenvectors, so are killed by $B$. Thus every middle vertex has a number of neighbors in that component proportional to $s$. The totals 14 and 700 give exactly $s/50$ neighbors. Every component therefore has size **$50t$**, and the component sizes come from partitions of 14, not arbitrary orders. Inside one component, each seed fiber partitions it into fifty $t$-sets. Across all seed fibers the pairs at distance beyond two in the component are partitioned into clique factors on those sets. Between components of sizes $50t$ and $50t'$ the analogous blocks partition all pairs into rectangles $K_{t,t'}$.

A pentagon seed is guaranteed, unlike a Hoffman–Singleton seed. Its middle layer is a 55-sheet cover of the complementary pentagon, hence cycles of lengths $5\ell$ for a partition of 55. A length-five middle cycle gives a Petersen extension. More generally for either a pentagon or Petersen seed, let $W_0$ denote the middle vectors centered within each fiber. The cross-equation makes $W_0\cap\ker B^T$ invariant under $X$. There, $X^2+X=56I$, impossible because the middle degrees are 2 and 6 respectively. Thus $B^T$ injects $W_0$ into the final shell and sends each eigenvalue $\lambda$ to $-1-\lambda$. Off this image and the constant, the final shell has only eigenvalues 7 and −8.

| Seed | Final-shell vertices and degree | Transferred dimensions | Additional 7 / −8 dimensions | Sum of middle and final pentagons |
|---|---|---:|---:|---:|
| Pentagon | 2,970, degree 52 | 270 | 1,454 / 1,245 | **36,383,149** |
| Petersen, if present | 2,700, degree 47 | 530 | 1,189 / 980 | **21,702,168** |

The multiplicities follow from trace; the displayed pentagon sums from the fifth trace, subtracting the known fiber quotient. The middle pentagon term for the first row is exactly the number of Petersen extensions of the chosen seed. These final shells are connected by the displayed spectra. Neither identity forces an extension or decides existence.

## Replay

Run the included scripts with the small chart inputs from the [Moore57 quadratic-tensor folder](../moore57_quadratic_tensors_20261125/) available on the Python path and in the working directory. Generate with `python grid_packet.py 1 > grid1.txt` (marker only) and `python grid_packet.py 2 > grid2.txt` (all references). Compile `latin_hall.cpp` with `g++ -O3 -std=c++17`; pass `0 28 57` as arguments to emit that vertex's candidate list on stderr. Use `fixed_rank.py`, `extra_routes.py` and `routes_and_moments.py` for the rank, extra equations, and exact moments.

The capped joint search uses `python star_cover.py grid1.txt chosen.dom 11 80000`.
