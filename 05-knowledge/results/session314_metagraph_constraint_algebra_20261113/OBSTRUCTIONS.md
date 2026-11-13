# Three graphs, three levels of structure

Use the [tiling conventions and closure law](../../../01-canon/theorems/THM-4729-closed-and-open-tilings-hamiltonian-cycle-bipartition-of-the-metagraph-lines-and-the-forced-arrangement.md). This note distinguishes the simple **support** of the actual merged line multigraph, a **score quotient** of its blue sector, and a **labelled cycle-space shadow**. Subgraph means not necessarily induced unless explicitly stated.

## Actual support witnesses

Write a tiling by a bitmask on the non-backbone pairs `(i,j)`, in increasing lex order, with `0<=i<j<n` and `j>i+1`. The fixed backbone is `0→1→...→n−1`, and bit 1 orients its pair backward. A node's label is the minimum mask among all Hamiltonian orderings of its tournament and the converse. As in THM-4729, complement all non-backbone bits for a line; call the line blue if its mask is fixed by mirror `(i,j)→(n−1−j,n−1−i)`. Discard loops and multiplicities for this table.

| Order | Nodes | All edges | Blue edges | Black edges |
|---|---:|---:|---:|---:|
| 4 | 3 | 2 | 1 | 1 |
| 5 | 10 | 18 | 8 | 10 |
| 6 | 34 | 173 | 26 | 151 |
| 7 | 272 | 6032 | 246 | 5836 |

The support can have an edge present in both colors, so the last two columns need not add to the third.

- The first actual (K_5) and (K_{3,3}) are at order 5. The (K_5) uses nodes `[2,4,5,6,13]`; the (K_{3,3}) uses `(2,5,22)|(6,12,13)`. The latter is not induced. The (K_5) is in the combined support, not wholly in either color at this order.
- The first actual Petersen **subgraph** is at order 6, with witness `[2,4,6,5,69,40,90,8,25,9]` in the lexicographic two-subset order on `{0,1,2,3,4}`: disjoint subsets require adjacency. The first **induced** Petersen is at order 7 in combined and black support: `[8193,83,876,297,146,55,80,17,48,18]` in that same order.
- A blue (K_{3,3}) first occurs at order 6: `(8,32,40)|(45,57,90)`. A blue (K_4) at the same order uses `[40,108,90,76]). No blue Petersen subgraph or blue (K_5) occurs through order 7; see the odd-order theorem below.

[Replay](replay_obstructions.py): for each unseen tiling, depth-first enumerate all Hamiltonian paths using its outgoing masks, relabel along each path, and include the mirrors. This produces the exact whole merged fibre, since every tiling in that fibre is one of the two types of path labelling. Then traverse one endpoint from each bit-complement pair. Clique search intersects common neighborhoods; (K_{3,3}) tests the common neighborhood of each triple; Petersen tests all injective assignments of vertices, also removing forbidden edges for the induced variant. All tests are finite-exact, not persistence claims at larger orders.

## Odd blue edges move one score particle

The earlier [THM-060](../../../01-canon/theorems/THM-060-bipartite-skeleton.md) proves odd-blue bipartiteness by triangle parity. A stronger score lift is available. Put (m=lfloor n/2floor). For a blue tiling, the mirror maps outdegree (s_i) to (n-1-s_i), so retain one score from each paired position, say the first (m) scores. The unpaired middle score when (n) is odd is (m). Center these by ((n-1)/2). Their sorted absolute values are an isomorphism invariant, with no dependence on the representative orientation within a pair.

On merged nodes, a complement line can be computed instead by reversing the backbone and taking no other changes: this is converse of the complement tiling. Reversing the path changes (s_0) to (s_0-1), and (s_{n-1}) to (s_{n-1}+1). Only the first pair is changed. At odd order its centered magnitude (|s_0-m|) changes by one, with all other magnitudes fixed. Thus the sorted magnitude partition changes by one Young-lattice box inside the (m)-by-(m) rectangle. The parity of its size bipartitions the odd blue graph, excluding all odd cycles, (K_5), and Petersen as blue subgraphs at all odd orders. A profile distance in the containing Young lattice is a lower bound on graph distance. At even order magnitudes are half-integers; the only way an edge stays in one score profile is a move from (+1/2) to (-1/2) at the initial endpoint, leaving magnitudes unchanged. Otherwise the same one-box rule applies after subtracting (1/2) from all magnitudes. In particular any odd blue cycle at even order needs a same-profile edge.

The possible odd score profiles alone are exactly the partitions dominated by ((m,m-1,ldots,1)): if the descending magnitudes are (a_1,ldots,a_m), Landau's inequalities for the (k) smallest scores say (sum_{ile k}a_ile km-k(k-1)/2) for (1le kle m), and the remaining ones follow by score symmetry. [Eplett's self-converse score theorem](https://doi.org/10.4153/CMB-1979-004-6) gives sufficiency. This is a vertex characterization for the **score** quotient only; individual nodes at one score can differ.

## Exact blue score-edge generating polynomial

Count blue tilings by their first-half *labelled* scores (s=(s_0,ldots,s_{m-1})), using one variable per paired index. The result is a Laurent-product expression with nonnegative integer coefficients:
[
P_n(z)=
 prod_{i=0}^{m-1} z_i^{n-1-i}(1+z_i^{-1})^{u_i}
 prod_{substack{0le i<j<m\j-ige2}}(1+z_j/z_i)
 prod_{0le i<j<m}(1+(z_i z_j)^{-1}).
]
For odd (n), (u_i=2) except (u_{m-1}=1). For even (n), (u_i=1) except (u_{m-1}=0). Start at the transitive score. A same-half mirror orbit moves one score from (i) to (j), except when it is a backbone edge; a cross-half orbit decreases both representatives by one. The fixed pair arc decreases its representative by one, unless that pair is the middle backbone edge at even (n). At odd (n) the orbit to the middle vertex gives a second such term, except at the innermost pair where it is backbone. This exhausts all free orbit bits. Equivalently insert a new outer pair variable (y) at (n+2): for (nge2), multiply (P_n(z)) by
[
 y^{n+1}left(prod_i z_iight)(1+y^{-1})^{u}
 prod_{ige1}(1+z_i/y)prod_{ige0}(1+(yz_i)^{-1}),
 qquad u=2	ext{ at odd }n,quad u=1	ext{ at even }n.
]
The old coordinates are reindexed one place inward.

For each monomial of coefficient (c) and exponent (s), deposit (c) on the unordered score-profile pair given by (s) and (s-e_0) (center, take magnitudes, sort). Divide every final bucket by two: each complement line is counted by its two tilings. This recovers exact blue score-quotient edge multiplicities, including pairs with equal profiles, not individual class adjacency. The [direct score replay](replay_score_quotient.py) found that the one-box rule after excluding nonstrong–nonstrong pairs is sufficient in the quotient through odd order 9 (profile/edge counts (2/1,5/4,15/17,48/71) at orders (3,5,7,9)). Sufficiency at all orders is not proved. Exact coefficients bypass that conjecture.

## (A_5) and Petersen live differently in the shadow

Use the [cycle-space shadow](../../../01-canon/theorems/THM-646-line-metagraph-score-complement-skeleton.md), where the line operation is addition of (K_n) for odd (n), and (K_n) minus the consecutive matching for even (n). At order five the translation is therefore (K_5) itself. The ten triangle vectors of (E(K_5)) are indexed by the two vertices *omitted* from their supports. For two distinct triangle vectors, their XOR is a four-cycle if the omitted pairs intersect, and a bow-tie of two triangles sharing one vertex if the omitted pairs are disjoint. Thus the bow-tie relation on the ten-point orbit is exactly Petersen, with the natural faithful (A_5) action. This complements the old [Petersen/Kneser and Walsh support dictionary](../sumner_t5_petersen_collatz_20261005.md); the bow-tie test is a relation of *triangle vectors*, not the translation edge relation. At even orders 6, 8 and 10 the translation graphs are respectively complete multipartite with 3, 4 and 5 pairs: at 6 it is the planar octahedron (neither Kuratowski graph as a subgraph); at 8 it contains (K_{3,3}) but no (K_5); at 10 it contains (K_5).

The actual order-five combined support has the same number of vertices as Petersen and contains a (K_5), yet its degree classes have sizes (3,2,3,2) for degrees (1,3,5,6), respectively. It has no automorphism of order five (such an automorphism would need an orbit of length five in a degree class), so simplicity of (A_5) implies every (A_5) action by automorphisms is trivial. The abstract (A_5) of the five-clique and the genuine (A_5) of the labelled triangle-vector shadow do not extend to an (A_5) of the actual order-five support. This is a concrete obstruction to reconstructing actual adjacency by carrying over a shadow symmetry.
