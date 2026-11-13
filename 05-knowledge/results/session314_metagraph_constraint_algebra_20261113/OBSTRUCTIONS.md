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
- A blue (K_{3,3}) first occurs at order 6: `(8,32,40)|(45,57,90)`. A blue (K_4) at the same order uses `[40,108,90,76]. No blue Petersen subgraph or blue (K_5) occurs through order 7; see the odd-order theorem below.

[Replay](replay_obstructions.py): for each unseen tiling, depth-first enumerate all Hamiltonian paths using its outgoing masks, relabel along each path, and include the mirrors. This produces the exact whole merged fibre, since every tiling in that fibre is one of the two types of path labelling. Then traverse one endpoint from each bit-complement pair. Clique search intersects common neighborhoods; (K_{3,3}) tests the common neighborhood of each triple; Petersen tests all injective assignments of vertices, also removing forbidden edges for the induced variant. All tests are finite-exact, not persistence claims at larger orders.

## Odd blue edges move one score particle

[THM-060](../../../01-canon/theorems/THM-060-bipartite-skeleton.md) already proves odd-blue bipartiteness by triangle parity. The following score lift gives more information. For a blue tiling, the mirror pairs scores $s_i$ and $n-1-s_i$. Retain one score from each of the $m=\lfloor n/2\rfloor$ pairs, using the first $m$ positions. Center by $(n-1)/2$, take absolute values and sort. The unpaired middle score at odd order is $m$. The profile can be read from the full score multiset, so it is invariant under isomorphism and converse.

Reverse the backbone, which gives the same line in the merged graph. Only its endpoint scores change, by $-1,+1$. At odd order the sorted magnitudes change by one Young-lattice box in the $m$ by $m$ rectangle. Their total parity recovers bipartiteness; their profile distance in the rectangle is a lower bound on blue distance. At even order subtract $1/2$ from the centered magnitudes; the rule is the same except that a move from centered score $+1/2$ to $-1/2$ at the initial endpoint leaves the profile unchanged. Any odd blue cycle at even order requires a same-profile edge.

At odd order the feasible profiles are exactly the partitions dominated by $(m,m-1,\ldots,1)$. For descending magnitudes $a_i$, Landau's inequalities on the lowest $k\le m$ scores become $\sum_{i\le k}a_i\le km-k(k-1)/2$; the rest follow by symmetry. [Eplett's self-converse score theorem](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/ED7BC2A1F6FF2D7F31547D771CE69B80/S000843950005743Xa.pdf/selfconverse_tournaments.pdf) gives sufficiency, and the inherited anti-Rédei law supplies a blue presentation of each self-converse type. This classifies the **score profiles**, not the classes at each profile.

## Exact blue score-edge generating polynomial

Keep the first-half labelled scores as exponents $s=(s_0,\ldots,s_{m-1})$. The polynomial is
$$
P_n(z)=
 \prod_{i=0}^{m-1} z_i^{n-1-i}(1+z_i^{-1})^{u_i}
 \prod_{\substack{0\le i<j<m\\j-i\ge2}}(1+z_j/z_i)
 \prod_{0\le i<j<m}(1+(z_i z_j)^{-1}).
$$
For odd $n$, $u_i=2$ except $u_{m-1}=1$; for even $n$, $u_i=1$ except $u_{m-1}=0$. Start at the transitive score. A same-half free mirror orbit moves one win from $i$ to $j$. A cross-half orbit decreases both representatives by one. A fixed pair arc decreases its representative by one except at the middle even backbone edge; the odd middle vertex gives a second such factor except at the innermost pair. These are all free mirror orbits.

Inserting a new outer pair with variable $y$ gives, for $n\ge2$, the update
$$
P_{n+2}(y,z)=
y^{n+1}\left(\prod_i z_i\right)P_n(z)(1+y^{-1})^u
\prod_{i\ge1}(1+z_i/y)\prod_{i\ge0}(1+(yz_i)^{-1}),
$$
with $u=2$ for odd $n$, $u=1$ for even $n$. Each product on the right ranges over the old indices meeting its condition. Each monomial deposits its coefficient on the unordered profile pair obtained from $s$ and $s-e_0$ by centering, taking magnitudes and sorting. Divide by two at the end to count lines rather than their two tilings. This recovers exact quotient multiplicities, including pairs with equal profiles, not the edges of the unquotiented graph.

The [direct score replay](replay_score_quotient.py) finds through odd order 9 that every feasible one-box pair occurs except when both profiles are nonstrong. The profile/edge counts are $2/1,5/4,15/17,48/71$ for $n=3,5,7,9$. That candidate sufficiency is not proved in general; the exact coefficient rule does not depend on it.

## (A_5) and Petersen live differently in the shadow

Use the [cycle-space shadow](../../../01-canon/theorems/THM-646-line-metagraph-score-complement-skeleton.md), where the line operation is addition of (K_n) for odd (n), and (K_n) minus the consecutive matching for even (n). At order five the translation is therefore (K_5) itself. The ten triangle vectors of (E(K_5)) are indexed by the two vertices *omitted* from their supports. For two distinct triangle vectors, their XOR is a four-cycle if the omitted pairs intersect, and a bow-tie of two triangles sharing one vertex if the omitted pairs are disjoint. Thus the bow-tie relation on the ten-point orbit is exactly Petersen, with the natural faithful (A_5) action. This complements the old [Petersen/Kneser and Walsh support dictionary](../sumner_t5_petersen_collatz_20261005.md); the bow-tie test is a relation of *triangle vectors*, not the translation edge relation. At even orders 6, 8 and 10 the translation graphs are respectively complete multipartite with 3, 4 and 5 pairs: at 6 it is the planar octahedron (neither Kuratowski graph as a subgraph); at 8 it contains (K_{3,3}) but no (K_5); at 10 it contains (K_5).

The actual order-five combined support has the same number of vertices as Petersen and contains a (K_5), yet its degree classes have sizes (3,2,3,2) for degrees (1,3,5,6), respectively. It has no automorphism of order five (such an automorphism would need an orbit of length five in a degree class), so simplicity of (A_5) implies every (A_5) action by automorphisms is trivial. The abstract (A_5) of the five-clique and the genuine (A_5) of the labelled triangle-vector shadow do not extend to an (A_5) of the actual order-five support. This is a concrete obstruction to reconstructing actual adjacency by carrying over a shadow symmetry.

## Missing simple roots and the 4-to-5 and 8-to-9 insertions

After dividing out the transitive monomial, the even half-score formula is precisely the product of $1+e^{-\alpha}$ over the **non-simple positive roots of type $B_m$**. The roots $e_i-e_j$ encode same-half score exchanges, $e_i+e_j$ encode cross-half choices, and $e_i$ encode the fixed arcs. The pinned backbone deletes the adjacent differences and the final short root, exactly the simple roots.

The center insertion at fixed $m$ gives
$$
P_{2m+1}(z)=P_{2m}(z)\prod_{i=0}^{m-1}(1+z_i).
$$
The old transitive scores each rise by one and the unary exponent each rises by one. There is a new free mirror orbit at every pair: at the innermost pair the old central backbone arc is freed and the new two center arcs are fixed. In centered coordinates the extra factor is $\prod_i(z_i^{1/2}+z_i^{-1/2})$, the spin-weight pattern. In particular $4\to5$ supplies two bits and $8\to9$ supplies four bits, multiplying the blue tiling totals by four and sixteen. This is a precise score transfer, not an eight-period theorem about whole metagraphs.

The [coefficient DP](replay_score_dp.py) also verifies the one-box support conjecture through odd orders 11 and 13, with respectively 161 profiles and 295 edges, then 557 profiles and 1219 edges. It enumerates score monomials, not the blue tiling cube. This is evidence for a score-quotient conjecture; it does not force edges between individual classes.
