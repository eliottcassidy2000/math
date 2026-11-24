# Last seven rows: the old Moore grid determines the missing incidence skeleton

**PROVED conditional.** This is an inverse for a partially filled edge-rooted grid. It does not construct the 49 known rows in a width-56 grid, or solve the column and matching problems in the remaining seven rows.

Use width q, with q rows and q columns. Each pair of different rows is joined by a perfect matching; the same is true of columns. Cells in the same row or column have no common grid neighbor. For cells in different rows and columns there is exactly one direct or two-step route in the grid. Thus for its adjacency B,

    B²+B = (q−1)I + (J_q−I_q)⊗(J_q−I_q).

Suppose s=q−r full rows and the matches between them are fixed. Let C be their known adjacency, and H the rectangular incidence from old cells to new cells in a completion. The old-old block yields

    HH^T = (q−1)I + (J_s−I_s)⊗(J_q−I_q) − C² − C.

Call this the residual. Its diagonal is r; its off-diagonal entries are all 0 or 1. The representation is jointly constrained: a column of H places one new cell's neighbors into ALL the known rows.

## Recover H before assigning a single new column

Temporarily name the unknown new cells as vertices of a hypergraph. Every old cell is a size-r hyperedge, containing its r new neighbors. The hypergraph is linear: no two hyperedges share more than one vertex (else the old cells have two new common neighbors). Its intersection graph is the off-diagonal support of the residual. A vertex of the hypergraph is incident to s old cells, so it corresponds to an s-clique (a star) of the intersection graph.

Assume r≥2 and s>r²−r+1. Take any adjacent pair of residual vertices e,f, whose hyperedges share vertex x. Form W={e,f} together with their common neighbors in the residual graph. It contains the whole star at x. There are at most (r−1)² outsiders: an outsider must contain a point of e outside x and a point of f outside x, and any such pair belongs to at most one hyperedge. An outsider meets at most r star members, since each of its r points occurs with x at most once. Consequently its degree within W is at most (r−1)²−1+r=r²−r. Every true star member has degree at least s−1 in W. The strict inequality gives an exact extraction rule: keep exactly the vertices of degree at least s−1 in W. Running it on residual edges recovers all the stars and therefore the columns of H, up to names. One can skip an edge once its star has assigned it.

This is not merely recognition of large cliques after a search: the degree filter recovers the complete star from any of its edges. In an arbitrary proposed partial grid, wrong star sizes, repeated covering of residual edges, or an old cell in the wrong number of stars are immediate rejections.

## Recover the new row classes

On the recovered new-cell vertices, connect two when they have a common OLD neighbor. No such pair lies within one new row. Between one new cell and any other new row, there are exactly s such contacts: take a two-step route via each of the s old rows, whose outputs must all be different. Thus in the complement graph each row is a q-clique, and every new vertex has exactly r neighbors in each other row.

Two vertices in the same row share at least q−2 complement neighbors in that row. Two in different rows have at most r² common complement neighbors, since in every row at least one of them has at most r neighbors there. If q−2>r² this distinguishes the row classes without labels. At q=56 both strict inequalities hold for 2≤r≤7. In particular **49 specified rows force the incidence skeleton for the other seven**. A synthetic width-56 replay over F7×F8 verifies the extraction at r=6,7; it is an incidence test, not a claimed Moore partial grid.

## Remaining selectors and joint constraints

Each recovered new cell gets one of its r absent column labels (no two of its s old neighbors may share a column). Each new row uses all q labels once. At each old cell, the r new neighbors must have its r missing labels. Given such labels, each pair of distinct new rows supplies at most r−1 admissible new matching partners per new cell: among r complement neighbors one has the same column, which cannot be adjacent. Let D be the new adjacency. The simultaneous tests are the original block equations

    (C+I)H + HD = old/new target,
    H^T H + D² + D = new/new target.

Here the targets are the corresponding blocks of the displayed grid polynomial. These equations, with the known residual and the required permutation matchings, are sufficient as well as necessary. The first stage has therefore recovered ALL incidence placements to the old grid and bounded the remaining matching domains by six in the seven-row case, not solved them independently row by row.

At r=2 the residual is the line graph of a connected bipartite graph on the new cells; recovering row classes by that bipartition works even at q=6, where the complement-count threshold is not strict. The two absent labels at each node give Boolean choices. Each old cell prescribes the unordered pair of labels of its two incident nodes, so one initial choice propagates across the connected bipartite graph, giving at most two global trials and at most one valid label assignment up to swapping the missing rows. All remaining cross guards can also be encoded as 2-SAT; see the [exact small replays](moore_grid_pipeline.md). For r=1 see the direct clique-component inverse there.

## The late-stage adjacency solver can omit all new/new quadratic tests

For the same width-56 regime r≤7, after fixing a viable column assignment the last displayed **new/new equation is redundant** if all old/new guards and the row/column matching conditions are enforced. This yields a fixed-label SAT model of width at most r−1 (and 2-SAT when r=3).

Indeed H^T H has diagonal row blocks sI; its off-diagonal row blocks are J−W_ij, where W_ij is the 0/1 matrix of the complement contacts and every row and column sum of W_ij is r. On the sum-zero space in each new row, the off-diagonal operator has norm at most r(r−1), by the block row norm bound. Since s>r(r−1), H^T H is positive on this space. On new-row constants H acts by summing those constants, so its kernel consists exactly of new-row constants with sum zero.

For the completed candidate B, exact one-per-other-row and one-per-other-column matching makes both row and column partitions equitable. Consequently B commutes with the full target Q=(q−1)I+(J_q−I_q)⊗(J_q−I_q): Q is a linear combination of I, J and the equivalence matrices of those two partitions. The defect F=B²+B−Q commutes with B as well. If the old residual and old/new guards hold, F is supported on its new/new block F_new, so the old/new block of [B,F] gives H F_new=0. Moreover F annihilates all vectors constant on each full row, since the action of B there is J_q−I_q and its quadratic is (q−1)J_q, just as Q. Thus the image of F_new lies in ker H while being orthogonal to it by symmetry, proving F_new=0.

Here is the resulting joint test, after labels:

1. A candidate new edge must be a complement contact in a different new row and have a different label. It also cannot use an old neighbor's column at either endpoint.
2. In each new target row, choose exactly one such edge for every source cell; enforce also one per still-missing column. A block for a row has at most r−1 candidate edges. A column may appear in at most one candidate per new row, so a missing-column block likewise has at most r−1 candidates.
3. For an old cell z and new cell x, require the sum of chosen edges from x to the members of z's incidence hyperedge to equal the old/new target entry minus ((C+I)H)_(z,x). There are at most r−1 eligible terms, since that hyperedge has one member per new row and no edge can stay in x's row. The right-hand side is 0 or 1 unless this candidate is immediately rejected.

The same undirected edge variable occurs in all these tests. To encode a required 1 use its coverage clause plus pairwise exclusions, and for 0 use unit exclusions. This verifies the completed candidate with no quadratic tests solely among new cells. At r=3, once the labels have been fixed, every clause has at most two literals, so standard 2-SAT decides the remaining matching phase. This does not claim the preceding search for those column assignments is 2-SAT at r=3.
