# The two sharp 54-sheet prism shells, with the distinguished blocks fitted

2026-11-26. Follows the [prime-24 and hexagon-prism packet](../prime24_boundary_and_hexagon_prism_20261126/README.md). **Proved:** the two local sign optima and their connected constructions. **Finite-exact:** both types of distinguished-block placement for each shell, and the complete optimal path-cover families on the four prime windows. **Open:** the remaining ordinary incidence design, the antineighborhood adjacency, and the degree-57 graph.

## 1. Two local extrema are attainable

At 54 sheets, the most four-cycles in one triangle transport is 13. When the two triangles have equal type their signs agree, so an odd number of the three square signs differ from a perfect pairing's sign. Such a square must change at least four sheets: three sheets in a longer cycle would leave 51 to pair. With all three squares made of swaps, the triangles have opposite signs, so the combined number of four-cycles is at most 25 (26 would make both types 4^13 2).

Both limits occur, with connected shells of 324 vertices and 486 edges:

| Case | Two triangle loop types | Three squares |
|---|---|---|
| Maximal equal triangles | 4^13 2 and 4^13 2 | 2^27, 2^27, 4 2^25 |
| All squares paired | 4^13 2 and 4^12 6 | 2^27, 2^27, 2^27 |

Use top parents (0,2,4), bottom parents (3,5,1), and identity vertical edges. Let top edge permutations in cyclic order be I,I,T. The corresponding bottom edges are R0,R1,T R2 (rightmost permutation applied first). So the triangle transports are T and U=T R2 R1 R0, and the square transports have the types of R0,R1,R2. A nonbacktracking closed walk of length three or four in the prism is a face, so excluding fixed points on these five loops excludes triangles and squares in the cover.

On six sheets, take T=[1,2,3,0,5,4], R0=R1=[1,0,3,2,5,4], and use R2=[1,3,4,5,2,0] for equal triangles or R2=[2,4,0,5,1,3] for all paired squares. Add a four-sheet block with T=[2,0,3,1], R0=R1=R2=[3,2,1,0], obtained from the Hoffman–Singleton hexagon (0,1,2,3,28,25) under the indexing of the [preceding replay](../prime24_boundary_and_hexagon_prism_20261126/prime_contacts.py). Immediately after appending, conjugate R2 in the equal case or R0 in the paired case by the transposition of sheet 0 with the new first sheet. Repeat twelve times.

For the splice proof, write U=K Rj, with K=T,j=2 for the equal case and K=T R2 R1,j=0 for the paired case. If s exchanges the first sheets, U changes to (K s K^-1) U s. Within each seed, K sends its first sheet to its predecessor under U, as a direct check of the displayed permutations shows. Let the two cycles at the exchanged sheets be (x0,x1,...,x_(l-1)) and (y0,y1,...,y_(r-1)). They become (x0,y1,...,y_(r-2),x_(l-1)) and (y0,x1,...,x_(l-2),y_(r-1)), of lengths r and l. In particular, the predecessor of x0 does not change, so induction applies. The changed square is conjugated and its type is unchanged. The six-sheet seeds are transitive under T and the three Rj. In each append, the new block is connected by T, the exchanged Rj edges incident at sheet 0 now reach the first new sheet, and sheet 0 reaches the new block along an exchanged edge. Hence the old components remain connected via the new block. The same proof gives both extrema with any number of sheets congruent to 6 modulo 4.

The ordinary missing cells in all fifteen pair boards follow by radius-two projection:

| Shell | Six adjacent-parent boards | Six distance-two-parent boards | Three opposite-parent boards |
|---|---|---|---|
| Equal triangles | Four with 4^27, two with 8 4^25 | Six with 8^13 4 | One matching of 54 edges per board |
| Paired squares | Six with 4^27 | Three with 8^13 4, three with 8^12 12 | One matching of 54 edges per board |

The numbers are component lengths of the graph of missing cells, not of a graph on U. A union of two disjoint matchings has components twice the lengths of the relative permutation's cycles.

## 2. All distinguished blocks can be placed jointly for each of the two shells

Use one column for each of six parents. Each has 54 ordinary symbols and one distinguished symbol, denoted 54. Its distinguished symbol represents the same 55-vertex block as that of its opposite parent. A prospective member of U is a six-tuple recording its unique block in each column.

The two replay inputs from §1 each admit **both** overlap choices for the three distinguished blocks:

- **Concurrent:** one tuple in all three large blocks and 54 tuples exclusive to each of the three, for 163 boundary tuples.
- **Separated:** three tuples in exactly two large blocks (one for each pair) and 53 tuples exclusive to each of the three, for 162 boundary tuples.

These are genuine simultaneous six-coordinate witnesses. In an exclusive group, the ordinary columns are permutations of their necessary labels, minus the labels already used by a double membership in the separated case. The replay searches by permuting entries inside individual columns and then verifies the six-column output directly: all pairs between a large block and a column other than its own or its opposite occur once, every large block covers 55 tuples, no ordinary pair is used twice, and no missing cell is used. The labels are in the graph of the displayed prism shell, not in a reshuffled collection of abstract margins.

There are 2,754 or 2,755 all-ordinary tuples left to place, respectively. These now have an exact clique-decomposition formulation: between different ordinary columns draw an edge if the pair is allowed by the shell and not already used by a boundary tuple. Each of the fifteen bipartite pieces has exactly 2,754 or 2,755 edges. Every ordinary label has degree 51 to **each** other column in the concurrent case. In the separated case just one label per column has degree 52 to every other column (the one appearing in a double membership); all others still have degree 51. Every remaining tuple must be a clique on six vertices with one vertex in each column, and every leftover edge must be covered exactly once. This would complete the incidence of U with the shell, but still would not determine the edges inside U. Neither whole problem is solved here.

## 3. Exactly which prime-sum path contacts are forced

The first 24 odd primes, with neither, either, or both of 2 and 101 added, have exactly three optimal spanning square-only linear forests. The only choice is which two of 3–13, 3–61 and 3–97 to keep. All other selected edges of the original forest stay fixed across the four windows. The new forced edges are 2–7 and 2–79 when adding 2; 43–101 when adding 101.

| Window | Square edges in optimum | Forest components | Nonsingleton paths | Optimal spanning orders up to reversal |
|---|---:|---:|---:|---:|
| P | 17 | 7 | 3 | 60,480 |
| P + 2 | 19 | 6 | 3 | 8,640 |
| P + 101 | 18 | 7 | 4 | 120,960 |
| P + 2 + 101 | 20 | 6 | 4 | 17,280 |

Deleting nonsquare edges of any optimal spanning order recovers a unique optimal forest. Conversely, any ordering and orientation of the forest paths is optimal because no endpoints in different paths have a square contact (that would increase the forest). Thus each forest contributes k! 2^(r-1) orders up to reversal, with k components and r nonsingleton paths. Adding 101 doubles the number of optimal orders without changing the nonsquare join count. The contact 2–7, **forced in every optimal order**, is forbidden in the unique best permutation on all 26 labels.

Run **node prism_boundary.js**, or add **--tuples** for the four explicit distinguished-block incidence certificates. The exact forest enumeration is include/exclude with degree and cycle checks; the boundary routine's result does not rely on optimality of its search.
