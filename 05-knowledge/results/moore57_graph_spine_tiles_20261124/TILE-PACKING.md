# A Hamilton spine turns Moore girth into an exact tile packing

Work with an **undirected graph spine** P=(v_0,...,v_(n−1)): every consecutive pair is an edge of G. The associated tournament orients all those spine arcs forward. On the free tiles, the sparse convention points a graph chord backward and all graph nonedges forward. Its native complement points graph chords forward and graph nonedges backward. The complement is one metagraph step in the same spine, not a change of graph.

| degree d | vertices n | all graph edges | off-spine graph edges | graph nonedges |
|---:|---:|---:|---:|---:|
| 2 | 5 | 5 | 1 | 5 |
| 3 | 10 | 15 | 6 | 30 |
| 7 | 50 | 175 | 126 | 1050 |
| 57, conditional | 3250 | 92625 | 89376 | 5187000 |

Indeed b=nd/2−(n−1). The numbers 6 and 126 are the free graph tiles. In the dense convention those tiles are **forward** and, including the spine, there are 15 and 175 forward arcs. To make all 15 or 175 graph edges backward, choose the fixed spine in the **complement graph** instead; then none of those graph edges has been reserved. The known Moore graphs have such complement spines, and the degree-57 complement would as well: its min degree exceeds n/2. This is a different choice of address.

## The plus and the fork

Let P also denote the undirected path adjacency and X the symmetric 0/1 matrix of chosen graph chords. Then A=P+X and the Moore condition is

    X²+PX+XP+X = (d−1)I+J−P²−P.

The left side makes a local picture on the upper triangular board. A chord cannot have span 2 or 3 (it would close a triangle or square in G); it has span at least 4. At every cell (i,j) of span at least 3, **exactly one** is supplied by:

- the tile itself or a one-cell shift of one endpoint of a chord along the spine;
- an unordered pair of chord edges meeting at a third vertex.

In particular each chord covers a five-cell plus (clipped at the first or last path vertex). These pluses are disjoint; e.g. two marked centers at Manhattan distance at most two cannot both occur. The uncovered cells are covered by the forks, neither twice nor by a plus. Span one already has its fixed graph edge; span two has its fixed midpoint, so no chord fork can land there. Conversely the prescribed chord degrees (d−1 at the endpoints, d−2 elsewhere), the gap rules, and this disjoint covering are sufficient by the matrix equation. With the degree profile the area is fixed, so avoiding collisions and keeping forks out of forbidden cells already forces coverage.

The plus area is 5b−2(d−1) (there is one missing arm for each incidence at a path endpoint). The fork area is 2 binom(d−1,2)+(n−2) binom(d−2,2). These sum to binom(n−2,2), the number of cells of span at least 3. In Petersen they give **26 plus cells and just two fork cells**. Its chord graph is exactly two 3-vertex paths (centered at the graph spine endpoints) plus two disjoint edges: the interior vertices have chord-degree one, and Petersen's endpoints cannot themselves be adjacent since it has no Hamiltonian cycle. In Hoffman–Singleton the ledger is **618 plus cells and 510 fork cells**; at degree 57 it would be 446,768 plus cells and 4,826,360 fork cells. It is a placement rule, not just an arithmetic check on the number of squares.

## Blue and black

The Moore graph's girth is undirected. The tournaments obtained from these tiles do have directed triangles. In fact a backward graph chord (i,j) in the sparse convention makes a directed triangle with i+2: i→i+2 is forward since that pair is not a graph edge, and i+2→j is forward since an i+2,j graph edge would give the 4-cycle through i+1 and i. The dense member has exactly (d−1)² more directed triangles, independent of the spine; see [the count and Hamilton-path bounds](moore-spine-chords.md).

No even-order direct graph spine can be blue. A fixed-free involution of a diameter-two girth-five graph would have to pair every vertex with an adjacent vertex: a nonadjacent pair has a unique common neighbor, fixed by the involution. But an edge between distinct adjacent swapped pairs makes a 4-cycle together with its image and the two pair-edges, contradicting connectedness outside a single pair. Reflection in an even spine is fixed-free. If its free chord pattern were invariant, the reflection would preserve the undirected graph (it already preserves the spine). Both members of every such directly marked native pair are therefore black. The 5-cycle's single closing chord is reflection-fixed and gives a blue pair. This obstruction is about the direct encoding, not the [previous compact blue constructions](../moore57_compact_blue_routes_20261124/).

For Petersen the direct obstruction goes further: the two normalized chord masks are mirror partners, and each tournament class has no blue incidence even after changing its Hamilton root. The replay enumerates 240 **oriented graph** Hamilton paths and gets just two masks, each 120 times, but each resulting sparse **tournament** has 5213 directed Hamilton paths. The distinction is which paths are legal words in the graph, and which are legal after completing every pair to a tournament.

In the unrestricted native metagraph both direct Petersen sparse classes are exactly **four moves** from transitive; their score-debt lower bound is 4 and the replay supplies literal directed reroots. They are different converse classes and exactly **two black moves** apart. Neither has any incident blue edge. By contrast the earlier compact blue Petersen uses a different rule, on the same number of vertices, and is seven moves from transitive. These are distinct size-10 nodes, not competing path lengths for one node.
