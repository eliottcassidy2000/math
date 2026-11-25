# Nested Moore cores and a rigid extension of the fourteen-site chart

25 November 2026. Continues the [one-square Latin obstruction](LOCAL-LATIN.md), [the three-star construction](THREE-STARS.md), [the actual fourteen neighborhoods](PLANTED.md), and [the separate 400-coclique trade star](../moore_midpoint_20261125/PROJECTOR-CODE-AND-TRADE-STARS.md). The partial construction survives. The new constraints say which placements it cannot have in a full completion; they do not decide existence.

## 1. Petersen is a pair of matched pentagons, and more than that

In the affine Hoffman–Singleton graph, let \(A_i\) be a pentagon in cyclic difference order 1 and \(B_j\) a pentagon in cyclic difference order 2. The contact \(y=x+ij\) is a perfect matching. Shift its \(B_j\) coordinate to make the matching the identity: the two five-cycles with different steps and these five matching edges are the generalized Petersen graph \(G(5,2)\), i.e. Petersen. This gives 25 immediate copies. The full affine replay enumerates the 1260 pentagons and every pair whose cross edges are exactly a matching, giving **525 distinct Petersen supports** (each with six decompositions and each pentagon in five copies). This is consistent with the classical [pentagon-edge completion in Hoffman–Singleton](https://www.degruyterbrill.com/document/doi/10.1515/advgeom-2025-0027/pdf). A general degree-57 graph is not known to contain even one Petersen copy; [Godsil's smaller-core question](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v2i1f1/pdf/) is a different problem from the guaranteed [Sylvester-like substructure](https://doi.org/10.1016/j.ejc.2018.02.018).

There is an intrinsic shell calculation for any induced degree-\(h\) Moore graph inside a degree-\(k\) one. Call its \(m=h^2+1\) vertices the **core**. Every vertex outside touches at most one core vertex: any two core vertices already have their unique route of length at most two inside. Put the \(k-h\) neighbours outside a core vertex \(i\) into its port class \(P_i\). The remaining vertices form \(R\). Two port classes have no edges if their core vertices are adjacent; otherwise they have a perfect matching, forced by the unique route from a port to the other core vertex. Every vertex of \(R\) has one port in each color, forced by its routes to all core vertices. Counting gives
\[
 |R|=(k-h)(k+h-h^2-1).
\]
A port has one core neighbour, \(h(h-1)\) port neighbours, and \(k-1-h(h-1)\) residual neighbours. A residual vertex has \(m\) port neighbours and \(k-m\) residual neighbours, if \(R\) is nonempty.

For Petersen inside Hoffman–Singleton, \((h,k)=(3,7)\): **all 40 remaining vertices are ports** in ten groups of four; there is no residual. This is an exact, intrinsic form of the embedding independent of the affine picture. For Hoffman–Singleton inside degree 57, \((h,k)=(7,57)\) instead gives \(50+2500+700\) with quotient
\[
 \begin{pmatrix}7&50&0\\1&42&14\\0&50&7\end{pmatrix}.
\]
The residual graph \(R\) is 7-regular of girth at least five. Its sites have 50-coordinate signatures, each selecting one of 50 symbols in the corresponding port class; every symbol in each color appears exactly 14 times. Two signatures agree in zero positions when their sites already have an edge or a common neighbour in \(R\), and in exactly one position otherwise. For each pair of colors the port-only routes account for 36 slots per initial symbol (36 two-step matchings if the core pair is adjacent; 35 and the direct matching if not). The remaining **14 slots** must be supplied by the same 700 signatures. This 14 is a parameter of a hypothetical embedded core, not the fourteen sites of the partial fourth branch.

These signatures include a useful extreme-case check. If \(D\) is the residual adjacency and \(M\) its residual-to-port matrix, then \(MM^{\mathsf T}=56I+J-D-D^2\). The residual kernel dimension is the number of components minus one, since \(56-\lambda-\lambda^2=(7-\lambda)(\lambda+8)\) for the nonconstant modes. There are at most 14 components by the radius-two bound of 50 for each 7-regular component of girth at least five. At equality all 14 are Hoffman–Singleton. Since a port cannot touch two vertices of the same such component, each of the 2500 ports then sees one vertex in each of **15 pairwise nonadjacent copies** (the original and the 14). For any two of those copies each vertex pair has a unique common port. The port data become an exact orthogonal array \(OA(2500,15,50,2)\), or 13 mutually orthogonal Latin squares of order 50, with the remaining port adjacency still to be filled. This is necessary for that extreme case only.

## 2. Connecting the coclique swaps to the nested chart

The coclique and rooted Latin pictures enforce the same kind of coverage: in the former, a pointed parallel class covers each coclique point once, and in the latter a direct matching and all two-step transports cover each destination slot once. The marked branch labels and the marked block labels matter in both pictures; an aggregate modular checksum does not ensure the needed common placements. They are different parametrizations, not an equivalence between every Latin chart and a projective completion.

There is a numerical invariant if both pictures coexist. In Hoffman–Singleton every 15-coclique meets **every** Petersen core in exactly three vertices. If it contains \(a\) core sites and \(b\) ports, summing its three contacts over the Petersen core gives \(3a+b=3(10-a)\); since \(a+b=15\), \(a=3\). In degree 57, if a 400-coclique and a Hoffman–Singleton core coexist, write its counts as \((a,b,c)\) in the three layers above. Summing its eight contacts over the core gives
\[
     7a+b=8(50-a), \qquad (a,b,c)=(a,400-15a,14a).
\]
The same relation applies to all 400 leaves of a projective trade star, not just to its center. An implementation simultaneously tracking such a core and a coclique therefore needs these same physical counts before worrying about calibrating the swaps. In the extreme 15-core array case, the 400-coclique must meet all 15 cores in the same number: its port count \(b\) is common to all 15 equations \(7a_i+b=8(50-a_i)\).

## 3. A single saturated Latin board detects a nested six-slot copy

Put the prescribed root inside a hypothetical Hoffman–Singleton core. Its seven internal neighbours are branch heads, and six sheets per such head belong to the core. Every between-head matching restricts to a bijection of these six slots. In each of the other 50 heads, the branch from a port of the root splits its 56 slots into \(42\) ports of nonadjacent core vertices and \(14\) residual vertices. Thus the core imposes a **different cut** on the existing physical branch files, not just a formal containment of order-6 Latin labels inside order 56.

For two heads in the core, normalize their direct matching as identity. Returning from a third head to the first along the second and then the direct inverse gives a triangle transport. On the six slots of a Hoffman–Singleton core, **every** triangle transport is three disjoint transpositions: the replay checks all \(50\binom73=1750\) rooted triples in the actual affine graph (also the involution phenomenon in [Smith–Montemanni](https://www.mdpi.com/2227-7390/14/23/4318)). Thus an exterior head is eligible as a third core head only if its physical triangle permutation on 56 slots has at least three transpositions. Ultimately five distinct additional heads must share *one invariant six-set* on the same starting fiber. Equivalently, the single saturated two-head board needs a six-row, six-column subsquare through the direct row, with the five other rows of type \(2^3\) on those columns. This is stronger than and independent of merely checking modular rank.

On the patched [fourteen-site placement files](PLANTED.md), a four-head selection containing all three old completed heads has no invariant six-set: add any of the 54 exterior heads, and the smallest orbit of its three normalized triangle maps is at least 14. More strongly, even allowing a core that uses just two completed heads, this is the exact candidate list for the other five:

| Completed pair | Extra heads with at least three transpositions | Common six-set from five of them? |
|---|---|---|
| \(56,0\) | 14, 19, 20, 29, 37, 45 | No: all six five-head choices act transitively on 56 sheets |
| \(56,28\) | 3, 10, 19, 48 | No: fewer than five |
| \(0,28\) | 5, 6, 10, 14, 20, 45 | No: all six five-head choices act transitively |

The obstruction uses only maps already fixed by the three stars. **In any completion keeping them, a Hoffman–Singleton copy containing the chart root uses at most one of these three branches** as a core branch. It does not exclude a copy based elsewhere or using different branches.

## 4. Remove the last symmetry without moving any edges

The exact 13,231-edge fourteen-neighborhood fixture has **one remaining involution** even after the cyclic symmetry is broken. Ordinary color refinement run to stability gives 56 singletons and 1597 pairs. Every automorphism is encoded by one flip bit per two-cell. Between two paired cells the bipartite edges are empty, complete, or a matching, by equitability. Only the matching constrains the flags: their flip bits must agree, regardless of whether the matching is parallel or crossed. Singleton incidences and within-cell edges are unchanged. The graph connecting matched two-cells is connected, so its only nontrivial automorphism swaps every pair. It fixes the root and exchanges branch heads 27 and 49.

The single extra contact **\((1,1)\)--\((2,7)\)** is legal: the old graph has no route of length at most three between its endpoints. With the contact added, color refinement is discrete, so the partial graph is rigid even unmarked. No old edge has to move. A larger verified patch retaining that contact seats **51 of 53** new neighbours at \((1,1)\); only branches **13 and 33** remain. It retains every original edge and all fourteen complete neighborhoods, has 13,282 edges, gives the new source degree **55**, and is again completely rigid. The verifier inserts each edge against the current graph and checks for old routes of length at most three, so it detects short cycles that would use two or more new edges as well.

An exact repair search shows that inserting the last two contacts **while keeping the marker** requires rearranging at least fourteen of the new placement rows of this 51-contact reference. This is a local lower bound, not a proof that no 15th neighborhood exists. Candidate tuples avoid target slots already occupied by the original fourteen sources, and must cover all three central slot complements with distinct coordinates. Starting with the two holes, try a candidate in an open row, mark it chosen, and open each still-held reference row whose coordinates it displaces. Charge each discarded reference once. This branching is exhaustive within the budget: follow a hypothetical solution's choices for whichever row is open; any conflicting held row must have changed in that solution. A solution modifying extra disconnected rows only makes the budget larger. The replay exhausts budget 13, with 24,571,716 nodes at that budget. The rest of the fourth star is not completed.

Run in this directory:
~~~sh
python3 nested_core_check.py
c++ -O3 -std=c++17 rigid_patch_cover.cpp -o /tmp/rigid_patch_cover
python3 nested_core_check.py --board | /tmp/rigid_patch_cover 13
~~~
The final command is expected to exit 1 because every tested repair budget fails. The script leaves intact the distinct question of simultaneously seating all remaining sources.
