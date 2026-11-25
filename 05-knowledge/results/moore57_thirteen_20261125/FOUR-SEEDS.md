# Thirteen Hoffman–Singleton extensions and the exact cost of a pentagon defect

For two nonadjacent vertices of a diameter-two Moore graph, call their unique common neighbor their midpoint. The closure of a set repeatedly includes midpoints of its pairs. This is two-neighbor bootstrap percolation. Ibrahim et al. prove that in a hypothetical degree-57 Moore graph the minimum number of generators is 3 or 4; they give the equivalent extremal Petersen-count test. We use their theorem to count how often a fourth seed succeeds and to recognize a three-seed generator from the previously constructed pentagon coordinates. [Ibrahim et al., Propositions 4.15–4.16 and Theorem 4.18](https://ajc.maths.uq.edu.au/pdf/93/ajc_v93_p060.pdf)

## How many choices does a labelled triple leave?

Fix any independent triple lacking a common neighbor. Either it generates all 3,250 vertices or it generates an induced Petersen graph \(P\) on ten vertices, by the cited closure theorem. If it generates \(P\), write \(b(P)\) for the number of Hoffman–Singleton subgraphs containing \(P\).

**PROVED.** One has \(0\le b(P)\le13\). Exactly \(3240-40b(P)\) choices of a fourth vertex outside \(P\) generate the entire Moore graph. In particular there are at least **2,720** such choices. For each specified \(p\in P\), exactly \(54-4b(P)\) of its 54 exterior neighbors work. This number is always at least two. Choosing \(p\) outside the starting triple gives four independent generators.

Indeed, any induced Moore subgraph of diameter two is closed, since it already supplies the unique ambient common neighbor of each nonadjacent pair. Intersections of closed sets are closed. A closed set containing a Petersen graph can only be Petersen, Hoffman–Singleton, or the whole degree-57 graph. Consequently two distinct Hoffman–Singleton extensions intersect precisely in \(P\). At any fixed \(p\in P\), they consume disjoint sets of \(7-3=4\) among its \(57-3=54\) remaining neighbors. This yields the bound 13, the exact count of unused neighbors, and (since the extensions have disjoint 40-vertex exteriors) the count of all successful fourth choices. A site outside \(P\) adjacent to \(p\) cannot touch any other site of \(P\): an adjacent pair would give a triangle, and a nonadjacent pair already has its midpoint inside \(P\). The independent choice follows.

## Reading three-generator density directly from pentagons

For an induced pentagon \(C\), its five exterior neighborhoods of size 55 are linked by matchings in pentagram order. Let \(f_C\) count the fixed sheets of the permutation obtained after following all five matchings. This is also the maximum number of missing rows insertable into its five-column partial orthogonal array. [Pentagon matching and array completion theorem](../moore_mobius_prime_packets_20261124/PENTAGON.md)

**PROVED.** If \(T\) is the number of unordered three-vertex generating sets in the hypothetical graph, then
\[
                   3T=5\sum_{C}(55-f_C).
\]
In particular four generators are necessary **if and only if every pentagon has identity monodromy**, equivalently every one of its local five-column arrays admits completion by all 55 missing rows. One nonclosing sheet anywhere certifies the existence of a three-vertex generator. Array completion here is local to each pentagon.

A fixed sheet and the original pentagon give an induced Petersen: each exterior site has its matching spoke and its two pentagram neighbors, and all other adjacencies among those ten sites are excluded by short cycles. Conversely a Petersen through \(C\) must select exactly one sheet closing after five steps. Each Petersen contains twelve pentagons, so \(\sum_C f_C=12P\) for the number \(P\) of Petersen subgraphs. Each nonbacktracking three-edge path belongs to exactly one pentagon, giving \(C_5=3250\cdot57\cdot56^2/10=58{,}094{,}400\). The type-II triple census and its closure theorem give \(T=5{,}325{,}320{,}000-20P\), since each Petersen accounts for twenty such triples and every remaining triple of that type generates. No other kind of triple generates: one or two edges close inside a pentagon or a path, and a triple sharing a neighbor closes inside its star. Eliminating \(P\) gives the displayed identity. The [ten-site port test](TRIPLE-PORT.md) detects the same alternative from a single attempted Petersen closure.

This relates an integer count of successful initial states to a directly measurable placement loss. Applying the midpoint operation on a known adjacency relation reconstructs vertices. Specifying the unknown ambient midpoint operation from these initial labels alone remains a construction problem.

Combining the ten-site count with the pentagon count gives a direct conversion between three observable defects. If \(e\) counts induced Petersen-minus-one-edge sets, then
\[
               4e=3T=5\sum_C(55-f_C).
\]
Indeed the ten-site port formula gives \(e=15(266266000-P)\) while the triple count above gives \(T=20(266266000-P)\). Thus counting one missing Petersen edge locally and counting nonclosing sheets around all pentagons recover exactly the same global obstruction. [Ten-site multiplicity proof](TRIPLE-PORT.md)
