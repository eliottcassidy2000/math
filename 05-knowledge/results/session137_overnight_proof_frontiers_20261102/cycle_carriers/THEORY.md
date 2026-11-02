# Natural OCF carriers: sharp packing cost, exact small images, and link constraints

**Status:** PROVED consequences of cited Moon cycle minima and inherited SCC/OCF theorems; FINITE-EXACT image classification for components of at most 15 cycle-vertices. Full Hamiltonian-path spectrum completeness remains OPEN. No priority claim.

## Findings first

1. A connected natural carrier with packing number k has at least `f(k)=floor((3k-1)^2/4)` cycle-vertices. This bound is attained for every k. The exact extension to r components is the balanced sum of f. This is a source-realizability constraint, not an inequality for arbitrary graphs.
2. Every connected carrier with m cycle-vertices has a strong source of order at most `1+floor(sqrt(4m+1))`. This gives a finite, exact recognition algorithm, component by component. The bound is sharp infinitely often.
3. There are exactly **29 connected natural carriers of order at most 15**, listed below. Arbitrary disjoint unions of these give exactly the natural carriers whose components have at most 15 vertices. The reduction makes this an all-source-order statement, not a small-order extrapolation.
4. For an odd cycle C, its closed-antineighborhood is exactly the natural carrier of the tournament left after deleting C. In particular a packing-height-two carrier's disjointness graph cannot have degrees 3, 10, or 14. More generally `3 alpha_3 >= 2 N_3+2 N_10+N_14`, where N_d counts odd cycles disjoint from exactly d other odd cycles.
5. H=31 has minimum total cycle count 9, attained on seven tournament vertices. Its minimum tournament order is six, and at that order at least 13 odd cycles are necessary. These two resources cannot be minimized simultaneously.

## 0. Inheritance and scope

Fixed source repository pin: `bb7cb2d1f61a14632389814469317131d897b640` (Session136). No source-repository mutation, fetch, or publication by this worker.

* OCF is inherited from **THM-002**, `01-canon/theorems/THM-002-ocf.md`: `H(T)=I(Omega(T),2)`.
* **Session81**, `05-knowledge/results/session81_counterexample_families_20261028/CONFLICT_IMAGES.md`, already proves that natural-carrier components correspond exactly to nontrivial tournament SCCs; each component is K1, K2, or contains K4; forests are exactly disjoint unions of K1 and K2. These are not new here.
* The classical input is **Moon (1966), Corollary 2.1**: a strong n-tournament has at least `n-l+1` distinct directed l-cycles for every `3<=l<=n`. This follows alternatively from Moon's Theorem 2 on induced strong subtournaments and a Hamiltonian cycle in each. The primary PDF was read on 2026-11-02 at https://www.cambridge.org/core/services/aop-cambridge-core/content/view/5D7CCBE66AB07BB3F90F771D1CD3CD08/S0008439500052875a.pdf/on_subtournaments_of_a_tournament.pdf . The in-repo intake is `05-knowledge/reference/CORE-PAPERS-TOURNAMENT-STRONG-EARS-2026-08-25.md`.
* The consecutive-edge tournament used for sharpness is Moon's extremal construction in the same paper, and the inherited tribonacci partner of **THM-4731**, `01-canon/theorems/THM-4731-tiling-line-metagraph-self-lines-by-twisted-burnside-sigma-twins-depth-and-the-tribonacci-partner.md`, and `05-knowledge/results/defect_boundary_20261020/grammar.md`. Its Hamiltonian-count sequence is not new. The exact packing-cost extremal use and carrier-budget classification are the present consequences.
* **THM-1370**, `01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md`, proves the global omissions 7 and 21. We recover them through the bounded carrier table, not claim the omissions or their priority anew.
* Session136 `gluing/PROOF.md` supplies exact module contexts, the same-profile-hostile pair masks40/76, and `w(31)=6`. Session136 `ear_geometry/positive_cut_geometry.md` supplies the positive path-cut kernel. The relation to these is described in Section 6.

Cycle-vertices always mean **literal directed simple odd cycles of length at least three, modulo cyclic rotation only**. Different cycles on the same vertex support remain different vertices. Packing number means independence number of Omega, equivalently maximum number of disjoint directed triangles: a directed cycle's support contains a directed triangle.

Concept board: source SCCs / cycle budget / independent packings / natural links / universal cycles / context loss. Canonical hostiles: K3 and K10 are not natural carriers; K7-minus-one-edge has an attained scalar evaluation but is not natural; same-order K7 sources masks40/76 have different module-context responses.

## 1. A sharp finite source bound and exact image recognizer

Let S be strong of order n>=3 and let m be its number of odd cycles. Summing Moon's lower bounds over odd l gives

    m >= sum_(odd l=3..n) (n-l+1) = floor((n-1)^2/4).       (1)

Since an integer square is 0 or 1 modulo4, (1) is equivalent to

    n <= 1+floor(sqrt(4m+1)).                              (2)

The floor/square-root is integer arithmetic. This is not a bound on tournament order before deleting invisible singleton SCCs.

**Recognition theorem.** A finite simple graph G is a natural Omega carrier iff, for every connected component F of order m, some strong tournament S of order

    3 <= |S| <= 1+floor(sqrt(4m+1))

has Omega(S) isomorphic to F. Necessity uses the inherited SCC/component theorem and (2). For sufficiency, choose such a source for each component and put the sources in dominance order. No directed cycle crosses these sources, so Omega is their disjoint union. For empty G use a one-vertex transitive tournament.

Thus finite enumeration of strong sources, literal cycles, and exact graph isomorphism is a terminating necessary-and-sufficient recognition algorithm. No polynomial-time claim is made. It considers at most `2^(2m+O(sqrt(m)))` labelled source tournaments for one m-vertex component, before cycle enumeration and graph isomorphism. A source with all singleton SCCs removed has at most 3m vertices in total, since (2)<=3m for m>=1. The global 3m bound is sharp: an edgeless m-vertex carrier requires m triangle SCCs.

For sharpness of (2), let M_n have vertices0,...,n-1 with forward consecutive arcs `i->i+1`, and all other arcs backward. It is strong. Starting a directed cycle at its least vertex a forces successive vertices a+1,a+2,...,b: from the current vertex the only possible unvisited larger out-neighbor is its successor, and an already visited smaller vertex cannot be revisited except to close at a. Thus every directed cycle is exactly a consecutive interval of length at least3, traversed in increasing order and closed by its backward long edge. Consequently M_n has exactly `floor((n-1)^2/4)` odd cycles. The source-order bound is attained for every n>=3 at these budgets.

## 2. Sharp cycle cost of independent packings

If Omega(S) is connected and has independence number k>=1, S contains k disjoint odd cycles, so n>=3k. Equation (1) gives

    m >= f(k),  f(k)=floor((3k-1)^2/4).                    (3)

This is sharp for every k: M_(3k) contains the k disjoint consecutive triangles, has no packing larger than k by order, and has exactly f(k) odd cycles. In particular the minimum connected cycle budgets for packing1,2,3,4,5 are

    1, 6, 16, 30, 49.

A connected carrier with at most 15 vertices necessarily has packing height at most 2. The abstract graph condition k<=m misses this quadratic cost.

More generally fix total packing number k and number r of nonempty connected components, with 1<=r<=k. Write k=rq+s, 0<=s<r. Then the exact minimum total cycle count is

    (r-s)f(q)+s f(q+1).                                 (4)

Indeed packing numbers k_i add over components and each component costs at least f(k_i). The successive increments of f are strictly increasing: `f(2t)=9t^2-3t`, `f(2t+1)=9t^2+6t+1`, so balancing two arguments differing by at least 2 reduces their sum. The minimum occurs at q and q+1. Dominance sums of r copies of the appropriate M_(3k_i) attain it. For r=k this reduces to m=k, realized by disjoint triangle SCCs; for r=1 it is (3).

The complete polynomial for the sharp source is also explicit. If P_n(x)=I(Omega(M_n),x), interval-packings give

    sum_(n>=0) P_n(x) z^n = (1-z^2)/(1-z-z^2+(1-x)z^3),  (5)
    [x^j]P_n = [z^(n-3j)] (1-z)^(-j-1)(1-z^2)^(-j).      (6)

A packing decomposes its line into uncovered single vertices and selected odd intervals; (5) is the sequence generating function. Equivalently distribute total extra even interval lengths and j+1 gaps to get (6). In particular the maximum packing of M_(3k) is unique. At x=2 this recovers the **inherited** tribonacci H sequence 1,1,1,3,5,9,17,31,57,105,...; it is not a full-spectrum constructor.

## 3. Universal cycles and source-preserving links

For a cycle C inside a strong source S,

    C is universal in Omega(S) iff S-V(C) is transitive.  (7)

A nontransitive complement contains a directed triangle disjoint from C; conversely any odd cycle disjoint from C makes the complement nontransitive. Thus universal cycle-vertices are precisely odd decycling cycles.

Cycles whose complements have at most 2 vertices are universal. Moon therefore gives at least 1 universal vertex at n=3, at least 2 for even n>=4 (the (n-1)-cycles), and at least4 for odd n>=5 (the n-cycle plus at least three (n-2)-cycles). Every nontrivial connected carrier has at least two universal vertices. These bounds are conservative: M_n has four universal vertices for even n>=6 and five for odd n>=7, but no all-order stronger minimum is claimed here.

For every odd cycle C in any T, there is a canonical literal-cycle identity

    Omega(T)-N[C] = Omega(T-V(C)),                       (8)

where N[C] is its closed graph neighborhood. For an independent packing P, deleting the union of its closed neighborhoods likewise gives Omega of the remaining tournament support. Thus natural carriers are closed under disjoint union and these antineighborhood operations. They are not closed under arbitrary vertex deletion or universal-vertex addition: deleting one cycle-vertex from realized K4 gives unrealizable K3; adjoining a universal vertex to realized K2 also gives K3.

Let D be the complement of Omega, so its edges mean disjoint cycles. For each cycle C, D[N_D(C)] is the complement of the natural carrier in (8). Its edges count packings of two cycles disjoint from C. Consequently

    sum_C e(D[N_D(C)]) = 3 alpha_3(Omega).               (9)

The finite classification below implies that a natural carrier with 3, 10, or 14 vertices has at least 2,2,1 independent pairs, respectively. For disconnected carriers a split of m vertices already creates at least m-1 independent cross-pairs, so the stated minima follow from the connected table (and the matching classification for m=3). Writing N_d for the number of C with degree_D(C)=d gives

    3 alpha_3 >= 2N_3+2N_10+N_14.                       (10)

In particular a tournament with packing number at most 2 has no odd cycle disjoint from exactly3,10,or14 odd cycles. Equivalently its disjointness graph has no such degrees. The degree 3 case uses only the inherited K3/forest obstruction; degree 10/14 use the present complete budget classification. No prohibition of induced graph shapes in an arbitrary subset of cycles follows.

## 4. Complete connected carrier list through 15

Write `K_m \ F` for deletion of the edges of F from K_m, with all additional vertices universal. The complete list is:

| m | F (deleted-edge graphs, isolates omitted) |
|---|---|
|1|empty|
|2|empty|
|3|none|
|4|empty|
|5|empty|
|6|empty; K2|
|7|empty|
|8|empty|
|9|empty; K2; P4|
|10|2K2; P3|
|11|empty; K2|
|12|empty; K2; 2K2; A|
|13|empty; K2; P3|
|14|K2; 2K2; P4; 4K2; B|
|15|empty; C|

Here the three exceptional deleted-edge graphs have these explicit edges:

    A on0..7: 07,17,26,36,45,47,56,67;
    B on0..5: 05,15,24,25,34,35;
    C on0..5: 05,15,25,34,45.

All unspecified vertices of K_m are isolates in F. A has 8 edges, B has 6, C has 5. The graph in each row has independence polynomial `1+mx+e(F)x^2`; no larger packing occurs by (3). There are29 graph types in total.

**Why exhaustive over all tournament orders.** For m<=15, (2) gives n<=8. The verifier builds complete isomorphism-class banks by augmentation, with counts1,1,2,4,12,56,456,6880 at orders 1,...,8. Every (n+1)-tournament has an n-vertex deletion, so augmentation completeness is inductive. Exact canonicalization sorts outdegrees and enumerates all permutations inside equal-degree classes; any isomorphism preserves those degrees. For each strong representative it lists literal directed cycles, builds the overlap graph, and canonicalizes the nonisolated part of its complement. The checks cover all 6,404 strong classes through 8. All candidates with m<=15 appear in the table, and each type retains an explicit source adjacency witness and literal cycle list. The known class counts are checks, not the proof of augmentation coverage.

No scalar polynomial or abstract graph is treated as automatically realizable. Each positive entry has an actual tournament, every absent entry inherits the rigorous source-order bound, and arbitrary unions have the explicit dominance-source construction.

**Universal/link tests are insufficient.** `K7\K2` contains K4, has five universal vertices, and every nonempty closed-antineighborhood is K1 (a natural carrier). Yet the complete table excludes it. Its OCF evaluation is19, which is attained by realized K9. Thus all those structural tests, and even an attained scalar evaluation, fail to certify natural realizability. Likewise adding a universal vertex to realized `K6\K2` gives this unrealizable graph: the outputs17 and19 are both attainable, but that graph operation is not source-preserving.

## 5. Exact cycle-budget compiler and the 31 tradeoff

Every natural carrier is the disjoint union of the listed component types when its total number of cycle-vertices is at most 15. Therefore multiplying their polynomials with additive cost m is a complete bounded-budget constructor, not a depth-limited heuristic. A component expression compiles to a dominance sum of its stored tournament sources. The verifier finds 300 distinct attained H-values at total cycle cost at most 15, and reports every attained H at each cost, and a source expression for each minimum-cost value. Distinct carrier graphs or polynomials can have the same H; the compiler does not identify them as source-equivalent.

Define q(h)=min_T |Omega(T)| over tournaments with H(T)=h, if any. Whenever the compiler first attains h at cost<=15, it has found q(h) globally: any smaller-cost witness would also be in the classified universe. Selected exact results are

    h:    1  3  5  9 11 13 15 17 19 23 25 27 29 31 35 37 41 53
    q(h): 0  1  2  2  5  6  3  6  9  9  4  3 10  9 13 14 14 14.

In particular, 35 cannot be obtained via a factor 7; its minimum-budget constructor is the actual connected carrier `K13\P3`.

In particular, 31 is prime, so every realization has one nontrivial SCC carrying31, with any other SCCs singleton. The table shows the least possible cycle count is9, achieved by M7 with carrier `K9\P4`. The complete n<=5 bank has no31, while n=6 does; the order 6 bank has minimum odd-cycle count13 among its31 witnesses, attained by `K13\K2`. Hence minimum order 6 and minimum cycle budget 9 cannot be attained simultaneously. This complements Session136's `w(31)=6`; it does not identify substitution width with either other resource.

Finally a hypothetical H=7 or21 has total cycle count at most(H-1)/2 by OCF, hence at most 10. The bounded compiler omits both, recovering the inherited all-order omissions without importing Moon-Busch's strong Hamiltonian-path minimum. This is an alternative finite carrier certificate, not a new omission or spectrum-completeness result.

## 6. Genuine connections and destroyed information

**Module contexts.** Omega is not a substitution state, even when source order is retained. Inherited five-vertex pair-bit masks40 and76 both have Omega=K7, but their profiles are `(15,39,33,10,1)` and `(15,45,35,10,1)`. The actual homogeneous context C3[T,1,1] returns123 and135, with carrier polynomials `1+45x+8x^2` and `1+47x+10x^2`, independently computed here. Source->carrier preserves H and all odd-cycle packings but destroys path-cover data. Session136's profile, not the conflict graph alone, is the needed context sidecar. Even the empty carriers of all transitive L_n separate under that context: the outputs are2^n+1.

**Positive ears.** In the Session136 ear formula, deleting the apex x from an odd cycle through x gives its even-cardinality path P. The contribution weight `H(T-V(P))` is exactly the activity2 evaluation of that cycle's natural link (8) in the child. Thus the positive path-cut kernel is a pointed aggregation of natural-link values. Forgetting which cycle goes through x and which ordered endpoints remain destroys the grouping needed for the cut kernel. The link theorem is compatible with that kernel, but does not reconstruct it from an unpointed Omega graph.

**Remaining frontier.** The finite recognition theorem is exact but not a practical uniform classification. The 29-type table does not extend beyond component budget 15. The stronger universal-count minimum at arbitrary order is unproved here. No unrestricted graph-join constructor, bounded packing-height full-spectrum theorem, or all-odd-except7/21 realization theorem is established.
