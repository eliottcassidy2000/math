# Simple-cycle census and length extension in the degree-57 Moore problem

2026-11-28. All statements about a complete degree-57 Moore graph are conditional. The known diameter-two Moore cases have degrees 2, 3 and 7; 57 is the fourth spectral possibility, not a known graph. See [the cyclic-derangement paper](https://www.mdpi.com/2075-1680/15/5/332). Cycle counts below are for unoriented simple cycles, not closed walks, and do not require the cycles to be induced.

## Exact small lengths and certified bounds at longer lengths

| Length | Exact count or approximate lower bound | Approximate upper bound if not exact |
|---:|---:|---:|
| 5 | 58,094,400 | exact |
| 6 | 2,662,660,000 | exact |
| 7 | 123,243,120,000 | exact |
| 8 | 6,042,906,870,000 | exact |
| 9 | 300,896,555,960,000 | exact |
| 10 | 15,140,589,832,368,000 | exact |
| 11 | 769,345,344,658,800,000 | exact |
| 12 | 3.94045 × 10^19 | 3.95941 × 10^19 |
| 21 | 1.17860 × 10^35 | 1.22500 × 10^35 |
| 63 | 5.16049 × 10^107 | 1.08341 × 10^108 |
| 81 | 2.51113 × 10^138 | 2.47210 × 10^139 |
| 100 | 1.03910 × 10^169 | 3.28952 × 10^172 |
| 170 | 6.80057 × 10^280 | 4.56913 × 10^294 |
| 500 | 9.90610 × 10^783 | 1.24170 × 10^871 |
| 1430 | 2.09950 × 10^1736 | 9.84224 × 10^2491 |
| 3250 | no positive lower bound here | 4.83476 × 10^4367 |

The exact integer program covers all lengths 5 through 3250. This gives a positive lower bound at **every integer length 5 through 1430**. An unweighted graph has integer edge-count lengths. If every edge is instead declared to have the same positive rational length u, the entry for k gives the count at length ku. Arbitrary edge-specific rational weights would be a different problem.

## 1. Nonbacktracking trace and exact collision mass

Write n=d²+1, q=d−1, and let B be the oriented-edge matrix whose transition is e→f when the head of e is the tail of f and f is not the reverse of e. [Ihara–Bass](https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.FSTTCS.2017.46) gives

$$
\det(I-zB)=(1-z^2)^{n(d-2)/2}\prod_{\lambda\in\operatorname{spec}A}(1-\lambda z+q z^2).
$$

Thus the contribution of an adjacency eigenvalue λ to tr(B^k) is the power sum for the two roots of X²−λX+q, computed by the recurrence p_0=2, p_1=λ, p_k=λp_{k−1}−qp_{k−2}. The 57 case uses eigenvalues 57,7,−8 with multiplicities 1,1729,1520. The trace counts closed words with no immediate reversal even at the wraparound; each simple k-cycle contributes 2k words. A repeat splits a word into two nonempty linearly nonbacktracking returns, each of length at least five. So through length 9 the trace divided by 2k is exact.

For longer lengths we can compute not just the trace but the sum of the **number of equal-position pairs** over all such words. Fix a vertex v and neighbors x,y. Let H_l(x,y) count l-edge nonbacktracking returns at v with first edge v→x and last edge y→v, without imposing a wraparound condition. There are just two channels, h_l on x=y and g_l on x≠y, even without any transitivity hypothesis. Explicitly, put Q(z)=(I−zA+qz²I)^−1; then

$$
H_l(x,y)=[z^{l-2}]\frac{Q_{xy}-2zQ_{vx}+z^2Q_{vv}-\mathbf1_{x=y}}{1-z^2}.
$$

For a derivation let R reverse oriented edges, S mark their heads (edge-by-vertex), and T mark their tails (vertex-by-edge). Then B=ST−R, TS=A, TRS=dI. Apply the matrix inverse lemma to I+zR−zST, using (I+zR)^−1=(I−zR)/(1−z²). Its (v→x,y→v) entry equals z times the displayed generating function before coefficient extraction. The adjacency equation reduces every Q coefficient to a linear combination of I,A,J. In particular its (x,y) entry is uniform on x≠y in the neighborhood of v, since distinct neighbors are nonadjacent.

Write L_l=d(h_l+qg_l) for the total linearly nonbacktracking returns at v. For a specified pair of equal positions with gaps a and b, the number of cyclically nonbacktracking concatenations per v is

$$
D(a,b)=(1-2/d)L_a L_b+d(h_ah_b+qg_ag_b).
$$

Without either junction rule there are L_a L_b choices. Either forbidden immediate reversal costs L_a L_b/d, by uniformity of each first/last neighbor marginal. Their intersection counts a pair with reversed first/last neighbor data, giving the term added back. Consequently the exact total collision mass is

$$
\mathcal D_k=\frac{nk}{2}\sum_{a=5}^{k-5}D(a,k-a).
$$

The factor k/2 sums over cyclic positions and counts each unordered equal-position pair once. Every nonsimple word has at least one repeat pair. Also each vertex appears at most p=floor(k/5) times, by spacing of its occurrences. Thus the maximum possible number M of equal-position pairs is bounded by floor(k/p)·binom(p,2)+binom(k mod p,2): pack multiplicities up to p to maximize the convex sum. For traces W_k this gives exact-integer bounds

$$
\left\lceil\frac{W_k-\mathcal D_k}{2k}\right\rceil
\le C_k\le
\left\lfloor\frac{W_k-\lceil\mathcal D_k/M\rceil}{2k}\right\rfloor .
$$

Use zero for a negative lower bound. This method alone gives the large improvement at length 63; its lower bound eventually goes to zero, so we use the constructive extension below.

## 2. Exact corrections at 10 and 11

Two distinct pentagons meet either in a path on 1, 2, or 3 vertices, or not at all. For any two common vertices the shorter connecting path (of length at most two) must be identical in both cycles: another path of length at most two would violate the girth. This convexity forces the listed intersections. Let s_i count unordered pairs sharing exactly an i-vertex path. There are q² pentagons through every edge, q through every two-edge path, and d q²/2 through every vertex. Therefore

$$
s_3=\frac{ndq}{2}\binom q2,\qquad
s_2+2s_3=\frac{nd}{2}\binom{q^2}2,\qquad
s_1+2s_2+3s_3=n\binom{dq^2/2}2.
$$

A bad length-10 word splits at an equal vertex into two pentagons. If they share no edge, the two directions of each loop can be chosen independently. If they share an edge, they must follow the shared path in the same direction to avoid a reversal at the junctions. All common vertices then give splits, and for an i-vertex overlap the proportion of starting positions admitting a split at offset 5 is 2i/10. Counting weighted marked splits gives 40s_1+20(s_2+s_3) bad words from distinct pentagons. Going twice around the same pentagon contributes 10C_5 more.

For length 11 the two pieces are a pentagon and a hexagon. Both are induced. A pair of vertices in different components of their intersection must be opposite in the hexagon: otherwise their identical shorter path forces them into the same component. Thus the intersection is a path on 1 through 4 vertices, or exactly two isolated vertices opposite on the hexagon. Write t_i for mixed pairs with i-vertex common path, and u for mixed pairs with two isolated common vertices. Between the endpoints of any simple three-edge path, there are q three-edge paths in total; they are internally disjoint since any shared interior vertex forces the other interior vertex by uniqueness of the length-two route. Thus a specified three-edge path occurs in one pentagon and q−1 hexagons. A specified two-edge path occurs in q pentagons and q(q−1) hexagons. Counting by paths and vertices gives

$$
\begin{aligned}
t_4&=(ndq^2/2)(q-1),\\
t_3+2t_4&=(ndq/2)q^2(q-1),\\
t_2+2t_3+3t_4&=(nd/2)q^4(q-1),\\
u+t_4&=(ndq/2)q\binom q2,\\
t_1+2t_2+3t_3+4t_4+2u&=n(dq^2/2)^2(q-1).
\end{aligned}
$$

The fourth equation marks a common vertex pair opposite on the hexagon, with q choices of pentagon and binom(q,2) hexagons for each nonadjacent vertex pair. For every valid concatenation, each shared vertex supplies one start with the five-edge segment first. For the disconnected intersection this is because the other shared vertex is at offset 2 or 3 on the pentagon and offset 3 on the hexagon. For a shared path, the directions must agree as above; the two directed offsets along that path then differ by zero or one. Counting such marked starts with weight 11/i per split when the intersection has i vertices gives 44(t_1+u)+22(t_2+t_3+t_4). This weighting does not assume that all splits decompose into the same two geometric cycles. Subtract these corrections from the traces before dividing by 20 or 22. The replay checks the two-channel formula and exact small census independently on Petersen and Hoffman–Singleton.

## 3. An injection-with-fibers that extends every length through 1430

Take a simple k-cycle X and one of its k two-edge windows u–v–w. The leaf contacts between N(u)−{v} and N(w)−{v} are a perfect matching of size q, by the [pentagon-board argument](../finite_zeta_and_moore_pentagon_boards_20261128). Each match is an alternate u-to-w path of length three. At most d_X(u)−1+d_X(w)−1 of the matches use a blocked vertex already on X, where d_X denotes degree in the graph induced on its vertex set. A remaining match replaces the window and constructs a simple (k+1)-cycle. Summing over windows gives at least k(q+2)−2·sum_{x in X} d_X(x) distinct marked expansions. Given a resulting cycle, there are at most k+1 possible three-edge windows to reverse, and for each the original two-edge path is unique.

The upper centered eigenvalue is 7, not 8: decomposing the subset indicator along the constant line gives mean induced degree at most 7+(57−7)k/3250=7+k/65. Hence

$$
C_{k+1}\ \ge\ \left\lceil C_k\frac{k(44-2k/65)}{k+1}\right\rceil.
$$

The code rounds the induced degree sum down to an even integer first. The factor is positive for k<1430. The procedure proves the contiguous range with no Hamiltonicity or long-cycle assumption.

For a complementary upper bound at large k, every directed Hamilton ordering of a k-subset gives a permutation supported on that induced adjacency; each unoriented cycle gives two permutations. [Brégman's bound](https://doi.org/10.1006/JCTA.1996.2727) bounds that permanent by the product of (j!)^(1/j) for its row sums j (zero rows give zero). The code uses integer bisection for rigorous rational upper approximations to these values, and tests that each affine majorant dominates them at every integer j=1,...,57. AM–GM and the same induced degree bound give an upper bound after multiplying by binom(3250,k)/2. Take the minimum with the trace. This matters especially near spanning lengths, without asserting that spanning cycles exist.

## 4. Fibonacci–triangular numbers and the two spectral discriminants

In the [previous branch-ring transfer](../moore57_fibonacci_triangular_branch_cycles_20261128) the restriction to nonadjacent edit positions introduced the Fibonacci recurrence, and the empty-or-two-edit schedules numbered T_{m−2}. [Luo Ming](https://www.fq.math.ca/Scanned/27-2/ming.pdf) proves that the positive triangular Fibonacci numbers are exactly 1,3,21,55. Here their occurrence is structural in two different statistics, not one degree-selection law:

| Degree d | Unordered pairs of incident edges at a vertex | Width d−2 of an off-diagonal row in the two-root grid |
|---:|---:|---:|
| 2 | **1** | 0 |
| 3 | **3** | **1** |
| 7 | **21** | 5 |
| 57 | 1596 | **55** |

The first statistic is binom(d,2). Asking for 55 in the *same* column would require degree 11, not 57. In degree 57 the 55 also counts the forward nonbacktracking choices at the end of a fixed three-edge path *other than* the unique next step that closes a pentagon with the starting vertex. It is exactly the centered transfer magnitude q−1 from the branch-ring calculation.

For a precise two-root interpretation, an all-detour ring of m branches has m hubs at every third position. Fix two candidate roots r,s. Adjacent roots cannot share such a cycle: consecutive hubs from the two phases and the edge rs would make a square. For nonadjacent r,s let t be their common neighbor. Pair a_i in N(r)−{t} with b_i in N(s)−{t} by the unique cross matching. For i≠j let w_ij be the common neighbor of b_i and a_j. It cannot be either root or in either root neighborhood; two distinct cells cannot coincide since that would give two common neighbors for two sites in one root neighborhood. The unfilled positions of the graph are exactly N(t)−{r,s}, also of size d−2: these sites cannot be grid cells by uniqueness of the common-neighbor routes to r or s, and vertex counts agree. The common 3m-cycles are precisely concatenations of a_i,b_i,w_ij,a_j around a cyclic order of m different labels. The common hub t is excluded because a consecutive s-hub would form a triangle with s. The pair of roots chooses the orientation r-hub then s-hub; divide by m rotations, not by 2m. Thus every nonadjacent root pair shares (d−1)_m/m cycles. The repo's one-root sum minus these pair overlaps gives a corrected all-detour-only bound above 2.9748×10^73 for length 63. This isolates the structured family even though the collision method counts many more cycles overall. Also, a permutation of 56 sheets **cannot have exactly 55 fixed sheets**: moving one forces a second to move. This removes the hypothetical f=55 insertion example in the preceding branch-ring note.

The adjacency eigenvalue labels 57,7,−8 collide on reduction exactly at 2,3,5,13, all [Ogg primes](https://uva.theopenscholar.com/files/ken-ono/files/1411.5354.pdf). Nonbacktracking introduces a different set of discriminants: λ²−4q for these same three labels gives 55², −175 and −160. Their prime support is 2,5,7,11, also all Ogg primes. These are collisions *within* the respective quadratic factors, not an assertion that the full nonbacktracking spectrum or the 15-prime list is explained by those factors. In particular the Boolean mod-2 certificate from the [shared-edge session](../moore57_shared_edge_boards_20261128) is exact because of the degree budget; recovering integer cycle counts uses the stronger integer trace and collision calculations.

## Replay

Run python3 count.py for the rounded summary and control enumeration (the Hoffman–Singleton census uses a rooted exhaustive search and may take a minute). Run python3 count.py bounds.csv --bounds-only to generate the exact integer table without the controls. Run python3 overlap.py for a direct Hoffman–Singleton check of the two-root overlap at length 9 and the degree-57 arithmetic.
