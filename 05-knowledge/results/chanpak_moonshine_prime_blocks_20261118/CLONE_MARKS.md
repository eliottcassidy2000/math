# Clone tests for a switching frame: one mark, many marks, and shared growth

The [same-lattice switches](FRAME_REPAIR.md) preserve the two-slot height while changing the Hamiltonian count. The [ordinal-gap test](GAP_INVERSE.md) only checked one special placement. Here the marks can be anywhere in a fixed quotient, even on the switched side. This gives exact algebraic criteria for retaining the height under enlargement. It does not yet give a qualifying seed or place a prescribed prime.

## Slot vectors and a finite test at a single mark

In a labeled tournament A replace vertex x by a nonempty tournament Q, making every arc between Q and an exterior vertex uniform, with its original direction at x. If A' is a cut switch of A, give all of Q the original cut membership of x. The resulting pair is still related by a cut, hence has the same unmarked switching lattice.

In an ordered cover by q optional paths, compress each maximal Q-run to x. For each compressed word use all exterior vertices once and let r count occurrences of x, with no two consecutive within a slot. Write c_r(A,x,q) for the number of valid words. A word is filled by an ordered cover of Q by r nonempty paths. Thus, with f_Q(j) counting j optional slots,
```
f_{A[Q/x]}(q) = sum_r c_r(A,x,q) N_Q(r),
N_Q(r) = sum_{j=1}^r (-1)^(r-j) binom(r,j) f_Q(j).
```
There is at most one first run per slot; every later run needs an intervening exterior path from an out-neighbor to an in-neighbor, disjoint from the separators already used. In particular r is at most q + min(outdegree(x), indegree(x)). Use the larger cutoff R from the two frames for a comparison.

The coefficients in the optional-slot basis follow by an invertible triangular change of basis. If the height difference is nonzero formally, let its largest nonzero coefficient be lambda_J. A neutral nonempty Q requires J >= 2 and
```
(|Q|+J-1)/(J-1) <= (sum_{j<J} |lambda_j|)/|lambda_J|.
```
Indeed splitting one slot at a selected position offers |Q|+J-1 choices, with at most J-1 adjacent merges recovering the original cover; hence f_Q(J)/f_Q(J-1) is at least the left side. All smaller evaluations are at most f_Q(J-1). If J=1 the difference cannot vanish. There is consequently unbounded neutral growth at a fixed mark **if and only if** the whole two-slot vector agrees. Literal transitive replacements of widths 1 through R suffice to certify agreement, since f_TT_k(j)=j^k and the evaluation matrix is an invertible Vandermonde matrix times column factors j. A nonzero one-slot vector difference makes the Hamiltonian count change at all sufficiently large transitive widths; the change grows unboundedly if its highest surviving base exceeds 1.

All 19 marks on the two elementary pairs fail already at width two. The three pairs in the H=1127,1139,1151 lattice triple also have no qualifying mark. The background exhaustive pass through quotient order seven found no mark on any switching pair with equal two-slot vectors and differing one-slot vectors. For order eight, the included search normalizes the marked vertex to a source and uses the 456 archived order-seven exterior representatives. Every cut through an eight-vertex labeled quotient and marked vertex appears after this normalization, up to relabeling. For each exterior, bucket all its switched states by the complete height evaluation vector at widths 1 through 5, then compare the one-slot vectors in the buckets. This tests through quotient order eight, **not** through expanded order eight. There were 13,023 additional equal-height-vector entries behind first representatives, all with equal ordinary vectors. A separate order-nine random scan of 1,400 exteriors found no hit (22,720 additional entries), without any claim of exhaustiveness at order nine.

## Simultaneous substitution

Mark s distinct quotient vertices and substitute independent nonempty blocks Q_i. Compress each monochromatic maximal run, where the color is the marked vertex whose module supplied it. Because every between-color arc is uniform, the same proof gives a tensor:
```
f_expanded(q) = sum_{r_1,...,r_s} c_{r_1,...,r_s} product_i N_{Q_i}(r_i).
```
For the run cutoff of each marked color, count possible intervening **quotient letters with multiplicity**, not just distinct exterior vertices: a neighboring marked module can be revisited. A finite bound follows from a simpler fact. For distinct colors i,j their order is uniform, so a return to color i must use a closed directed walk of length at least 3 on quotient colors. This alone does **not** bound the number of runs when several modules can be revisited. In particular the all-marked directed triangle has the legal word (abc)^k, for every k: the tensor has infinite support. A fixed Vandermonde certificate therefore requires an actual finite-run hypothesis.

The exact finite-run condition is that the marked subtournament in **each** quotient being compared is transitive. Conversely, if it contains a directed triangle x→y→z→x, take a Hamiltonian path through that quotient and replace the occurrence x by (xyz)^k x. This preserves both exterior adjacencies, uses every ordinary vertex once, and makes the number of marked runs unbounded. For a transitive marked subtournament, delete the ordinary exterior vertices from each compressed word by cutting there. Each remaining uninterrupted stretch consists of strictly increasing marked vertices, with no repeats. Thus every color occurs at most q plus the number of exterior vertices. Under this hypothesis, or any separately justified uniform bound, the optional-slot tensor is finite; use its actual cutoff R_i at each mark. A rectangular grid of independent transitive widths 1 through R_i certifies the whole tensor, by a product of invertible Vandermonde matrices. An alternative is one transitive ray with widths w_i k where w_i=(1+floor(log_2 max R_i))^(i-1), for max R_i >= 2. The bases product_i j_i^w_i are all distinct: for every prime, the exponent is a base-(1+floor(log_2 max R_i)) digit encoding of the original exponents. Distinct bases reduce the ray to one Vandermonde system.

## Shared identical blocks: a stronger, physical test

Still assume finite run. When the same Q is installed at all marks, the tensor becomes a homogeneous polynomial F(z_1,...,z_R) of degree s, with z_j=f_Q(j). This groups terms by the **multiset** of slot readings. Setting Q=TT_k groups them further, by the product of their slot indices: for example the monomials z_2 z_6 and z_3 z_4 agree on every transitive width. A zero on the equal transitive ray does not imply neutrality for a shared nontransitive block.

In fact there is a finite physical test for the full shared-block polynomial, and the same block can be used at every mark. Let P_d be the d-vertex transitive tournament with its longest arc reversed. Each directed odd cycle uses that edge, and its other vertices form an increasing subset between the endpoints. No two such cycles can coexist. By the odd-packing formula, for d>=3,
```
f_{P_d}(j) = j^d + j((j+1)^(d-2) - (j-1)^(d-2)).
```
Take any distinct nonnegative exponent vectors a,b indexing monomials. Their evaluations at P_d eventually differ. To see this, put m=d-2 and take the log of the ratio. It has an affine part in m, from the leading terms j(j+1)^m. The remaining log terms use
```
1 + j (j/(j+1))^m - ((j-1)/(j+1))^m.
```
If the affine part is nonzero it eventually dominates. If it is identically zero, take the largest j with a_j != b_j. Its correction is asymptotic to (a_j-b_j) j (j/(j+1))^m; all corrections from lower indices are little-o of it. Either way, equality is impossible for sufficiently large d, also along even d alone.

For the finitely many distinct monomials of interest, find one even d by exact integer comparisons until all the monomial values are distinct. Call them v_a. Use Q equal to the ordinal join of k copies of P_d. Then f_Q(j)=f_{P_d}(j)^k and F evaluates to sum_a F_a v_a^k. Testing k=1 through the number of monomials decides whether F vanishes identically, by Vandermonde. This is a terminating procedure, not a heuristic genericity claim. It gives an all-order dichotomy for shared growth: either the shared two-slot polynomial is zero and *every* identical replacement is neutral, or there is a single explicit nontransitive joining ray along which neutrality eventually fails. Likewise, a nonzero shared one-slot polynomial certifies a genuine move on a suitable such ray, choosing d to separate both supports together.

## If the marks have a directed cycle, use a rational series instead

The infinite run tensor still has a finite description. Fix the order and slot assignment of all ordinary vertices. Between two consecutive ordinary vertices (or a slot boundary), any nonempty stretch of marked colors contributes
`l^t X (I-A_M X)^{-1} r`, where A_M is the marked adjacency, X is diagonal with run variables t_i, and l,r are the allowed entrance and exit indicator vectors. A missing exterior endpoint gives the all-ones vector on that side. Add 1 for an allowed empty stretch. The formal matrix inverse sums all possible lengths; each term records exactly one run at each color occurrence. Multiply these factors across stretches and sum over the finitely many exterior arrangements. Finally apply inclusion–exclusion, zeroing subsets of t_i, to keep monomials containing every mark. For a concrete tuple of modules only terms with r_i <= |Q_i| contribute. This also yields a **finite physical equality test even with cyclic marks**, by switching the evaluation basis. At independent transitive widths $k_i$, the nonempty-slot reads are $N_{TT_{k_i}}(r_i)=r_i!\,{k_i\brace r_i}$, where the braces count set partitions. They are zero for $k_i<r_i$ and nonzero on the diagonal. Thus the evaluations on the finite simplex $k_i\ge1$, $\sum_i k_i\le D$ determine every coefficient with $\sum_i r_i\le D$ by triangular inversion. This is different from the finite-support Vandermonde test.

Take $g=$ (number of ordinary vertices) $+q$, the number of stretches in every exterior skeleton. Before extracting positive exponents, a common denominator is $\det(I-A_M X)^g$. The numerator and denominator both have total degree at most $sg$: for each stretch the inverse's adjugate has degree at most $s-1$, the factor $X$ adds at most one, and adding an empty stretch uses the determinant itself. After extracting positive exponents, use the product of these denominators over all subsets of the $s$ colors. At subset $U$ the degree is at most $g|U|$. For a difference between two frames, use the product for both; the numerator and denominator have degree at most $D=2g\sum_{U}|U|=g s 2^s$, and the denominator has constant term one. A nonzero rational difference therefore has at least one nonzero series coefficient of total degree at most $D$. We conclude:
 
**For any fixed pair of quotient frames, agreement at all independent transitive widths with $\sum_i k_i\le gs2^s$ certifies agreement for arbitrary independent nonempty blocks, including with cyclic marks.** This remains a large test; it certifies universality, not equality for a special shared cyclic family. A failed test gives an actual transitive witness, not merely a signed formal obstruction.

This identifies exactly what information is lost by the transitive test and what brings it back. It also tells us when **not** to use this finite algebra: if either quotient's marked vertices admit a directed cycle, there are arbitrarily many runs, so the finite-support shared-block argument does not apply without more work. Finding a changing frame with a nonzero one-slot motion and neutral two-slot response, and then choosing its absolute height from the prescribed prime, remain separate open steps.

### Reproduction

Run `python3 clone_replay.py` alongside `replay.py` and `frame_replay.py`. The order-eight census uses `g++ -O3 -std=c++17 clone_find.cpp -o clone_find` and `./clone_find 8 456 input < n7_representatives.txt` with the [archived input](../actual_hamilton_multiplicities_20261021/session24-multiplicities/n7_representatives.txt).
