# Regular cross-block enumeration for two cyclic tournament cells

## Scope and definitions

Let p be an odd prime, m=(p−1)/2, and let R_p have vertex set F_p and arcs i→i+d for 1≤d≤m. Consider tournaments consisting of two induced R_p cells X,Y, with a binary cross-matrix B: B_ij=1 means X_i→Y_j. Suppose every row and column of B has sum k. This note completely counts this restricted family. It does not count all tournaments or prove a universal metagraph wiring rule.

The two cells have vertex scores m+k and m+p−k. For 0≤k≤m these differ, so the cells are intrinsic score classes. Within R_p, the arc i→i+d has exactly m−d common outneighbors. Thus the short arcs i→i+1 are intrinsic, and Aut(R_p)=C_p. Isomorphism within a fixed k-layer is exactly independent cyclic rotation of rows and columns.

Swapping the two cells replaces B by 1−B^T and k by p−k. Consequently the unordered two-cell family is counted exactly once by k=0,1,…,m. The k=0 layer has one ordinary class, one self-converse class, and one merged class. Every nonuniform layer is strongly connected: each cell is strongly connected and there are cross-arcs in both directions. The uniform layer alone is not strongly connected.

## Exact all-prime coefficient formulas

Define

A_(p,k) = [x_1^k⋯x_p^k] e_k(x_1,…,x_p)^p,

S_(p,k) = [x_1^k⋯x_p^k] ∏_i(1+x_i) ∏_(i<j)(1+x_i x_j).

A counts all binary matrices with every row and column sum k. S counts the symmetric ones. The diagonal factors in S encode cross-edges X_iY_i, not tournament loops.

For 1≤k≤p−1, the numbers of ordinary classes N, self-converse ordinary classes F, and classes merged with their converses M are

N_(p,k) = (A_(p,k)+(p−1)^2 binom(p,k))/p²,

F_(p,k) = (S_(p,k)+(p−1) binom(m, floor(k/2)))/p,

M_(p,k) = (N_(p,k)+F_(p,k))/2.

### Ordinary Burnside proof

Use the action (g_(a,b)B)_ij=B_(i+a,j+b). Identity fixes A matrices. If exactly one of a,b is nonzero, all rows or all columns are identical; regularity with 0<k<p is then impossible. If both are nonzero, the grid has p translation orbits, each meeting each row and column once. Choosing k of these p orbits gives binom(p,k) matrices. Averaging over p² group elements proves N.

### Converse Burnside proof

Converse, followed by cell exchange and internal negation, is f(B)_ij=B_(-j,-i). It normalizes the rotation group. Fixed rotation-orbits under f are counted by averaging the fixed sets of the p² maps g_(a,b)f.

Their fixed equation is B_ij=B_(-j-b,-i-a). If a=b, relabel the column j as −j−a; the condition becomes symmetry and contributes S. There are p such maps.

If a≠b, squaring this transformation gives translation by (a−b,b−a). Since p is prime, invariance forces B_ij=c_(i+j). The remaining condition is c_t=c_(-t-a-b). This reflection of F_p has one fixed point and m two-cycles. An invariant k-subset selects floor(k/2) two-cycles, and contains the fixed point precisely when k is odd. Hence each of the p(p−1) remaining maps fixes binom(m,floor(k/2)) matrices. Dividing by p² proves F. Pairing the nonfixed ordinary classes proves M.

The k=1 specialization has A=p! and S equal to the number I_p of involutions in S_p, recovering the matching-cell formulas. The complete two-cell totals are 1+Σ_(k=1)^m N, 1+Σ F, and 1+Σ M.

## Exact finite values

| p | k | A | S | N | F | M |
|---|---|---:|---:|---:|---:|---:|
|3|1|6|4|2|2|2|
|5|1|120|26|8|6|7|
|5|2|2,040|112|88|24|56|
|7|1|5,040|232|108|34|71|
|7|2|3,110,940|6,912|63,504|990|32,247|
|7|3|68,938,800|35,150|1,406,940|5,024|705,982|

Including the uniform layer gives:

| cells | ordinary | self-converse | converse-merged |
|---|---:|---:|---:|
|two R_3 cells|3|3|3|
|two R_5 cells|97|31|64|
|two R_7 cells|1,470,553|6,049|738,301|

The two-R_5 family therefore has 96 strongly connected ordinary classes and 63 strongly connected merged classes. These integers arise from the specified quotient enumeration; no identification with unrelated occurrences of 63, 64, or 88 is claimed.

## Independent degree-two decomposition

A degree-two binary bipartite graph is a disjoint union of even cycles. For p=5 there are only two possibilities: one 10-cycle, or a 4-cycle and a 6-cycle. The first contributes 5!4!/2=1440. The second contributes binom(5,2)^2·3!2!/2=600. Thus A_(5,2)=2040.

Every component has two alternating choices in an ordered pair of disjoint perfect matchings. The 1440 connected matrices consequently occur twice, and the 600 disconnected matrices four times, in the ordered matching-pair universe. This explains an actual 2-versus-4 fiber mechanism, rather than a numerical analogy. Independently enumerating that universe gives 88 rotation orbits: eight of size 5 and eighty of size 25. Exactly 24 are self-converse, leaving 56 merged classes.

A symmetric degree-two matrix is a graph with optional diagonal marks counting one toward degree. Its components are cycles of length at least three, or paths with at least two vertices and diagonal marks at both ends. For five vertices the contributions are: 5-cycle, 12; marked 5-path, 60; triangle plus marked edge, 10; marked 3-path plus marked edge, 30. Their sum is S_(5,2)=112.

## Degree-two generating functions for arbitrary size

Let A_(n,2) and S_(n,2) denote the same matrix counts for any n. Include the empty matrix with count one. Bipartite cycle decomposition gives the two-sort exponential generating function

Σ_(n≥0) A_(n,2) z^n/(n!)²
 = exp(Σ_(j≥2) z^j/(2j))
 = exp(−z/2)(1−z)^(-1/2).

Indeed, a connected component with j labels on each side has j!(j−1)!/2 realizations. The symmetric component decomposition gives

Σ_(n≥0) S_(n,2) z^n/n!
 = exp(Σ_(j≥3) z^j/(2j) + Σ_(j≥2) z^j/2)
 = (1−z)^(-1/2) exp(−z/2−z²/4+z²/(2(1−z))).

The cycle and marked-path counts on j labels are (j−1)!/2 and j!/2 respectively. These are component-counting identities, with no external novelty claim.

For reproduction, put a_n=A_(n,2)/(n!)². Then a_0=1, a_1=0 and n a_n=(1/2)Σ_(j=0)^(n−2) a_j. For S, put c_2=1 and c_j=(j+1)(j−1)!/2 for j≥3. Then S_0=1, S_1=0 and S_n=Σ_(j=2)^n binom(n−1,j−1)c_j S_(n−j).

## Verification and provenance

`check_regular_blocks.py` multiplies the labelled truncated coefficient polynomials literally, without sorting capacity states. It independently enumerates all unions of two disjoint perfect matchings at p=5, verifies complement-transpose between k=2 and k=3, computes quotient and converse counts, and compares the degree-two generating-function recurrences to the coefficient calculations. Run with Python normally and with `-O`; explicit RuntimeError checks remain active.

A separately implemented check at `../check_regular_block_coefficients.py` uses sorted remaining-column capacities for A and symmetric-row degree recursion for S. Its p=3,5,7 outputs agree with this packet. The p=7 values above are coefficient computations, not an enumeration of 68,938,800 matrices. All-prime statements follow from the proofs, not extrapolation from finite controls.
