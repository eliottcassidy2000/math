# Matrix products, odd-pair gates, and compatible Moore and Pollock lifts

This extends the [Collatz carry packet](../collatz_mod12_quadruples_and_divisor_grades_20261126/README.md), [weighted expanding passages](../collatz_weighted_passages_first_crossing_20261126/README.md), [Moore route colors](../moore_route_colors_tournament_ports_20261126/README.md), and [Pollock's centered moment filter](../pollock_joint_local_moments_and_selector_20261126/README.md). The common rule is to distinguish a matrix product from the constraints on the states that realize it. The constraints are different in the three problems.

## 1. Collatz: two odd steps are the threshold for an expanding run

Represent the halved map by projective matrices, acting on \((x,1)^\top\):
\[
E=\begin{pmatrix}1&0\\0&2\end{pmatrix},\qquad
O=\begin{pmatrix}3&1\\0&2\end{pmatrix}.
\]
The next matrix is chosen by the parity of the dehomogenized first coordinate; an arbitrary product only acts on the integer Collatz orbit if its parity word is legal. The odd pair is
\[
O^2=\begin{pmatrix}9&5\\0&4\end{pmatrix},\qquad
T^2(x)+1=\tfrac94(x+1)\quad(x\equiv3\bmod4).
\]
The contrasting \(EO\) pair for time word \(10\) scales \(x-1\) by \(3/4\). A run with no adjacent odd pair parses into single zeros and \(10\)'s, with possibly a last unmatched 1. Each completed piece reduces \(x-1\) by at least the factor \(3/4\) while the state is above 1. In particular a positive orbit must reach 1 or hit an \(11\) within logarithmically many such pieces. This states where possible growth enters the word, not that every odd pair leads to net growth.

**Exact return law.** Start at positive odd \(x\), run \(r\ge1\) consecutive odd steps, and then \(s\ge1\) consecutive even steps up to the next odd \(x'\). (The odd run is finite: \(x+1\) cannot be divisible by every power of 2.) The entire block is \(E^s O^r\), and
\[
x'=\frac{3^r(x+1)-2^r}{2^{r+s}},\qquad
x'+1=\left\lceil\frac{3^r}{2^{r+s}}(x+1)\right\rceil .
\]
The ceiling identity uses the legality of the word: the unrounded argument is exactly \(x'+2^{-s}\). As a result, **\(x'>x\) if and only if \(3^r>2^{r+s}\)**, with no extra low-height exception from the carry. If the slope is below 1, the ceiling is at most \(x+1\); if above 1, the ceiling is strictly greater. A one-odd block never expands, and every \(r\ge2\) block with \(s=1\) expands. Extra even steps may erase the expansion. The replay independently tests 10,000 complete runs. This is an exact complement to prior shifted-height bounds, not a Collatz convergence argument.

## 2. Moore: a Gram matrix of the missing branch contacts

Let \(d\) be the putative degree of a Moore graph of diameter two, \(b=d-1\). The usual root chart has \(d\) level-two branches of \(b\) leaves each. There are no leaf-to-leaf edges within a branch. Between distinct branches, the leaves match by a permutation: the unique path of length two from a leaf to the other branch's parent must use just one leaf in that other branch. Write \(X_{ij}=X_{ji}^{\top}\) for those permutation matrices, with rows in \(i\) and columns in \(j\).

Suppose all pairwise matchings among \(q=d-k\) old branches are known, and \(k\) whole branches remain. The off-diagonal old-route residual is
\[
R_{ij}=J_b-X_{ij}-\sum_{\ell\ {\rm old},\ \ell\ne i,j}X_{i\ell}X_{\ell j}.
\]
A negative entry rules out completion. If no entry is negative, \(R_{ij}\) is a 0–1 \(k\)-regular bipartite adjacency matrix, since every used route is a permutation. Write \(C\) for the old leaf adjacency matrix, \(D=I_q\otimes J_b\). The single residual matrix with these off-diagonal blocks and \(kI_b\) on the diagonal is
\[
K=(d-1)I-C^2-C+J-D.
\]
Every completion factors \(K=QQ^\top\), where \(Q\) concatenates the old-to-new matching stars. In particular \(K\) is positive semidefinite with \(\operatorname{rank}_{\mathbb Q}K\le1+k(b-1)\). The same rank upper bound holds mod every prime. The loss of \(k-1\) dimensions is because column blocks for different new branches all have the same column sum vector, the all-ones vector.

On the old within-branch zero-sum space, the polynomial is simply \((d-1)I-C^2-C\); this subspace is invariant under \(C\) because its blocks are permutations or zero. Therefore at degree 57 the restricted \(C\) has spectrum in \([-8,7]\), with at least \((57-2k)55\) eigenvalues in the set \(\{-8,7\}\), counted with multiplicity, whenever \(q\ge k\). Here the constant all-ones direction accounts for one dimension of the Gram rank. In particular:
- **One branch missing**: \(\operatorname{rank}K\le56\), and at least 3,025 restricted eigenvalues are endpoints.
- **Two branches missing**: \(\operatorname{rank}K\le111\), and at least 2,915 are endpoints.

These tests apply once all mutual old matchings are fixed. They do not assume the degree-57 graph exists or that the current partially filled branch bank has already reached that stage.

### One missing branch: an exact completion criterion

When \(k=1\), each nonnegative residual is a permutation. The remaining star exists if and only if these residuals satisfy
\[
R_{ij}R_{j\ell}=R_{i\ell}
\]
on triples of old branches. Indeed the missing contact blocks \(B_i\) must have \(R_{ij}=B_iB_j^\top\); fixing an anchor and its new-branch labels reconstructs every \(B_i\). Given this factorization, old-to-new routes are automatically disjoint: a collision of the direct route with a route through old branch \(j\) contradicts disjointness of \(R_{ij}\) from the direct old edge; a collision through \(j\) and \(\ell\) contradicts disjointness of \(R_{j\ell}\) from the route through \(i\). These routes then cover \(b\) slots. Routes between two leaves in the new branch are already blocked from sharing an old leaf by the matching property; their parent supplies the unique common neighbor. Hence every needed leaf route is correct, and the root and parent routes follow from the full star and matching chart.

Equivalently, with nonnegative residuals, \(K\succeq0\) alone suffices at \(k=1\): a Gram representation assigns orthonormal frames to each branch, and a permutation cross-block forces the paired unit vectors to coincide. Also, \(\operatorname{rank}K\le b\) alone suffices here, by the Schur complement at any anchor block \(I_b\).

Pairwise coverage does not imply the product law. A small adversarial bank has three old branches of width 3 with \(X_{01}=X_{02}=I\), \(X_{12}=P\) a 3-cycle. Every old pair has disjoint occupied routes and one residual permutation, but \(R_{01}=P,\ R_{12}=R_{02}=P^\top\) violate the product law.

### Two missing branches: a synchronous split and one final matching

For two missing branches \(U,V\), split every \(R_{ij}\) coherently as \(B_i^U(B_j^U)^\top+B_i^V(B_j^V)^\top\), with \(B_i^U,B_i^V\) permutations. To enumerate splits, anchor both stars at a common old branch with independent new coordinate labels. Each residual to that anchor is a union of even cycles in a bipartite graph, so its splits are its two alternating edge colorings independently on each component. The proposed split at every other pair is then forced and checked.

Given a coherent split, for old \(i\) let \(D_i^U\) be the residual to \(U\) after taking the direct edge and all intermediate old branches. These routes are disjoint by the same old-pair argument as at \(k=1\), and \(D_i^U\) is a permutation. The final matching from \(V\) to \(U\) exists if and only if \((B_i^V)^\top D_i^U\) is the **same** permutation \(Y\) for all \(i\). To see sufficiency, a collision in the route to \(V\) via \(U\) with a direct route to \(V\) would, read backwards through \(Y\), give a collision at \(U\) with its direct edge. A collision with a route to \(V\) through old \(j\) would instead give a collision from \(j\) to \(U\) between routes through \(i\) and \(V\). Between \(U,V\), two routes through distinct old branches would create duplicate contributions to one old-pair residual; overlap with the direct matching would contradict the old-to-\(U\) disjointness. All surviving route families have the required size. Thus only one of the two new-destination matching tests has to be run explicitly.

The replay removes every one or two branches at every root of the Petersen and Hoffman–Singleton graphs and independently checks these reconstruction conditions: respectively 60 and 1,400 cases.

## 3. Pollock: the matching Gram condition has a genuine collision

At fixed center \(c\) write \(x_i=c+y_i\), \(\sum_i y_i=0\), \(Q=\sum y_i^2/2,\ B=\sum y_i^3/6,\ p_4=\sum y_i^4\). Here the representation \(N=5F(c)+cQ+B\) uses \(F(x)=(x^3-x)/6\). Let \(V\) be the \(5\times3\) matrix of rows \((1,y_i,y_i^2)\). Its Gram matrix is
\[
V^\top V=\begin{pmatrix}5&0&2Q\\0&2Q&6B\\2Q&6B&p_4\end{pmatrix}.
\]
Cauchy–Binet gives the determinant as the sum over triples of squared Vandermonde products. It vanishes exactly on at most two distinct offsets, the rank-deficient stratum in the previous joint moment test. For integer offsets with at least three distinct values its minimum positive value is **12**: every triple of distinct integers contributes at least 4, and five slots with at least three distinct values give at least three such triples. Equality requires three consecutive values with multiplicities \(3,1,1\) in some order; centering forces \((-1,0,0,0,1)\). Thus the determinant distinguishes a regular stratum, but not individual solutions.

For a concrete collision the centered multisets
\[
(-17,-5,-4,12,14),\qquad (-16,-10,3,7,16)
\]
have **the same Gram matrix** with \((Q,B,p_4)=(335,-105,143554)\), and determinant 178158400. They yield distinct positive representations of every \(5F(c)+335c-105\) for \(c\ge18\); at \(c=18\), both give 10770. Newton identities show exactly what the matrix omitted: after fixing these first four power sums, the monic root polynomial has only its constant term free,
\[
z^5-Qz^3-2Bz^2+(Q^2/2-p_4/4)z-e_5.
\]
The two constant terms are 57120 and \(-53760\). This gives a family of integer collisions inside the regular moment stratum, not just local congruences. It also contrasts the two inverse problems: the last Moore star has orthonormal permutation blocks, which turn the Gram data into a rigid gluing; the Pollock Vandermonde rows retain no corresponding equality of paired unit vectors and do not give uniqueness.

## 4. The tournament boundary

A tournament adjacency matrix already encodes the oriented graph, just as \(X_{ij}\) encodes a Moore matching. But powers of an ordinary matrix count walks and forget revisits. For tournament Hamiltonian paths one must lift the state, e.g. mark each visited vertex with a variable \(x_v\) with \(x_v^2=0\), and extract the top squarefree coefficient. A Collatz matrix product forgets legality unless one retains the parity cylinder, a Moore residual forgets which missing branch supplies a contact unless one synchronizes the split, and a Pollock Gram matrix can forget the fifth root coefficient even after integrality. The transferable technique is to identify the missing state and the corresponding compatibility test, not to identify their numerical invariants.

Run `python3 replay.py` to check both Moore controls, the adversarial bank, the return law and the Pollock examples.
