# Closing the order-five decagon: rigid seed and exact completion test

Continuation of [the star/Hoffman–Singleton comparison](README.md). A degree-57 graph is still hypothetical. This note isolates the subcase in which the ten independent orbits of the tight \(a_1=325\) trace case form one whole connected component of the 700-site residual. It does **not** require the other components to close until explicitly stated.

## 1. The closed decagon has no internal phase freedom

Let \(F\) be the fixed Hoffman–Singleton graph, \(C\) the 2500 sites adjacent to it, and \(Z\) the remaining 700 sites. Write \(R=A_Z\) and \(T=A_{C,Z}\), so \(Z\) is 7-regular and
\[
R^2+R+T^\top T=56I+J.
\]
On any zero-sum combination of component indicators the right side minus \(R^2+R\) vanishes, so \(T\) annihilates it. As every \(C\)-row has 14 neighbors in \(Z\), its number of neighbors in each component \(Z_i\) is exactly \(|Z_i|/50\). Hence all component orders are divisible by 50. If the remaining simple edges of the ten-orbit alternating decagon stay inside those orbits, it is one 50-site component \(F'\). Girth and the degree-seven Moore bound identify it as Hoffman–Singleton.

There is a stronger *marked* rigidity statement, with a replayable finite proof. Label the doubled ring \(0,1,\ldots,9\) in order, with alternating half-differences \(\pm1,\pm2\), and gauge nine ring centers to zero. Its unweighted orbit matrix \(Q\) has 2 on the ring, 0 on the diagonal, and three simple edges at each site. No simple edge has circular distance 3, since a length-three ring segment would give a square in the lift. The choices are therefore at distances 2, 4 and 5. Of 814 labelled simple degree completions, only two satisfy
\[
Q^2+Q=6I+5J,
\]
and they differ by one ring rotation. In one, the even sites use their distance-four pentagon, the odd sites use their distance-two pentagon, and the five diameters supply the third simple edge at each site. The other swaps the parity.

For each allowed \(Q\), introduce the fifteen simple voltages and the one remaining center \(h\). Let \(W_{ij}\) be the **sum** of the phase shifts from \(i\) to \(j\) in \(\mathbb F_5\) (twice the center on a double). For every ordered endpoint pair the phase sum of all walks of length one or two must be zero: distinct quotient sites have exactly one route at each of five shifts, and the diagonal has zero-sum shifts as well. Thus
\[
QW+WQ+W=0\pmod5.
\]
The 45 off-diagonal equations have rank **16** on the 16 variables in both quotients. The forced solution has \(h=0\) and every simple shift zero. Both alternating color choices pass the full coefficient-by-coefficient lift test by direct computation. In particular **the four nonzero fluxes cannot close**, even though they pass the sparse-ring girth check before the simple edges are supplied. The closed 50-site block is completely specified once the ring is numbered and colored. The executable test counts all degree completions, solves the linear system, and checks the required five coefficients individually.

## 2. Two Hoffman–Singleton blocks give a grid, not independent attachment choices

If \(F'\) closes, each site in \(C\) has exactly one neighbor in \(F'\). For \((f,h)\in F\times F'\), these nonadjacent sites have exactly one common neighbor; it must be in \(C\). Thus \(C=F\times F'\), where the automorphism acts by \((f,h)\mapsto(f,t(h))\). The \(42\) neighbors of \((f,h)\) **within \(C\)** must use every nonneighbor of \(f\) other than \(f\) exactly once in the first coordinate, and every nonneighbor of \(h\) other than \(h\) exactly once in the second. This follows by taking the unique common neighbor with each vertex in the corresponding block. The remaining 13 neighbors in \(Z\setminus F'\) are not assigned by this chart.

In the extreme case where all of \(Z\) consists of fourteen size-50 components, the fifteen copies (including \(F\)) are pairwise anticomplete. Mačaj and Širáň [already identified this potential configuration](https://www.sciencedirect.com/science/article/pii/S0024379509003735) in their concluding remarks. The attachment array on \(C\) is then \(\operatorname{OA}(2500,15,50,2)\): every pair of symbols in two distinct columns singles out the unique common neighbor. Each closed \(C\)-neighborhood of size 43 projects injectively in every column, with the seven missing symbols precisely the neighbors of the row's symbol in that copy. The identical reasoning applies to five pentagons and 25 rows in the degree-seven model.

In the order-five tight case, all fourteen residual copies are invariant: at least thirteen are needed to support the 130 cyclic orbits, and the invariant component count is 14 modulo 5. Within an invariant size-50 copy, the trace on its 28-dimensional \(2\)-eigenspace is \(a/5-2\) by the projector \((15I+5R_i-J)/25\); here \(a\) counts its vertices moved to an adjacent image. This trace is congruent to 28 modulo 5, and \(0\le a\le25\) because \(t,t^2\) have equal traces and no orbit has more than one step class. So \(a=0\) or \(25\). This gives one independent-orbit decagon block and thirteen all-cyclic blocks.

An all-cyclic block also has no internal search. Its five pentagons of each step class have one simple matching between every opposite pair: equal-class edges make a square, and each site needs five external neighbors. On this \(K_{5,5}\), the phase matrix has each pair of distinct rows differ by a permutation of \(\mathbb F_5\), or there would be a square. Gauge its first row and column to zero and relabel columns to make its second row \(0,1,2,3,4\). The replay enumerates 6 remaining matrices, all row permutations of \(ij\). This is the usual affine Hoffman–Singleton formula. Hence if all fourteen components close, **both kinds of internal block can be fixed in advance**; the still-open decisions are contacts and the graph on \(C\). In the quotient array, each of the 50 fixed symbols indexes ten \(C/t\) rows, and in every other column those ten rows use the ten symbol-orbits once. For two nonfixed columns each symbol-orbit pair uses five row-orbits, whose relative phases run through the five residues. These are constraints on the array, not a proof of its existence with the required \(C\)-edges.

## 3. A completion certificate without any automorphism

For an induced Hoffman–Singleton \(F\), let \(P=A_F\), \(B=A_{C\cup Z}\), and let \(N\) be the 3200-by-50 matrix of attachments to \(F\). Take \(N\) to have disjoint columns of size 50, and \(B\) any simple symmetric graph with degrees 56 on \(C\), 57 on \(Z\), satisfying \(BN=J-N(I+P)\). These two linear constraints imply the correct two upper blocks of \(A^2+A=56I+J\). The last defect is
\[
D=B^2+B-56I-J+NN^\top.
\]
Let \(V\) be the span of the columns of \(N\) and the all-ones vector. On the centered columns \(B\) acts as \(-I-P\). On the two constant sectors its matrix in value coordinates on \(C,Z\) is \(\left(\begin{smallmatrix}42&14\\50&7\end{smallmatrix}\right)\): the \(C\)-to-\(C\) and \(Z\)-to-\(C\) counts follow by summing \(BN\), and degrees give the rest. Direct substitution in each sector shows \(D|_V=0\). On \(V^\perp\), the remaining desired equation is \(B^2+B=56I\). If completed, the spectrum is exactly
\[
(49\pm5\sqrt{161})/2,\qquad
(-3)^{28},\quad 2^{21},\quad 7^{1679},\quad (-8)^{1470}.
\]
The last two multiplicities follow from 3149 remaining dimensions and trace \(-7\).

More strongly, for any candidate \(B\) satisfying the linear constraints above,
\[
\|D\|_{\rm F}^2
=\operatorname{tr} B^4+2\operatorname{tr} B^3-20\,048\,700
=12 C_3(B)+8 C_4(B).
\]
For the first equality expand \((B^2+B-56I)^2\) on \(V^\perp\): \(V\) is fixed, and the first three moments (orders 0 through 2) on all of \(B\) are fixed by its dimensions and degrees. Subtract the formal target moments of orders 3 and 4 from the displayed spectrum. For the second equality \(\operatorname{tr}B^3=6C_3\); the nonsquare contribution to \(\operatorname{tr}B^4\) is \(2\sum d_v^2-\sum d_v=20\,048\,700\). Thus **checking just internal triangles and squares suffices**, once the attachments and degrees are exact. This completion criterion does not assume order five, \(F'\), or an orthogonal array.

### The phase-sum test extends to the entire order-five branch

In the [640-orbit certificate](README.md#exact-simultaneous-test-without-flattening-the-voltages), write each existing voltage as an integer residue in \(\mathbb F_5\). Let \(Q=M(1)\) and \(W\) be the matrix of first moments (total phase per entry). The equation \(QW+WQ+W=0\) still holds, since the sum of the five shifts in \(U\) is zero mod 5. It is a linear prefilter for **any** proposed orbit support, with changes of phase origins in its kernel. In characteristic 5 put \(z=\exp u=1+u+u^2/2!+u^3/3!+u^4/4!\) in \(\mathbb F_5[u]/u^5\), and let \(M=\sum_{i=0}^4 M_i u^i/i!\), so \(M_i\) records the \(i\)th raw phase moment and \(M_0=Q\). The change of variables is invertible since the first derivative is 1, and \(U=(z-1)^4=u^4\). For \(r=1,2,3,4\) the equations are
\[
M_r+QM_r+M_rQ+\sum_{i=1}^{r-1}\binom ri M_iM_{r-i}
= \mathbf 1_{r=4}\,4!\,(J-S^\top S).
\]
With integral augmentation, fixed/free equations, degrees, and the support rules imposed, these mod-five equations recover the **exact** coefficient constraints. In an off-diagonal entry the count is nonnegative, total 0 or 5, with desired coefficients 0 or 1 each. On the diagonal the zero-shift count is already fixed by reciprocal shifts and the free degree, leaving total 0 or 4 across nonzero phases. Thus none of the modular certificates can hide an excess multiple of five. The closed decagon is a case where the first moment alone fixes the gauge class. Outside the closed subcase no such rank claim is made.
