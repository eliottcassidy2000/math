# Two root contacts produce a physical 2/3 clock and an exact unit window

Status: **proved**. The two-contact formula holds for arbitrary underlying tournaments. The sparsity result applies only to chains with the specified first-anchor outlet, not to arbitrary cores. The unit window is not universally placed.

## An arbitrary-core two-contact formula

Write H(T) for the number of directed Hamiltonian paths (and H(empty)=1). For distinct ports r,u of T, append a transitive block of m ordinary vertices, then b then v, initially with all new arcs forward in this ordinal extension. Reverse only the arcs b→u and v→r, meaning make those the backward arcs. Take m≥1. For a nonempty directed path P, write |P| for its number of vertices. Define
\[
Q_x=\sum_{P_x} H(T\setminus P_x),\quad
A_{ru}=\sum_{P_r\perp P_u} H(T\setminus(P_r\cup P_u)),
\]
\[
B_{ru}=\sum_{P_r\perp P_u} H(T\setminus(P_r\cup P_u))
\left((-1)^{|P_r|}+(-1)^{|P_u|}\right),\quad
C_{ru}=\sum_{\substack{P:r\leadsto u\\ |P|\ {\rm odd}}}H(T\setminus P).
\]
Here each P_x starts at x; perpendicular means vertex-disjoint. A and B are symmetric in r,u; C is directional. The physical count of the extended tournament is
\[
\boxed{N_m=H(T)+2^m(2Q_r+Q_u)+2\,3^mA_{ru}+B_{ru}-2C_{ru}.}
\]

Use the [Grinberg–Stanley odd-cycle packing formula](https://arxiv.org/abs/2307.05569), which assigns weight 2 to every selected odd directed cycle. If neither backward arc is used, H(T) remains. If only b→u is used, choose a rooted path in T and a subset of the m ordinary vertices to close the cycle, giving weight 2^m per rooted path and residual packing. An odd cycle cannot close directly from the singleton u using the missing u→b, since the purported cycle would be even. If only v→r is used, b is an extra ordinary candidate, and the initial baseline is 2^{m+1}Q_r; but the deleted forward arrow u→b must not be used. It occurs exactly when the old path goes from r to u in odd order, and finishes by u→b→v without visiting the transitive block. Subtract 2C_{ru}.

When both reverse arcs are used, two disjoint rooted paths in T are selected. The ordinary vertices are assigned to the first path closure, the second path closure, or neither. If the backward arcs belong to distinct cycles, their even-order path parity conditions select weighted count
\[
3^m+(-1)^{|P_r|}+(-1)^{|P_u|}+(-1)^{|P_r|+|P_u|+m}.
\]
This is the two-variable parity filter applied to the three choices per vertex, including the weight 4 for two cycles. If the two arcs are on the same cycle, its route from v→r follows P_r, its first chosen ordinary set, b→u, P_u, the second ordinary set and then closes at v. Its weighted count is
\[
3^m-(-1)^{|P_r|+|P_u|+m}.
\]
The old paths are disjoint, so P_r cannot end at u and P_u cannot end at r: no deleted-edge correction applies to this coupled case. In the disjoint case a singleton P_u cannot close directly at b under the odd parity constraint. Adding the two expressions cancels the mixed parity term. This gives the claimed formula after summing all rooted paths and packings in their complement.

**Physical unit windows at every height.** For T=C3, each Q_x=3; each A=2, B=−2; and for r fixed, C=0 at its outneighbor, C=1 at its inneighbor. Thus the two choices of u realize
\[
9\,2^m+4\,3^m-1,\quad 9\,2^m+4\,3^m+1\qquad(m\ge1).
\]
The starting pair is 29,31. This uses the *same core* with different contact orientations, preserves both exponential coefficients, and changes exactly the absolute term. More generally whenever Q_r=Q_u, exchanging the roles of r and u leaves the 2^m and 3^m coefficients and B untouched; the only possible discrepancy is 2(C_{ru}−C_{ur}) with the sign depending on which graph is subtracted. This is a constructive test for larger cores, not a universal existence theorem.

## Exact rooted chain grammar

For the backward chains of [BACKWARD_CHAINS.md](BACKWARD_CHAINS.md), let r=a0. Append v beating only r and receiving all other arcs. Denote the chain count by H, the count after the outlet by N=H+2E. Then
\[
N=H+F-1,\qquad
\begin{pmatrix}F'\\U'\\R'\end{pmatrix}
=M_m\begin{pmatrix}F\\U\\R\end{pmatrix},\qquad
M_m=\begin{pmatrix}2\,2^m-1&2\,3^m&2\,3^m+2\\2^m&0&0\\0&3^m&3^m\end{pmatrix},
\]
starting at (1,0,0) and reading gaps left to right.

**Transfer derivation.** Define Z(t) as the weighted sum of all nonempty paths starting at r, where the path carries weight t to the power of its number of vertices and the complement contributes H. First regard its complement as a packing of odd cycles. Selected reverse arcs partition into runs, either belonging to the rooted path (P-runs) or to a packing cycle (C-runs). Runs with adjacent selected reverse arcs share an anchor and hence belong to a single path or cycle. A P-run has all its anchors consecutive in descending order on the rooted path. Distinct P-runs and singleton P anchors occur along that path in increasing order of their disjoint anchor intervals: between runs only increasing arcs are available. An ordinary vertex inside a P-run can be placed before/after that run or left unused, weight 1+2t. Inside a C-run on r consecutive gaps and M ordinary vertices, each such vertex can be assigned to the cycle, to the path, or omitted. The parity filter and weight 2 give contribution \((2+t)^M+(-1)^r t^M\), in two separate transfer states. If a gap is not part of either kind of run, any ordinary vertex there can be on the rooted path or unused. The factor is A=(1+t)^m, except if both endpoints are singleton P anchors: at least one ordinary vertex must be used between them because the direct forward arc is reversed, making the factor A−1. These descriptions are reversible, as all unspecified forward contacts are present and the ordinal cut rule identifies every reverse run's compulsory anchors.

Thus after each anchor, track counts X (unselected singleton anchor), Y (singleton P anchor), P (endpoint of P-run), and C+,C− (endpoint of C-run with each parity summand). An open run is eligible both to finish or to continue at each vertex but must have at least one edge. With \(L=X+Y+P+C_++C_-\), update on the next gap m:
\[
X'=A L,\quad Y'=t(AL-Y),\quad P'=t(1+2t)^m(P+tX),
\]
\[
C_+'=(2+t)^m(C_++X),\qquad C_-'=-t^m(C_-+X).
\]
For example starting a new run at the preceding anchor uses X as the count of admissible configurations ending in an unselected singleton there, before activating that anchor. At t=−1, one has \(P'+C_-'=-(-1)^m(P+C_-)\), so this sum stays zero from the initial state (0,t,0,0,0); also L'=L, giving Z(−1)=−1 at all orders. At t=1, P=C+ from initialization and identical updates. Put \(F=L,U=X,R=P=C_+\) to get the three-dimensional formula above. Odd-cycle packing applied to the added outlet gives \(2E=Z(1)+Z(-1)=F-1\), since the cycle through v consists precisely of v and an even-order rooted path in the old tournament, closing through an old last vertex distinct from r. The first-anchor orientation ensures the signed sum remains exactly −1; this is not asserted at a general root or in a general tournament.

**Full state remembers every gap.** Since F is odd, \(v_2(U')=m\) uniquely decodes the last gap of any nonempty word. From U' decode F, then recover
\[
R=(F'-(2\,2^m-1)F-2R')/2,\quad U=R'/3^m-R.
\]
The empty word is characterized by U=0. Thus the rooted state remembers the ordered word even where the unrooted arithmetic loses it under rank-one reset. The two C3 prefix encodings (1) and (0,0) have the same H=3 but states (3,2,0) and (3,1,1). After the same positive gap m, their unrooted state (3,3·2^m) agrees, and their rooted states differ by exactly (2,0,0), which gives the same unit window. It does not replicate unmodified at an arbitrary interior cut: \(M_0^2-M_1\) changes U+R by 2R before the common gap, and the unrooted B2/B0 replacement also reads the already nonzero selected contact.

## Sparse support even after adding the outlet

The number of rooted-gap **words** with N≤X is O(X^{24/25}), hence also bounds the distinct path counts. This is stronger than a small-window experiment. In particular changing every gap and this outlet cannot place a universal window.

Let \((F_r,U_r,R_r)=M_0^r(1,0,0)\). Group the letters into blocks of r≥0 zeros followed by one positive gap m, then a final zero tail. All matrix entries are nonnegative. Each block multiplies the first coordinate by at least
\[
a_{r,m}=(2\,2^m-1)F_r+2\,3^mU_r+(2\,3^m+2)R_r.
\]
The final r-zero tail multiplies it by at least F_r. The lower bound \(F_r\ge2^{\lfloor r/2\rfloor}\) follows by positivity and \(F_2=3\). At s=24/25, the exact rational certificate in the replay gives \(\sum_{r\ge0,m\ge1}a_{r,m}^{-s}<99/100\). For terms outside r≤39 and m≤19, use \(a_{r,m}\ge 2^{m+\lfloor r/2\rfloor}\) and \(2^{-24/25}<3/5\); the same geometric bounds as for the unrooted edge transfer apply. The final-tail partition sum converges, and the concatenation of finitely many blocks gives a geometric series because the block mass is less than 1. Since N≥F, the result follows from finiteness of the Dirichlet series over words.

The obstruction leaves the arbitrary-core formula open as a separate positive strategy. The cutoff is not a property of using two outlets; it is a consequence of the word grammar's restricted prefix and the growth of the 2/3 clock below any fixed target.

## Change the reservoir without pretending the coefficients are independent

Let an arbitrary tournament D be an ordinal tail of the marked core T, before the two-contact construction above. Write \(P_j(D)\) for the number of ordered decompositions of its vertices into j labeled directed paths, permitting empties. The formula becomes exactly
\[
N_m(T\to D)=
(H+B-2C)H(D)+2^m(2Q_r+Q_u)P_2(D)+2\,3^m A\,P_3(D),
\]
where all capital coordinates on the right other than D's slot counts belong to T.

Indeed each rooted path starting in T can only append an optional directed segment in D, so one such path contributes a factor \(P_2\) and two contribute \(P_3\). A rooted path to the other port cannot leave T, giving a factor \(H(D)\) in the direction-sensitive term. For B, consider the square-zero ring \(\mathbf Z[x_v:v\in V(D)]/(x_v^2)\) and
\[
h_D(x)=\sum_{S\subseteq V(D)} H(D[S])\prod_{v\in S}x_v.
\]
By odd-cycle packing applied in every subset, it factors as \(\prod_v(1+x_v)\prod_{C\;{\rm odd}}(1+2x_C)\). Every factor is supported on an odd-size vertex set, so \(h_D(-x)=h_D(x)^{-1}\) in that ring. One sign-weighted path slot and two unsigned slots consequently have total generating function \(h_D(-x)h_D(x)^2=h_D(x)\), giving B's claimed factor.

For the fixed triangle core with alternate ports this gives the unrestricted literal family
\[
9\,2^m P_2(D)+4\,3^m P_3(D)\pm H(D).
\]
The clock pairs can move with D; the offset in a matched D is of size H(D), unlike the original unit pair. A transitive D with \(\ell\) vertices merely increases m by \(\ell\); nontransitive D may change the two clock coefficients at different rates.

## A factorial rigidity theorem for retaining clocks while changing the base

Suppose two same-order n>0 reservoirs D,D' preserve \(P_j\) for \(j=2,\ldots,t\), with t≥2. Then
\[
H(D)-H(D')\equiv0\pmod{(t-1)!(t+1)!}.
\]
In particular, preserving the 2- and 3-clock coefficients forces H congruent modulo 48. The argument needs only the odd-cycle structure, not any census. In the square-zero factorization above, raise h to a positive integer power x. A vertex factor becomes \(1+x x_v\), and a cycle factor becomes \(1+2x x_C\), so the polynomial \(P_D(x)\) has leading term \(x^n\), and all lower coefficients even, all exponents positive and congruent to n modulo 2. For \(\epsilon=1\) when n is odd and 2 when n is even, the difference of the two slot polynomials is \(2x^\epsilon f(x^2)\), with integer polynomial f. The matched readings force divisibility by the monic product \(\prod_{j=2}^t (z-j^2)\) in \(\mathbf Z[z]\). At z=1 its magnitude is \((t-1)!(t+1)!/2\), giving the formula including its factor two. If n<2t+1 for n odd or n<2t+2 for n even, f is identically zero by degree. This proves the divisibility, not the existence of reservoir pairs attaining equality.

Thus even allowing arbitrary same-order reservoirs cannot yield a unit step **by shifting H alone at fixed slots**. In the two-port construction a unit step can instead come from changing the core's directional odd-path correction, as the triangle demonstrates, or from changing one or both clock coefficients with the prescribed absolute target. This is a structural obstruction at the exact lost interface, not an exclusion from the full tournament spectrum.
