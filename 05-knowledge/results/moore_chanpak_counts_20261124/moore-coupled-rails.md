# Coupled two-port unit rails seeded by a compressed pentagon

This extends the three single-port rails in [the compact observer note](moore-compact-count-rails.md). The same two seven-vertex pentagon observer charts also supply two growing **equal-width two-port** \(+2\) rails. One of these is genuinely coupled: its literal reversal has unit magnitude if and only if the transitive widths are equal.

For two marked vertices \(p,q\) and module tournaments \(X,Y\), a Hamilton path is obtained by first writing a compressed word with the exterior vertices appearing once and the port symbols recurring. Each consecutive arc in the word must be legal in the uninflated chart. Adjacent repeats of one port cannot occur, since they would be just one visit. If the ports occur \(a,b\) times, fill them by ordered path covers counted by \(P_a(X)P_b(Y)\). Between successive copies of either port there must be a genuine exterior vertex: the other port alone cannot support a return because the direction between the two ports is fixed. Thus each visit count is at most the exterior order plus one, and a finite subset recurrence gives the coefficient table \(b_{a,b}\).

For equal transitive widths \(k\ge1\), expand each ordered-cover count
\[
 P_a(\operatorname{Tr}_k)=\sum_{\ell=1}^{a}(-1)^{a-\ell}\binom{a}{\ell}\ell^k .
\]
Then group by the integer base \(\ell s\) on the diagonal of the two width variables. Distinct positive integer exponentials are linearly independent on positive integer \(k\) by Vandermonde. Thus equality of all nonconstant base coefficients and a difference of two in the constant term is an exact test for an all-width seam, not unproved finite extrapolation.

## Complete finite result in the specified architecture

Use either compressed pentagon chart and at most one preparatory reversal, require strong hosts and growth, mark any two distinct vertices, and make one exterior reversal not incident to either port. The [exact replay](moore-coupled-rails.py) checks all such choices. No direct growing equal-width unit seam exists from either unedited chart. There are 12 marked choices after one preparation (8 from \(H=147\), 4 from \(H=135\)), yielding just two lower curves:
\[
\begin{aligned}
 D_k&=16\cdot9^k-20\cdot6^k+33\cdot4^k-29\cdot2^k+11,\\
 E_k&=66\cdot9^k-180\cdot6^k+187\cdot4^k+64\cdot3^k-196\cdot2^k+71 .
\end{aligned}
\]
Each upper curve is exactly two more, by one specified arc reversal on the same \(2k+5\) vertices. The former occurs twice (both from \(H=147\)) and the latter ten times.

For \(D\), start at the \(H=147\) rows \((106,100,25,82,35,76,20)\), reverse \((0,4)\), mark vertices \((3,4)\) as \(X,Y\), and reverse the exterior \((5,6)\) for the upper curve. The two-variable difference is
\[
\begin{aligned}
 \Delta={}&2P_1(X)P_1(Y)
 +17(P_2(X)P_1(Y)-P_1(X)P_2(Y))\\
 &+8(P_3(X)P_1(Y)-P_1(X)P_3(Y))
 +4(P_3(X)P_2(Y)-P_2(X)P_3(Y)).
\end{aligned}
\]
Write \(h=P_1\), \(Q=P_2+2h\), and \(R=P_3+3P_2+3h\); these are the one-, two-, and three-slot counts with empty slots permitted. The same formula simplifies to
\[
 \Delta=2h_Xh_Y+5(Q_Xh_Y-h_XQ_Y)+4(R_XQ_Y-Q_XR_Y).
\]
For transitive widths \(a,b\), \(h=1\), \(Q=2^a,2^b\) and \(R/Q=(3/2)^a,(3/2)^b\). Both ratios increase strictly. At \(a>b\), \(\Delta\ge12\); at \(a<b\), \(\Delta\le-8\). Thus \(|\Delta|=2\) iff \(a=b\), and then \(\Delta=2\). The antisymmetric terms also cancel if \(X,Y\) have identical profiles, but then the unit becomes \(2h_X^2\), not a unit for equal nontransitive fillers. With different arbitrary fillers the normalized equation
\[
 \frac{\Delta}{h_Xh_Y}
 =2+5(q_X-q_Y)+4q_Xq_Y(r_X-r_Y),\qquad
 q_X=Q_X/h_X,\quad r_X=R_X/Q_X
\]
is an exact signed cancellation condition. In particular increasing both profile ratios cannot give \(+2\) with a nontransitive filler.

For \(E\), take the same \(H=147\) rows, prepare by reversing \((1,2)\), mark \((0,6)\), and use final flip \((2,4)\). Its slot difference is simply \(2P_1(X)P_1(Y)\), so it is also a unit on independent transitive widths, not just equal ones.

## CRT improvement and the finite scope

Together with the three single-port rails, the lower residues modulo \(3\cdot5\cdot7\cdot11\cdot13=15015\) have per-rail cardinalities \(60,48,44,48,48\) under a common period 60. They form 243 distinct lower residues. Their translation by \(+2\) intersects them in exactly two phases, for 484 distinct endpoint residues in all.

The \(D\) example proves that adaptive coupling is a genuine way around the independent-width interpolation restriction: a nonconstant antisymmetric tensor can vanish on the matching-width diagonal. It does **not** defeat the sparsity of finitely many bounded-order charts, give arbitrary nontransitive unit fillers, or show that a degree-57 Moore graph exists. A possible role for degree57 would be to provide many mutually compatible, independently steerable coupled port profiles, rather than merely a larger fixed chart.
