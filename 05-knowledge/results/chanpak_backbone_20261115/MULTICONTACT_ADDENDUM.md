# Multi-contact transfer: three exact neighbors, cut coordinates, and a growing third digit

All counts below come from literal tournament orientations. The results construct families of windows, not windows at every prescribed absolute target.

## Three consecutive odd counts at every height

In the two-contact construction of [ROOTED_CLOCK.md](ROOTED_CLOCK.md), write the count as \(L2^m+M3^m+D\). The triangle with a fixed outer root and either inner port gives \((L,M,D)=(9,4,\pm1)\). A transitive triple with the outer root its source and inner port its sink also has \(L=9,M=4\), and gives \(D=-3\): \(Q_r=4,Q_u=1,A=2,B=-2,C=1,H=1\). Thus
\[
9\,2^m+4\,3^m+\{-3,-1,+1\}\qquad(m\ge1)
\]
is a three-count window, starting at 27, 29, 31. Against an arbitrary shared ordinal tail R, the same cores keep clocks \(9P_2(R),4P_3(R)\) and give offsets \(\{-3,-1,+1\}H(R)\), so the resulting window is consecutive only if \(H(R)=1\).

## Sign-free cut interpretation

In the square-zero ring let \(X=\operatorname{diag}(x_v)\), \(A\) be the 0/1 tournament adjacency, and \(M=X(I-AX)^{-1}\) be the matrix of path generating functions with length measured by vertices. Write \(f=M\mathbf1,\ g^t=\mathbf1^tM,\ h=1+\mathbf1^tM\mathbf1\). Matrix inverses exist by nilpotence; adjacent repeated vertices contribute zero, so the \(-I\) in \(A^t=J-I-A\) may be dropped in this path readout. Sherman–Morrison for \(I+AX-JX\), followed by \(x\mapsto-x\), yields
\[
h(x)\big(M(x)+M^t(-x)\big)=f(x)g^t(x).
\]
In fact \(X(I+AX)^{-1}=-M(-x)\) and its rank-one correction has denominator \(h(-x)\), giving the displayed formula. Summing against the all-ones vector also gives \(h(x)f(-x)=-g(x)\), using \(h(x)h(-x)=1\). This last identity is the previously recorded rooted slot reciprocity; the matrix version retains both endpoints.

Let \(J_{ru}\) count all nonempty directed \(r\to u\) paths, weighted by \(H\) of their complement. Let \(K_{ru}\) count covers by two nonempty directed paths, one starting at r and the other ending at u. Extracting full-set coefficients, \(2C_{ru}=J_{ru}+J_{ur}-K_{ur}\), and \(B=-K_{ur}-K_{ru}\). Hence
\[
D=H-K_{ru}-J_{ru}-J_{ur}.
\]
In particular for \(Q_r=Q_u\), swapping the two root roles holds both clock coefficients fixed and changes the offset by the directional difference of the full two-path covers. Neither weighted reachability total changes. This is a direct count condition for selecting a cut, rather than a test on sign-weighted data.

## Many contacts: exact first two scales

Append \(m\ge1\) ordinary transitive vertices and then ordered new vertices \(b_1,\ldots,b_k\), with \(k\ge2\). All new edges point forward except \(b_i\to r_i\), where \(r_1,\ldots,r_k\) are distinct old roots. For \(S\subseteq[k]\), let \(A_S\) be the sum over tuples of nonempty vertex-disjoint directed paths \(P_i\) starting at \(r_i,\ i\in S\), each weighted by \(H\) of the uncovered old vertex set. The Hamiltonian path count has expansion
\[
k! A_{[k]}(k+1)^m+(k-1)!\sum_{j=1}^k(k+1-j)A_{[k]\setminus j}\,k^m+O_{T,k}((k-1)^m).
\]

Use [odd directed cycle packings](https://arxiv.org/abs/2307.05569) and expand each odd-cycle filter as 1 minus its even/odd size sign, so selected odd cycles receive weight 2. A contact is active if its backward arc is in the packing. Once active rooted paths are fixed, each ordinary vertex is unused or assigned to one active path, always in ascending order. The paths connect active new apices by a permutation \(\pi\); any inactive apex used in a cycle belongs to a suffix going to a later active apex. For fixed apex assignments the signed filters contribute exponential terms with integer bases of absolute value at most one plus the number of available active path colors. Deleted forward arcs affect configurations only by forbidding a color on all ordinary vertices: \(r_j\to b_j\) could be used only if inactive \(b_j\) is the first apex on a suffix from a path ending at \(r_j\), and that suffix color receives no ordinary vertices. Inclusion-exclusion for these events reduces the number of available colors.

With all k active there are no inactive apex assignments or deleted-arrow problems: a would-be deleted edge from an active path to its own starting apex closes an even 2-cycle, while every other active root is already reserved. Expand parity filters by the subset \(I\subseteq[k]\) formed by the selected permutation cycles. Permutations preserve I; the factor from signs of the selected cycles includes the alternating sum \(\sum_{\sigma\in S_I}(-1)^{\#\text{cycles}(\sigma)}\), zero for \(|I|\ge2\) by permutation signs. For \(I=\varnothing\) the contribution is \(k!(k+1)^m\). Singleton I contributes \((k-1)!(k-1)^m\sum_i(-1)^{|P_i|}\). These are the only all-active terms. With \(k-1\) active contacts and unique inactive \(b_j\), ignoring deleted edges the leading contribution for each permutation counts \(1+k-j\) options: omit \(b_j\), or assign it to one of the slots whose terminal apex is later than j. A nonempty parity selection has base of absolute value at most \(k-2\), and an edge correction costs one color. Thus that case contributes the indicated \(k^m\) coefficient. Fewer contacts cannot contribute at these scales.

## A growing physical third digit

Now assume every old vertex is a port, \(k\ge3\), and the old tournament is doubly regular with common indegree \(d=(k-1)/2\) and common pairwise in-neighbor count \(c=(k-3)/4\). Vary only the correspondence between its vertices and the ordered b's. The third coefficient, at base \(k-1\), is
\[
D_k+(k-2)!\,b(\rho),
\]
where \(D_k\) depends only on k and \(b(\rho)\) is the number of old backward arcs in the port order \(\rho\). The residual is \(O_{T,k}((k-2)^m)\). The first two coefficients are independent of \(\rho\), since \(A_{[k]}=1\) and \(A_{[k]\setminus j}=1+d\).

Here are all sources of the third coefficient:

1. All k contacts contribute \(-k!\), since all k old paths are singletons. With k−2 contacts and two inactive indices \(i<j\), each of those apices has a fixed number of legal terminal apices depending only on i,j, and omission is allowed. Their choices are independent even if both go on one path, since they occur there in increasing order. The associated old path-tuple count is \(d^2+2d-c\), independent of the old pair: neither omitted vertex is used (1), exactly one is used (\(2d-1\)), both are used in distinct paths (\(d(d-1)-c\)), or both are used in one path (d). All parity filters or deleted-arrow corrections in this k−2-contact case are lower order.
2. With k−1 contacts the only contribution at base \(k-1\) is deletion of the forward arrow into the single inactive \(b_j\). Its bad use requires an old arc \(r_i\to r_j\), as the active path at i must end at \(r_j\) and has exactly two vertices. It also requires that the i-slot receive no ordinary vertices and terminate at an active index greater than j. The unselected parity term contributes \((k-1)^m\) for each such permutation. A selected-cycle parity term can contribute at the same scale only when that cycle uses just the forbidden-color slot, i.e. \(\pi(i)=i\). Then the actual cycle there is even (two old vertices, \(b_j\), \(b_i\)), so its filter kills the leading term. Thus the count subtracted for this old arc is, at the third scale, \((k-2)!((k-j)-[i>j])\). Summing over all old arcs gives a constant depending on d,k plus \((k-2)!b(\rho)\). No other cases have this scale.

Every integer value of b between its extrema occurs: permutations are connected by adjacent transpositions, and each such transposition flips b by one. The range is at least \(2\sqrt{k(k-1)/24}\). Indeed with uniform random order each arc contributes variance 1/4, disjoint pairs are independent, and adjacent pairs have covariance ±1/12 according as their arrows at the shared vertex point the same way or opposite ways. At a regular vertex the sum of these covariance signs over unordered adjacency pairs is \(-(k-1)/2\), so the total variance is \(k(k-1)/24\). Order reversal makes the max and min equidistant from the mean, giving the stated variance bound.

This creates a contiguous growing alphabet in the third *coefficient*, not a contiguous interval of actual odd counts. For fixed k, raising m separates the digit since the residual is strictly lower exponential order, but it does not remove that residual. An exact-target theorem must still link the chosen digit to a physical lower-scale repair.

For the three-cycle the formula specializes further, with \(b=1,2\), to exact counts
\[
6\,4^m+24\,3^m+(12+b)\,2^m-7\qquad(m\ge1).
\]
The two lower allowed bases are \(+1,-1\) (base zero vanishes at positive m); the two residual evaluations at m=1,2 by literal subset path DP determine their coefficients. The replay includes checks through m=8 for both port orders.
