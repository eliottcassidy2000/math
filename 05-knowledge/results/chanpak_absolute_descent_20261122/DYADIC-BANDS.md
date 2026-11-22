# Absolute dyadic bands and a fixed-backedge compiler

This packet produces exact strong counts with explicit order bounds, including a literal unit edge at every sufficiently high dyadic height on one rail. The constructions are not a proof that all sufficiently large odd integers occur. A counting argument shows that no bounded number of exceptional backedges can suffice on its own.

## 1. One physical edge gives \(2^m+3\leftrightarrow2^m+5\)

Use the seven-vertex tournament with outgoing masks
\[
 (4,9,34,21,71,27,47)
\]
on vertices \(0,\ldots,6\). The upper quotient reverses only the pair \((3,5)\). Both quotients are strong and have the same port \(z=6\); in both, its only inneighbor is vertex 4. Substitute any nonempty tournament \(U\) homogeneously for \(z\). Every return between two successive visits to \(U\) must use an exterior arc from an outneighbor of \(z\) to an inneighbor of \(z\), hence it must pass 4. There can therefore be at most two visits.

Let \(Z_U(q)\) count ordered \(q\)-slot path covers of \(U\), allowing empty slots, and let \(h=H(U)\). One-visit skeletons are ordinary Hamilton paths of the quotient; two visits are filled in \(Z_U(2)-2h\) ways. The quotient counts and skeleton counts are exactly
\[
 (c_1,c_2)=(67,32),\qquad (69,32).
\]
They can be independently read by literal path counts at transitive port widths one and two: these are \(67,131\) and \(69,133\). Hence the pair of counts for the same incoming \(U\) is
\[
 \bigl(32Z_U(2)+3h,\ 32Z_U(2)+5h\bigr).                    \tag{1}
\]
The six old exterior vertices are not replaced, so the difference is still exactly the one labelled reversal \((3,5)\), for every \(U\). Strongness lifts from the quotient: paths between quotient blocks lift, while a directed cycle through the mark connects any two vertices within the replaced block.

In particular set \(U\) to the transitive tournament of width \(m-5\). Then \(Z_U(2)=2^{m-5}\), \(h=1\), and for every \(m\ge6\) this gives a literal strong pair
\[
 (2^m+3,\ 2^m+5)
\]
on exactly \(m+1\) vertices. The height address \(m\) is an exact target selector, not a congruence. It does not claim a unit seam elsewhere in the same dyadic interval.

## 2. Five independent multiplicative offsets at every dyadic height

The following small seed masks give strongly connected tournaments with first-slot and two-slot profiles \((h_i,2^{q_i})\).

| \(h_i\) | \(q_i\) | outgoing row masks |
|---:|---:|---|
| 5 | 5 | (2,8,3,5) |
| 19 | 8 | (2,16,11,19,37,15) |
| 217 | 12 | (2,24,19,133,201,31,47,103) |
| 223 | 12 | (2,136,11,145,71,31,47,117) |
| 599 | 13 | (28,73,146,148,98,15,45,115) |

For an ordinal sum \(A\to B\), each slot splits uniquely into its \(A\)-part and \(B\)-part in order. Thus both \(H\) and \(Z(2)\) multiply; adding an ordinal singleton preserves \(H\) and doubles \(Z(2)\). For nonnegative integers \(e_i\), set
\[
 D=\prod h_i^{e_i},\qquad b=\sum q_ie_i.
\]
For any \(m\ge\max(1,b)\), take \(U\) to be the ordinal sum of these copies and \(m-b\) singletons. Then \(H(U)=D\), \(Z_U(2)=2^m\). Put \(U\) into one vertex of the directed triangle, leaving the other two singleton. Its skeleton profile is \((3,1)\), so the resulting strong count is
\[
 D+2^m,
\]
on exactly \(m+2-\sum_i(q_i-|A_i|)e_i\le m+2\) vertices. This extends the nested long-backedge chart beyond a single residual component.

The seed counts are multiplicatively independent (each has a private prime divisor). The number of nonnegative vectors \((e_1,\ldots,e_5)\) with \(\sum q_ie_i\le m\) is \(\Theta(m^5)\), by containment between fixed dilates of a five-dimensional positive box. The resulting offsets are all distinct, and
\[
 D\le 2^{\rho m},\qquad \rho=\max_i(\log_2 h_i)/q_i<0.71.
\]
Thus this gives \(\Theta(m^5)\) distinct absolute count values in the thin band \([2^m+1,2^m+2^{0.71m}]\), each on at most \(m+2\) vertices. To invert a specified band target, take \(m=\lfloor\log_2 N\rfloor\), subtract \(2^m\), and check whether the residue is a product of the five counts with the indicated weighted exponent budget. This is an exact factor test in a fixed prime set. It does not assert that every integer in the band passes.

If instead \(m\ge\max(6,5+b)\), give the same blocks \(m-5-b\) ordinal singletons and use (1). The resulting physical arc joins
\[
 (2^m+3D,\ 2^m+5D),
\]
on \(m+1-\sum_i(q_i-|A_i|)e_i\) vertices. Its gap is \(2D\), not a unit seam unless \(D=1\). This is a literal summand/multiplicand rule with a fixed quotient and variable directed components.

### A sixteen-factor certified band

The same construction admits the following eleven additional strong seeds, independently checked by the replay. Together with the five above, they give a sixteen-parameter band with \(\Theta(m^{16})\) distinct counts and offsets at most \(2^{0.73m}\), still on at most \(m+2\) vertices. To prove multiplicative independence first eliminate these new seeds by their private primes \((37,73,139,23,761,773,101,163,103,179,389)\), in the row order below, and then use independence of the original five. The bound on the weighted logarithms is exact using \(h_i^{100}<2^{73q_i}\).

| \(h_i\) | \(q_i\) | outgoing row masks |
|---:|---:|---|
| 259 | 13 | (256,257,259,71,47,335,407,319,24) |
| 803 | 14 | (256,1,27,147,99,79,143,311,126) |
| 695 | 14 | (256,257,43,307,327,147,303,95,164) |
| 713 | 14 | (256,261,289,7,15,411,63,95,216) |
| 761 | 14 | (256,5,49,39,331,339,143,319,78) |
| 773 | 14 | (256,269,297,193,271,283,55,375,72) |
| 707 | 14 | (256,273,267,451,141,287,311,103,144) |
| 815 | 14 | (256,21,25,67,41,79,407,319,62) |
| 1957 | 15 | (258,36,401,295,75,405,47,91,210) |
| 1969 | 15 | (258,132,353,23,39,331,155,61,218) |
| 1945 | 15 | (260,69,112,327,43,139,177,287,118) |

### Why the unit rail does not iterate as a dyadic-slot seed

For the literal \(3\leftrightarrow5\) quotients of §1 with transitive insertion width \(w\), their two-slot counts are respectively
\[
 Z_-(2)=8[11(3^{w+1}-1)+18\cdot 2^w],\qquad
 Z_+(2)=8[11(3^{w+1}-1)+19\cdot 2^w].
\]
Indeed in two ordered cover slots there can be at most three total runs of the mark: at most two first runs, and each return needs the same exterior vertex 4. Thus as functions of \(w\ge1\) the counts are linear combinations of \(1,2^w,3^w\), by the onto-allocation formula for each fixed number of runs. The first three literal widths have respective counts \(992,2864,8192\) and \(1008,2896,8256\), determining the displayed polynomials.

For \(w\ge4\), \(v_2(3^{w+1}-1)<w\), by LTE (the valuation is 1 for odd \(w+1\) and \(2+v_2(w+1)\) otherwise). Therefore the bracket's 2-adic valuation is exactly that of its first summand, and the positive bracket is strictly larger than that power of two. Neither profile is dyadic in those widths. The only dyadic profile at smaller widths is the lower one at \(w=3\), giving the seed \(h=259\) above. The non-dyadic \(3^w\) channel thus gives a concrete obstruction to recursively feeding this unit-height rail back into the same dyadic compiler, even though the one-slot count has a pure dyadic leading term.

## 3. Compiling an arbitrary fixed backedge pattern

Start with an ordered transitive tournament and reverse a prescribed finite list of arcs. Mark every endpoint of a reversed arc; these endpoints come in a fixed linear order. The clean vertices in the gaps between marks have no reversed incidences, and their gap widths vary. There may be zero widths. Only gaps covered by some reversed arc matter; clean vertices in uncovered gaps are ordinal singleton factors and can be erased.

**Exact compiler.** For any fixed ordered pattern of \(r\) backwards arcs, the Hamilton count is a finite integer sum
\[
 H(w_1,\dots,w_s)=\sum_\alpha c_\alpha
                 \prod_{j=1}^s q_{\alpha,j}^{\,w_j},\qquad
 q_{\alpha,j}\in\{0,1,\dots,r+1\},
\]
where \(0^0=1\). The table depends on the pattern, not on the widths. To construct it, enumerate the ordered list of actual descents of a spanning path, i.e. the backward arcs it uses and their order. Between descents the path is an ascending run. Assign the finitely many marked endpoints to these runs, enforcing the prescribed starts and ends. Each clean vertex is then independently assigned to any run whose endpoint bounds permit it. There is one remaining validity condition: an unused backward arc cannot appear as an ascending adjacency of its two endpoints in a run. For each such event possible given the marked assignment, perform inclusion-exclusion, forbidding all intermediate clean vertices in that run when the event is selected. Each gap then has an integral number \(q_{\alpha,j}\) of permitted runs, constant across that gap. Sorting inside each run is forced, giving the displayed factorization. Because descents and marked assignments are recovered uniquely from the path, there is no multiplicity in this compilation.

For example, with four ordered marks \(0<1<2<3\) and clean widths \(x,y,z\), the nested reversals \(3\to0,2\to1\) yield, at \(y\ge1\),
\[
 H=1+2^y+2^{x+z+1}(3^y+2^{y+1}-1).
\]
Setting \(y=2\) gives \(2^m+5\) on \(m+1\) vertices for \(m=x+z+5\). The crossing reversals \(2\to0,3\to1\), at \(y=0,\ x,z\ge1\), give
\[
 H=2(2^x+1)(2^z+1)-1.
\]
These formulas are checked by the accompanying compiler and independent subset-endpoint DP, including zero-width controls.

One further affine unit seam comes from this crossing pattern at \(x=1+w,\ z=3,\ y=0\). The upper seed at \(w=0\) has rows \((246,252,120,241,224,192,128,4)\); reverse only \((5,7)\) for the lower seed. Insert \(w\) clean vertices between 0 and 1 in both. The exact heights are
\[
 (36\cdot2^w+15,\ 36\cdot2^w+17)
\]
on \(w+8\) vertices. Both literal hosts remain strong because the consecutive forward spine is untouched and the two original backward intervals cover it and overlap. With at most three backward arcs the dependence on \(w\) is in the span of \(0^w,1^w,2^w,3^w,4^w\). The independent DP evaluations for \(w=0,1,2,3,4\) are respectively \(51,87,159,303,591\) and \(53,89,161,305,593\), proving the formula by the invertible five-column Vandermonde test.

## 4. A quantitative obstruction to any bounded-backedge key

Let \(L=\lfloor\log_2 X\rfloor\). If a reversed arc spans \(u\ge1\) clean vertices, every odd-size subset of those clean vertices, taken in increasing order between its endpoints, yields a directed odd cycle closed by that arc. There are \(2^{u-1}\) choices. The positive odd-cycle-packing formula of [Grinberg–Stanley](https://arxiv.org/abs/2307.05569) counts each such single-cycle family with weight 2, in addition to the empty family. Therefore \(H\ge1+2^u\). In particular every retained gap in a witness with \(H\le X\) has width at most \(L\). With at most \(2r\) marked endpoints, deleting uncovered ordinal vertices yields an equivalent count on at most \(2r+(2r-1)L\) vertices. This is a **witness order upper bound conditional on the chart**, not a claim that every target is in it.

There are at most \(O_r((1+L)^{2r-1})\) distinct heights below \(X\) with any such bounded-backedge presentation: there are finitely many ordered endpoint patterns, and each has at most \(2r-1\) widths in \([0,L]\). Consequently for fixed \(r\) one can decide exact target membership in this entire set of charts in time polynomial in the bitlength, by compiling all endpoint patterns and checking the bounded widths. Neither the compiler nor the order bound can make this architecture cofinite.

For varying \(r\le L\), count ordered backedge sets by
\(\binom{\binom{2r}{2}}{\le r}\le (r+1)(2er)^r\), with an additional factor \(2r\) for the endpoint count, and bound widths by \((L+1)^{2r-1}\). The base-two logarithm of the count is at most \(r(\log_2 r+2\log_2 L+O(1))+O(\log r)\). Thus for any fixed \(c<1/3\), all but \(o(X)\) odd integers below \(X\) have **no** presentation with \(r\le c\log_2 X/\log_2\log_2 X\) backwards arcs; any possible larger witness must genuinely enlarge its exceptional directed pattern. This does not rule out a fixed *local* edge rule repeated a growing number of times.
