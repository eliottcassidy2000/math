# One sum, three blocks: the rooted split normal form

Fix a strong source with retained root and initial state \((H_0,h,L_0)\). Here \(h\) is the Hamilton count with the root deleted, and \(L_0\) is the rooted relaxed two-piece count. A seam adding \(r\ge2\) vertices (block exponent \(r-1\)) acts as
\[
(H,h,L)\longmapsto (H+2^{r-1}L,\ h,\ 2^r L-h).
\]
Use the [breakpoint expansion](../session325_absolute_cut_overlap_20261114/READY_OVERLAPS.md): for total added size \(R\), record the completed block positions \(B=\{s_1,\dots,s_k=R\}\), consecutive gaps at least two. Then
\[
2(H_{\rm out}-H_0)=L_0\sum_{s\in B}2^s-h\sum_{s<t\in B}2^{t-s},\qquad
L_{\rm out}=2^R L_0-h\sum_{s\in B}2^{R-s}.
\]

## A sharp obstruction to independent additive switches

Separately changing an interior breakpoint in each of two disjoint windows looks like two independent size-preserving repairs. It never is from one fixed source. Let the two replacement choices in the left window have breakpoint sets \(A_0,A_1\), and similarly \(B_0,B_1\) in the right window; all four whole words must be legal, and the other breakpoints are fixed. Write \(H_{ij}\) for the outputs.

The mixed difference is exactly
\[
2(H_{11}-H_{10}-H_{01}+H_{00})
=-h\left(\sum_{a\in A_1}2^{-a}-\sum_{a\in A_0}2^{-a}\right)
    \left(\sum_{b\in B_1}2^{b}-\sum_{b\in B_0}2^{b}\right).
\]
Indeed internal-window terms and all terms involving only one varied window cancel; the only cross terms are the products \(2^{-a}2^b\). Distinct finite sets of powers of two have distinct sums, even with negative exponents after clearing denominators. Thus both parentheses are nonzero when both choices are genuine. In particular just adding optional breakpoints \(u<v\) gives the nonzero defect \(-h\,2^{v-u}\) in twice the count. This interaction is the exact obstruction to an additive cube in this fixed-source grammar; it contrasts with physically independent square-sum moves and does not rule out changing roots.

## A section that retains the whole split

At fixed total size \(R\) and fixed deletion count \(h\), the single reading \(L_{\rm out}\) recovers the entire breakpoint word and the initial \(L_0\), even if the source was not supplied. Since \(h\) is odd, reduce
\[
b=\sum_{s\in B}2^{R-s}\equiv-h^{-1}L_{\rm out}\pmod{2^R}.
\]
The true binary word \(b\) lies in \([0,2^R)\) and is therefore uniquely determined. Its one-bits are the breakpoints read backwards, and \((L_{\rm out}+h b)/2^R=L_0\). Given also \(H_{\rm out}\) the initial \(H_0\) follows by subtraction. In particular different decompositions of the same total size and same deletion count cannot have the same full seam state, even if the source changes. This is an exact observation about what a section retains; it is not reconstruction of the underlying tournament.

## One block versus exactly three is rigid

Put \(x=L_0/h\). Compare the single block of size \(R\) with a triple having breakpoints \(a,b,R\), all three gaps at least two. Let \(d=b-a\) and \(u=R-a-b\), \(v=d-a\). The counts agree exactly when
\[
 x=2^u+\frac{2^v}{1+2^d}.
\]
This follows by cancelling the common final term \(L_0 2^R\) and using
\[
\frac{2^{b-a}+2^{R-a}+2^{R-b}}{2^a+2^b}
=2^{R-a-b}+\frac{2^{d-a}}{1+2^d}.
\]
The reduced odd part of the denominator of \(x\) must be exactly \(1+2^d\), since a dyadic term cannot cancel any of that factor against the power of two in the numerator. So \(d\) is forced, or there is no solution. After multiplying by \(1+2^d\), the binary expansion is
\[
(1+2^d)x=2^u+2^{u+d}+2^v.
\]
The legal final gap \(R-b=u+a\ge2\) gives \(v\le u+d-2\). Hence \(u+d\) is the highest bit; there is no carry into it, since \(d\ge2\) and \(v\) is at least two positions below it. Thus \(u\) is forced by the highest bit, then \(v\) by subtracting the other two terms. This determines \(a,b,R\) uniquely, subject to the stated gap tests. Per source there is **at most one** one-block versus three-block collision. The older one-versus-two identity can instead recur indefinitely when \(x\) is dyadic; the two phenomena even require disjoint ratio classes.

The first concrete control in the packet uses an order-five source with rows \((2,4,9,19,7)\), root \(0\), and \((H_0,h,L_0)=(9,5,11)\). The block sizes \((7)\) and \((2,2,3)\) each give \(713=9+704\) at order twelve, with final \(L\) values 1403 and 1203 respectively. Only the three-block route uses all genuine Fermat blocks. Another same-size realization is ordinal multiplication \(713=23\cdot31\) on two strong six-cores. The replay checks all three literal counts and the different port readings. None of this makes those three constructions interchangeable.
