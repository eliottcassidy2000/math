# Selecting the opposite end of a line at a prescribed path count

Fix a source height \(h\), the number of directed Hamilton paths. A tiling presents a Hamilton path in the order \(0\to1\to\cdots\to n-1\). Its line partner retains the consecutive arcs and reverses every other arc. Repeated copies of the same rooted strong component can sit at different distances along the path. Those distances are the variables in the inverse below.

**PROVED.** At every fixed attained \(h\), a finite list of tensors with three entries per axis determines every integer height opposite \(h\) on a single tiling line, at all orders. Each axis satisfies the tribonacci recurrence. In particular, for any nonconstant rational polynomial \(P\) and integer \(b\ge2\), only **finitely many** exponents \(p\ge0\) can give \(P(b^p)\) as an opposite height of a fixed finite set of source heights. The proof of finiteness uses the finite-rank unit-equation theorem. Testing a *prescribed* integer target is effective by the elementary recurrence and its growth bound.

## The finite geometric data

The strong components of the source are intervals of its path. Heights multiply over this ordered list. Every nonsingleton strong component has height at least three and bounded order at fixed \(h\), by Busch's lower bound. Keep the ordered word of nonsingleton components, each with its **literal anchored tiling**. This leaves a finite collection of words. A word of \(k\) components of total order \(s\) has \(d=k+1\) available gaps for strings of singleton vertices. Write their lengths \(r_0,\ldots,r_k\). Every such placement gives a valid source of the required height; the word and gaps are recovered uniquely from its tiling. Reflection reverses the word, reflects the individual panels, and reverses the gap vector. A line is blue precisely when this full record is fixed. These facts combine the prior interval and path-fiber theorems. [THM-4953, tiling line cut sets](../../../01-canon/theorems/THM-4953-every-metagraph-line-joins-a-strong-closed-end-to-an-open-end-whose-cut-set-is-counted-exactly-blue-lines-join-only-eternally-blue-values-fibres-1-3-5-carry-1-n-2-5n-15-lines.md) · [THM-796, tiling/line/node incidence](../../../01-canon/theorems/THM-796-three-sorted-recursive-tiling-line-node-incidence.md)

Let \(F(\mathbf r)\) be the partner's height for this word. For **each** gap, with all others fixed, including at zero,
\[
 F(\ldots,r+3,\ldots)=F(\ldots,r+2,\ldots)+F(\ldots,r+1,\ldots)+F(\ldots,r,\ldots).
\]
To prove this for an interior gap, write its left and right portions in the partner as \(L,R\) with last left vertex \(u\) and first right vertex \(v\). All cross edges point right to left except \(u\to v\). Append the \(r\) gap vertices to the left. In its converse, each new vertex has its sole exit to the preceding vertex. For a rooted tournament let \(z=(H,e,B)\), where \(e\) counts full paths ending at its root and \(B\) counts ordered covers by two paths with the second starting at its root, the first possibly empty. The root moves to the new vertex by the transformation
\[
 Mz=(H-e+B,\ H-e,\ H+B),\qquad M^3=M^2+M+I .
\]
At the cut the partner height is \(H(L)H(R)-e_u(L^{\mathrm{op}})e_v(R)+B_u(L^{\mathrm{op}})B_v(R)\). This follows by subtracting the old full paths using \(v\to u\) and adding the ordered covers concatenated at \(u\to v\). Thus the final height is a fixed linear readout of \(M^r z(L^{\mathrm{op}})\). An outside gap is a one-sided tower, using converse as needed. With an empty word both sides of each nonconsecutive source arc are transitive in path order before taking the partner, giving the classical tribonacci board with counts \(1,1,1\) at orders \(0,1,2\). [THM-4865, sole-exit and two-root join identities](../../../01-canon/theorems/THM-4865-the-sole-exit-tower-is-tribonacci-h-k-equals-h-k-1-plus-h-k-2-plus-h-k-3-for-every-rooted-base-the-lineage-c2-law-at-any-vertex-and-the-three-root-response-numerations-attain-every-odd-number-in-69-702189.md)

For an explicit selector, let \(w(0),w(1),w(2)\) be the three unit coordinate vectors, and propagate \(w\) by tribonacci. Index a tensor by \(\mathbf j\in\{0,1,2\}^d\). Its entries \(A_{\mathbf j}=F(\mathbf j)\) come from actual boards of order at most \(s+2d\). Then
\[
          F(\mathbf r)=\sum_{\mathbf j} A_{\mathbf j}\prod_i w(r_i)_{j_i}.
\]
With \(\lambda>1\) the real root of \(x^3-x^2-x-1\), induction gives \(\lambda^{r-2}\le\sum_jw(r)_j\le\lambda^r\). Consequently
\[
 (\min A)\lambda^{n-s-2d}\le F(\mathbf r)\le(\max A)\lambda^{n-s}.
\]
For a target \(N\) this confines \(n\) to a bounded-length interval around \(\log_\lambda N\), uniformly in \(N\) at fixed \(h\). At a candidate order there are \(O_h((1+\log N)^k)\) placements per word. For a blue line the gaps must form a palindrome, giving at most \(O_h((1+\log N)^{\lceil d/2\rceil-1})\). Evaluating the tensor recovers an actual board for every positive match. The calculation allows more than one source isomorphism class or presentation.

## Counting physical placements at each size and modular target

Fix a modulus \(q\ge2\). The companion matrix of tribonacci has determinant one. Let \(L\) be its finite order modulo \(q\). The vectors \(w(r)\) repeat modulo \(q\) with period \(L\) from the initial position. For a word as above and a target residue \(c\), the exact series over its source tilings, weighted by order, is
\[
 z^s\frac{\displaystyle\sum_{0\le j_1,\ldots,j_d<L\, :\,F(\mathbf j)\equiv c\pmod q}
                   z^{\sum_i j_i}}{(1-z^L)^d}.
\]
Subtract the empty tournament term if excluded. Every summand independently extends a physical gap by a nonnegative multiple of \(L\). For blue lines use only a reflection-fixed word, identify opposite gap variables, and give each paired variable weight two in numerator and denominator. Thus the blue denominator is \((1-z^{2L})^{\lfloor d/2\rfloor}(1-z^L)^{d\bmod2}\). Subtraction gives black incidences. Summing the finite bank proves that the counts by size are eventually quasipolynomials of degrees at most \(k\) in total and \(\lceil d/2\rceil-1\) in blue for this word. These count **anchored source incidences**; a line whose two ends both meet the source criterion contributes twice. At sufficiently large order the opposite end exceeds \(h\) uniformly.

There is also an exact rational series for the sum of the *integer* opposite heights over all placements. Replace each basis vector coordinate in the tensor by its ordinary generating series:
\[
 z^s\sum_{\mathbf j} A_{\mathbf j}\prod_i W_{j_i}(z),\qquad
 (W_0,W_1,W_2)=\frac{(1-z-z^2,\ z-z^2,\ z^2)}{1-z-z^2-z^3}.
\]
Both series retain the order of the physical panels. Averaging over abstract permutations of the panels would lose source classes.

## Why a bounded source cannot follow Euclid exponents indefinitely

The three distinct roots \(\alpha_i\) of the tribonacci polynomial are algebraic units. Expanding each recurrence in its eigenbasis expresses \(F(\mathbf r)\) as a sum of at most \(3^d\) fixed nonzero algebraic coefficients times monomials \(\prod_i\alpha_{j_i}^{r_i}\); discard zero coefficients. Work in the splitting field and consider \(F(\mathbf r)-P(b^p)=0\). Select a minimal zero subsum containing its highest-degree term \(c b^{Dp}\), where \(D=\deg P\). Such a subsum exists by partitioning a zero sum into minimal zero subsums. Divide by that term. Every remaining term belongs to a fixed finite-rank multiplicative group, generated by \(b\), the three roots and finitely many coefficients. For each of finitely many index subsets, the unit-equation theorem gives finitely many tuples of term ratios. [Evertse, Schlickewei and Schmidt, *Linear equations in variables which lie in a multiplicative group*](https://annals.math.princeton.edu/2002/155-3/p04)

At a prime ideal over any rational prime dividing \(b\), every tribonacci monomial has valuation zero. Hence each of those ratios has valuation equal to a fixed constant minus \(Dp\,v(b)\), or, for a polynomial term of degree \(j<D\), a fixed constant minus \((D-j)p\,v(b)\). There is at least one ratio and either kind determines \(p\). This proves the claimed finiteness; it does not locate the final exponent. For \(P(X)=(X^2-X)/2\pm1\) this applies to every formal Euclid neighbor, including the perfect-number subsequence. A bounded source bank can realize only finitely many of them by one line. Other target constructions and paths using several lines are outside this restriction.

There is also an elementary obstruction at a real scale. Normalize by \(\lambda^n\). Each vector \(w(r)/\lambda^r\) has a limit at infinity since the other two characteristic roots have modulus \(\lambda^{-1/2}\). The normalized height therefore extends continuously to \((\mathbb N_0\cup\{\infty\})^d\), a countable compact space. Its image is countable, compact and bounded away from zero. The logarithms modulo one in base \(\lambda\) form a countable compact subset of the circle, over all words of any fixed finite bank. Some open interval of phases is missed at every order. The unit-equation argument further prevents infinitely many *exact* exponential coincidences even at exceptional phases.

## A sharp inverse when the source has three paths

The sole strong panel is the directed triangle, with one presentation. Put \(r,s\) singletons on either side. The nine partner seeds, in increasing gap order, are
\[
 A=\begin{pmatrix}1&5&9\\5&13&23\\9&23&41\end{pmatrix},\qquad
 F(r,s)=w(r)^T A w(s).
\]
For \(m=r+s\ge1\) the least partner is \(a_m=F(0,m)\) and the greatest is \(b_m=F(1,m-1)\). Successive bands are disjoint: \(b_m<a_{m+1}\). Thus there is at most one possible order for any prescribed opposite height. The roots of this constraint are the physical boundary and the one-site interior placement.

For the lower bound the difference \(F(r,s)-a_{r+s}\) satisfies tribonacci in each variable. On \(r,s\in\{1,2,3\}\) its rows are \((4,8,12),(8,12,20),(12,20,32)\), so it is strictly positive in the interior. For the upper bound the difference \(b_{r+s}-F(r,s)\) vanishes when either coordinate is one. On \(r,s\in\{2,3,4\}\) its rows are \((0,4,4),(4,12,16),(4,16,20)\), so the only further interior equality is \((2,2)\). The boundary values satisfy \(a_m\le b_m\) by the lower bound. Finally \(a_{m+1}-b_m\) has positive first three entries for \(m=1,2,3\) and follows tribonacci. All partner counts at this order have the **same mod-four class**, equal to \(a_m\pmod4\): check the nine seeds for \(F(r,s)-a_{r+s}\) and propagate. The residue is 3 at \(m\equiv3\pmod4\) and 1 otherwise.

Only \(r=s\) is blue. The generating function for the opposite count along those placements is
\[
 \sum_{r\ge0}F(r,r)z^r=
 \frac{1+11z+12z^2+2z^3-z^4-z^5}{1-2z-3z^2-6z^3+z^4+z^6}.
\]
Indeed the six potential exponential bases are the three squares and three pair products of the roots, with cubics \(x^3-3x^2-x-1\) and \(x^3+x^2+x-1\). The first six values \(1,13,41,129,457,1533\) give the numerator. At order 25, placements \((8,14)\) and \((11,11)\) give the same opposite height 2,295,169, respectively on a black and a blue line. The partner is the full tribonacci board with the arc at distance two over the triangle reversed. At the first placement its outdegree multiset has three 9s, at the second three 12s; the remaining interior degrees away from the defect occur at most once. Even merging with converse keeps the two nodes distinct because 9 reflects to 15, not 12.

## Exact finite search at exponent 61

**FINITE-EXACT.** For sources of every height at most 23, every anchor and every order, the exponents \(3\le p\le61\) with an opposite height equal to either \(E_p-1\) or \(E_p+1\) are exactly \(3,4,5,7\).

| Exponent | Sources at \(E_p-1\) | Sources at \(E_p+1\) |
|---:|---|---|
|3|5, 9, 13, 15, 17|3, 5, 9, 11, 13, 15, 17, 19, 23|
|4|5, 15, 17, 19|9, 17, 19|
|5|13, 15, 17, 19|9, 15, 17, 23|
|7|15|9, 15|

At exponent 61 the nearest partner below the midpoint is at offset \(-2,213,982,461,498,407,132,682,959,208,946,669\). The nearest above has offset \(+1,485,690,430,693,459,244,528,868,155,770,979\). Both use height 23 and order 137. The lower witness has a six-vertex panel with backward tiles \((0,3),(0,4),(1,3),(1,4),(1,5),(2,5)\) and gaps \((4,127)\). The upper has tiles \((0,2),(0,3),(0,4),(0,5),(1,3),(1,4),(1,5),(2,4),(3,5)\) and gaps \((128,3)\). Each pair of heights is joined by a tiling-complement line. A one-arc tournament reversal is a different operation.

The [replay](line_targets.py) enumerates every anchored strong panel at heights at most 23. Orders through six suffice since the strong minimum at seven is 25 and increases thereafter. Any source at these heights has at most two nonsingleton components. The search ranges over every allowed ordered product and exact gap lengths. Since \(\lambda>9/5\) and \(d\le3\), total gap length at least 150 already exceeds \(2E_{61}\) in opposite height; the reported witnesses lie on both sides within the searched range. Run `python3 -B line_targets.py` from this directory. This is a decision and exclusion for the specified lines. Ordinary realizability of the exponent-61 upper target, a single-arc crossing there, and the full Chan–Pak spectrum remain open. [Existing strong order bounds](../path_spectrum_two_angles_metagraph_20261106/README.md) · [Exponent-61 frontier](../chanpak_odd_square_spines_20261122/README.md)
