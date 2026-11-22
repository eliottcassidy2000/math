# A compatible second count does not always exist in the backchain

**Proved.** The whole [fixed-backedge chain](../chanpak_absolute_descent_20261122/BACKCHAIN.md), with arbitrary gap widths and arbitrary runs of zero gaps, has only \(O(X^{23/24})\) distinct ordinary heights at most \(X\). The same upper bound holds even after taking all finite ordinal products of these heights. Thus a density-one set of odd targets has no compatible second count in this architecture, already without any restriction on the endpoint order. This neither excludes another tournament witness nor refutes Chan–Pak.

## Reset is a genuine multiplicative decomposition

Use the nonnegative state \((D,U,V)\), height \(F=D+U+V\), initial state \((0,0,1)\), and updates
\[
 Z(D,U,V)=(D+2U,V,D+U),\qquad
 P_m(D,U,V)=(2^m(D+2U+V),0,D+U+V)\quad(m\ge1).
\]
After a positive gap \(U=0\), \(V\le F/3\), since the update has \(D'\ge2V'\). If another positive gap immediately follows, its new state is
\[
 P_m(D,0,V)=F(2^m,0,1)=F\,P_m(0,0,1).
\]
Every suffix then applies linear maps. Thus every adjacent-positive boundary splits off the earlier height as an exact factor of the final height. This is a count decomposition at a shared spine endpoint, not an assertion that the tournament is an ordinal sum. A leading zero can be normalized: \(Z(0,0,1)=(0,1,0)\); the next \(Z\) is the same state as the fresh \(P_1\), and a next \(P_m\) is the same state as the fresh \(P_{m+1}\). If there is no next step the height is 1.

Terminal zero runs can be removed **without changing the height or order** for nontrivial outputs. On any state the final scalar after \(ZZ\) is \(3D+5U+3V\), the same as after \(P_1\); replace the final two zeros with that one positive gap of width one. A lone final \(Z\) after \(P_m\) can be absorbed by replacing \(P_m\) with \(P_{m+1}\): both give \(F+2^{m+1}(F+U)\) from the same preceding state. These are scalar identities, not equalities of full state. Each replacement adds the same number of vertices as the original suffix.

Call a word an atom when it starts and ends with a positive gap and has at least one zero between any two positives. Let \(A\) be its height set. The exact scalar language is the multiplicative monoid \(\langle A\rangle\): normalization and adjacent-positive resets give the decomposition, and any finite product is realized by concatenating the atom words. The empty product is the height 1.

## A real-place pressure bound, not a congruence obstruction

At a positive boundary, divide the state by \(F\). It lies on the closed segment from \(a=(1,0,0)\) to \(b=(2/3,0,1/3)\). A block of \(j\ge1\) zeros followed by a positive \(m\ge1\) multiplies height by at least
\[
 r_{j,m}=\min_{t\in\{a,b\}}\bigl[f_j(t)+2^m(f_j(t)+u_j(t))\bigr],
 \qquad f_j(t)=\boldsymbol 1^\top Z^jt,\quad u_j(t)=(Z^jt)_2 ,
\]
where the second coordinate is \(U\). The minimum over the segment lies at an endpoint because this expression is linear. For \(s=23/24\), set \(R_{j,m}=r_{j,m}^{-s}\). The following strict contraction suffices:
\[
 \sum_{j,m\ge1} R_{j,m}<0.999<1.
\]
Here are explicit bounds, with a standard-library integer replay. The entrywise inequality \(Z^3\ge2I\) implies \(R_{j+3,m}\le(3/5)R_{j,m}\), since \(2^{-s}<3/5\). For \(m\ge1\), increasing \(m\) multiplies each endpoint growth by at least \(5/3\): the powered term is at least twice the unpowered term. Hence \(R_{j,m+1}\le(2/3)R_{j,m}\), since \((5/3)^{-s}<2/3\). If
\[
 B_j=\sum_{m=1}^{20} R_{j,m}+2R_{j,20},
\]
the full sum is bounded above by \(\sum_{j=1}^{20}B_j+\tfrac32(B_{18}+B_{19}+B_{20})\). Each finite entry can be rounded upward to a multiple \(k/10^6\) by an integer comparison: if \(r_{j,m}=c/3\), test \(k^{24}c^{23}\ge (10^6)^{24}3^{23}\). The resulting rational upper bound is less than \(999/1000\); it needs neither a floating-point calculation nor truncation without a tail estimate.

A first positive gap contributes a summable series of weights \((1+2^m)^{-s}\). Each additional positive, with its nonempty preceding zero run, has total relative weight bounded by the contraction above, independently of the current state. Consequently the sum of \(F^{-s}\) over all atom presentations is finite, and so is \(\sum_{a\in A}a^{-s}\). This already proves \(O(X^s)\) for atom outputs.

Finally all atoms have height at least 3. The Euler product \(\prod_{a\in A}(1-a^{-s})^{-1}\) converges because \(\sum_a a^{-s}<\infty\). It counts each factor multiset, so the sum of \(n^{-s}\) over the ordinal-product outputs \(n\) is bounded by this product even when distinct multisets coincide. This proves \(O(X^s)\) for the full architecture. In particular neither extra zero gaps, freely selecting the terminal register, nor taking arbitrary ordinal products can make the backchain cofinite.

## Three extra families of genuine parts still do not suffice

In the [generalized shared-gap chain](SHARED-CONTINUANTS.md), a nonempty gap with profile \((h,Q)\) instead updates to \((Q(F+U),0,hF)\). Enlarge the palette of transitive gaps by allowing any of the profiles \((5,32),(11,100),(19,256)\), with any number of ordinal singleton padding vertices in each gap. Their strong seed row masks are \((2,8,3,5)\), \((18,4,9,19,6)\), and \((2,16,11,19,37,15)\), respectively. The replay recounts their profiles. Thus the four lanes of allowed profiles are
\[
 (1,2\cdot 2^k),\ (5,32\cdot2^k),\ (11,100\cdot2^k),\
 (19,256\cdot2^k)\quad(k\ge0).
\]
Padding preserves \(h\) and doubles \(Q\). The pressure proof applies unchanged with each step growth replaced by the endpoint minimum of \(h f_j+Q(f_j+u_j)\); the positive-boundary cone is preserved because \(Q\ge2h\) on these lanes. The adjacent positive resets are still homogeneous. A terminal \(ZZ\) can still be replaced by the transitive width-one gap, and a lone terminal \(Z\) after a nonempty module is absorbed by doubling its \(Q\) via padding.

With \(s=199/200\), over all four lanes the same two geometric tail bounds and an exact integer rounding for \(j=1,\ldots,20,\ k=0,\ldots,19\) give total pressure below one. Thus even this palette, with infinitely many choices of physical gap modules, only produces \(O(X^{199/200})\) heights, also after arbitrary ordinal multiplication. This is a test for specific new parts, not a universal obstruction for unrestricted modules.

## All gaps assembled from strong pieces through order five still miss most primes

Now permit an arbitrary ordinal sum of strong tournaments of order at most five inside each gap, not just one nontrivial seed. A prime final height forces every used gap to obey \(\gcd(h,Q)=1\): after that positive step any common factor divides every state coordinate and stays present, while the output already exceeds that factor. Also the word must be a single atom, since adjacent nonempty gaps factor the height into two integers greater than one. The same terminal normalization works: padding doubles \(Q\).

Complete strong profiles through order five are
\[
 (3,12),(5,32),(9,84),(11,100),(13,104),(15,108),(15,120).
\]
Only \((5,32)\) and \((11,100)\) have coprime components. Moreover they cannot be used together in the same gap, since \(5\mid100\). The available prime-bearing gap profiles are therefore either nonempty transitive, or \((5^i,32^i2^k)\), or \((11^i,100^i2^k)\), with \(i\ge1,\ k\ge0\). At fixed zero length and padding, raising \(i\) by one multiplies the lower endpoint growth by at least 23 in the first lane and 70 in the second. Indeed \(Q(f+u)\ge2hf\), so weighted growth is at least \((5+2\cdot32)/3=23\), or \((11+2\cdot100)/3>70\). At \(s=199/200\), the negative powers of these bounds are below \(1/20\) and \(1/60\), respectively. Thus multiply each first-power lane's pressure bound by \(20/19\) or \(60/59\). The replay checks that adding these to the transitive lane still gives total pressure below one. It follows that this larger grammar hits only \(O(X^{199/200})\) prime heights up to \(X\), whereas there are on the order of \(X/\log X\) primes. This is a compositional obstruction based on actual strong subparts, not merely a modulus test on the total count.
