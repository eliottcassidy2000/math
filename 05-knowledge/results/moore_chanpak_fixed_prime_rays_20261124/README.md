# The coupled pentagon: fixed prime support and exact padding rays

**Scope.** This is about the specific coupled-pentagon seam from the [nontransitive-filler session](../moore_chanpak_nontransitive_seams_20261124/), not all tournament reversals. The two fillers have profiles \((h,Q,R)=(Z_1,Z_2,Z_3)\), with \(Z_k\) counting ordered \(k\)-path covers including empty paths. The signed gain is
\[
 \Delta=2h_Xh_Y+5Q_Xh_Y-5h_XQ_Y+4R_XQ_Y-4Q_XR_Y .
\]
Ordered joins multiply each profile coordinate. In particular, adjoining \(a\) ordinal singleton vertices changes \((h,Q,R)\) to \((h,2^aQ,3^aR)\), irrespective of where the singletons sit among strongly connected components.

## Fixed prime alphabet theorem

Fix a finite prime set \(S\) and a nonzero integer \(d\). Restrict to fillers for which all six positive profile coordinates have prime factors in \(S\). Apart from finitely many pairs, any pair with \(\Delta=d\) satisfies
\[
 d=2h_Xh_Y,\qquad
 \theta_X=\theta_Y,\qquad \theta_T=(4R_T-5h_T)/Q_T.
\]
Moreover all pairs with these two equalities come from finitely many base profile pairs by simultaneous ordinal padding, except for a finite exceptional set; the base pairs are effectively computable from \(d>0\), without any prime restriction. Precise meaning: each ray has profiles
\[
 (h_X,2^m Q_X,3^m R_X),\quad (h_Y,2^m Q_Y,3^m R_Y)
 \qquad(m=0,1,\ldots)
\]
with the same \(m\). A ray predicts profile values, not unique isomorphism classes of the padded tournaments. If one of \(2,3\) is absent from \(S\), even these rays contribute at most finitely many pairs in the prime-restricted chart.

For \(d=-2\), there are only finitely many solutions. For \(d=2\), the only infinite family that can remain is two equal-length transitive fillers. Thus **any infinite family of distinct unit pairs involving a nontransitive filler must use unbounded prime factors in its slot profiles**. In particular, choosing both fillers as arbitrary ordinal joins of a fixed finite catalog, even with unlimited singleton padding, cannot do it. This does not require the new primes to occur in \(h\) rather than \(Q\) or \(R\), nor does it cover substitutions that introduce new primes. Even on a ray from two nontransitive cores the gain can stay constant but need not be a unit: two equally padded directed triangles have gain \(2\cdot3^2=18\).

### Proof of finiteness away from the resonance

Append the term \(-d\) to the five displayed gain terms. Choose a minimal zero subsum **containing** \(-d\). None of the other selected terms has a nonempty vanishing subsum: otherwise it could be deleted without deleting \(-d\). Divide by \(d\). This is a nondegenerate linear equation in a fixed finite-rank multiplicative group, namely a Cartesian power of the positive \(S\)-units. The coefficients are fixed rationals, so their prime factors need not be in \(S\). The [Evertse–Schlickewei–Schmidt theorem](https://annals.math.princeton.edu/2002/155-3/p04) gives finitely many tuples for each choice of terms.

If the only selected gain term is \(2h_Xh_Y\), then \(d=2h_Xh_Y\); the remaining gain is \(Q_XQ_Y(\theta_X-\theta_Y)\), so \(\theta_X=\theta_Y\). Otherwise at least one bounded product includes \(Q\) or \(R\) from one filler. Since \(Z_k\ge k^n\) for order \(n\) (every coloring has an ordered path cover), this bounds that filler's order. A fixed filler has only finitely many partners for a fixed \(d\), by the following estimate.

Indeed for a filler \(Y\) of order \(n\),
\[
 Q_Y\ge(n+1)h_Y,\qquad R_Y\ge(1+n/2)Q_Y .
\]
For the first bound split a single Hamilton path at any prefix boundary. For the second split the first path of an ordered two-cover, whose mean length is \(n/2\) by exchanging the two colors. Both maps are injective, because concatenation recovers the split path. The identity \(\Delta=2h_Xh_Y+Q_XQ_Y(\theta_X-\theta_Y)\) then gives
\[
 2n+4\le \theta_X+
 \frac{5+(|d|+2h_X)/Q_X}{n+1}.
\]
The right side decreases in \(n\), giving an explicit finite order bound with \(X\) fixed. The inequality also holds when the names are exchanged, because the estimate uses absolute values. The \(S\)-unit step itself gives finiteness, not an effective order bound for all off-ray pairs.

## Effective description when both ordinary heights are bounded

This part does not require an \(S\)-unit theorem and includes \(d=0\). Fix a common height bound \(B\). Remove all singleton strong components from either filler. The remaining core has order at most \(B\): a strong tournament of order \(j\ge3\) has at least \(j\) Hamilton paths, by cutting its Hamilton cycle at all \(j\) locations, and the height is the product over its ordered components. A product of two or more integers each at least three is at least their sum; with one component the bound is direct. Thus enumerate a finite set of possible cores. Include the empty core with profile \((1,1,1)\), and remember that a resulting filler must be nonempty.

For fixed cores \(X,Y\), let \(a,b\) count the removed singleton components to restore. Set \(A=R_X/Q_X\), \(B'=R_Y/Q_Y\), and \(r=3/2\). Dividing the gain equation by \(4Q_XQ_Y2^{a+b}\) gives
\[
 |A r^a-B' r^b|\le K,\qquad
 K=\frac{|d-2h_Xh_Y|+5Q_Xh_Y+5h_XQ_Y}{4Q_XQ_Y}.
\]
If \(a\ge b\), this implies \(A r^{a-b}\le B'+K\). If \(b\ge a\), it implies \(B'r^{b-a}\le A+K\). Hence an exact rational loop bounds the finite set of possible differences \(\delta=a-b\).

For each \(\delta\), put \(a_0=\max(\delta,0)\), \(b_0=\max(-\delta,0)\), and \((a,b)=(a_0+m,b_0+m)\). Direct substitution gives
\[
 d-2h_Xh_Y=L\,2^m+M\,6^m
\]
with
\[
 L=5(2^{a_0}Q_X h_Y-h_X 2^{b_0}Q_Y),\qquad
 M=4(3^{a_0}R_X2^{b_0}Q_Y-2^{a_0}Q_X3^{b_0}R_Y).
\]
If the left side is nonzero, its 2-adic valuation bounds \(m\); test all the resulting nonnegative integers. If it is zero and \(M\ne0\), test the one possible power \(3^m=-L/M\). If it is zero and \(M=0\), there are no solutions unless \(L=0\), in which case every \(m\) with nonempty fillers works.

The ray conditions \(L=M=0\) say exactly that after offset \((a_0,b_0)\) the two fillers share \(Q/h\) and \(R/h\). These normalized readings, not just a coincidence of \(\theta\) at one width, are what preserve the resonance under indefinite simultaneous growth. When \(d=2h_Xh_Y>0\), both heights are automatically bounded by \(d/2\), so the effective classification applies to the exceptional possibility in the \(S\)-unit argument. When \(d=2\), both heights are one, hence both tournaments are transitive. Their \(Q/h=2^n\) forces equality of lengths on any ray; in fact the strict transitive width comparator from the original session gives the same equality for every transitive unit pair.

**Consequence for the Chan–Pak search.** A fixed finite stock of filler profiles plus multiplication is now obstructed even if both ports grow and all matching is adaptive. For this seam, a nontransitive unit mechanism must generate new profile primes and cannot keep both ordinary filler heights bounded. This is a constraint on a proposed exact selector, not a solution to Chan–Pak or a restriction on all possible tournament charts.
