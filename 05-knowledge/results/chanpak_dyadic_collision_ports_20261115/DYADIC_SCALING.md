# Nearby after opposite dyadic scaling

There is an analogous power bound after arbitrary opposite dyadic scaling of the axes. This treats, for instance, pairs whose original long coordinate differs by a factor of two. It is still an estimate on a restricted collision set.

Let \(F,\tau,\delta,W\) be as in [the local bound](COLLISIONS.md). In the central shell consider pairs with
$$
 |x-\lambda x'|\le W,\qquad |\lambda y-y'|\le W
 \quad\text{for some power of two }\lambda\ge1.
$$
Exclude the exact transfer pairs \(x=\lambda x'\), \(y'=\lambda(y+1)-1\), whose total energy is already counted in [the orbit note](MASK_ZERO_ORBITS.md). Then the ordered collision count is
$$
 O_\delta\left( (\sqrt T W^\tau+W^4+T^{1-\delta})\log^4(2T)\right).
$$
Thus these pairs contribute \(o(T)\) when \(W\le T^{1/(2\tau)-\epsilon}\), apart from the known orbit itself. The same estimate follows with \(\lambda\) in the other direction by exchanging primed and unprimed addresses. Both transposed versions satisfy the same estimate, as proved in the last section.

## Normalize the axes and classify singular lines

First suppose \(\lambda\le W\). Write \(x=\lambda X+U,y'=\lambda Y+V\) with \(0\le U,V<\lambda\), and \(Z=W/\lambda\ge1\). Re-use the common-prefix splits of [the local theorem](COLLISIONS.md) on the pairs \((X,x')\), \((y,Y)\):
$$
 (X,x')=(Au+B,Au+B'),\qquad (y,Y)=(Cv+D,Cv+D').
$$
Put \(h=B'-B,k=D'-D\), and let \(d_x,d_y\) be the reversed differences as defined there. The slopes of \(F(\lambda X+U,y)-F(x',\lambda Y+V)\) as an affine form in \(u,v\), on fixed prefix parities, are
$$
 \alpha=A[-100\lambda k+100(\lambda-V-1)-11r_\lambda(U)]+11d_x,
$$
$$
 \beta=C[-100\lambda h+100U-11R_\lambda(V)]+11d_y.
$$
Here \(r_\lambda,R_\lambda\) are reversal and complementary reversal. When \(\lambda=1\) the corresponding bounded parity terms have their literal zero-width interpretation.

If \(\alpha=0\), reduction modulo \(A\) gives \(d_x=0\), implying \(h=0\); the remaining bracket equation gives \(k=0\), since \(100(\lambda-V-1)-11r_\lambda(U)\) lies strictly between \(-100\lambda\) and \(100\lambda\). Similarly \(\beta=0\) implies both normalized shifts are zero. If both vanish, the [parallel-cylinder classification](MASK_ZERO_ORBITS.md) forces the exact transfer, excluded here. If just one vanishes, there are at most \(\lambda\) possible residual pairs \(U,V\) for its equation, because either suffix determines the other. The full equality, viewed in \(X=x',Y=y\), fixes the coordinate whose slope is nonzero up to four choices (split the actual parities). If \(y\) is fixed, then \(x'\ll T/(\lambda T^\delta)\) from \(y'\ge\lambda y\); if \(x'\) is fixed use \(x\ge\lambda x'\) for the same bound on \(y\). Therefore the singular-line count is \(O(T^{1-\delta})\) at each width.

## Nonzero slopes obey the same moment bound

The two normalized products are comparable to \(T/\lambda\), first with possible unit shifts \(U/\lambda,V/\lambda\), which do not affect fixed-factor comparisons. Their longer compared coordinate is \(\gg\sqrt{T/\lambda}\); its displacement \(O(Z)\) is \(o(\sqrt{T/\lambda})\), so it is comparable in the two pairs. The product annulus then makes the shorter coordinate comparable as well. Each coordinate is bounded below by a constant times \(T^\delta\), since one compared coordinate in each axis was unscaled and central.

Take dyadic scale parameters \(L_{\rm short},L_{\rm long}\) for this normalized rectangle and \(\rho=L_{\rm long}/L_{\rm short}\ge1\). Dividing the collision identity by \(\lambda\), the lower-order terms are \(O(L_{\rm short}+L_{\rm long}+1)\) because \(0\le p(z)\le z\) and \(U,V<\lambda\). Thus the normalized shorter displacement is \(O(1+Z/\rho)\). The longer displacement is at most \(Z+1\) by hypothesis. As before take power-of-two base radices of those orders, and classify the carry types. In a fixed pair of carry types and fixed \(U,V\), the slopes still determine the two displacement and reversal signatures uniquely: the affine coefficients are spaced by \(100\lambda A\) or \(100\lambda C\), larger than the full variation of the corresponding reversal terms.

For \(1\le\rho\le Z\), redo the primitive-slope summation from [the local proof](COLLISIONS.md) after scaling the allowable slope magnitudes by \(\lambda\): the mixed coefficients can now range over \(O(\lambda AZ)\), \(O(\lambda CZ/\rho)\) when the first axis is shorter. The upper bound \(O(Z/M)\) for the number of possible common divisors becomes \(O(\lambda Z/M)\); nothing else in the calculation changes. This adds a factor of \(\lambda^{1/Q}\) to the Hölder estimate. If the second axis is shorter, exchange the two coordinate labels in this estimate. Sum over \(U,V\) (at most \(\lambda^2\) choices), carry types and normalized dyadic rectangles. The nonconstant lattice term is at most
$$
 O\!\left(\lambda^{2+1/Q}\sqrt{T/\lambda}\, Z^\tau\log(2T)\right)
   = O\!\left(\sqrt T W^\tau \lambda^{3/2+1/Q-\tau}\log(2T)\right),
$$
with power of \(\lambda\) strictly negative. The additive line term is \(O(\lambda^2 Z^4\log^3(2T))=O(W^4\lambda^{-2}\log^3(2T))\).

If \(\rho>Z\), take the shorter radix \(O(1)\). Simply count line solutions by first fixing the shorter common prefix, giving \(O(1+L_{\rm short})\) choices for each of the \(O(Z^2)\) suffix choices per carry pair. Including \(\lambda^2\) residual pairs this is
$$
 O(W^2(1+L_{\rm short})\log^3(2T))
 \ll (W^2+\sqrt T W^{3/2})\log^3(2T),
$$
since \(L_{\rm short}\ll\sqrt{T/(\lambda Z)}=\sqrt{T/W}\). Summing all cases over the permitted dyadic widths proves the estimate.


## The transposed comparisons have no exceptional two-slope orbit

The estimate also holds for \(|x-\lambda y'|,|\lambda y-x'|\le W\) and for \(|\lambda x-y'|,|y-\lambda x'|\le W\). Use the same normalized comparisons \((X,x_1)\), \((y_0,Y)\) with low differences \(h=x_1-X,k=Y-y_0\) and prefix radices \(A,C\). The algebra differs only by its fixed slope offsets:
- For \(F(\lambda X+U,y_0)-F(\lambda Y+V,x_1)\) the slopes are
  \(\alpha=A[-100\lambda k+100(\lambda-V)-11(1+r_\lambda(U))]+11d_x\)
  and
  \(\beta=C[-100\lambda h+100(U-\lambda)+11(1+r_\lambda(V))]+11d_y\).
- For \(F(X,\lambda y_0+U)-F(Y,\lambda x_1+V)\) the slopes are
  \(\alpha=A[-100\lambda k+100U+89-11R_\lambda(V)]+11d_x\)
  and
  \(\beta=C[-100\lambda h-100V-89+11R_\lambda(U)]+11d_y\).

In both lists the bracket's fixed offset has absolute value less than \(100\lambda\). Thus one zero slope still forces both normalized shifts to vanish, and there are at most \(\lambda\) residual combinations solving its zero equation, by determining the ordinary suffix from the one under reversal. In the two-zero case the equations are exactly the transposed matches ruled out in [the orbit classification](MASK_ZERO_ORBITS.md). There is no exceptional orbit of a transposed comparison. Counting the remaining single-zero lines works exactly as before, since the two products are at least \(\lambda X y_0\) and both compared unscaled axes are central. The entire nonzero-slope lattice proof is unchanged.


## Larger widths require fewer residual suffixes

The restriction \(\lambda\le W\) is unnecessary. If \(\lambda>W\), each compared residual suffix belongs to \([0,W]\cup[\lambda-W,\lambda)\), since the comparison is within \(W\) of a multiple of \(\lambda\). There are \(O(W^2)\) residual pairs. In normalized coordinates each displacement is at most one; a base radix of two on each axis suffices. The normalized longer coordinate is still large: one unscaled coordinate in each pair is at least \(T^\delta\), and the near comparisons imply this for the scaled-down partners up to a unit. In particular the possible widths are \(O(T)\), and the number of dyadic widths is \(O(\log(2T))\). For each residual pair and carry type with both slopes nonzero, fixing the shorter common prefix leaves at most \(O(1+\sqrt{T/\lambda})\) points, as the normalized products are comparable to \(T/\lambda\). Hence the nonzero-slope contribution at one such width is at most
$$
 O\bigl(W^2(1+\sqrt{T/\lambda})\log^3(2T)\bigr)
 \ \ll (W^2+\sqrt T W^{3/2})\log^3(2T).
$$
Single-zero degeneracies have at most \(O(W)\) residual choices by the same coefficient equations (one suffix is uniquely determined by the other), so their bound is \(O(WT^{1-\delta}/\lambda)\), again \(O(T^{1-\delta})\). The previous two-zero classification still applies. This proves the same estimate over all positive dyadic scales, with transposition also allowed.

