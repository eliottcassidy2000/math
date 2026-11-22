# A target-selected lattice search for two Farey flanks

**FINITE-EXACT.** In the canonical tagged tower with left tag mask 12, no pair of binary flank addresses at spine lengths **20, 24, 28 or 32** realizes the lower count \(N=2^{60}(2^{61}-1)-1\). The other orientation always adds two. Each exclusion covers independent arbitrary words on both sides, including either or both empty. All lengths divisible by four from 36 onward were already excluded. [Inherited tower definition and large-length bound](FERMAT-RECEIVER-SEARCH.md)

## Removing a congruence before choosing either flank

The lower chart from the seven slot readings is
\[
 N=k+t x+c y+qxy+r(x\bmod2)+e(y\bmod2)+h\bigl(f(x)+f(y)\bigr).
\]
Write \(g=\gcd(h,c,q,e)\). Regardless of the second flank, \(N\equiv k+tx+r(x\bmod2)\pmod g\). For this target every admissible left word at length 12 is \(141\) or \(170\) modulo 290. For every other tower length from 4 through 32 it is \(0\) or \(1\) modulo 10. Here \(g=145\) at length 12 and 5 otherwise.

Thus replace \(x\) by \(2gX+b\) for each allowed \(b\). More generally write the two addresses as \(s_x X+o_x,s_y Y+o_y\), with one scale 1 and zero offset. Expansion leaves positive coefficients \((Q,C,T,K)\) on \(XY,Y,X,1\). An even scale absorbs a fixed parity correction into \(K\). The unknown nonlinear readings still have coefficient \(h\). Choose the smaller of \(X,Y\) first. Its bound is \(\lfloor\sqrt{(N-K)/Q}\rfloor\), and the second orientation excludes the diagonal to count a solution only once. After normalization these bounds at lengths 20, 24, 28 and 32 are approximately 61 billion, 6.8 billion, 753 million and 84 million for each surviving residue. Both raw words are reconstructed before testing the chart.

## Why thin hyperbola strips can be decoded together

Put \(D=(N-K)/Q,\ \alpha=C/Q,\ \beta=T/Q\). Both offsets lie between 0 and 1 in the normalized charts. Dropping the nonlinear terms and the remaining parity corrections gives the exact real upper boundary
\[
 Y_0(X)=\frac{D+\alpha\beta}{X+\alpha}-\beta.
\]
For \(X\ge L\ge1\), \(Y\le D/X\) and \(0\le f(z)\le z\). Consequently a solution lies in \(Y_0(X)-G\le Y\le Y_0(X)\), where
\[
 G=\frac{h}{Q}\left(s_x+s_yD/L^2\right)+
       \frac{h(o_x+o_y)+r+e}{QL}.
\]
Parity coefficients already absorbed in \(K\) can be dropped from this bound.

On a block \(X_0,\ldots,X_0+W-1\), put \(\lambda=-Y_0'(X_0)\). Pick a reduced rational \(a/b\) near \(\lambda\), with \(b\le W\). For any integer anchor \(y_*\) near \(Y_0(X_0)\), set \(F=Y_0(X_0)-y_*\). Every solution has an *integer* line index
\[
 z=b(Y-y_*)+a d,\qquad d=X-X_0,
\]
in the interval
\[
 b(F-G)+\min(0,(a-b\lambda)W)
 \le z \le b(F+B)+\max(0,(a-b\lambda)W),\qquad
 B=(D+1)W^2/L^3 .
\]
Indeed \(Y_0\) is convex with second derivative at most \(2(D+1)/L^3\). On any surviving line, \(ad\equiv z\bmod b\) fixes the progression of offsets and the line determines \(Y\) exactly. The code takes \(W\) of the order \(L/D^{1/3}\) and sweeps adjacent fractions of the Farey sequence of order \(W\) as the slope decreases. Its neighboring fractions also give \(a^{-1}\bmod b\) without another Euclidean algorithm. The outer receiver uses a determinant of tagged path counts. The second Farey choice uses curvature at the prescribed target. Any reduced rational satisfies the interval proof.

Near the axes the code uses congruences first. A divisor of \(h\) removes the two nonlinear readings. At fixed address residue and both parities this gives one arithmetic progression for the opposite address or eliminates that branch. For each surviving outer word the interval above rejects an empty progression near \(Y_0\). For longer opposite words the exact integer stage splits off \(b\) low bits, squeezes the remaining quotient by \(0\le f(z)\le z\), and intersects it with the necessary slope congruence modulo \(h\). It uses the [low-word concatenation identity](../chanpak_e61_fibonacci_target_20261122/SUFFIX-INVERSE.md). The fixed-address step always tests the literal exact height of a candidate.

The floating-point operations enclose the geometric strips. All reported heights use unsigned 128-bit integer arithmetic. For the vectorized row test the absolute rounding allowance is \(2^{-47}(D/(mL)+2)\) for a progression of spacing \(m\). The omitted shift \(\alpha\beta/(X+\alpha)\) gets an additional \(2/L\) in unscaled coordinates. For the slanted lines the anchor residual numerator \(N-K-TX_0-(QX_0+C)y_*\) is formed as a signed integer before division. With IEEE binary64 unit rounding error \(2^{-53}\), converting and dividing the positive chart entries costs at most six units of relative error. The slope \((y_*+F+\beta)/(X_0+\alpha)\) uses another two additions and a division. The possible cancellation in \(a-b\lambda\) is charged in **absolute** error proportional to \(a+b\lambda\). The implementation uses twice the positive bounds on curvature and thickness, and expands each line interval by \(2^{-46}[1+b(1+|F|+B'+G')+W(a+b\lambda)]\) using the enlarged bounds \(B',G'\). Thus neither rounding nor an inaccurate initial choice of the nearest rational can discard a solution. The searched normalized small coordinates are below \(2^{52}\). The tangent stage caps \(D/L\) at \(3\cdot10^{15}\).

## Reproduction

From this directory:

~~~sh
python3 -B hyperbola_controls.py > /tmp/farey-e61-charts
c++ -O3 -std=c++17 farey_hyperbola.cpp -o /tmp/farey-e61-hyperbola
python3 -B hyperbola_controls.py --check /tmp/farey-e61-hyperbola
for n in 20 24 28 32; do /tmp/farey-e61-hyperbola "$n" < /tmp/farey-e61-charts; done
~~~

The controls compare the complete inverse with direct pair enumeration on small planted targets and test long nonconstant flanks in both orientations, including a pair in the tangent regime. For another lower height pass its decimal integer after the spine length. Extra arguments select direction (0 or 1), inclusive normalized small-coordinate range, and the allowed physical left residue. Each returned pair is checked against the integer chart. [Exact inverse](farey_hyperbola.cpp) · [Generator and controls](hyperbola_controls.py)
