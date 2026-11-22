# Absolute height in a prescribed residue cell of two Farey words

The previously proved Gamma correction for one unrestricted two-word bank gives a more sensitive statistic after fixing the low digits of both words. At spine length 12 the relative proportion of E61-compatible pairs modulo 17 converges to \(49/34\) **from below** when pairs are counted up to an increasing absolute height. The first correction is determined explicitly by the two residues of the addresses, including the dyadic digits that a prescribed residue forces. At the actual 121-bit scale the entire canonical tower can cover **fewer than one in 13,258** odd lower targets. These statements concern ordered address pairs and lower endpoints of the marked flip. [Earlier unconditioned Gamma count](../chanpak_viability_shape_20261115/LOCAL_DENSITY.md) · [Height-ordered leading local densities](ABSOLUTE-WORDS.md)

## Two-term law in every physical residue cell

**PROVED.** Fix a chart of the form
\[
 H(x,y)=k+qxy+tx+cy+r[x\ {\rm odd}]+e[y\ {\rm odd}]+h(f(x)+f(y))
\]
with \(x,y\ge0\), positive \(q,c,t,h\), nonnegative \(r,e\), and \(c+h,t+h,r+e\le q\). All eight small tower charts satisfy these inequalities. Put \(a=c/q,\ b=t/q,\ d=h/q,\ S=(Z-k)/q\). Fix \(m\ge1\) and residues \(0\le i,j<m\). Let \(B=2^{v_2(m)}\), and let \(\operatorname{rev}(j)\) reverse the low \(\log_2 B\) bits of \(j\), including leading zeros. For \(B=1\) set the reversal to zero. Define
\[
 \lambda_j=(B-1-\operatorname{rev}(j))/B,\qquad
 A=(i+a+d\lambda_j)/m,\quad C=(j+b+d\lambda_i)/m,\quad D=d/(mB).
\]
Then, as \(Z\to\infty\), the number of address pairs in this cell of height at most \(Z\) is
\[
 \frac{S}{m^2}\left(\log\frac{S}{m^2}+K_{ij}\right)+o(S),\qquad
 K_{ij}=-1+
 \frac{\log\Gamma(A)-\log\Gamma(A+D)+\log\Gamma(C)-\log\Gamma(C+D)}{D}.
\]
In particular the dyadic part of the right residue shifts the vertical parameter \(A\), while the dyadic part of the left residue shifts the horizontal parameter \(C\).

Write \(x=m u+i,\ y=m v+j,\ s=S/m^2\). The whole feasible domain contains the divisor region \((u+1)(v+1)\le s\) and lies in \(uv\le s\). Indeed
\[
 qxy\le H-k\le q(x+1)(y+1),
\]
and \(x+1\le m(u+1)\) on this cell. The extra points relative to the inner divisor region have asymptotic mass concentrated at the two boundaries. More precisely, on \(J<u\le v\) their number is bounded by the sum over \(J<u\le\sqrt s\) of \(s/u-s/(u+1)+2\), which is \(O(s/(J+1)+\sqrt s)\). Exchanging the coordinates proves the same bound on the other half.

On a fixed row \(u\), the inherited binary reading satisfies \(f(y)=y(1-\rho_2(y))+O(\log(y+2))\), where \(\rho_2\) reflects the binary digits after the point. Since \(m/B\) is odd, it permutes residues modulo every power of two. Accordingly \(1-\rho_2(mv+j)\) has limiting distribution \(\lambda_j+B^{-1}U\), with \(U\) uniform on \([0,1]\). This statement also holds over any interval of \(v\) of length tending to infinity. Truncating \(U\) to its first \(w\) digits gives \(2^w\) arithmetic progressions; bounding the slope by the two endpoints on each progression gives matching Riemann sums as \(w\to\infty\). The logarithmic reading error and fixed row offset are \(o(s)\). The row count divided by \(s\) therefore tends to
\[
 \int_0^1\frac{dz}{u+A+Dz}=\frac1D\log\frac{u+A+D}{u+A}.
\]
The inner divisor region contributes \(1/(u+1)\) in the same units. The column calculation exchanges \(A,C\). The overlap of the first \(J\) rows and columns has bounded size for fixed \(J\). Let \(s\to\infty\), then \(J\to\infty\), using the tail bound above. This gives the sum of the two absolutely summable boundary corrections. [Radical-inverse height estimate](../chanpak_radical_inverse_height_20261118/README.md)

The divisor region has \(s\log s+(2\gamma-1)s+o(s)\) points by splitting the hyperbola at \(\sqrt{s}\). The Gamma product and \(H_l-\log l\to\gamma\) give
\[
 \sum_{u\ge0}\left(\frac1D\log\frac{u+A+D}{u+A}-\frac1{u+1}\right)
 =\frac{\log\Gamma(A)-\log\Gamma(A+D)}D-\gamma.
\]
The two boundary terms cancel the two occurrences of Euler's constant in the divisor law. This proves the cell formula.

## What the singular prime 17 measures below the leading term

At \(n=12\) the reading coefficient is divisible by 17, so modulo 17 the lower count is constant on each residue cell modulo 34. Let \(\mathcal C\) be the cells matching \(2^{60}(2^{61}-1)-1\bmod17\). Direct reduction of the seven integers gives 98 cells among 1,156, hence the leading density \(98/1156=(49/34)/17\). Put \(\bar K_{\mathcal C}=|\mathcal C|^{-1}\sum_{\mathcal C}K_{ij}-2\log34\) and use \(K_0\) for the cell \(m=1\). The proportion of actual word pairs up to height \(Z\) meeting that fixed residue is
\[
 \frac{49}{578}\left(1+\frac{\bar K_{\mathcal C}-K_0}{\log S}
                       +o(1/\log S)\right).
\]
Numerically \(\bar K_{\mathcal C}=1.187465\ldots\) and \(K_0=3.007980\ldots\). The sign also follows without numerical Gamma evaluation. For \(v,D>0\) the average of \(-\psi\) over \([v,v+D]\) lies between its values at the two endpoints because \(\psi\) is increasing. For \(z>0\), bound \(\psi(z)\) by the first 256 terms of
\[
 \psi(z)=-\gamma+\sum_{l\ge0}\left(\frac{1}{l+1}-\frac{1}{l+z}\right).
\]
Here \(0<z<1\) for every compared cell endpoint, so the remaining tail lies in \([-(1-z)/(255+z),0]\) by comparison to the integral of \(t^{-2}\). Rational interval arithmetic in the replay bounds the difference between \(-1.9\) and \(-1.7\). The common \(2\gamma-1\) cancels in this comparison. These are asymptotic second terms and give no certified numerical error at the specific E61 height. [Local 17 calculation](LOCAL-WORDS.md)

## A finite shortage at the perfect-number scale

**PROVED by integer bounds.** Among the \(2^{119}\) positive odd integers in \([2^{120},2^{121})\), at most **50,128,304,294,132,497,289,253,221,463,461** are lower endpoints of this canonical tower, with arbitrary addresses and every positive spine length divisible by four.

In particular some consecutive block of **13,258 odd integers** in this 121-bit band contains no lower endpoint at all. If there are at most \(U\) occupied places among \(M=2^{119}\), the remaining places form at most \(U+1\) consecutive gaps. Their maximum length is at least \(\lceil(M-U)/(U+1)\rceil=13{,}258\) for the displayed bound. Every upper endpoint is two above its lower one, so at least 13,257 consecutive odd integers in this block miss both endpoints of every tower flip. The block's location within the band remains undetermined.

For an upper height \(Z=2^{121}-1\), the two physical axes allow at most \(1+\lfloor(Z-k)/t\rfloor+\lfloor(Z-k)/c\rfloor\) presentations. The positive interior allows at most the ordinary divisor sum \(\mathcal D(s)=\sum_{i=1}^s\lfloor s/i\rfloor\) for \(s=\lfloor(Z-k)/q\rfloor\). Subtract at least \(\mathcal D(u)\) pairs below the band, where \(u=\max(0,\lfloor(2^{120}-1-k)/q\rfloor)\), by the included region \((x+1)(y+1)\le u\). The elementary bounds \(v(\log v-1)\le\mathcal D(v)\le v(1+\log v)\) for \(v\ge1\) follow by bounding each floor and using integral bounds for the harmonic sum. In the replay, \(0.69<\log2<0.70\) and integer bit lengths suffice. Sum through length 68. The bare height at 72 already exceeds the band, and the tower is nested. Counting every presentation gives the asserted upper bound on its distinct endpoints. [Nested tower](FERMAT-RECEIVER-SEARCH.md) · [Exact bound and phase calculation](farey_census.py)

~~~sh
python3 -B farey_census.py
~~~
