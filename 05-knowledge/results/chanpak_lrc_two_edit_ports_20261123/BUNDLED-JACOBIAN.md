# Four bundled orientation ports in planar constant-Jacobian charts

**PROVED** as a reduction to established degree theorems, with a direct star inverse below. This extends the [two-edit planar charts](UNRESTRICTED-JACOBIAN.md) to a control variable repeated on arbitrarily many arcs. Polynomial invertibility does not make a Boolean output admissible.

Let \(F=(P,Q)\) be two expected Hamilton-path readings (possibly on different hosts), with independent orientation weights polynomial in \(x,y\). The probability for the reverse orientation is the complementary polynomial. Suppose every weight depending on \(x\) is affine in \(x\) with coefficients in \(k[y]\), over a field \(k\) of characteristic zero. In each host write \(G_x\) for the undirected graph of these arcs. The product corresponding to a path has \(x\)-degree at most the number of its edges in \(G_x\). Such edges form a linear forest; if the matching number of \(G_x\) is at most four, every such forest has at most eight edges (alternate edges on each component path give a matching). This includes all arc bundles incident to at most four vertices, since any matching then uses distinct covered vertices. No restriction on \(y\) or its overlap with \(x\) is needed.

**Partial-degree conclusion.** If \(J(F)\) is a nonzero constant and \(\max(\deg_x P,\deg_x Q)\le8\), \(F\) is a polynomial automorphism. A coordinate independent of \(x\) is affine nonconstant in \(y\), and the other is affine in \(x\) with constant slope, by the determinant equation. Otherwise use [Moskowicz, Theorems 1.1 and 2.7](https://arxiv.org/pdf/1810.08202): Theorem 2.7 accepts a positive partial degree \(r\) if the gcd of \(r\) and the degree in \(y\) of its leading \(x\)-coefficient is 1, 4, or prime. Up to 8, the only partial degrees not automatically handled are 6 and 8.

For equal positive degrees \(m=n\ge2\), the leading coefficients \(a,b\in k[y]\) obey \(a b'-a'b=0\) from the coefficient of \(x^{2m-1}\) in \(J(F)\). Thus \(b/a\in k\), and an invertible linear output subtraction lowers one degree. The only possible exceptional pair left is \(6,8\) (in either order).

More generally, for distinct positive \(x\)-degrees \(m,n\) with \(m+n>1\), comparison at \(x^{m+n-1}\) gives \(m a b'-n a'b=0\). If \(g=\gcd(m,n)\), comparison of leading degrees in \(y\) gives \(\deg a=(m/g)c,\ \deg b=(n/g)c\) for an integer \(c\ge0\), including the case of both coefficients constant. Precompose by \(x\mapsto x+y^N\) with \(N\) large enough that the leading \(x\)-terms dominate all other terms in total degree. The new total degrees have gcd
\[
 \gcd(mN+\deg a,\ nN+\deg b)=gN+c .
\]
For \(g=1\) choose this number prime. For \(g=2\) choose it an odd prime if \(c\) is odd, or twice a prime otherwise; there are arbitrarily large choices in either case. The total-degree theorem now implies invertibility, including for \(6,8\). This additional shear is why two extra complete port vertices can be admitted beyond the initial degree-four transfer.

The same reasoning with maximum partial degree 10 leaves, for this test, only the reduced pair \(6,9\) with \(\deg a=2c,\ \deg b=3c,\ 3\mid c\): the other distinct pairs among \(6,8,9,10\) have gcd at most 2, and if \(3\nmid c\), the individual leading-degree gcds for \(6,9\) already qualify. This is an explicit boundary of the test, not a noninvertibility claim.

## Direct inverse for one controlled star

A single controlled star gives degree at most two. Here the inverse can be written without importing a degree theorem. If neither degree is two, use the [affine-input argument](UNRESTRICTED-JACOBIAN.md#4-what-the-jacobian-analogy-yields-and-where-it-stops). Otherwise a linear output change, using the top determinant relation, puts the map in the form
\[
 F_1=a(y)x^2+b(y)x+c(y),\qquad F_2=d(y)x+e(y).
\]
If \(d=0\), the determinant cannot be a nonzero constant while \(a\ne0\), so \(d\ne0\). The \(x^2\) and \(x\) coefficients of the determinant give
\[
 a=k_0d^2,\qquad b/d-2k_0e=r_0\qquad(k_0,r_0\in k).
\]
Indeed those equations are \((a/d^2)'=0\), then \((b/d)'=2k_0e'\). Hence
\[
 F_1=k_0 F_2^2+r_0F_2+s(y),\qquad s=c-k_0e^2-r_0e .
\]
The remaining determinant is \(-ds'\). Therefore \(d\) is constant and \(s\) is affine with nonzero slope. To invert, recover \(y\) from \(F_1-k_0F_2^2-r_0F_2\), then \(x=(F_2-e(y))/d\). In particular the original vector of \(x^2\) coefficients must have been **constant in \(y\)**: the output change was linear and \(a=k_0d^2\) is constant. This is a direct obstruction when the genuine pairwise interaction of a star varies with the other control.

The Jacobian restrictions constrain the polynomial interpolation; they neither prove existence of an unrestricted prime-to-sevenfold square nor replace checking that an inverse output has literal Boolean controls.
