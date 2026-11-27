# Additive triples, Pythagorean spinors, and the half-length Fermat cycles

This extends the [abc radical-defect law](../abc_euler_product_phases_and_radical_defect_20261127/README.md), the [Fermat odd-pair calculation](../fermat_fields_moore_gram_and_pollock_lifts_20261126/README.md), and the [Berggren factor coordinates](../chanpak_berggren_divisor_twins_20261123/README.md). The cycle length here is \((p-1)/2\).

## 1. A Pythagorean triangle contains two linked primitive additive triples

Write the primitive triangle in odd-leg, even-leg, hypotenuse order as
\[
(A,B,C)=(m^2-n^2,\;2mn,\;m^2+n^2),
\qquad m>n>0,\quad (m,n)=1,\quad m+n\text{ odd}.
\]
The previous Berggren coordinates \(u=m-n,\ v=m+n\) give not just the squared additive equation \(A^2+B^2=C^2\), but also the **two** primitive linear equations
\[
u+2n=v,\qquad u+v=2m.
\]
Every positive primitive \(a+b=c\) with \(b\) even and \(a,c\) odd occurs in the first equation: take \(m=(c+a)/2,\ n=(c-a)/2\). Conversely it determines those spinor coordinates uniquely. The union of prime supports of the two primitive equations is exactly the support of \(AB=2mn(m-n)(m+n)\); the hypotenuse \(C\) adds new primes only, all \(1\bmod4\) since for an odd prime dividing \(C\), \(m/n\) is a square root of \(-1\). At the root \((3,4,5)\) the two linear equations are \(1+2=3\) and \(1+3=4\), while the squared equation is \(9+16=25\).

This is an explicit gluing distinction: the first \(a+b=c\) supplies three divisibility lines, a second equation sharing the same \(m,n\) leaves four distinct lines at an odd prime, and the norm \(m^2+n^2\) adds either two lines or none. It is not a product of independent densities for the two additive equations.

## 2. The same half-length cycle as a rotation of the reduced Pythagorean circle

For \(p\equiv1\bmod4\) choose \(i\in\mathbb F_p\) with \(i^2=-1\). The affine conic \(X^2+Y^2=1\) is then identified with \(\mathbb F_p^\times\) by \(z=X+iY\), with \(z^{-1}=X-iY\). In the primitive Pythagorean parametrization, for \(C\not\equiv0\bmod p\),
\[
z=\frac{A+iB}{C}
  =\frac{m+in}{m-in}.
\]
This is the finite-field counterpart of the ordinary doubled angle \(2\arctan(n/m)\). The two omitted projective slopes \(n/m=\pm i\) are exactly the hypotenuse-divisibility lines; they have zero norm and no affine circle phase.

For a Fermat prime \(p\ge17\), the previous result gives \(g=3/2\) of order \(p-1\) and \(\lambda=g^2=9/4\) of order \(L=(p-1)/2\). Identify the centered nonzero odd-pair coordinate \(x+1\) with \(z\). Then the two modular Collatz orbits transfer exactly to two cycles of the Pythagorean circle under multiplication by \(\lambda\), one on each square class. In spinor slope \(t=n/m\), this action is the projective transformation
\[
t\longmapsto\frac{t+t_0}{1-t t_0},\qquad
t_0=\frac{\lambda-1}{i(\lambda+1)}.
\]
This follows from the Cayley expression \(z=(1+it)/(1-it)\), and includes a possible passage through the slope at infinity. For \(p=17,\ i=4\), the reduced \((3,4,5)\) has \(z=14\) and goes around the nonsquare cycle of length 8. If \(p\equiv3\bmod4\), the circle instead is the norm-one group in \(\mathbb F_{p^2}\) of order \(p+1\); there is no split-field identification giving these \((p-1)/2\) cycles. None of these finite rotations is a positive-integer Collatz cycle or a move up the Berggren tree.

### The two phases give different answers on one full cycle

Let \(z\) lie in the coset \(u\langle\lambda\rangle\). A multiplicative character \(\rho\) satisfies
\[
\sum_{j=0}^{L-1}\rho(u\lambda^j)=
\begin{cases}
L\rho(u)&\rho(\lambda)=1,\\
0&\text{otherwise}.
\end{cases}
\]
Thus the only multiplicative frequencies left by averaging around the orbit are the trivial and quadratic characters. This is the finite-field analogue of a logarithmic phase: multiplication becomes addition of phase.

The ordinary additive character \(e_p(kz)=\exp(2\pi i kz/p)\) instead gives, for \(k\ne0\),
\[
\sum_{j=0}^{L-1}e_p(ku\lambda^j)
=\frac{-1+\chi(ku)\tau_p}{2},\qquad
\tau_p=\sum_{z\ne0}\chi(z)e_p(z),
\]
where \(\chi\) is the quadratic character. Insert the coset indicator \((1+\chi(z/u))/2\); the two terms are ordinary additive orthogonality on the units and the quadratic Gauss sum. For these \(p\equiv1\bmod4\), the [classical Gauss evaluation](https://math.indianapolis.iu.edu/~bnica/gauss_sums.pdf) gives \(\tau_p=\sqrt p\). Already at 17 the two sums are \((-1+\sqrt{17})/2\) and \((-1-\sqrt{17})/2\), instead of the vanishing from primitive multiplicative frequencies. The additive phase distinguishes the two half-cycles even after averaging over every step.

There is a third readout if one insists on observing the actual normalized legs \(X,Y\) rather than the split coordinate \(z\): a linear additive phase \(e_p(aX+bY)\) becomes \(e_p(\alpha z+\beta z^{-1})\), where \(\alpha=a/2+b/(2i)\) and \(\beta=a/2-b/(2i)\). The half-orbit average is then the corresponding half-sum of the Kloosterman sum and its quadratic-character-twisted version, with the latter signed by \(\chi(u)\). One cannot silently replace an additive phase on the legs by an additive phase on the group coordinate.

## 3. A second exact Euler product: radical defect in the primitive triangles

Sample the Euclid parameters uniformly from \(1\le n<m\le X\), subject to coprimality and opposite parity. Let \(D_P=ABC/\operatorname{rad}(ABC)\). For odd \(p\), put \(r_p=5+\chi_4(p)\), where \(\chi_4(p)=1\) for \(p\equiv1\bmod4\) and \(-1\) for \(p\equiv3\bmod4\). Then for \(k\ge1\),
\[
\lim_{X\to\infty}\Pr(v_p(ABC)\ge k)=\frac{r_p}{p^{k-1}(p+1)}.
\]
The four linear factors \(m,n,m-n,m+n\) have distinct primitive projective roots at every odd prime. For \(p\equiv1\bmod4\), the quadratic \(m^2+n^2\) has two more distinct roots, each lifting uniquely to every \(p^k\); otherwise it contributes none. Each line has proportion \(p^{-k}(1-p^{-1})/(1-p^{-2})\) among primitive residue pairs, as in the abc calculation. Thus at split primes the density of high powers is exactly twice the abc *local* density; at inert primes it is \(4/3\) times that density. At 2 the even spinor has a geometric valuation with mean 2, starting at 1, and \(v_2(D_P)\) is exactly that valuation.

For clarity, the product over primes here is a true limit, including the hypotenuse. For a fixed finite set of primes, primitive opposite-parity pairs obey the local densities independently by CRT and the common-prime sieve. For a uniform large-prime bound on the norm, each root modulo \(p^k\) defines a lattice \(m\equiv r n\bmod p^k\); distinct lattice points have Euclidean separation at least \(p^{k/2}\), since the squared norm of their nonzero difference is a positive multiple of \(p^k\). Disk packing in the box bounds the count by \(O((1+X/p^{k/2})^2)=O(X^2/p^k)\) whenever \(p^k\le2X^2\). Beyond that cutoff there are no positive norm hits. The linear factors have the analogous \(O(X^2/p^k)\) bound by slicing. Uniform domination of \(\sum_p\sum_{k\ge2}\log(p)/p^k\) now justifies both the defect distribution and the finite-box expectation limit. In particular
\[
\lim\Pr(D_P=2)
=\tfrac12\prod_{p>2}\left(1-\frac{5+\chi_4(p)}{p(p+1)}\right)
\approx0.21447,
\qquad
\lim\mathbb E\log D_P
=2\log2+\sum_{p>2}\frac{(5+\chi_4(p))\log p}{p^2-1}
\approx2.9623.
\]
This is the previous logarithmic loss calculation with the split/inert character retained instead of averaging it away. For the squared additive triple \((A^2,B^2,C^2)\), the radical stays the same; its quality tends to \(2/3\) in this parameter box, rather than \(1/3\) for an ordinary primitive additive box. Indeed, away from the arbitrarily small boundary strips in \(m/X,n/X\), the three base sides have logarithmic size \(2\log X+O(1)\), and the logarithmic defect is bounded in mean. This is a typical-value theorem, not a bound on exceptional \(abc\) triples.

Run `node replay.js` to check the two complete orbits at 17, 257 and 65537, the local prime-power counts, and the finite triangle boxes.
