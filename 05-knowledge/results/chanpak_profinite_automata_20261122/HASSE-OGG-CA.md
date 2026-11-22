# The Ogg test as a central read of a linear cellular automaton

For a prime \(p\ge5\), put \(m=(p-1)/2\), and work in
characteristic \(p\). Let \(S\) translate a configuration by one site,
and start with a single unit at site zero. The linear cellular
automaton
\[
T_\lambda=(1+S)(1+\lambda S^{-1})
          =(1+\lambda)I+S+\lambda S^{-1}
\]
has the following central column:
\[
c_n(\lambda)=[X^0]\big((1+X)(1+\lambda/X)\big)^n
             =\sum_{j=0}^n \binom nj^2\lambda^j.
\]
Its row at time \(m\) thus tests the supersingularity of the Legendre
curve \(y^2=x(x-1)(x-\lambda)\): the coefficient of \(x^{p-1}\)
in \((x(x-1)(x-\lambda))^m\) is
\((-1)^m c_m(\lambda)\). This coefficient is the Hasse invariant
in this coordinate (the Cartier operator on \(dx/y\) takes the
\(p\)-th root of this coefficient). Every elliptic curve in
characteristic \(p\ge5\) has a Legendre model over the algebraic closure,
by sending two of the three distinct finite 2-torsion abscissas to 0 and 1.

Write \(H_p=c_m\), \(A=256(1-\lambda+\lambda^2)^3\), and
\(B=\lambda^2(1-\lambda)^2\), so the Legendre \(j\)-invariant is
\(J=A/B\). Then the exact test is
\[
p\text{ is an Ogg prime}\quad\Longleftrightarrow\quad
H_p\mid A^p-A B^{p-1}\quad\text{in }\mathbf F_p[\lambda].
\]
Indeed \(H_p(0)=1\), while \(H_p(1)=\binom{p-1}{m}\ne0\) by
Vandermonde. The polynomial \(H_p\) is squarefree: coefficient
comparison using \((j+1)^2 a_{j+1}=(m-j)^2 a_j\) proves
\[
\lambda(1-\lambda)H_p''+(1-2\lambda)H_p'-\tfrac14 H_p=0.
\]
At any multiple root \(\lambda_0\notin\{0,1\}\) of order
\(2\le d\le m<p\), the first term has order exactly \(d-2\)
with nonzero leading coefficient, and the other terms have
higher order; this is impossible. Consequently the divisibility
criterion tests exactly \(J(\lambda)^p=J(\lambda)\) on all Hasse roots.
It does not test that the Legendre *parameters* themselves are rational.

The same automaton has an exact digit recursion, even for symbolic
\(\lambda\). For \(0\le r<p\),
\[
c_{pn+r}(\lambda)=c_n(\lambda^p)c_r(\lambda).
\]
Indeed the support of the \(r\)-th power lies in \([-r,r]\), while
the \(pn\)-th power has only exponents divisible by \(p\), by
Frobenius, after replacing \(X,\lambda\) by \(X^p,\lambda^p\).
Their product can reach zero only from two zero exponents. Reading
digits \(n_i\) least significant first gives
\(c_n(\lambda)=\prod_i c_{n_i}(\lambda^{p^i})\). For
\(\lambda\in\mathbf F_{p^f}\), an automatic readout therefore uses
only a multiplicative accumulator and position modulo \(f\).
This is a separate CA in each characteristic, not a single apparatus
recognizing the integer prime list in an integer digit base.

References: [Ogg, original 1974/75 paper](https://numdam.org/item/SDPP_1974-1975__16_1_A4_0/);
[modern arithmetic revisit](https://celebratio.org/Ogg_AP/article/1157/);
[Furstenberg, algebraic series as rational diagonals](https://doi.org/10.1016/0021-8693(67)90061-0);
[Rowland–Yassawi, automatic columns of linear CA with memory](https://doi.org/10.1016/j.aam.2014.10.002).
The repository already checks the Hasse polynomial and connects Ogg
class numbers to the marked Paley half-row in
`05-knowledge/results/ogg_triangular_chowla_20260927.md`; its
`chanpak_midpoint_ogg_chambers_20261118/README.md` separately uses a
finite Ogg-prime bank for arithmetic chamber factorization. Neither
uses this CA cut to obtain a tournament Hamilton count.

The standalone finite replay `check-hasse-ogg-ca.py` certifies the
divisibility test on primes 5 through 127, and separately verifies
the literal CA layer, combinatorial formula and independent floating-free
digit evaluation at representative finite field parameters.
The finite range is a test of the bridge, not a substitute for Ogg's
all-prime theorem or the separate characteristic-2/3 cases.
