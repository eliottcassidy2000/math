# Infinitely many locally admissible exponents at a finite square-sum star

A signed Collatz square-sum edge with center \(y\), sign \(e\) and exact division exponent \(j\) has square root \(m\) of the sum of its endpoints. The equation to invert is
\[
             3m^2=2^j y+3y-e .
\]
Here \(y\equiv1\pmod3\), \(y\equiv-e\pmod4\), \(m\) is even and \((-1)^j=e\). Once \(y,e,j\) are given, a positive root determines the unique leaf. [Signed-edge classification and the Pell construction](SIGNED-SQUARES.md)

## Three different questions about finiteness

**PROVED using the integral-point theorem for cubics.** Every center has finite degree, even when every division exponent is allowed. Infinitely many centers have two incoming edges with the same two prescribed exponents 4 and 8. For any three *fixed* distinct exponents, only finitely many centers admit all three at once.

For a fixed center and sign, write \(j=3q+d\), \(0\le d<3\). Put \(X=2^q\), \(Y=3m\). Then
\[
                 Y^2=3\cdot2^d y X^3+3(3y-e).
\]
Both coefficients are nonzero. The cubic therefore has three distinct complex roots, and Siegel's theorem gives finitely many integral points. Each of the three curves accounts for finitely many exponents. These conclusions are effective in principle by the bounds for hyperelliptic equations. [Baker, *Bounds for the solutions of the hyperelliptic equation*](https://doi.org/10.1017/S0305004100044418) · [Mason, statement of the integral-point theorem](https://www.cambridge.org/core/journals/mathematical-proceedings-of-the-cambridge-philosophical-society/article/abs/hyperelliptic-equation-over-function-fields/257B19914D45C19F80DE807A7239347D)

For three prescribed exponents put \(A_i=2^{j_i}+3\). Any common center uses a common sign. Multiply \(3m_i^2=A_i y-e\) to obtain
\[
               (9m_1m_2m_3)^2=3\prod_{i=1}^3(A_i y-e).
\]
The three roots \(e/A_i\) are distinct. This again gives finitely many integral centers. For two prescribed exponents the corresponding curve is a conic. If \(AB\) is nonsquare, every positive seed yields infinitely many centers by the Pell action. If \(AB\) is square, write \(A=d a^2, B=d b^2\) with \(d\) squarefree. The difference of squares \((a v-bu)(a v+bu)=e(B-A)/(3d)\) is fixed and nonzero. If the right side is integral it has only finitely many factor pairs. These statements do not give a uniform bound on the degrees of all centers.

## Exactly what the ternary residues allow

Fix \(y\equiv1\pmod3\) and an exponent \(j_0\ge2\) with \((-1)^{j_0}=e\). For \(j=j_0+2t\) put
\[
              F(t)=(2^{j_0}4^t y+3y-e)/3 .
\]
For every \(k\ge1\), the values of \(F(t)\) run once through **every** residue mod \(3^k\) as \(t\) runs mod \(3^k\). Indeed \(4\) generates \(1+3\mathbb Z\) modulo \(3^{k+1}\), by \(v_3(4^u-1)=1+v_3(u)\), and \(2^{j_0}y\) is a unit. The exact proportion of exponent classes for which \(F(t)\) is a square modulo \(3^k\) is
\[
       3^{-k}\left(1+\sum_{i=0}^{\lfloor(k-1)/2\rfloor}3^{k-2i-1}\right)
                    \longrightarrow \frac38 .
\]
Count zero once. At valuation \(2i<k\), the unit must be one of the \(3^{k-2i-1}\) unit squares at the remaining precision. In particular, at every center any three consecutive same-parity exponents fail jointly already at precision \(3^2\) on the original undivided equation: after division by three their proposed squares exhaust the three classes modulo three.

## A fixed center defeats every finite battery of congruences

**PROVED.** At \(y=79\) with the positive sign, for every positive modulus \(L\) and every lower bound \(J\), there are infinitely many even exponents \(j>J\) admitting even square roots of \((79\cdot2^j+236)/3\) modulo \(L\). Yet only finitely many exponents admit an actual integer square root. This is a fixed-center, unbounded-exponent obstruction to using modular selection for absolute equality.

The seed is \(1685\to79\) at \(j=6\), with \(1685+79=42^2\). Write \(L=2^b q\) with \(q\) odd. Taking \(j\equiv6\pmod{\operatorname{ord}_{3q}(2)}\) fixes the quotient at \(42^2\pmod q\) and keeps \(j\) even. For the binary part and any such \(j\ge6\), the quotient is four times the odd integer
\[
                  (79\cdot 2^{j-2}+59)/3\equiv1\pmod8 .
\]
Every odd unit congruent to 1 modulo 8 has a square root modulo every power of two. Doubling that root gives an even root of the original quotient. CRT combines it with 42 modulo \(q\); when \(b=0\), impose evenness separately. This gives arbitrarily late exponent candidates simultaneously at any finite battery by taking its least common multiple. The finite-degree theorem removes all sufficiently large exponents in the integers. The quantifiers matter: the exponent can vary when the modulus changes.

The tournament analogue inherited from the Farey spines admits exact target congruences at all finite precisions while the absolute crossing remains open. Here the same logical gap has an example in which finite integral support can be proved. [Farey local solvability and the absolute search](../chanpak_odd_square_spines_20261122/README.md)
