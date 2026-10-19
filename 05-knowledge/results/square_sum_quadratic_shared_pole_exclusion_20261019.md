# No degree-two contact map with exactly one shared pole pair

**Status: PROVED by exact coefficient identities, independently checked against symbolic elimination.** This closes the remaining shared-input-pole case. It makes no claim here about two pole pairs disjoint from the input pair.

Put

\[
f(t)=\frac{4t(t^2-1)}{(t^2+1)^2},\qquad
k_\eta(t)=\frac{t^2+2\eta t-1}{t^2+1}.
\]

## Theorem

There are no rational functions `h,y in Q(t)` with `deg h=2` satisfying

\[
y^2=1+\varepsilon f(t)+\eta f(h(t)),\qquad
\varepsilon,\eta\in\{1,-1\},
\]

if `k_eta(h(t))` has four distinct simple poles, exactly two of them `i,-i`, and its other two poles lie in `Q(i)`.

Together with the previously established Gaussian splitting theorem, this excludes the complete degree-two case in which exactly one simple pole pair is shared with the input.

## 1. Divisors force a short identity

Write `D=t²+1`, `T=t(t²−1)`, and `z=k_eta(h)`. Let the new conjugate pole pair have monic polynomial

\[
q=(t-a)^2+b^2,\qquad a,b\in\mathbb Q,\quad b\ne0,
\quad q\ne D.
\]

Set `r=y−z`, `s=y+z`, so `rs=epsilon f`. At the shared pair both factors have simple poles: neither can have a pole of order greater than one, and their product has a double pole. At the new pair exactly one has a simple pole and the other a simple zero, because their product is finite and nonzero there. The choice is constant on the conjugate pair because both functions have rational coefficients.

After swapping the names of `r,s` if necessary, and accordingly replacing `z` by `−z`, the entire divisor of `r` consists of zeros at the new pair and poles at the shared pair. Thus, for some nonzero rational `c`,

\[
r=\frac{cq}{D},\qquad
s=\frac{4\varepsilon T}{cqD},\qquad
z=\frac{4\varepsilon T-c^2q^2}{2cqD}.                 \tag{1}
\]

There are no missing zeros or poles at infinity: both `q,D` are monic quadratics, and the displayed `r` has nonzero finite value there.

## 2. Every eligible composition has a Gaussian quadratic form

Write `h=P/Q` with coprime rational polynomials of maximum degree two. One of `P+iQ` and `P−iQ` has root `i`. Choose that factor, whose other root is `a+bi`, changing the sign of `b` if necessary. Up to a nonzero constant it equals

\[
L=(t-i)(t-a-bi)=X+iY,
\]

where

\[
X=t^2-at-b,\qquad Y=-(b+1)t+a,\qquad
X^2+Y^2=Dq.                                         \tag{2}
\]

The usual expression for `k_eta(P/Q)` in terms of `(P+iQ)²` and its conjugate gives rational `U,V` such that

\[
z=\frac{U(X^2-Y^2)-2VXY}{Dq},\qquad U^2+V^2=2.       \tag{3}
\]

The optional sign change of `z` in Section 1 simply changes both `U,V`. Explicitly, if `P+iQ=C L`, then the complex coefficient is `(1−i eta)C/conjugate(C)`, which has norm two; choosing the conjugate factor gives the conjugate version of the same identity.

Set `S=c²` and `d=2cV`. The leading coefficient of (1)=(3) gives `2cU=−S`. Equating numerators therefore yields

\[
E(t):=S(q^2-X^2+Y^2)-2dXY-4\varepsilon T=0.           \tag{4}
\]

In particular `S!=0`. The norm constraint `S²+d²=8S` follows as well, but is not needed below.

## 3. Three two-by-two determinants exclude rational parameters

The four coefficients of (4), after removing harmless nonzero numerical factors, are

\[
\begin{aligned}
F_3&=aS-(b+1)d+2\varepsilon,\\
F_2&=(5a^2+3b^2+4b+1)S-2a(b+2)d,\\
F_1&=(2a^3+2ab^2+2ab+a)S+(-a^2+b^2+b)d-2\varepsilon,\\
F_0&=(a^4+2a^2b^2+a^2+b^4-b^2)S+2ab d.
\end{aligned}                                       \tag{5}
\]

All four vanish. Consequently the three homogeneous equations `F2=0`, `F3+F1=0`, and `F0=0` have a nonzero solution `(S,d)`. Their coefficient matrix is

\[
M=\begin{pmatrix}
5a^2+3b^2+4b+1&-2a(b+2)\\
2a(a^2+b^2+b+1)&-a^2+b^2-1\\
a^4+2a^2b^2+a^2+b^4-b^2&2ab
\end{pmatrix}.
\]

Put `H=a²+(b+1)²`. The three row-pair determinants are exactly

\[
\begin{aligned}
\det M_{12}&=H\bigl((4b+3)a^2+3b^2-2b-1\bigr),\\
\det M_{13}&=2aH\bigl((b+2)a^2+b^3+b\bigr),\\
\det M_{23}&=H\bigl(a^4+(2b+1)a^2-b^2(b-1)^2\bigr).
\end{aligned}                                       \tag{6}
\]

Here `H!=0`: for rational real `a,b`, its only zero is `a=0,b=−1`, which would give `q=D`, excluded by hypothesis.

If `a=0`, the first and third equations in (6) force both `(3b+1)(b−1)=0` and `b²(b−1)²=0`; hence `b=1`. But then `F2=8S`, contradicting `S!=0`.

If `a!=0`, the first two equations in (6) imply

\[
(4b+3)a^2+3b^2-2b-1=0,
\quad
(b+2)a^2+b^3+b=0.
\]

Eliminating `a²` without dividing by any parameter gives

\[
\begin{aligned}
0&=(4b+3)(b^3+b)-(b+2)(3b^2-2b-1)\\
 &=2(2b^4+4b+1).                                    \tag{7}
\end{aligned}
\]

The rational-root theorem restricts a rational root of `2b⁴+4b+1` to `±1,±1/2`; none is a root. This is the desired contradiction. QED.

## 4. Exact reproducibility and scope

`verify_quadratic_shared_pole_exclusion_20261019.py` checks (2), the coefficient identities (5), all three determinant factorizations (6), and elimination (7) over `Q`. It also checks the four rational-root candidates. The verifier uses only sparse polynomial arithmetic over Python’s standard-library `Fraction`, with no solver or third-party algebra package. The optional discovery script `explore_shared_sympy_crosscheck.py` separately confirms that the five coefficient/norm equations have Groebner basis `[1]`; that computation is unnecessary for the elementary proof.

The proof itself is elementary and needs neither a bounded coefficient search nor the Groebner computation. It handles both signs simultaneously and both orientations of every Gaussian pole pair. The distinct-pole hypothesis is essential to its stated scope; repeated poles are handled separately in the low-degree obstruction work.

All verifier conditions are explicit exceptions, rather than removable assertions. Normal and `python -O` runs give byte-identical output, recorded by the published `quadratic_shared_pole_controls.json` receipt and matching optimized-mode replay.
