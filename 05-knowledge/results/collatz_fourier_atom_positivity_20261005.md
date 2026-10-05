# Positivity as a Fourier-coefficient problem, with the atom retained

2026-10-05. **PROVED, scoped:** lossless continuous encoding of the injection
measure, effective finite approximation, a quantitative refinement criterion,
non-Hoelder behavior and a modular-transformation obstruction. **FINITE-EXACT:**
the accompanying controls. **OPEN:** positive weight at every integer and
universal Collatz. No automorphic correspondence for the Collatz operator is
claimed.

## 1. Inheritance and the research target

The anchor is positivity of the critical flow at every integer. The niche
is lossless passage between atomic measures and continuous characters. The
wildcard is the level-22 oldform representation and its positive graph lift.
The board is **individual atom / inverse section / residue representation /
Fourier coefficient / positive cone / refinement floor**.

Incoming commit `bde61bf16`, [adaptive mixture flow](collatz_adaptive_mixture_flow_20261005.md),
sharpens the previous mass bound and identifies the exact injection measure.
Its weight on a rooted source with base-step count T and total sibling depth K is

\[
 W(n)=\frac{2K!(T+1)!}{(T+K+2)!};
 \qquad W(n)=0\text{ off the rooted component}.                \tag{1}
\]

It proves \(\sum W\le16/3\), and on the positive odd nonroot sources
\(W-\mathcal KW=\lambda\), where

\[
 \lambda(n)=\begin{cases}W(n),&3\mid n,\\0,&3\nmid n,\end{cases}
 \qquad \sum_n\lambda(n)=1.                                 \tag{2}
\]

The subsequent incoming commit `598879da4` improves this W-mass bound to
less than41/10 using a checked cost-eight kernel, and constructs a singular
mixture Z with nonroot mass less than23/8. Both injection measures still
have total mass one and the same rooted support. The formulas here concern
W and lambda in (1)-(2); the stronger bounds do not alter that choice or
establish its missing atoms. Equation (6) below uses the simpler analytic
bound, so it does not depend on that new finite kernel.

This is the closest proved mechanism. The hostile is a missing atom despite
positive mass in every finite residue class. The corrected near miss is
mistaking continuity of a transformed measure for positivity of its individual
coefficients. The least-used sidecar is a lower bound attached to one fixed
integer as the representation is refined.

The inherited section sends every odd 3-unit n to an odd multiple of three
z with U(z)=n and z<(64/3)n<22n. The companion
[inverse sections](inverse_predecessor_sections_20261005.md) proves sharp
versions, including the size-induction boundary and the failure of a fixed
size window to capture a fixed share of incoming weight. Thus positivity
of (2) at every odd multiple of three would settle the entire target.
Its unit total mass alone does not establish that support statement.

Two policy cards guide the transfer: **the same representation is not the
same carrier**, and **audit the positive kernel before lifting through a
flat square**, in [META-PATTERNS](../../00-navigation/META-PATTERNS.md).

## 2. A genuine continuous-discrete equivalence

Index the injection sources by n=6m+3, m>=0, and put
\(p_m=\lambda(6m+3)\). Define

\[
 F(z)=\sum_{m\ge0}p_m z^m\quad(|z|\le1),\qquad
 \Phi(t)=F(e^{2\pi it}).                                    \tag{3}
\]

**F1.** The series converges absolutely and uniformly on the closed disk.
It is holomorphic inside and continuous on the circle, with |F|<=1 and
Phi(0)=1. It retains every atom exactly:

\[
 p_m=\int_0^1\Phi(t)e^{-2\pi imt}\,dt.                       \tag{4}
\]

Uniform convergence and elementary character orthogonality prove these
claims. It also gives a positive-definite kernel, since for any finite
complex scalars c_j and real t_j,

\[
 \sum_{i,j}c_i\overline{c_j}\Phi(t_i-t_j)
 =\sum_m p_m\left|\sum_i c_i e^{2\pi imt_i}\right|^2\ge0.    \tag{5}
\]

This is an actual bridge, not a comparison of shapes: atoms map to characters,
and integration recovers the atoms. The complete transform loses nothing.
However, nonnegative coefficients and positive coefficients at *every* index
are different predicates. The target is precisely p_m>0 for all m.

### Effective finite approximation

The improved incoming bound also strengthens the earlier numerical compiler.
Let S(n)=4n+1. Keep only rooted inverse-base depth T<N, individual inverse
sibling depths k<=B, and final sibling depth j<=B. Give these finitely many
sources their full weights (1). Call the resulting submeasure W_(N,B).
For N>=2, B>=1, and rational 0<epsilon<1,

\[
\begin{aligned}
 \|W-W_{N,B}\|_1\le{}&4\epsilon+
 4(1-\epsilon^2/3)^{N-2}\\
 &+\frac6B+\frac6{B+1}+\frac{12}{B+2}.                      \tag{6}
\end{aligned}
\]

Here is a proof retaining the parameter endpoints. At a fixed 0<r<1,
write D=1+r+r^2 and kappa=(1+r)/D. Incoming P1 proves the rooted base mass
bound P(r)=2+2r-r^2<=3. Its first and second inverse generations have masses
A and at most A2, with

\[
 1+A+\frac{A_2}{1-\kappa}=P(r),\qquad
 \frac{2A_2}{1-\kappa}\le4.
\]

After mixing with density 2(1-r), the T>=N tail is at most
\(4\int_0^1\kappa^{N-2}dr\). Since kappa<=1-r^2/3, splitting at epsilon
gives the first line of (6). Removing inverse depths above B perturbs the
base operator by norm at most r^(B+1). The resolvent estimate therefore
gives base error at most 3D r^(B-1). Sibling extension and mixing give
\(6\int_0^1D r^{B-1}dr\). Finally the outer sibling tail costs at most
\(6\int_0^1r^{B+1}dr\). These are exactly the three rational terms in (6).

First reduce epsilon, then increase N and B. Restriction to multiples of
three cannot increase this error. Consequently the finite Fourier polynomials
from this compiler approximate (3) uniformly with a proved effective error,
without assuming universal convergence. They approximate a measure that might
still have zero atoms. This does not turn approximation into a support decision.

## 3. The exact finite-representation test is a uniform atom floor

For a modulus q define the residue mass

\[
 p^{(q)}(a)=\sum_{m\equiv a\pmod q}p_m
 =\frac1q\sum_{\ell=0}^{q-1}\Phi(\ell/q)e^{-2\pi ia\ell/q}. \tag{7}
\]

This is finite Fourier inversion, justified by the same absolutely convergent
sum. For q>m,

\[
 0\le p^{(q)}(m)-p_m
 \le R(q):=\sum_{j\ge q}p_j\longrightarrow0.                \tag{8}
\]

Thus the full numerical tower of residue representations recovers every atom;
one finite representation does not. In particular

\[
 p_m>0\quad\Longleftrightarrow\quad
 \inf_{q>m}p^{(q)}(m)>0.                                   \tag{9}
\]

It suffices to take any unbounded sequence of moduli, with the infimum over
that sequence. The [representation package](collatz_representation_positivity_20261005.md)
proves the corresponding statement for the full W measure, its exact affine
and Fourier actions, and the information lost by omitting valuation colors.

**A precise sufficient certificate.** If a spectral or arithmetic argument
gives a rigorous lower bound L for p^(q)(m), and a proved tail bound R for
R(q), then

\[
                        p_m\ge L-R.                         \tag{10}
\]

An inequality L>R proves the selected atom positive. The desired new input
would establish such a bound without first assuming that source's ROOT word.
Computing both sides solely from a bank of completed ROOT words does not
automatically provide new sources; the independent lower-bound step matters.

### Positive at every finite level can still vanish at one integer

Delete atom1 from the geometric probability 2^(-m-1), and normalize:

\[
 \widetilde p_1=0,\qquad
 \widetilde p_m=\tfrac43\,2^{-m-1}\quad(m\ne1).
\]

Every residue class of every modulus has positive mass. Nevertheless for q>=2,

\[
 \widetilde p^{(q)}(1)=\frac1{3(2^q-1)}\longrightarrow0.      \tag{11}
\]

Its generating function is even the rational analytic function
\(\tfrac43[(2-z)^{-1}-z/4]\). The circulant Fourier matrices at every
finite modulus are strictly positive definite, because their eigenvalues
are q times these positive residue masses. Smoothness and finite spectral
positivity together still miss the omitted atom. This is a toy measure,
not an assertion that the actual Collatz atom at9 vanishes.

## 4. The actual continuous encoding has a strong regularity boundary

The earlier root ray is S^j(1), with W-weight 2/[(j+1)(j+2)]. In the odd
index (n-1)/2 its exponents are

\[
              \frac{2(4^j-1)}3=0,2,10,42,170,\ldots.        \tag{12}
\]

Thus the user's sequence enters the continuous encoding through the actual
Fourier indices of a certified family, with its measure retained.

For the injection measure itself, put k=9t+1, t>=1, and

\[
 n_k=S^k(1),\quad z_k=\frac{2n_k-1}3,\quad
 m_t=\frac{z_k-3}6=\frac{16(4^{9t}-1)}{27}.
\]

Here n_k=5 mod9, z_k is an odd multiple of three, and
z_k -> n_k -> 1 is an actual two-step certificate. The inverse-section
package gives the exact atom

\[
 p_{m_t}=\lambda(z_k)=\frac4{(k+1)(k+2)(k+3)}.               \tag{13}
\]

**F2.** Phi is not Hoelder at zero of any positive exponent. At
h_t=1/(2m_t), all summands in the following real difference are nonnegative,
and the displayed family alone gives

\[
 1-\operatorname{Re}\Phi(h_t)
 =\sum_m p_m(1-\cos(2\pi m h_t))\ge2p_{m_t}.                \tag{14}
\]

The right side decays as t^(-3), while h_t decays exponentially. Hence
(14)/h_t^alpha is unbounded for every alpha>0. The full W Fourier transform
has the same obstruction from (12), already with decay j^(-2).

Consequently this exact circle bridge is continuous but supplies no positive
Hoelder regularity or finite positive size moment. The fast finite mean
stopping time under lambda is compatible with enormous two-step sources.
These are different observables.

Also, a nonzero summable weight on positive integers cannot itself be
continuous in a p-adic topology at a point with positive weight: a distinct
integer sequence tending p-adically to that point escapes to infinity in
ordinary size, along which the summable weights tend to zero. Transforming
the measure as in (3) avoids this false pointwise continuity assertion.

## 5. How the level-22 bridge can be used honestly

The companion [oldform and eight-state lift](collatz_level22_oldform_bridge_20261005.md)
recovers a concrete answer to 22=2*11. The level-11 eta form f and f(2tau)
span the weight-two cusp space at level22, an established oldform example
in [Stein, Example9.6](https://wstein.org/books/modform/stein-modform.pdf).
In that basis its U2 matrix is

\[
 A=\begin{pmatrix}-2&1\\-2&0\end{pmatrix},\qquad A^4=-4I.
\]

Modulo3 it is the golden operator of order8. Over the reals, A preserves
no nonzero pointed cone. The package constructs an integral surjective
intertwiner from a **positive weighted eight-cycle** to A. This is a genuine
graph/representation quotient, but its signed projection loses positivity.
The eight states retain phase that the two-dimensional quotient compresses.
Neither the local Galois/Hecke match nor this intertwiner supplies a map from
actual Collatz sources preserving the ROOT predicate.

There is a further elementary boundary for identifying (3) directly with
a modular form. **F3.** A nonzero holomorphic function bounded throughout
the upper half plane cannot obey a positive-weight modular transformation
law on Gamma0(N), with a multiplier of absolute value one. Indeed

\[
 \gamma_t=\begin{pmatrix}1&0\\Nt&1\end{pmatrix}\in\Gamma_0(N),
 \qquad |H(\gamma_t\tau)|=|Nt\tau+1|^k|H(\tau)|.
\]

For k>0 and a point with H(tau)!=0 the right side is unbounded as t tends
to infinity. In contrast H(tau)=F(exp(2pi i tau)) has absolute value at
most1. Thus this particular probability generating function cannot itself
be a nonzero positive-weight modular form. This does not exclude other
automorphic constructions, nonlinear encodings, or an appropriately typed
intertwiner; it rejects the most direct equality of the two functions.

The bounded predecessor constant and the oldform level are therefore two
separate proved roles of22. The latter provides an exact arithmetic phase
module and graph lift. To advance positivity, a further map must preserve
the nonnegative mass and individual source, or prove the refinement floor
in (9)-(10). A correspondence preserving only spectra or traces is insufficient.

## 6. Reproduction and scope

[Program](../../04-computation/experiments/collatz_fourier_atom_positivity_20261005.py)
and [output](collatz_fourier_atom_positivity_20261005.out):

```text
python -B 04-computation/experiments/collatz_fourier_atom_positivity_20261005.py
python -B -O 04-computation/experiments/collatz_fourier_atom_positivity_20261005.py
```

The 915 exact checks include Fourier inversion in rational cyclotomic
coordinates at primes2,3,5,7,11,13 for512 certified injection atoms;
all toy residue projections through modulus100; one hundred root-ray
indices; forty two-step lacunary injection sources; and the rational
constants in (6). The script imports only the existing certified finite-box
generator, whose checks belong to its earlier package. There are no floating
Fourier or eigenvalue approximations. The modular and regularity obstructions
are proved above; their finite controls do not replace their proofs.

The output sought next is an independently grounded lower bound for a
particular unresolved injection atom, stable through refinement. This gives
a precise role to new representation-theoretic information while keeping
universal positivity open.
