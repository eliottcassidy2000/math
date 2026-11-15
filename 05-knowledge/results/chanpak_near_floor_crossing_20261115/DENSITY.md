# What the uniform carrier Fourier profile does not prove

The order-five mask-zero carrier has no local obstruction and its Fourier
peaks are now controlled to leading relative scale on the entire circle.
One still needs a genuinely exact-scale argument for density one. This note
proves that the present type of analytic result cannot establish that
conclusion on its own: a density-zero selection of *actual count fibers*
retains the same scaled uniform profile and local laws, even with diverging
average multiplicity. It also computes the exact second-moment constant a
collision proof must match for this one carrier.

## A density-zero model using only physical presentations

Fix a form $F$ in the physical finite bank and $0<\delta<1/2$. At
$T=2^k$ let $D=\lceil T^\delta\rceil$, and let $R_T(n)$ count
presentations $x,y\ge D$ at $F(x,y)=n\in (T/2,T]$. Then
$$
  \sum_n R_T(n)\ll T\log T,\qquad \max_n R_T(n)\ll T/D.
$$
For the second estimate split into $x\ge y$ and $y>x$, fix the
larger coordinate, and use strict monotonicity in the smaller one from
the [sharp-cut packet](../chanpak_true_cut_fourier_20261115/README.md).
The larger coordinate has $O(T/D)$ possibilities by $Qxy\le T$.
Both bounds also hold for a fixed finite union of forms.

Fix $0<\beta<1$. For each odd $n\in (2^{k-1},2^k]$, independently
keep its whole fiber with probability $p_k=(k+1)^{-\beta}$; let $B_n$
indicate selection. This defines a single random subset of the odd
integers. Its retained presentations on a shell are
$\widetilde R_T(n)=B_n R_T(n)$. Almost surely, for all sufficiently
large $k$, uniformly in $\alpha\in\mathbb R/\mathbb Z$ and
$W\in[T/2,T]$,
$$
 \left|\sum_{T/2<n\le W}(\widetilde R_T(n)-p_k R_T(n))e(\alpha n)\right|
  \ll T^{1-\delta/2}\log T.                                   \tag{1}
$$
Indeed the sum of the squared summand ranges is bounded by
$\sum R_T(n)^2 \ll (T/D)T\log T$. Hoeffding on the real and imaginary
parts gives failure probability $O(T^{-A})$ at any one frequency and
integer sublevel, for an arbitrary fixed $A$ by choosing the constant
in (1). Use a circle mesh of spacing at most $T^{-3}$ and at most
$T+1$ distinct sublevels. The derivative in frequency is bounded by
$O(T^2\log T)$, so mesh interpolation costs only $O(T^{-1}\log T)$.
Taking $A>6$ and summing over dyadic shells, Borel–Cantelli gives (1).
The same event gives (1) for any finite bank, by applying the finite
union bound over its forms. Chebyshev and Borel–Cantelli show that the
number of selected odd integers in a shell differs from $p_k$ times
their number by $O(T^{3/4})$ eventually almost surely. Summing shells
shows that the selected set has natural density zero.

Consequently there exists an unweighted set of genuine physical
presentations having all of the following:

* Its support has density zero because the whole integer fiber is
  retained or removed together.
* On dyadic shells its total presentation mass is
  $\asymp T(\log T)^{1-\beta}$ for any one carrier with positive
  leading mass. In particular mean presentations per target diverge.
* Its entire [uniform Fourier profile](../chanpak_voltage_profile_20261115/README.md)
  is that of the original bank multiplied by $p_k$, up to a power
  error. The same holds for every shell subinterval.
* Residue-class and logarithmic-window laws scale by $p_k$. A
  residue restriction is an average over rational phases, so applying
  (1) to it incurs no extra factor in the modulus. Thus all the
  growing-modulus, genuinely positive local predictions in the
  original logarithmic range persist.

An even stronger analytic foil is the nonnegative sequence
$\widehat R_T=p_k^{-1}\widetilde R_T$: it has the *unscaled* mass,
uniform Fourier main term and exact local factors up to a power error,
but is still supported on a zero-density set. These are weights, not
new unweighted physical presentations. This argument does **not**
show that the full actual carrier is sparse. It shows why uniform
relative Fourier estimates, positivity of local factors, and
diverging average multiplicity do not by themselves force density one.

## The exact off-diagonal test for the mask-zero carrier

Now use the fixed carrier $Q=100, a=11, t=89, f=11, s=e=0, n_0=11$.
For clarity use $n$ for the output count and $n_0$ for the constant.
On $T/2<n\le T$, $n$ odd, put
$$
 \ell_T=\frac{1-2\delta}{Q}\log T,\qquad
 m_T(n)=2\ell_T c(n),\qquad
 c(n)=\begin{cases}21/11&11\mid n,\\10/11&11\nmid n.\end{cases}
$$
The mask-zero local law in the sharp-cut packet gives, for each odd
class $b\bmod22$, $\sigma_{11}(b)=c(b)/11$. The hyperbolic area
derivative on the shell is $\ell_T+O(1)$; hence
$$
 \sum_{n\equiv b\;(22)}R_T(n)
     =\frac{T}{2}\frac{c(b)}{11}\ell_T+o(T\log T).
$$
This agrees to the same error with the sum of $m_T$ in that class.
Consequently $\sum R_Tm_T=\sum m_T^2+o$T\log^2T$$. The average
of $c^2$ over the eleven odd residue classes is
$$
 \frac{(21/11)^2+10(10/11)^2}{11}=\frac{131}{121}.
$$
All sums in the following display run over odd counts in the shell:
$$
 \sum (R_T-m_T)^2
    = \sum R_T(R_T-1)
      -\frac{131}{121}T\ell_T^2+o(T\log^2 T).                  \tag{2}
$$
In particular the right off-diagonal leading coefficient cannot be
smaller than $131/121$ after normalization by $T\ell_T^2$.
A matching upper bound with $o$T\log^2T$$ error would prove density
one for this one physical carrier. Indeed each unsupported odd
target contributes at least $(20\ell_T/11)^2$ to the left in (2).
The same proof at dyadic endpoints suffices for natural density.
This criterion is sufficient, not asserted necessary.

Here is the equivalent row-intersection problem. For each larger
coordinate $r$, split into the two orientations $x=r\ge y$ and
$y=r>x$. In each such row the shorter-coordinate function is
strictly increasing and its images form a set $A_i$, rather than a
multiset. Then $R_T(n)=\sum_i 1_{A_i}(n)$, and
$$
 \sum_n R_T(n)(R_T(n)-1)=\sum_{i\ne j}|A_i\cap A_j|.
$$
The issue is therefore not the fine distribution *inside* one row,
but the exact intersection mass *between* rows, summed over different
shorter-coordinate transitions. Any fixed bounded-degree relation
among presentations (such as one specific transpose alias, or a
fixed finite bank of checked bijective cylinder aliases) contributes
only $O(T\log T)$. Such aliases should be accounted for in a
structural de-duplication, but removing only them cannot close the
order $T\log^2T$ second-moment obligation. Analogously, Parseval
reduces (2) to the integrated squared Fourier residual after the
periodic main term; relative pointwise power-log savings are far
above this scale. The exact collision bound remains open.

The [replay](verify_density.py) computes the ordered pair intersections and
independent exact residue-centered variance accounting at four dyadic
sizes. These are validation controls for the identity, not evidence that
the required asymptotic is attained.
