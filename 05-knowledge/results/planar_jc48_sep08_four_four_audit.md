# Independent audit of the distinguished 4+4 quartic exclusion

**Status: full independent analytic/source audit PASS; all 96 normal,
optimized and frozen gates agree.** This accepts
[the distinguished four-four theorem](planar_jc48_sep08_four_four.md)
for its actual fixed DG section space and unrestricted polynomial mate.
The conclusion is not a rational-mate exclusion. No arbitrary pair of
boundary points is normalized to zero and infinity.

## 1. Complete entry and the actual inversion

The complete binary octic is `alpha*x^4`, `alpha!=0`: as a degree-eight
section on the graph, it has multiplicity four at the distinguished
finite point zero and multiplicity four at infinity, with no other
zeros. It is not merely a local fourth-order vanishing hypothesis.

The full octic restriction map from the fifteen-dimensional `(4,2)`
box has rank nine. Its kernel is the defining section of the graph
times the full six-dimensional `(2,1)` box. Since `w=x^2*t` is
global and `w^2` has boundary numerator `x^4`, every permitted
approximate root is precisely

    H=alpha*w^2+K, K in L1,
    L in L1,
    L1=span{1,t,xt,w,v,h},
    v=x+x^3*t, h=x^2+x^4*t.

This gives the complete twelve lower coefficients in the primary.
The source verifies both the restriction rank and a unit determinant
for this actual L1 basis; all numerator boxes are retained. The proof
places no bound on the degree of a proposed polynomial mate. Constant
`L` is correctly separated using the nonunit factor `2H` in its
Jacobian; `H` is nonconstant because `alpha!=0`.

The inversion used to study infinity is the genuine automorphism of
the ambient product `P1 x P1`, taking `(x,z)` to `(1/x,1/z)` and
preserving the graph. In its affine chart,

    u=1/x, T=-x^2-x^4*t,
    dx wedge dt=u^2 du wedge dT.

The old six basis elements become exactly

    1, -h_new, -v_new, -1-w_new, -uT, -T.

Moreover `w_old^2=(1+w_new)^2`. Thus the full distinguished 4+4
section family is preserved with changed coefficients and the same
nonzero alpha. The original source-volume weight is `u^2`. On the
new second chart that weight cancels the usual `r^2`, giving
`dr wedge db`. These are literal identities. The proof does not
transport polynomiality of an arbitrary mate to a different affine
source, or apply an unweighted local conclusion to the weighted form.

The proof uses the currently proved common-root and shared-root local
lemmas, together with the global section spaces. The separate
infinity-octuple result is motivational context only; no conclusion
from it is needed as a dependency here.

## 2. Complete regularity at the weighted infinity point

In the new finite coordinates the original form is

    eta=u^2*s^2 du/E_s,
    E=N^2+s^3 M-cs^4,

with the complete N and M displayed in the primary. I checked each
case and its local sheet count. In particular a shared-root normal-unit
lemma is not being silently applied outside its domain: the M-unit
case is also regular, as follows either from the proved M-unit analysis
or directly as below.

If `kt!=0`, `N=0` has an exact analytic centre `s=psi(u)` of
order four. At that centre `ell=ord(M-c psi)` lies generically in
`0,...,4`, including ell zero when M is a unit. The two fibre
determinations satisfy `ord s=4`, `ord N=6+ell/2`, and `N_s`
is a unit. Its product with N dominates the other terms of `E_s`.
The unweighted differential has exponent `2-ell/2` in u, hence is
regular on either actual normalization; multiplying by `u^2` improves
it. The exact Weierstrass degree is two, so there are no missing
branches in this case.

If `kt=0`, `lt!=0`, the local degree is three. Its boundary order
is `m=4`, so the exceptional balanced case `m=3j` cannot occur.
For normal order `j=1`, the low face has one nonzero simple root at
`s` of order two, and two further determinations cancel N at order
three. Both portions are unweighted regular. The high determinations
have unweighted exponent `1/2` in u and combine into a branch with
`u=tau^2`, where the differential has order two. For `j>=2`,
including an identically zero first normal coefficient, the simple
cubic face at `s` of order `8/3` is `alpha^2+lt Z^3`.
On its normalization `u=tau^3`, `s=tau^8` times a unit, `E_s`
has order sixteen and `s^2 du` has order eighteen. These are regular
unweighted and hence also weighted. One-plus-two or the three cubic
determinations exhaust degree three.

Now assume `kt=lt=0`. The exact specialization is
`E(0,s)=(k0^2+l0-c)s^4`, so the generic local degree is four.

- If `kxt!=0`, the two low simple nonzero roots at `s=uZ` have
  unweighted logarithmic forms, made regular of order one by `u^2`.
  The exact centre for the other determinations has order three,
  `N_s` of order one, and `1<=ell=ord(M-cs)<=3` generically.
  The weighted exponent is `(5-ell)/2`. It is positive on the
  actual normalization, including the ramified ell-two case.
- If `kxt=0`, `lxt!=0`, there is one such low simple root and
  the remaining cubic face is `alpha^2+lxt Z^3` with `s` of
  order `7/3`. On the actual normalization `u=tau^3`,
  `s=tau^7` times a unit, the orders of `E_s` and the full
  numerator `u^2*s^2 du` are seventeen and twenty-two. The
  relative differential therefore has order five.
- If `kxt=lxt=0`, the face at `s=u^2 Z` is

      (alpha+k2 Z+k0 Z^2)^2+l2 Z^3+(l0-c)Z^4.

  It has generic degree four and nonzero constant term alpha squared.
  Any nonzero multiple root must satisfy `ZP0'-4P0=0`, whose
  constant term is `-4alpha^2`; hence only finitely many fibre
  values are exceptional. All four roots are generically simple and
  nonzero. Each gives a smooth branch with `ord_u E_s=6`, exactly
  matching `ord_u(u^2 s^2)=6`, so its differential is regular.

The high/low determinations in each row exhaust the declared local
Weierstrass degree, rather than being a list of selected branches.
All higher coefficients of N and M occur in the exact equations and
are either retained in the face or shown to have strictly greater
order. Thus every normalized branch at the original infinity root is
regular for every coefficient choice under the stated generic-fibre
qualification.

## 3. Finite reductions and compact exactness

At the original finite point, `lt!=0` gives the unweighted M-unit
m-four regularity just discussed; `kt!=0` gives normal-unit regularity.
The infinity point is always weighted-regular by Section 2. There are
no further boundary roots. The relative form is regular inside W on
a generic fibre and nonzero on the original source part of a generic
component. Such a component cannot lie entirely on the added divisor:
if F is constant there, that is at most one exceptional fibre value;
otherwise that divisor is not a fibre component. Therefore the form
is not identically zero on a compact generic component.

A rational mate restricts meromorphically to that compact normalization
and would have derivative equal to the relative form. A regular exact
meromorphic differential on a compact curve vanishes, a contradiction.
This forces `lt=kt=0`. The argument is componentwise and does not
assume generic irreducibility or global regularity of the proposed mate
on the enlarged surface.

The finite low nonzero tangent roots for `kxt!=0` or `lxt!=0`
give nonzero unweighted logarithmic residues. One such actual point
already excludes a rational primitive; residues at other points cannot
cancel that local obstruction. Thus `kxt=lxt=0`. At this stage the
remaining four local branches correspond exactly to the four generic
roots of the quartic `f(w)=c` below; generic c also excludes zero
roots so the `w=1/Z` relation is valid.

## 4. Full residue identity and its coefficient consequences

In the original rational chart `w=x^2t`,

    F=f(w)+x g(w)+O(x^2),
    A=alpha*w^2+k2*w+k0,
    f=A^2+l2*w+l0,
    g=(1+w)(2k3*A+l3),
    omega=x^(-2) dx wedge dw.

At a simple generic root `w0`, the actual fibre displacement is
`w-w0=-g(w0)x/f'(w0)+O(x^2)`. Consequently

    Res eta=(g'f'-f''g)(w0)/f'(w0)^3.

I independently checked that this includes the induced displacement of
the denominator, rather than differentiating at a fixed approximate
root. Higher source rows cannot alter this residue.

Exactness for a putative fixed rational mate holds for generic original
fibre values. As those values vary, the roots cover all but finitely
many w, so the polynomial `g'f'-f''g` must vanish identically.
Equivalently `(g/f')'=0`, giving `g=C f'` with `C in C`.
The nonzero leading coefficient alpha makes `f'` cubic.

Independent coefficient comparison gives exactly:

- When `k3=0`, the degree comparison forces `C=0` and then `l3=0`.
  The original source derivatives of F both vanish along `x=0`.
  This excludes polynomial mates, without claiming a rational
  critical-point obstruction.
- When `k3!=0`, the leading coefficient fixes `C=k3/(2alpha)`.
  The successive remaining coefficient differences are
  `k3(2alpha-k2)`, then `l3`, then `-k3*l2/(2alpha)`.
  Thus `k2=2alpha`, `l3=l2=0`. Nonconstant L now means `l4!=0`.

No coefficient is omitted from these comparisons. Residue compatibility
is only a necessary condition; the proof correctly pays the remaining
polynomial obstruction separately.

## 5. Actual-source critical points exhaust the residue-compatible class

Let `q=1+x^2t`, `c0=k0-alpha`. The surviving exact function is

    H=c0+q(alpha*q+k3*x+k4*x^2),
    F=H^2+l4*x^2*q+l0,
    alpha*k3*l4!=0.

On `x!=0`, `(x,q)` is a regular coordinate chart on the original
source, with determinant `x^2` and inverse `t=(q-1)/x^2`.
Every critical point constructed there is therefore an actual affine
critical point of the original polynomial.

At `q=0`, one derivative vanishes identically and the other is

    x[2c0*k3+(2c0*k4+l4)x].

For `c0!=0`, `2c0*k4+l4!=0`, its root
`x=-2c0*k3/(2c0*k4+l4)` is nonzero. It yields the required actual
critical point. This point lies on a legitimate source curve; `q=0`
is not removed from the original affine plane.

For the remaining cases the choice `q=-k3*x/(4alpha)` is structural.
Independently eliminating the two derivatives gives

    x*F_x-2q*F_q=-2q*H*(4alpha*q+k3*x).

On that selected line, with `Cx=k3+2k4*x`, the two derivatives are
respectively `2q(H Cx+l4*x)` and `x(H Cx+l4*x)`. Thus the primary's
one-variable equations indeed force both derivatives to vanish.

If `c0=0`, the equation is

    k3*x*(3k3+4k4*x)*(k3+2k4*x)=16alpha*l4.

It is cubic with nonzero leading coefficient when `k4!=0`, and remains
nonconstant linear when `k4=0` because its linear coefficient is
`3k3^3`. The right side is nonzero. The fundamental theorem of algebra
therefore supplies a root; every root has nonzero x and nonzero displayed
linear factors. The resulting q is nonzero, and the common derivative
factor is exactly zero.

If `c0!=0`, `l4=-2c0*k4`, then `k4!=0`. The equation is

    x^2*(3k3+4k4*x)*(k3+2k4*x)=16alpha*c0.

Its leading coefficient is `8k4^2!=0`; the right side is nonzero.
Again a complex root exists and avoids all forbidden factors. Direct
substitution makes the same common derivative factor vanish. These two
exceptional cases together with the q-zero case exhaust the complete
coefficient space. No numerical search, assumption of a real root,
or unstated division by a possibly zero coefficient enters the proof.

Every resulting source critical point contradicts a polynomial
nonzero constant Jacobian. This proves the stated no-polynomial-mate
conclusion, with unrestricted mate degree.

## 6. Hostiles, exact controls and acceptance

The residue-compatible control `alpha=k3=l4=1`, `k2=2`, `k0=1`,
`k4=0` has the actual critical point
`x=16/3`, `q=-4/3`, `t=-21/256`. Both original source derivatives
vanish there. Thus the residue identity does not falsely imply a mate.

The stronger scope control is also correct: with `H=w^2`, `L=w`,

    F=w^4+w, G=1/[x(4w^3+1)], J(F,G)=1.

It is inside this exact global 4+4 class and has a rational mate with
poles. Its source criticality only excludes polynomial mates. Constant-L
composites give another boundary check. Neither control is discarded
by an unstated generic-coefficient hypothesis in the main theorem.

I read the complete standalone source and independently ran

    python3 -B 04-computation/planar_jc48_sep08_four_four.py
    python3 -B -O 04-computation/planar_jc48_sep08_four_four.py

Both completed with **96 always-active gates** and matched the frozen
405 bytes exactly. The producer imports no inherited mathematical
implementation. The generic branch and root-existence claims are paid
by the analytic arguments above, not the source's named sample fibres
or arithmetic order checks. A separate symbolic reconstruction checked
the critical-locus elimination, both exceptional polynomial equations
and all four coefficient comparisons in the residue identity.

Frozen pins:

- [Source](../../04-computation/planar_jc48_sep08_four_four.py),
  10,853 bytes, SHA256
  `f03a0843dd2521a35b0df92481a05db8af0c8130e210a5a5fb315afb9d5e8089`.
- [Output](planar_jc48_sep08_four_four.out), 405 bytes, SHA256
  `a3c7c6823ea672b09a459803e75111f1392e810d13d8897d10882724e92ccb90`.
- Semantic digest
  `e6de99bdcde1e30ca8f3d177c22ec9b24f73de73ade3eaff6a5e47f1f7e49eae`.

No mathematical or source correction was needed. The candidate is
accepted for owner-controlled promotion, with the distinguished
finite-zero/infinity boundary pair and polynomial-only conclusion
preserved. Arbitrary pairs of boundary points, other partitions and
the general quartic Keller problem remain outside its scope.
