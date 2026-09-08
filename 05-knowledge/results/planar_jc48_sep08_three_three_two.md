# Every (3,3,2) quartic boundary placement

**Status: PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED.**
This is a complete boundary stratum on the fixed DG surface, not a
reduction of arbitrary Keller maps to that surface. JC(2) remains open.

## 1. Statement and inheritance

Work over `C` with

    W=(P1_x x P1_z) minus S, S={z=x^2}, t=1/(z-x^2),
    x=1/r, t=-r^2-r^4 b_D,
    omega=dx wedge dt=r^2 dr wedge db_D, D={r=0}.

Let `H in L2`, `L in L1`, `F=H^2+L`. Assume that the leading
binary octic of `H` has exactly three distinct zeros with
multiplicities three, three and two.

**Theorem.** A rational mate of `F` forces the double
point to be the original infinity point and `L` to be constant.
Every such `F` has no polynomial mate of any degree. There are
actual constant-`L` rational mates with the double point at
infinity and every finite triple-point position.

The closest proved mechanism is the constant-`D` reduction by
primitive pole capacity in
[the fourfold-infinity (4,3,1) proof](planar_jc48_sep08_four_three_one_infinity.md).
The essential local suppliers are
[M-unit Newton/Morse forms](planar_jc48_sep08_quartic_common_root.md),
[complete shared-root jets and weighted multiplicity two](planar_jc48_sep08_shared_roots.md),
[leading inverse exactness](planar_jc48_sep08_leading_exactness.md),
and [the complete global filtration](planar_jc48_sep08_dg_quadratic.md).
The nonsquare-leading noncomposition argument is the one independently
audited in [the all-finite (4,3,1) proof](planar_jc48_sep08_four_three_one.md).

The live concepts are weighted boundary forms, primitive capacity,
complete global rows, Riemann--Roch membership and the generic fibre
constant. The source-to-target map restricts the actual volume to the
original generic fibre; it preserves a mate's exact differential.
Local pole data alone lose global sections and ordinary affine points,
so the complete primitive-space and chart checks in Section 5 are
necessary sidecars. The sharp rational family below is the canonical
hostile. A second hostile is the single exceptional level on which the
final bracket becomes constant; specializing to that level destroys
the generic-fibre obstruction.

## 2. The complete coefficient universe and all locations

If the double point is finite, leading exactness requires
`du/sqrt(N)` exact on its actual quadratic algebraic field.
At a double zero, `sqrt(N)=u*unit`, so this differential has a
nonzero residue on each normalization branch. This forbids every
rational mate. It also covers a triple point at infinity.

It remains to put the double point at the original infinity.
Actual source/surface scalings normalize the two finite points'
nonzero separation and the leading scalar. Write `u=x-p` with
arbitrary original finite position `p`. Naming the active triple
first is merely a choice of the two actual points before this
normalization. This ordinary polynomial coordinate translation
is not an automorphism of the compactification.

Now `N=u^3(u-1)^3`, and the entire global section space is

    H=N t^2+P t+Q,      P=sum_(i=0)^4 A_i u^i,
    Q=(A4-1)u^2+(A3-2p A4+4p+3)u+C0,
    L=M t+R,           M=sum_(i=0)^4 B_i u^i,
    R=B4 u^2+(B3-2p B4)u+E0.                         (1)

For completeness, in original `x=u+p` coordinates the numerator
rows of `H` are `Q`, `P-2Qx^2`, `N-Px^2+Qx^4`, each of degree
at most four. The coefficients of `u^6,u^5` in the last row give
exactly the two displayed coefficients of `Q`; all lower rows
then satisfy the degree bound. Conversely the complete section
box forces `deg Q<=2` and `deg P<=4`. The linear numerator rows
`R,M-Rx^2`, each of degree at most two, give exactly `R` in (1).
Thus these are necessary and sufficient global rows, not a chosen
subfamily. The companion verifies them with symbolic `p`.

## 3. Weighted infinity and the one-active-triple entry

At the original infinity the leading section has order two:

    r^8 N(1/r-p)=r^2(1-pr)^3(1-(p+1)r)^3.

The entire relative form is regular on every normalized generic
branch there, for all coefficients in (1). A normal unit is
regular by the shared-root lemma. A lower-section unit is
regular by the M-unit lemma, since `m=2` cannot be its balanced
multiple `3j`. In the remaining shared case with vanishing normal
unit the full tangent equation has the form

    (a r^2+b r s+h s^2)^2+s^3(l r+e s)-zeta s^4,
    a!=0.

Putting `s=r Z` gives four simple nonzero roots generically.
For `P0=(a+bZ+hZ^2)^2+lZ^3+eZ^4`, the repeated-root eliminant
`Z P0'-4P0` has constant `-4a^2`, so only finitely many fibre
values could produce a multiple root. These four branches
exhaust the local Weierstrass degree. Their unweighted form has
order minus one; the actual multiplier `r^2` makes its order one.
This is precisely where a finite-point residue obstruction cannot
be transported to infinity.

At either finite triple, an M-unit gives regular forms or forbidden
nonzero logarithms. A shared normal unit is regular. Otherwise
the finite first-jet obstruction forces

    P=P'=M=M'=0 at that triple.                      (2)

Call this last possibility active. At an active triple, the
substitution `u=tau^2,s=tau^3 Z` has four nonzero simple
determinations paired into two actual normalized branches.
Each relative differential has order minus two; its primitive
can have only a simple pole there. Thus one active triple gives
total primitive pole capacity two, independent of coefficient
degenerations. The leading tangent quartic has nonzero constant
and generic fibre coefficient, which pays that exhaustion.

Assume a rational mate and `L` nonconstant. There must be an
active triple: otherwise a primitive is holomorphic, hence
constant, on every compact normalized generic component, contrary
to its nonzero differential on the original source.

There cannot be two active triples. Their four prescribed zeros
and `deg M<=4` give `M=lambda u^2(u-1)^2`. The complete inverse
coefficient `T2=-M/(4N)` and normalized field trace require
`M du/N` rational-exact. But this is
`lambda du/[u(u-1)]`, with residues `-lambda,lambda`. Hence
`M=0`, which in the complete global rows also makes `L` constant.

Name the unique active triple `u=0`. The total primitive pole
degree is at most two, also on any one component. If `F|D` were
nonconstant, a generic transverse point on the original divisor
`D` would have differential order two, from the actual volume
factor `r^2`. Its primitive would have local degree three,
exceeding its pole degree. Therefore `F|D` is constant. This
entry needs no irreducibility assumption.

## 4. The full residual family and its genus

The slopes of `H|D,L|D` in `b_D` are `2-A4,-B4`.
The square of a nonconstant affine function cannot be canceled
by an affine function. Thus constant `F|D` gives `A4=2,B4=0`.
Combined with the active jets, (1) becomes

    P=2u^4+a u^3+b u^2,    Q=u^2+(a+3)u+c,
    M=u^2(lambda u+mu),    R=lambda u+e.             (3)

All `p` remain included: their disappearance in (3) is the result
of intersecting the complete original global rows. It is not an
assumed surface translation. The residues of `M du/N` are
`-mu,mu`, so `mu=0`. Nonconstant `L` means `lambda!=0`.
Absorb the additive `e` into the formal fibre constant `zeta`.

At `u=1`, the lower section is now a unit. If `P(1)!=0`, the
normal-unit triple has one relative-form zero of order four,
forcing primitive local degree five, again exceeding capacity
two. Hence

    b=-a-2.                                         (4)

For `P(1)=0`, the M-unit triple has unit relative forms in every
case admissible under a mate. Its three Puiseux determinations
have `s~(u-1)^2`; the leading cubic is
`(a0+b0 Z)^2+m0 Z^3`, with `a0,m0!=0`. Simple roots give
order zero. At a repeated root the inherited Morse analysis
has contact one or two: contact one gives a ramified unit form,
whereas contact two gives a forbidden nonzero logarithm. The
generic fibre coefficient pays the upper bound two. Thus every
admissible case contributes no zeros or poles.

Geometric integrality must be paid before counting the canonical
divisor. A polynomial decomposition of `F` has outer degree two
or four, by its degree four in `t`. A quadratic outer polynomial
can be completed to a square; uniqueness of the square prefix
in characteristic zero then forces `L` constant. An outer
quartic with an inner polynomial linear in `t` would make `N^2`
a fourth power, hence make `N` a polynomial square. Its two
odd triple multiplicities forbid that. There is consequently no
nontrivial polynomial composition. The classical closed-polynomial
criterion, with the precise source and hypotheses supplied in
the cited all-finite (4,3,1) proof, gives geometric integrality
of the generic fibre. Work over the algebraic closure of `C(zeta)`.

At infinity constant `D` gives four simple branches `s=rZ` as
in Section 3; their form orders are one. The two active points
over `u=0` have orders minus two each. The other triple gives
order zero, and the constant divisor `D` misses the generic
fibre. On the original smooth affine source the volume and
generic fibre are regular and nonzero, so there are no further
zeros or poles. Thus the canonical degree is `4-4=0`, and the
compact normalized generic fibre has genus one.

## 5. Complete primitive space and direct contradiction

Let `E` be the sum of the two points over the active triple.
A mate's primitive belongs to `L(E)`, and genus one gives
`dim L(E)=deg E=2`. The complete basis is

    1, J=u(u-1)^2 t+u.                              (5)

Here are all its regularity checks. At either active point,
`u~tau^2,t~tau^-3`, so `J` has exactly one simple pole.
At the inactive triple `t~(u-1)^-2`, canceled by its double
factor. It has no affine denominator. It is an actual global
section in `L1`: the linear numerator rows in Section 2 hold
with `M_J=u(u-1)^2,R_J=u`. In the original second chart it is
a polynomial in `r,b_D`, has value `2p+2` on `D`, and on each
deleted infinity branch `b_D=1/(rZ+O(r^2))` it tends to

    2p+2-1/Z.

The tangent roots are nonzero, so these are finite. This checks
the deleted boundary point as well as the retained divisor.
Finally `J` is nonconstant on the generic fibre: it has actual
simple poles, or alternatively its degree one in `t` cannot
vanish modulo the irreducible degree-four fibre equation.
Riemann--Roch therefore proves completeness of (5), not merely
membership of two trial functions.

A putative mate must be `A(zeta)J+B(zeta)` on that fibre,
with coefficients allowed even in the algebraic closure of
`C(zeta)`. Its Jacobian is `A(zeta)J_(u,t)(F,J)`. Reducing
this bracket modulo `F-zeta`, the coefficient of `t^2` is

    2c u^3(u-1)^4.

Since the bracket must equal one, `A!=0`, and hence `c=0`.
The entire remaining bracket is then

    2zeta-(2zeta+lambda)u.                          (6)

It is nonconstant over the formal generic fibre field, because
`2zeta+lambda!=0`. No scalar `A(zeta)` makes (6) equal one.
This contradiction excludes every rational mate when `L` is
nonconstant. Notice the hostile specialization `zeta=-lambda/2`:
then (6) is the nonzero constant `-lambda`. A single special
fibre must not replace the original generic one.

If `L` is constant instead, `J(F,G)=2H J(H,G)` cannot be a
nonzero constant for polynomial `G`, since `H` is nonconstant.
Together with Section 2, this proves the polynomial exclusion
at every boundary placement.

## 6. Exact rational boundary, verification and scope

Put

    A0=u(u-1), Z0=A0 t+1,
    H=A0 Z0^2+c,
    G_H=(2u-1)/(A0 Z0).

For every `p`, the original global coefficient constraints hold,
the leading polynomial is `u^3(u-1)^3`, and `J(H,G_H)=1`.
Therefore `H^2+e` has rational mate `G_H/(2H)`. The denominator
is essential. This proves the claimed rational sharpness without
producing a polynomial mate or an arbitrary Keller envelope.

The [exact companion](../../04-computation/planar_jc48_sep08_three_three_two.py)
checks the entire symbolic global row space, its necessary high
coefficient equations, both-active and single-active residues,
the actual infinity limits, every coefficient of the generic
bracket, the special-fibre hostile and the all-p rational family.
No arbitrary coefficient box is enumerated. Reproduce with

```sh
python3 -B 04-computation/planar_jc48_sep08_three_three_two.py
python3 -B -O 04-computation/planar_jc48_sep08_three_three_two.py
```

Both normal and optimized runs pass **56 gates** and reproduce
the same [268-byte frozen output](planar_jc48_sep08_three_three_two.out).
Source: 5,623 bytes, SHA256
`eabf9adaeeeefe68c543fab37963af1c14d75d9e39f3ee34b964963617de3ea0`.
Output SHA256:
`37f4bd8ad522d01ccaf93a20531025e0a7b3e3645b1e91286eb6503bc458c03e`.
Semantic gate SHA256:
`facd85633ee6f36fb3d7851442e670d799ea602633b65301389f16ed1403dd11`.

The [independent full audit](planar_jc48_sep08_three_three_two_audit.md)
accepts all original global rows, every normalized branch and complete
primitive-space membership. It includes the repaired Morse contact-one
case and the special-fibre hostile; both56-gate replays agree.
