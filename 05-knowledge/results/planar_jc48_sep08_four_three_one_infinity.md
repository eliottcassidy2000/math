# The remaining (4,3,1) location: fourfold infinity

**Status: PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED.**
This concerns complete global square-prefix sections on the fixed DG
surface. No arbitrary quartic, projective coordinate transport, or JC(2)
closure is claimed.

## 1. Precise statement and inherited mechanisms

Use

    W=(P1_x x P1_z) minus S,  S={z=x^2},
    t=1/(z-x^2),  x=1/r,  t=-r^2-r^4b_D,
    omega=dx wedge dt=r^2 dr wedge db_D,  D={r=0}.

Let `H in L_2`, `L in L_1`, `F=H^2+L`. Suppose the boundary
octic of `H` has a fourfold zero at infinity and two distinct
finite zeros of multiplicities three and one.

**Theorem.** A rational mate of `F` forces `L` to be
constant. Consequently `F` has no polynomial mate of any degree.
The rational boundary is sharp: Section 6 gives actual rational
mates with constant `L` in this same class, for every finite
triple-point position.

The closest mechanisms are the complete coefficient filtration,
[the weighted order-four infinity lemma](planar_jc48_sep08_four_four_transport.md),
[shared-root first jets](planar_jc48_sep08_shared_roots.md),
[M-unit local bounds](planar_jc48_sep08_quartic_common_root.md),
and [leading inverse-series exactness](planar_jc48_sep08_leading_exactness.md).
The needed second inverse coefficient is derived directly below,
as in [the all-finite companion](planar_jc48_sep08_four_three_one.md).
The source-to-target map retains the original fibre, the actual
volume, and all induced global coefficients. The constant-`D`
condition changes the degree of the lower coefficient, making
the final residue count stronger than a boundary multiplicity
count alone.

Live concepts are weighted infinity, same-component pole degree,
full global rows, all inverse coefficients, and field trace.
The corrected near miss is to analyze a residual quadratic field
before checking the residue of the already available `M/N`.
Section 5 preserves that independent trace obstruction but does
not treat it as necessary work for the shorter proof. The sharp
rational family in Section 6 prevents replacing the final
polynomial conclusion by a blanket rational exclusion.

## 2. Complete coefficients and the actual infinity weight

Normalize the nonzero separation and leading scalar by actual
source/surface scalings. Write `u=x-p`, retaining arbitrary `p`.
The leading polynomial is

    N=u^3(u-1).

The translation to `u` is only an ordinary polynomial coordinate
change. It is not an automorphism of `W`. The full global rows,
derived in the original `x` coordinates from
[the DG filtration](planar_jc48_sep08_dg_quadratic.md), are

    P=A4 u^4+A3 u^3+A2 u^2+A1 u+A0,
    Q=A4 u^2+(A3-2p A4)u+C0,
    M=B4 u^4+B3 u^3+B2 u^2+B1 u+B0,
    R=B4 u^2+(B3-2p B4)u+E0,
    H=N t^2+P t+Q,  L=M t+R.                         (1)

Let `s=z-x^2`, `calN=N+sP+s^2Q`, `calM=M+sR`. Thus
`H=calN/s^2`, `L=calM/s`; `N,M` always denote leading
polynomials, whereas normal derivatives concern `calN`.

At the original infinity point the entire relative differential
is regular on every normalized generic branch, for all coefficients
in (1). This is a local use of the audited infinity lemma, not
an inference from the earlier 4+4 global theorem. To make its
scope explicit, the actual surface involution gives

    r=1/x,  T=-x^2-x^4t,  omega=r^2 dr wedge dT.

Its transformed leading polynomial is
`r^8 N(1/r)=r^4(1-pr)^3(1-(p+1)r)`, order four times a unit.
Near that point write the entire section equation

    calN=r^4 a(r)+s B(r)+s^2 C(r,s),  a(0)=1,
    calM=D(r)+s E(r,s),
    Phi=calN^2+s^3 calM-zeta s^4,
    eta=(unit) r^2 s^2 dr/Phi_s.                     (2)

Here is the complete generic case split from the cited proof.
A normal unit gives two regular unweighted branches. An M-unit
with normal order positive has `m=4`, which is never balanced
with `3j`; all its branches are regular. The additional factor
`r^2` can only improve these cases.

Otherwise the local Weierstrass degree is four. If `ord B=1`,
there are two low determinations `s~r`, weighted form order one,
and two cancelling determinations `s~r^3`. The order `ell` of
`calM-zeta s` at their analytic centre is generically between
one and three; the `-zeta s` term pays the upper bound. Their
weighted form is `(unit)r^((5-ell)/2)dr`, regular on the actual
normalization. If `ord B>=2` and `ord D=1`, there is one low
branch and three determinations with `(ord r,ord s)=(3,7)`;
the latter have weighted differential order five. Finally, if
both first jets vanish, put `s=r^2Z`. The complete tangent
quartic is

    (1+B2 Z+C0 Z^2)^2+D2 Z^3+(E0-zeta)Z^4.

Its four roots are simple and nonzero generically: the repeated
root eliminant `ZP0'-4P0` has constant `-4`. On every branch
the weighted form has order zero. Higher terms of the unit
`a(r)` have higher weight in every face and do not change the
centre-order bound. These cases exhaust (2), including every
coefficient degeneration and the complete infinity chart.

## 3. Finite pole capacity forces the constant-D coefficient stratum

At the finite simple point every generic branch is regular:
use the M-unit lemma when it is not shared, and the shared
simple-root lemma otherwise. At the finite triple point, an
M-unit is regular or has only forbidden nonzero logarithms
in its balanced case. A normal unit is regular. If the triple
point is shared and its normal derivative vanishes, the finite
first-jet obstruction requires

    P(0)=P'(0)=M(0)=M'(0)=0.                          (3)

If the triple point is not active in this last sense, the
compact generic relative form has no possible primitive poles:
any logarithm already forbids a rational derivative, and all
remaining branches are regular by Section 2. A primitive would
be constant on each compact normalized component, contradicting
the nonzero relative form on the original smooth source part.
Thus a hypothetical rational mate must satisfy (3).

At an active triple point, write `calN=a u^3+s O(u^2)+s^2 C`
and `M=O(u^2)`. The rescaling `u=tau^2`, `s=tau^3 Z` gives
the leading quartic

    (a+C(0)Z^2)^2+(R(0)-zeta)Z^4.

It has four simple nonzero determinations, paired into two
actual normalized branches, and each differential has order
`-2`. The local Weierstrass degree four pays exhaustion.
Its total possible primitive pole degree is exactly two;
there are no possible primitive poles at the other boundary
points by Section 2 and the simple-root lemma.

If `F|D` were nonconstant, a generic finite point on `D` would
be a transverse point of the original fibre. The factor `r^2`
in the actual source volume makes `eta` vanish to order two
there, forcing primitive local degree three. That exceeds even
the total possible pole degree two, hence also the degree on
the same component. No geometric irreducibility is assumed.

Therefore `F|D` is constant. From the original full rows,
`H|D` is affine in `b_D` with slope `-A4`, and `L|D` is
affine with slope `-B4`. The leading quadratic term of a
nonconstant square cannot be cancelled by an affine function.
Consequently

    A4=B4=0.                                        (4)

This is the essential global reduction: it gives `deg M<=3`.
It is not obtained by moving the infinity point or translating
the surface.

## 4. The second inverse coefficient finishes the class

For a rational mate, the complete inverse expansion
`F(x,T)=v^-4` has `T2=-M/(4N)`. Its coefficient in the mate
satisfies `g6'=T2/2`. Normalized field trace to `C(x)` therefore
makes `M dx/N` rational-exact. The direct coefficient derivation
is the one given in Section 2 of the all-finite companion: the
`v^3` equation is `4a T2+M/a=0`, `a^2=N`. This consequence
does not require geometric connectedness of any special fibre.

By (3)--(4), `M=u^2(B3 u+B2)`. At the simple root `u=1`,
rational exactness forces `B2=-B3`. Hence

    M=lambda u^2(u-1),  M/N=lambda/u.

The nonzero residue at the triple point forces `lambda=0`.
Equivalently, exactness requires a triple zero there and a
further simple zero at `u=1`, too many for `deg M<=3`.
Thus `M=0`, and full globality in (1) makes `L` constant.
For a polynomial mate, the nonconstant factor `2H` in
`J(H^2+constant,G)` then gives the final contradiction.

## 5. Independent residual-field obstruction and its lost-data contract

This section is not needed after the last residue in Section 4.
It records a second, exact mechanism for the intermediate
residual with `lambda!=0`, if one keeps only the simple-point
residue. The full coefficients are

    H=u^3(u-1)t^2+u^2(a u+b)t+a u+c,
    Z=u^2(u-1)t+u,
    F=H^2+lambda Z.                                 (5)

The parameter `p` disappears because the actual complete
constant-`D` section intersection was calculated in (1)--(4).
This does not assert that source translation extends to `W`.

Put `h=H`, `z=Z`, and

    A=az+1-a-b+c-h,  B=(b-2)z+h-c.

The exact relation and Jacobian are

    A u^2+B u+z^2=0,
    y=2Au+B,  y^2=Delta=B^2-4Az^2,
    J_(u,t)(H,Z)=u y.                               (6)

Since `t=(z-u)/(u^2(u-1))`, these are actual inverse field
relations: `C(u,t)=C(h,z)(u)`. The extension has degree two
for every `a,b,c`. First, the coefficient of `t^2` in the
actual `J(H,Z)` is `-u^4(u-1)(2u-1)`, so the map is dominant
for every parameter and `h,z` are algebraically independent.
Then `Delta` is monic quadratic in `h`,
with discriminant

    disc_h Delta=16z^2(z-1)(z+a+b-1),

a nonzero polynomial for every parameter value. Thus it is
not a square in `C(z)(h)`: a rational square that is a monic
quadratic polynomial would be a polynomial square and have
zero discriminant. In particular this map is not declared
birational, and no fibre genus or connectedness is assumed.

For the derivation along the original generic fibre `F=zeta`,
use `z=(zeta-h^2)/lambda`. The relative form is

    eta=-dh/(lambda u y)
       =(1+B/y)dh/(2lambda z^2).

Trace through the actual quadratic field, with involution
`y -> -y`, gives

    Tr eta=dh/(lambda z^2)
          =lambda dh/(zeta-h^2)^2.                   (7)

Differentiation along fixed `zeta` commutes with field trace.
After the harmless extension of the constant field by
`rho^2=zeta`, (7) has residue `-lambda/(4rho^3)` at `h=rho`,
which is nonzero. It therefore cannot be the differential of
the rational function `Tr G`. This independently excludes every
nonconstant-lower-row member of (5). It preserves the original
derivation and all field degrees; retaining only its quadratic
equation without the source volume would not pay (7).

## 6. Sharp rational controls and all-location integration

For every finite triple point `p`, put

    u=x-p,
    H=u^3(u-1)t^2+c,
    G_H=-1/(u^2t).

Then `J(H,G_H)=1`. The original second chart is

    H=(1-pr)^3(1-(p+1)r)(1+r^2b_D)^2+c,

so `H` is global and has precisely the declared binary
partition `(4 at infinity, 3 at p, 1 at p+1)`. Therefore
`F=H^2+c0` has rational mate `-1/(2u^2tH)`. These are genuine
rational controls, not polynomial mates or Keller examples.

The exact producer checks full shifted global rows, the literal
weighted infinity transformation and complete local leading
equations, the counted-residue elimination, both directions of
the residual field map, its nonsquare discriminant and traced
residue, and the same-partition rational family. It does not
replace the analytic branch completeness or degree argument
with a finite scan. Reproduce with

```sh
python3 -B 04-computation/planar_jc48_sep08_four_three_one_infinity.py
python3 -B -O 04-computation/planar_jc48_sep08_four_three_one_infinity.py
```

Both replays pass **87 gates**, with exactly the same 399 bytes
as [the frozen output](planar_jc48_sep08_four_three_one_infinity.out).
The [source](../../04-computation/planar_jc48_sep08_four_three_one_infinity.py)
has 10,605 bytes and SHA256
`b4bf39bad5a819d2435379795d4bee97fbbd6c79e271b507192bd4d3217df77c`.
Output SHA256:
`92a0617b4594bc84c293a743ac79d6921310400924f06e8341c7894363e05268`.
Semantic gate digest:
`a747eb90da21add105fe24156322e989d0f302a74337e88862db4aee6ab5f057`.

Root independently reconstructed the quadratic trace and its
nonzero residue. The [independent complete audit](planar_jc48_sep08_four_three_one_infinity_audit.md)
accepts the full local/global proof, weighted infinity faces, counted
residues, exact field map, controls and both87-gate replays.

Combined with the separately audited all-finite companion and
the proved leading classification for a triple or simple point
at infinity, the result excludes polynomial mates for every
location of a binary `(4,3,1)` divisor. The all-finite companion
is now PROVED and independently audited in checkpoint be2e74d62f;
this all-location corollary is therefore fully discharged.
No claim for other binary partitions is made.
