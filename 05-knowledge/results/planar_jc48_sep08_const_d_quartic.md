# The fixed finite-six, constant-D quartic stratum has no rational mate

**Status: PROVED / FINITE-EXACT / INDEPENDENTLY AUDITED.**
This is a fixed-surface, fixed-boundary-point theorem. It does not
settle JC(2), move an arbitrary boundary point to zero, or exclude
rational mates outside the stated stratum.

## 1. The complete statement and inherited boundary

Use the same actual DG surface as the
[quartic boundary theorem](planar_jc48_sep08_quartic_boundary.md):

    X=P1_x x P1_z,       S={z=x^2},       W=X minus S,
    s=z-x^2,            t=1/s,
    r=1/x,              t=-r^2-r^4 b_D,
    D={r=0} in W,       omega=dx wedge dt=r^2 dr wedge db_D.

Let `H=N/s^2` and `L=M/s` be global functions in `L_2,L_1`, so
`N,M` are sections of `O(4,2),O(2,1)`. Put `F=H^2+L`.

**Theorem.** Suppose the complete binary boundary octic is
`N|S=a x^6` with `a!=0`, the value `M(0,0)=m0` is nonzero,
and `F|D` is constant. Then there is no rational function
`G in C(x,t)` satisfying `J_(x,t)(F,G)=lambda!=0`.

Here `N|S=a x^6` means multiplicity six at the specified point
`(x,z)=(0,0)` and multiplicity two at the specified infinity point.
There are no other boundary zeros. In particular `deg_t F=4`.
Translating `x` need not extend to this compactification, and is
not used to infer a statement about sixfold zeros at other points.

This is exactly the family

    N=a x^6+s[beta0+beta1 x+b x^2+c1 x^3+2a x^4]
                 +s^2[c0+c1 x+a x^2],
    M=m0+m1 x+d x^2+n1 x^3+s(n0+n1 x),                 (1)

where `a*m0!=0` and the other parameters are arbitrary complex
constants. Indeed the global filtration gives `deg H|D<=2` and
`deg L|D<=1`; constancy of `H|D^2+L|D` forces both restrictions
to be constant. Expanding the complete section boxes then gives
(1), with

    H|D=c0-b,       L|D=n0-d,
    F|D=(c0-b)^2+n0-d.                                (2)

For clarity, the section space for `N` has six parameters, including
`beta0,beta1`; the balanced subfamily has those two equal to zero.
The exact source independently checks the nine linear constraints
on all fifteen `O(4,2)` monomials and the resulting six-dimensional
basis. The five displayed `M` parameters give the full constant-D
`O(2,1)` space. No restriction on a proposed rational mate is made.

The closest proved mechanism is the componentwise primitive budget
in [the common-root theorem](planar_jc48_sep08_quartic_common_root.md),
using [the actual pole-degree gate](planar_jc48_sep08_pole_degree.md).
The common-root theorem does not subsume this stratum: the two
boundary forms share the infinity point. The concrete stopping
object in [the shared-root first-jet note, Section 4](planar_jc48_sep08_shared_roots.md)
lies here. The proof below does not depend on that note; its
relevant infinity calculation is supplied directly below.

The canonical hostile is the actual rational mate for `F=h^4+h`,
`h=x^2+x^4t`, whose boundary octic instead has multiplicity eight.
The corrected near miss is to call the original constant-D stopping
object residue-free: it is in the nonzero logarithmic layer.
The least-used sidecar is the order-one zero of the actual relative
form at every generic infinity branch. It replaces the unavailable
transverse-D order-two zero.

The live concepts are complete section spaces, critical-value
orders, actual infinity branches, residues, and the pole-constrained
function space on an elliptic curve. The maps below preserve the
literal differential and rational exactness; a genus label or a
pole count alone discards the final coefficient obstruction.

## 2. Four actual infinity zeros give a degree-two gate

At infinity set `q_inf=1/z`, `sigma=q_inf-r^2`. Define the actual
section numerators

    Ni=r^4 q_inf^2 N(1/r,1/q_inf),
    Mi=-r^2 q_inf M(1/r,1/q_inf).

When `beta0=beta1=0`, direct transport gives

    Ni=a r^2-c1 r sigma+(c0-b)sigma^2-b r^2 sigma,
    Mi=-n1 r-d r^2-m1 r^3-m0 r^4
                      +sigma(n0-d-m1 r-m0 r^2).        (3)

The two extra terms in (1) add
`-beta0 r^2(r^2+sigma)sigma-beta1 r(r^2+sigma)sigma`
to `Ni`. They do not change its quadratic part.

For a generic original fibre value `c`, the complete local equation
and relative form are

    Ei=Ni^2+sigma^3 Mi-c sigma^4=0,
    eta=omega/dF=(unit sign) r^2 sigma^2 dr/(Ei)_sigma.

The sign has no bearing on its order; the factor `r^2` does.
The substitution `sigma=r Z` has tangent quartic

    T_c(Z)=[a-c1 Z+(c0-b)Z^2]^2
                           -n1 Z^3+(n0-d-c)Z^4.       (4)

For generic `c` it has four simple nonzero roots. Its constant
term is `a^2!=0`, and its leading coefficient is `F|D-c`.
If `T0=T_c+c Z^4`, a repeated nonzero root must solve
`Z T0'-4T0=0`, whose constant term is `-4a^2`. This is a nonzero
polynomial independent of `c`, so only finitely many fibre values
can fail simplicity. Each simple root lifts to a genuine branch
`sigma=r Z(r)`, with `r` a parameter. Since
`Ei(0,sigma)=(F|D-c)sigma^4`, these four branches exhaust the
local Weierstrass degree, with no omitted branch. On each,

    ord_r (Ei)_sigma=3,       ord_r eta=2+2-3=1.       (5)

Consider each compact normalized generic component separately.
The function `x` is nonconstant on every component: at nonzero
fixed `x`, the leading `t` coefficient of `F` is `a^2 x^12`,
and at zero, `F=(beta0 t+c0)^2+m0 t+n0` still depends on `t`:
it is quadratic when `beta0!=0`, and otherwise linear with
nonzero coefficient `m0`.
Thus every component has a pole of `x`. It cannot meet the finite
part of `D`, since `F|D` is constant and `c` is generic. Every
component therefore contains one of the four points in (5).

If a rational mate existed, divide it by its nonzero Jacobian
constant. On a generic normalized component its restriction would
satisfy `dG=eta`. A pole of this restriction inside `W` is
impossible, since `eta` is regular on the smooth generic fibre
there. The finite set of bad fibre values, including vertical
denominator values of a rational `G`, can be discarded. A pole of
order `p` of `eta` allows a primitive pole only of order `p-1`;
nonzero residues instead preclude a primitive. At any point (5),
`G-G(P)` has local degree two. Consequently the pole degree of
that same component's primitive is at least two. This argument
does not assume generic irreducibility or sum local degrees over
different values of `G`.

## 3. The finite cubic face and all lower critical layers

At zero, `M` is a unit and the original fibre is

    E=N^2+s^3 M-cs^4=0,
    eta=s^2 dx/E_s.                                   (6)

The previously proved local M-unit classification applies with `m=6` and
`j=ord_x N_s(0,x)`. If `beta0` or `beta1` is nonzero then
`j=0` or `1`; the `m>3j` regime is regular on every actual
normalized branch. If both vanish and `b=0`, then `j>=3`,
possibly infinite, and the `m<3j` regime is again regular.
These cases have no finite primitive pole budget and contradict
Section 2.

It remains to assume `beta0=beta1=0`, `b!=0`. Put `s=x^4 Z`.
Here the Weierstrass degree in `s` is exactly three. (For the
already excluded `beta0!=0` case it was two.)
The complete equation divided by `x^12` is the polynomial

    Q=A(x,Z)^2+Z^3[m0+m1 x+d x^2+n1 x^3]
                       +x^4 Z^4(n0-c)+n1 x^5 Z^4,
    A=a+bZ+c1 xZ+x^2(2aZ+c0Z^2)+c1 x^3 Z^2+a x^4 Z^2.

Its initial cubic is `P(Z)=(a+bZ)^2+m0 Z^3`, with discriminant
`a^3 m0(4b^3-27a m0)`. If the face is simple, all three branches
are regular. The only remaining face is

    q=-3a/b!=0,      a=-bq/3,       m0=-4b^2/(9q),
    P(Z)=-4b^2/(9q) (Z-q)^2(Z-q/4).                   (7)

The simple branch at `q/4` is regular. At the double branch,
`Q_ZZ(0,q)=-2b^2/3` is nonzero. The critical centre `Z=z(x,c)`
exists by the implicit function theorem. Let `R=Q(x,z(x,c),c)`.
Exactly as in the inherited proof,

    partial_c R=-x^4 z(x,c)^4,

so generically `lambda=ord_x R` is one of `1,2,3,4`. Parametric
Morse normalization makes `eta` a unit times `dx/sqrt(-R)`.
For `lambda=1` it is regular on its ramified branch. For two it
has two nonzero simple poles. For three it has one double pole on
the actual normalization `x=tau^2`; its possible primitive budget
is one. For four it has two double poles, budget two. These are
normalized branches, not a count of Puiseux determinations.
For `lambda<=3`, total possible primitive pole degree is less
than the degree two required in Section 2. This excludes those
layers without any assertion about their component genera.

## 4. The order-four layer forces an even family

Requiring the first three coefficients of `R` to vanish is exactly

    m1=-4b c1/(3q),
    d=4(2b^2 q-3b c0 q-3c1^2)/(9q),
    n1=4c1(b^2 q-3b c0 q-c1^2)/(9bq).                 (8)

These conditions are successive necessary equations, not merely
a sufficient tuning. The critical centre through order three is

    z=q+z1 x+z2 x^2+z3 x^3+O(x^4),
    z1=-c1 q/b,
    z2=q(2b^2 q+3b c0 q+3c1^2)/(3b^2),
    z3=-c1 q(b^2 q+9b c0 q+3c1^2)/(3b^3).

Writing `R=r4 x^4+r5 x^5+O(x^6)` gives

    r4=q^3[4b^2q-24bc0q-27cq+36c0^2q+12c1^2+27n0q]/27,
    r5=-4c1 q^3[5b^2q-27bc0q-27cq+36c0^2q+12c1^2+27n0q]/(27b).

In particular `partial_c r4=-q^4!=0`, so no unanalysed generic
higher-order layer remains. Set `h=Q_ZZ(x,z)/2=h0+h1x+...`.
Then `h0=-b^2/3`, `h1=-2bc1/3`. The two actual branches have
`Z=z+e x^2+f x^3+...`, where

    e^2=-r4/h0,       f/e=-h1/(2h0)+r5/(2r4).

It follows directly from `eta=Z^2 dx/Q_Z` that its double-pole
coefficient is nonzero and its residue divided by that coefficient
is

    2z1/q-h1/(2h0)-r5/(2r4).                          (9)

Twice (9), multiplied by `r4`, is

    -2c1 q^3[2b^2q-18bc0q-27cq+36c0^2q+12c1^2+27n0q]/(27b).

Its derivative in `c` is `2c1 q^4/b`. A rational mate must have
zero residues for generic `c`; hence necessarily `c1=0`. Equations
(8) now say

    m1=n1=0,        d=8b^2/9-4bc0/3.                  (10)

The complete fibre and polar denominator are then even in `x`.
Each of the two branches over the double slope is an even series
in the actual parameter `x`, so the two double poles are genuinely
residue-free. This is the surviving second-kind layer; a residue
test alone no longer suffices.

## 5. A literal elliptic model and its only possible primitive

Assume (7), (10), with arbitrary `c0,n0` and nonzero `b,q`. Define

    e0=c0-b,        C=4b^2q/9,
    B=2a e0-C,      K0=e0^2+n0-d.

The following are birational coordinates on the original rational
source field, not a change of the fibre value:

    u=1/x,       v=x+x^3t,
    w=v-u/q,     z=vw,
    v=z/w,       u=q(z-w^2)/w.                         (11)

Literal substitution in `F` gives

    F=-3a^2 z^2+Bz+K0+(4a^2 z+C)w^2.

Thus over `C(c)`, and after taking its algebraic closure, the
generic curve has equation

    A(z)w^2=P2(z),
    A(z)=4a^2 z+C,       P2(z)=3a^2 z^2-Bz+c-K0.       (12)

For generic `c`, the two roots of `P2` are simple, nonzero, and
different from the root of `A`. Indeed the quadratic discriminant
has nonzero `c` slope `-12a^2`; the latter two forbidden values
also have nonzero `c` slopes. The rational function `P2/A` has
odd valuations at these two roots, the root of `A`, and infinity.
It is therefore not a square over the algebraically closed
constant field. This proves geometric irreducibility, while the
degree-two map to the `z` line, with exactly these four simple
branch points, gives genus one by Riemann--Hurwitz. This genus
claim is for this exact remaining family, not for arbitrary
parameter degenerations of the earlier strata.

The Jacobians of (11) give the actual differential

    omega=-q^2(z-w^2)/w^2 dz wedge dw,
    eta=q^2(z-w^2)/(2A w^3) dz
        =q^2(Az-P2)/(2A^2 w^3) dz.                   (13)

At either zero of `P2`, `w` is a local parameter and `z-z_i`
has order two. Since `z_i!=0`, (13) has exactly a double pole
there. At the root of `A`, use a parameter with `z-z_A` of
order two and `w` of order minus one: (13) is regular. At
infinity, `z` has order minus two and `w` minus one, and the
form is again regular. There are no other poles. Consequently
a rational primitive could have only simple poles at the two
points where `w=0`.

For completeness, its whole allowable function space can be found
without importing a special elliptic integration criterion. Under
the involution `w -> -w`, those two pole points are fixed. The even
part of an allowable function belongs to `Cbar(c)(z)`; every pole
at a branch point has even order and so cannot be a permitted
simple pole. It is constant. Write the odd part as `w S(z)`.
Then `S` has at most simple poles at the two roots of `P2`, must
vanish at the root of `A` to cancel the pole of `w`, and must
vanish at infinity for the same reason. These conditions force

    S=K A/P2,       G=K/w+constant,                    (14)

for a scalar `K` in the algebraically closed generic constant
field. This is also the degree-two Riemann--Roch space on the
elliptic curve. Allowing the enlarged constants field only
strengthens the impossibility below.

Differentiate (14) using (12). Exactness of (13) would require

    q^2(Az-P2)=-K(P2' A-P2 A').                        (15)

The two polynomial numerators in this identity are

    Az-P2=a^2 z^2+(B+C)z-(c-K0),
    P2' A-P2 A'=12a^4 z^2+6a^2 C z-BC-4a^2(c-K0).

Their leading coefficients force `K=-q^2/(12a^2)`, independently
of `c`. After that substitution the constant-coefficient
difference has derivative `-2q^2/3` in `c`, which is nonzero.
Thus (15) is impossible over the generic constants field.
The obstruction is the incompatible fibre-value coefficients
`-1` and `-1/3`; it does not depend on a numerical period test.

This excludes the last layer and proves the theorem.
A rational mate on the original surface would restrict to the
forbidden generic primitive. Before this last layer, the argument
was componentwise and made no irreducibility assumption; here
geometric irreducibility has been proved directly.

## 6. Hostiles, scope, and the resulting connection

With `b=1`, `m0=-1`, `a=-4/27` and all other coefficients zero,
the earlier constant-D stopping object has `R1=0` but `R2!=0`.
It has nonzero logarithmic residues; it is not the second-kind
example. Tuning instead `d=8/9`, with `c0=c1=n0=0`, reaches
the even order-four layer. Its critical coefficient is
`r4=(256/6561)(4/27-c)`, so generic fibres retain the two double
poles. They are residue-free but fail (15).

The named `c1=1` order-four tuning retains all lower critical
equations and has nonzero generic residues. Separate named
controls retain orders one and three. These are exact boundary
checks of the analytic classification, not a parameter census.

The actual outside-family example

    h=x^2+x^4t,       F=h^4+h,
    G=1/[3x^3(4h^3+1)]

satisfies `J(F,G)=1`. Its octic is `x^8`, not `a x^6`.
It prevents extending the conclusion to every shared-root
quartic. The separate actual pair `F=t^4`, `G=-x/(4t^3)` has
Jacobian one and constant `F|D=0`, but its boundary pattern and
`M=0` violate the present hypotheses. Thus constant-D alone is
not an obstruction. Likewise `F=x^2`,
`G=t/(2x)` shows why an affine critical point alone cannot prove
the rational-mate exclusion established here.

The source of the final connection is the *actual* translated
relative form on the original fibre. The target is its
pole-constrained meromorphic function space. The birational map
(11) preserves the field, original fibre parameter, and exact
symplectic coefficient. Recording only the genus would lose
(13); recording only the two-pole budget would lose (15).
The underused operation is to retain both those data while
turning the last necessary condition into a two-coefficient
identity. The next unaddressed strata include other boundary
multiplicity partitions and sixfold points away from this fixed
point. No normalization moving such points is assumed here.

## 7. Exact reproduction and audit boundary

From the repository root:

```sh
python3 -B 04-computation/planar_jc48_sep08_const_d_quartic.py
python3 -B -O 04-computation/planar_jc48_sep08_const_d_quartic.py
```

The source imports no inherited mathematical implementation. It
checks complete section-space dimensions, both literal charts,
the complete infinity tangent polynomial, all successive
critical-value equations and the residue sign, actual normalized
pole budgets, literal birational/Jacobian transport, generic
branch-point separations, the primitive coefficient obstruction,
and the named positive and hostile controls. The all-parameter
analytic proof and the completeness of the primitive space are
the arguments above, not an inference from finite controls.

Both executions pass **97 always-active exact gates**, with
byte-identical output. The frozen source is **13,598 bytes**,
SHA256 `cfc977adb5dbfb0b1eb3f0118669551bd40dde758d2bab112734b9bf4f5a2194`.
The [frozen output](planar_jc48_sep08_const_d_quartic.out) and both
independent-mode producer replays are **462 bytes**, SHA256
`7b38b00e6ac5c6282122cfafb3d3fcfbc99c05cb3bb1fd52d299e0d85919e92f`.
The semantic record has SHA256
`7602857f67214f96f1000d116cdf16c3d80036bb0fb83ea44daf64655de095d6`.

The source and output are frozen. The [independent audit](planar_jc48_sep08_const_d_quartic_audit.md)
passes the complete analytic/source proof and both97-gate frozen replays.
Root also independently read and accepted the full analytic argument.
The audit is11231 bytes, SHA-256
`1d4cb6a367c047ddda99ffedad80b6397766d6983d4d76f3945e9e73ccd42ac8`.
