---
id: THM-3785
title: "Linear higher-pole Russell pseudo-plane maximal polynomial observable"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  The
  third-order rational Keller seed F=x(t+c), P=1/(3x^3) has complete
  polynomial target-field intersection equal to the smooth Russell--Miyanishi
  pseudo-plane R^2E=Z^3-c^3R.  The source affine plane is a surjective
  nonfinite etale cubic atlas.  The target has only constant units and
  Picard group Z/3.  It has no birational or affine-carrier Darboux pair,
  and no Darboux pair supported on at most two Euler weights per output; any
  survivor is nonlinear and at least 2-by-3 multigraded, has target field
  degree at least three by THM-3794, and would pull back to a planar Keller
  map of field degree divisible by three and at least nine.  Arbitrary Darboux pairs remain
  open, so no JC(2) counterexample is claimed.
source: jc_quartic_c3_construct / higher-pole escape lane, 2026-08-23
audit: >
  INDEPENDENTLY HOSTILE-AUDITED by jc_zero_debt_lift and root, 2026-08-23.
  Both audits reconstructed the atlas/intersection, units and Picard group,
  homogeneous modules, every two-by-two seam, affine-carrier/genus-one
  gates, birational exit, and degree floor.  The exact companion verifies the rational
  Keller identity, both polynomial coordinate systems, hypersurface and
  field reconstruction, full Poisson packet, smooth and etale unit ideals,
  arm/off-arm inverse formulas, triple boundary valuation, 125 graded
  monomial profiles, homogeneous bracket formula, two-by-two endpoint-power
  and leading-degree laws, affine critical residual, genus-one differential
  orders, and rational-mate formula.  Normal, optimized, frozen, and
  independently replayed outputs byte-match.  An audit correction removed an unsupported
  identification with Miyanishi's exact type parameter r=2; the integral
  closure and its nonstandard deck character are now stated explicitly.
depends_on:
  - THM-3794-constant-unit-surfaces-have-no-quadratic-etale-plane-map
related:
  - THM-3779-three-component-tower-maximal-danielewski-polynomial-observable
  - THM-3782-simple-pole-spectral-danielewski-completion-and-target-field-gate
  - THM-3783-quadratic-tower-etale-surface-maximal-polynomial-observable
  - THM-3784-rational-keller-tower-different-codifferent-trace-duality
script: 04-computation/jc2_linear_higher_pole_russell_pseudoplane_thm3785.py
output: 05-knowledge/results/jc2_linear_higher_pole_russell_pseudoplane_thm3785.out
script_sha256: 86d58649f16086e3bdf099a0a827b2e0c2e0a19faf15ac2f749b45cfedcbb192
output_sha256: 5fd5632563574ba05197bf4eac38669ffa44084a86eae5e230c01e4f61a9e9fd
semantic_sha256: a0ff1432a148cfa5fe716ca86694feaaf37a06c169fa4dad6842ab26d1178976
hash_basis: raw LF bytes
---

# THM-3785 -- a third-order pole completes to a Russell pseudo-plane

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  This
theorem turns the first linear member of a higher-pole rational Keller family
into an exact smooth target rather than another failed target shear.  It does
not construct a planar Jacobian counterexample.  It identifies the complete
surface on which the first such counterexample in this family would have to
live and closes its birational, affine, and homogeneous entrances.

Let `k` be an algebraically closed field of characteristic zero, fix
`c in k*`, and put

```text
t=x^2(1+xy),                 q=t+c,
Q=q^3/3,                     F=xq,
P=Q/F^3=1/(3x^3).                                  (1)
```

Because `Q'=q^2`, direct differentiation gives

```text
J(F,t)=x^3q,                 J(F,P)=1.               (2)
```

Thus `(F,P)` is a rational Keller pair with a third-order axis pole.  The
point of the theorem is to compute every element of `k(F,P)` that is a
source polynomial.

## 1. The hidden exponent-(2,3) completion

Define three source polynomials

```text
R=x^3,
Z=F-R=x(c+x^3y),
E=3c^2y+3cx^3y^2+x^6y^3
 =[(c+x^3y)^3-c^3]/x^3.                              (3)
```

They satisfy

```text
R^2E=Z^3-c^3R.                                       (4)
```

Conversely,

```text
F=R+Z,                       P=1/(3R),                (5)
```

so `k(F,P)=k(R,Z)`.  Let

```text
B=k[r,z,e]/(r^2e-z^3+c^3r),              Y=Spec B.   (6)
```

The map from `(6)` to `k[x,y]` defined by `(3)` is injective: the
hypersurface ring is a two-dimensional domain and its image has fraction
field `k(R,Z)` of transcendence degree two.  Hence `(4)` is the complete
relation among the three observables.

The less symmetric coordinates produced by direct spectral interpolation
are

```text
Y_0=(F^3-3RF^2-c^3R)/R^2
   =3c^2y+3cx^3y^2-3cx+x^6y^3-3x^4y-2x^3.            (7)
```

They are triangularly equivalent to `(3)`:

```text
Z=F-R,                       E=Y_0+3F-R,              (8)
R^2Y_0=F^3-3RF^2-c^3R.                               (9)
```

Thus the apparently singular cubic spectral completion becomes the smooth
Russell--Miyanishi pseudo-plane `(6)` after retaining exactly one additional
axis jet.

## 2. A smooth surjective etale cubic atlas

The gradient of the equation in `(6)` is

```text
(2re+c^3, -3z^2, r^2).                                (10)
```

It has no common zero because `r=z=0` leaves the nonzero entry `c^3`.
Therefore `Y` is smooth and normal.  The source Jacobian minors are

```text
J(R,Z)=3R^2,
J(R,E)=9Z^2,
J(Z,E)=3c^3+6RE.                                     (11)
```

They are the signed gradient packet `(10)`, scaled by `3`, and generate the
unit ideal on `Y`.  Hence

```text
phi:A2_(x,y) -> Y,             (x,y) |-> (R,Z,E),     (12)
```

has full differential rank everywhere.

It is surjective on geometric points.  Given `(r,z,e) in Y` with `r!=0`,
choose any `xi` with `xi^3=r` and put

```text
x=xi,                         y=(z/xi-c)/r.           (13)
```

Then `Z=z`, while `(4)` gives the required value of `E`.  If `r=0`, equation
`(4)` forces `z=0`, and

```text
x=0,                          y=e/(3c^2)              (14)
```

is a lift.  Thus every off-arm point has three lifts and the arm

```text
L=V(r,z) ~= A1_e                                      (15)
```

has one lift.  In particular `phi` is quasi-finite.  Since both source and
target are smooth surfaces, quasi-finiteness supplies equality of the local
dimensions, and miracle flatness gives flatness.  Full differential rank then
makes `phi` etale.  Its generic degree is exactly three:

```text
k(x,y)=k(R,Z)(x),                 x^3=R,              (16)
```

and `R` is not a cube in `k(R,Z)` by its `R`-adic valuation.  The map is not
finite, since an etale finite map of degree three could not have the single
reduced boundary fibre `(14)`.  Two Kummer sheets escape at infinity along
`L`.

## 3. The exact maximal polynomial intersection

Inside `k(x,y)`, set

```text
A=k[x,y] intersection k(F,P).                         (17)
```

Equations `(3)--(6)` give `B subset A`.  Conversely, take `H in A` and a
height-one prime `q` of the normal ring `B`.  Surjectivity and etaleness of
`phi` supply a height-one prime `p` of `k[x,y]` over `q`.  The corresponding
DVR extension has ramification index one, so for `H in Frac B`,

```text
ord_p(H)=ord_q(H).                                    (18)
```

The left side is nonnegative because `H` is a source polynomial.  Hence
`ord_q(H)>=0` for every height-one prime of `B`.  The Krull intersection of
the height-one DVRs of the normal domain `B` proves

```text
k[x,y] intersection k(F,P)
 =k[R,Z,E]/(R^2E-Z^3+c^3R).                           (19)
```

This is the maximal polynomial target-field observable, not merely a finite
list of convenient cancellations.

## 4. Units, Picard torsion, and the cubic boundary

Localizing `(6)` at `r` eliminates `e`:

```text
B_r=k[r,r^(-1),z].                                    (20)
```

The unique height-one prime over `r=0` is `L=(r,z)`.  At its generic point,

```text
r(c^3+re)=z^3,                                        (21)
```

and the parenthesized factor is a unit.  Thus `z` is a uniformizer and

```text
ord_L(r)=3,                         div(r)=3L.         (22)
```

Nagata localization says that `Cl(B)` is generated by `[L]`, while `(22)`
gives `3[L]=0`.  The class is nonzero: if `L=div(h)`, then `h` is a unit in
`B_r`, hence `h=alpha r^n`; its divisor is `3nL`, never `L`.  Smoothness
identifies the class and Picard groups.  The same localization calculation
shows that a global unit can only have zero `r`-exponent.  Therefore

```text
B^*=k^*,                         Pic(Y)=Z/3.           (23)
```

The triple class is the geometric shadow of the Kummer equation `x^3=R`.
It survives even though `(12)` is everywhere etale because the other two
sheets are nonproper at the arm.

### 4.1 Exact modern pseudo-plane type and the deck-character sidecar

The projection `r:Y -> A1` has general fibre `A1_z` and the unique multiple
fibre `3L`: scheme-theoretically, setting `r=0` gives `z^3=0`.  Normalize

```text
u=r,                 v=e/c^3,                 w=z/c.                 (24a)
```

Then `(6)` is

```text
u(1+uv)=w^3.                                                   (24b)
```

This is exactly Dubouloz--Palka's pseudo-plane `S(3,3,1)` in
arXiv:1701.01425, equations (5.1)--(5.2).  Its finite universal cover and
deck action are

```text
T=k[X,V,W]/(X^3V-W^3+1),
mu_3:(X,V,W) |-> (zeta X,V,zeta^(-1)W),
(u,v,w)=(X^3,V,XW).                                      (24c)
```

In the source coordinates, `X=x`, `V=E/c^3`, and
`W=(c+x^3y)/c`.  Thus `(24c)` is also the integral closure of `B` in the
cubic source field.  The action is free: a nontrivial fixed point would have
`X=W=0`, contrary to its equation.  A weight-zero monomial reduces using
`X^3=u`, `XW=w`, and `W^3=1+uv`, proving directly that `T^(mu_3)=B`.

Miyanishi's arXiv:1504.07179, Lemma 2.5.2, gives the abstract affine-plane
Galois pseudo-cover mechanism for rational homology planes with one multiple
`A1` fibre.  Equations `(3),(12)--(16)` give the explicit one-branch affine
atlas associated to `(24c)`.

Theorem 2.5.10 of the same source excludes etale morphisms from a
pseudo-plane of its exact older type `(d,n,r)` to `A2` when `r!=2`, using

```text
Pic(Y)=Z/d,                     K_Y=(r-2)[L].          (24d)
```

Our symplectic form below trivializes `K_Y`, so that canonical-class
obstruction is silent here.  There is no contradiction with `(24d)`: after
rewriting Miyanishi's standard cover in the form `X^rV=W^d-1`, its deck
character is `a=d-1`, whereas `(24c)` has `a=1`.  For `d=3` these are
different pseudo-plane families.  The older `(d,n,r)` label suppresses this
character and must not be applied to `S(3,3,1)`.  The citation is context,
not a dependency; every intrinsic assertion here was proved directly.

## 5. Poisson packet and the exact symplectic form

The source bracket induces exactly `(11)`, equivalently

```text
{r,z}=3r^2,
{r,e}=9z^2,
{z,e}=3c^3+6re.                                      (24)
```

On `r!=0`, its inverse form is

```text
omega=dr wedge dz/(3r^2),              phi^*omega=dx wedge dy.         (25)
```

The bracket packet is a nowhere-vanishing global bivector on the smooth
surface, so its inverse extends `(25)` to a global regular two-form; the
displayed denominator does not signal a pole on the arm.

Unlike the exponent-one multiarm completion of THM-3779, this form is
exact.  Give `(r,z,e)` weights

```text
wt(r)=3,                       wt(z)=1,       wt(e)=-3.                 (26)
```

The relation has weight three and the Poisson bracket has degree `+2`.
For

```text
Eul=3r partial_r+z partial_z-3e partial_e,             (27)
```

one has `Lie_Eul(omega)=-2omega`.  Cartan's formula gives the global regular
primitive

```text
alpha=-(1/2)i_Eul(omega),                 d alpha=omega.               (28)
```

On the source this is

```text
phi^*alpha=-(3y dx+x dy)/2,              d(phi^*alpha)=dx wedge dy.   (29)
```

Exactness removes the quickest de Rham obstruction.  It is a construction
permission, not a construction.  The direct coordinates `(R,F,Y_0)` have

```text
{R,F}=3R^2,
{R,Y_0}=9F^2-18RF,
{F,Y_0}=3c^3+9F^2+6RY_0.                              (30)
```

The constant term in the last bracket is real, but Sections 6 and 7 show
why its finite polynomial correction cannot begin with an affine or a
single-weight answer.

## 6. Exact grading and homogeneous nonentry

Put

```text
w=c+x^3y,                      h=w^3-c^3.             (31)
```

Inside `k[x,x^(-1),w]`, the generators are

```text
R=x^3,                         Z=xw,          E=x^(-3)h.               (32)
```

For `a in Z`, let `rho(a) in {0,1,2}` be its residue modulo three and put

```text
m(a)=max(0,ceil(-a/3)).                                 (33)
```

The exact homogeneous module of weight `a` is

```text
B_a=x^a w^rho(a) h^m(a) k[w^3].                       (34)
```

Necessity follows by writing a generator monomial `R^iZ^jE^k`: its
`x`-exponent is `3i+j-3k`, so `j` has the required residue and a negative
exponent forces at least the power in `(33)`.  Conversely, the displayed
base profile is respectively a monomial in `R,Z` when `a>=0` and in `Z,E`
when `a<0`; multiplication by `w^3=c^3+RE` stays in `B`.

For homogeneous functions

```text
A=x^a f(w),                      C=x^b g(w),           (35)
```

direct differentiation gives

```text
{A,C}=x^(a+b+2)[a f g'-b f'g].                        (36)
```

A scalar requires `a+b=-2`.  If both weights are negative, `(34)` makes
both profiles divisible by `h`, so the bracket vanishes at every root of
`h`.  If `a>=0>b`, a nonzero value at `h=0` requires the negative profile to
have exactly one `h` factor and also requires `a!=0`.  With `b=-2-a`, this
leaves only

```text
(a,b)=(1,-3),                       or its swap.       (37)
```

In that case write

```text
f=wF(w^3),                         g=hG(w^3).          (38)
```

If `p=deg F` and `q=deg G`, the coefficient in `(36)` has degree

```text
3(p+q+1)>=3                                      in w, (39)
```

and leading coefficient

```text
3(3p+q+2) lc(F)lc(G)!=0.                             (40)
```

It cannot be a scalar.  Thus there is no homogeneous Darboux pair.  If one
output were homogeneous and the other had several weights, the unique
weight-zero bracket bucket would again be one forbidden homogeneous scalar
bracket.  Hence both outputs of any live pair must use at least two weights.

### 6.1 Complete two-by-two support nonentry

The same exact modules close every pair having at most two nonzero weights
per output.  Suppose the active supports are

```text
supp(A)={a,a+d},                 supp(C)={b,b+e},      d,e>0.          (41)
```

If `d!=e`, all four bracket buckets are distinct, so a scalar would be one
forbidden homogeneous bracket.  Hence `d=e`.  The scalar cannot lie at an
endpoint bucket, so it lies in the double middle bucket and

```text
a+b+d=-2,                                                (42)
{A_a,C_b}=0,                    {A_(a+d),C_(b+d)}=0.    (43)
```

After subtracting irrelevant constants, if a pure weight-zero component
disappears then the support shrinks to an already-excluded homogeneous-output
or smaller-support case.  Otherwise the low weights are both negative:
opposite-sign homogeneous functions cannot commute, while a weight-zero
profile commuting with a nonzero weight is constant.  Gaps `d=1,2` die
directly at a root of `h`; when `d=2`, the only residual endpoint has both
top weights zero and both middle scalar summands still vanish there.

Let `d>=3` and write

```text
A_0=-a>0,                 B_0=-b>0,       A_0+B_0=d+2,
C_0=a+d=B_0-2>0,          D_0=b+d=A_0-2>0.             (44)
```

Endpoint commutation and unique factorization give, after absorbing nonzero
scalars,

```text
f=p^(A_0/delta),          h=lambda p^(B_0/delta),
F=q^(C_0/epsilon),        H=mu q^(D_0/epsilon),         (45)

delta=gcd(A_0,B_0),       epsilon=gcd(C_0,D_0).
```

Here `f,h` are the low profiles and `F,H` the high profiles.  If
`A_0!=B_0`, the degrees of the two middle scalar summands differ by

```text
(B_0-A_0)[deg(p)/delta+deg(q)/epsilon],                (46)
```

and the higher summand has nonzero leading coefficient, so cancellation to
a scalar is impossible.  If `A_0=B_0`, the scalar bucket is either zero or a
nonzero scalar multiple of

```text
A_0 p q'+(A_0-2)p' q.                                 (47)
```

This has positive degree because the negative profile `p` is divisible by
`h=w^3-c^3`.  Thus no `2 by 2` pair exists.  Combined with the one-output
homogeneous observation, any live support has at least two weights in each
output and at least three in one output.

## 7. Complete affine-carrier classification

Let

```text
A=kappa+a r+b z+d e                                  (48)
```

be nonconstant.  The Poisson packet `(24)` gives the critical equations

```text
b r^2+3d z^2=0,
a r^2-d(2re+c^3)=0,
3a z^2+b(2re+c^3)=0.                                 (49)
```

If `d=0`, these have no solution exactly when `b!=0`.  If `d!=0`, scale
`d=1`.  Any critical point has `r!=0`; put `z=lambda r`.  The first two
equations give

```text
b+3lambda^2=0,                 e=(ar^2-c^3)/(2r),     (50)
```

and the surface equation becomes

```text
(a-2lambda^3)r^2+c^3=0.                              (51)
```

Unless `a=b=0`, one of the two choices of `lambda` from `(50)` has nonzero
coefficient in `(51)`, and algebraic closedness supplies a critical point.
When `a=b=0`, equations `(49)` together with `(6)` are incompatible, so `e`
is smooth.
Therefore the smooth affine-linear carriers, up to scale and translation,
are exactly

```text
z+lambda r,                         and e.             (52)
```

The first family has a rational mate:

```text
{z+lambda r,1/(3r)}=1.                              (53)
```

Moreover `k(r,z)=k(z+lambda r,r)`, so every rational mate is
`1/(3r)+H(z+lambda r)` with `H in k(z+lambda r)`.  None is polynomial in
`B`.  On the source,

```text
z+lambda r=x(c+lambda x^2+x^3y),                     (54)
```

whose zero fibre has the disjoint components `x=0` and
`c+lambda x^2+x^3y=0`.  Cancelling the third-order pole of `1/(3x^3)` on
the first component forces a third-order pole of `H` at the common target
value zero; that same pole is then introduced on the second component,
where `1/(3x^3)` is regular.  This is impossible.

The carrier `e` has no rational mate at all.  Over the generic value
`eta`, its projective fibre is the cyclic cubic

```text
z^3=eta r^2+c^3r.                                    (55)
```

It has three fully ramified points over the two finite roots and infinity,
so Riemann--Hurwitz gives genus one.  On this curve

```text
-dr/(9z^2)                                           (56)
```

is a nonzero holomorphic differential: it has order zero at both finite
branch points and at infinity.  If `{e,H}=1`, then
`dH=-dr/(9z^2)` on the generic fibre, impossible because a nonconstant
rational function has a pole and its differential cannot be holomorphic.
Scaling does not change the obstruction.  Thus no affine-linear carrier has
a polynomial Darboux mate.

## 8. The exact remaining Darboux floor

Now take `k=C` and suppose `A,C in B` satisfy `{A,C}=1`.  Put

```text
d=[Frac B:C(A,C)].                                    (57)
```

This is a finite positive integer.  It cannot equal one.  Indeed, a
birational etale morphism `Y -> A2_(A,C)` is an open immersion by Zariski
Main.  Since `B^*=C^*`, its complement cannot contain a divisor: the inverse
of an equation for that divisor would be a nonconstant unit.  The complement
is therefore finite.  Normal Hartogs extension gives the open set the same
global functions as `A2`, and affineness would force `Y=A2`, contradicting
`Pic(Y)=Z/3`.

THM-3794 also excludes `d=2`, because `B^*=C^*`.  Consequently

```text
d>=3.                                                 (58)
```

The source atlas has degree three by `(16)`, so the pullback polynomial
Keller pair would have function-field degree

```text
[C(x,y):C(A,C)]=3d>=9.                                (59)
```

Thus the exact first survivor is not another target shear or sparse
one-weight correction.  It must be a nonlinear, multigraded, nonbirational
Darboux pair on the smooth pseudo-plane `(6)`.  Such a pair would immediately
give a planar Jacobian counterexample of field degree divisible by three.
The theorem neither constructs nor excludes that final object.

The deterministic companion named in the metadata verifies every displayed
source identity and bracket, smoothness and the unit-minor ideal over
`Q(c)`, the off-arm and arm lifts, the triple valuation identity, 125
graded generator profiles, the exceptional homogeneous leading term, the
two-by-two power and degree identities, the affine critical residual, and
the genus-one differential orders.  Normal and optimized executions
byte-match the frozen transcript.  **QED.**
