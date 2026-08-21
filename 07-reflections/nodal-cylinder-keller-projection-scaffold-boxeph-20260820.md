# Nodal-cylinder projection scaffold: a counterexample cell past the normal-ideal gate

**Status: PROVED structural/width/period/first-conductor identities +
FINITE-EXACT two-scale controls + CITED degree and Newton transfer; OPEN
pointwise correction problem.  No Keller counterexample is constructed.**

The cusp-square packet of THM-3556 fails before integrability: its descended
normal multiplier does not belong to the image singular-minor ideal.  The
one small useful hostile is therefore an immersed nonnormal surface for which
that multiplier *does* belong to the ideal.  The nodal-cubic cylinder supplies
one, and it leaves polynomial-potential realizability as the exposed missing
sidecar.

## 1. A differential immersion with an automatic collision

On source coordinates `(u,t)`, put

```text
X=t^2-1,             Y=t(t^2-1),
Z(u,t)=(u,X,Y).                                           (1)
```

Its image is the hypersurface

```text
H=Y^2-X^2(X+1)=0.                                        (2)
```

The two tangent minors involving `u` are

```text
M_uX=2t,             M_uY=3t^2-1=3X+2,                  (3)
```

and they never vanish together, so `Z` has differential rank two everywhere.
Nevertheless

```text
Z(u,1)=Z(u,-1)=(u,0,0)                                  (4)
```

for every `u`.  Every polynomial projection of the image retains this
collision.  A unit-Jacobian projection would therefore be a genuine
noninjective planar Keller map.

### 1.1 The every-line gate forces a second collision

Gwozdziewicz's **CITED**
[injectivity-on-one-line theorem](https://arxiv.org/abs/alg-geom/9305008)
says that a planar Keller map over an algebraically closed characteristic-zero
field is an automorphism if its restriction to even one affine source line is
injective.  Hence a counterexample projection through `(1)` cannot map the
line `t=1` injectively.  If

```text
gamma(u)=F(u,1)=F(u,-1),                              (4a)
```

then `gamma` must itself be noninjective.  It is nevertheless immersed,
because a zero derivative would contradict the unit Jacobian.  Choosing
`u_0!=u_1` with `gamma(u_0)=gamma(u_1)` produces the four distinct preimages
`(u_i,+/-1)` of one target point.  Thus an admissible boundary is an immersed
but nonembedded polynomial curve.  More globally, a genuine counterexample
must be noninjective on *every* affine source line.

This already fixes the smallest boundary model.  If a polynomial curve
`gamma(z)=c_2z^2+c_1z+c_0` of degree at most two has a collision `r!=s`, then

```text
(r+s)c_2+c_1=0,
gamma'((r+s)/2)=0.                                    (4b)
```

Thus an immersed noninjective boundary has degree at least three.  The bound
is sharp: up to affine reparameterization and triangular target changes, the
minimal degree pair `(2,3)` is the nodal normalization

```text
gamma_0(u)=(u^2-1,u^3-u),
image: q^2=p^2(p+1).                                  (4c)
```

Indeed a quadratic first coordinate can be centered to `u^2`; subtracting a
multiple of it from the cubic second coordinate leaves `u^3+cu`, and
immersion forces `c!=0`, which scales to `(4c)`.  Consequently every
nodal-cylinder Keller projection has a fibre of size at least four, and its
least possible conductor restriction is exactly an immersed nodal cubic.

## 2. The normal-multiplier gate passes

The image normal and source tangent normal obey

```text
(H_u,H_X,H_Y)(Z)
  =X (0,-M_uY,M_uX).                                    (5)
```

Thus the Hodge multiplier is the descended function `X`.  Unlike the
cusp-square packet, it lies in the image Jacobian ideal.  Indeed, modulo `H`,

```text
X=-((3X+1)/2) H_X-(9/4)Y H_Y.                          (6)
```

Dividing the pulled-back identity by the common multiplier gives the explicit
unit tangent-minor certificate

```text
1=((3X+1)/2)M_uY-(9/4)Y M_uX.                          (7)
```

So arbitrary descended coefficients are already sufficient here.  This is a
strictly stronger positive control than a source-ring Bezout certificate.

The coefficient two-form can even be made closed without changing `(7)`:

```text
omega=-(9/4)Y du wedge dX
      +(3X+1)/2 du wedge dY
      +(15/4)u dX wedge dY.                            (8)
```

Direct differentiation gives `d omega=0`, and the last coefficient multiplies
the identically zero tangent minor `M_XY`.

Closedness is not polynomial-potential realizability.  The vector field dual
to `(8)`, after writing `V=3X+1` and multiplying by four, is

```text
D=15u partial_u-6V partial_V-9Y partial_Y.              (9)
```

If `omega=dA wedge dB`, then `A` and `B` are polynomial invariants of `D`.
Every invariant monomial `u^a V^b Y^c` satisfies

```text
5a=2b+3c.                                               (10)
```

At `u=0`, `(10)` forces `b=c=0`; hence every polynomial invariant restricts
to a constant there.  Two invariant gradients are therefore parallel on
`u=0`, so their cross product vanishes, while `(9)` is generically nonzero on
that plane.  The particular closed certificate `(8)` cannot equal
`dA wedge dB`.  Other certificates remain open; this negative control shows
exactly why “closed and decomposable” must not be conflated with two global
polynomial potentials.

## 3. Classical degree gates move the first cell to target cap 38

Give the target generators `(u,X,Y)` their pulled-back source degrees
`(1,2,3)`.  A target polynomial of total degree at most `E` pulls back to
source degree at most `3E`.  The cited Guccione--Guccione--Horruitiner--Valqui
bound therefore starts with `E>=36`; below source height `125`, target
reduction must give the degree pair `(72,108)`.

The raw bound is not sharp for this scaffold.  At source degree `108`, the
possible top binary monomials have `u`-exponent sets

```text
E=35: empty,
E=36: {0},
E=37: {0,1},
E=38: {0,1,2,3}.                                      (11)
```

At cap `36` the degree-`108` top is a scalar `t^108`, excluded by the
one-point-at-infinity gate.  At cap `37` it is
`t^107(a t+b u)`: for `b!=0` it is not a proper power, while for `b=0` it
again has one root.  The current proper-power and infinity-root gates thus
give the sharper scaffold-specific floor

```text
E>=38.                                                 (12)
```

At `E=38` the first admissible two-root power is

```text
K=t^35(t+lambda u),       Q_108=K^3,                  (13)
```

with `lambda!=0`.  The common-leading-base law then forces the reduced
degree-`72` top to be `P_72=K^2`.

The common base has a target-ring lift with the required top form:

```text
C=Y^12+lambda u X Y^11.                               (14)
```

Its pullback is `t^11(t^2-1)^12(t+lambda u)`, not literally `K`; its
ordinary top form is exactly `K`, which is the property used below.

Although the unreduced ambient cube `C^3` has cap `39`, the image relation
gives `X^3=Y^2-X^2`.  Therefore its exact representatives

```text
P_0=C^2
   =Y^24+2 lambda u X Y^23+lambda^2 u^2 X^2 Y^22,

Q_0=Y^36+3 lambda u X Y^35+3 lambda^2 u^2 X^2 Y^34
    +lambda^3 u^3(Y^35-X^2Y^33)
   =C^3                       modulo H                 (15)
```

both have target cap at most `38`.  Their pullbacks are literal square/cube
powers, so the highest-`u` resonance

```text
2a_2 b_3'-3a_2' b_3=0                                 (16)
```

holds identically.  Omitting the correction `-lambda^3u^3X^2Y^33` gives the
same ordinary top but a false nonzero width defect; reduction by the image
relation is load-bearing.

To meet the every-line gate on the two conductor lines, use the immersed
nodal boundary

```text
gamma(u)=(u^2-1,u(u^2-1)),
gamma(1)=gamma(-1)=(0,0).                             (17)
```

Its derivative `(2u,3u^2-1)` never vanishes.  The low transverse rows

```text
A_0=u^2-1+Y/2+P_0,
B_0=u(u^2-1)+(3/4)uY+Q_0                              (18)
```

do not disturb `(13)`, and exact differentiation gives

```text
Jac(A_0,B_0)|_(t=1)=Jac(A_0,B_0)|_(t=-1)=1.           (19)
```

All four points `(u,t)=(+/-1,+/-1)` map to `(0,0)`.  Thus `(18)` passes the
global degree architecture, common two-root leading base, automatic packet
collision, selected-line noninjectivity, four-point fibre, differential-rank,
normal-ideal, exact highest-row resonance, and collision first-jet gates.
Its full Jacobian is not constant, so it is a scaffold rather than a
counterexample.

There is also a decisive independent hostile.  On the horizontal line
`t=c`, write

```text
C=beta(c)+alpha(c)u,
alpha=lambda t^11(t^2-1)^12,   beta=[t(t^2-1)]^12.    (20)
```

Because `alpha` is nonconstant, some `c` has `alpha(c)^2=-1`.  At such a
point `beta(c)!=0`, and the first component in `(18)` becomes

```text
A_0|_(t=c)=2 alpha(c) beta(c) u+constant.              (21)
```

It is nonconstant affine-linear, so the restriction of the pair to this
source line is injective.  Gwozdziewicz's theorem therefore rules out the
displayed scaffold even before its nonconstant Jacobian is inspected.  More
generally, if a trial first component has horizontal `u`-degree at most two,
`a_2(t)u^2+a_1(t)u+a_0(t)`, the every-line gate forces
`rad(a_2)` to divide `a_1`: otherwise a root of `a_2` not shared by `a_1`
produces the same affine-linear hostile.  This necessary divisibility test is
cheap and survives arbitrary lower-row corrections.

### 3.1 Newton similarity forces the second scale into both components

The preceding hostile is the local shadow of a sharper global width law.
Lang's **CITED**
[*Newton polygons of Jacobian pairs*](https://doi.org/10.1016/0022-4049(91)90128-O)
proves that the two Newton polygons of a Keller pair are similar.  In the
reduced degree cell `(72,108)`, with the origin included in each polygon,

```text
Newt(B)=(3/2) Newt(A).                                  (22)
```

If `(m,n)=(deg_u A,deg_u B)`, taking the maximal `u` coordinate in `(22)`
gives

```text
n=(3/2)m,       hence (m,n)=(2r,3r).                   (23)
```

The minimal width `r=1` is impossible for this nodal source.  Write its
leading `u` coefficients as `a_2(t),b_3(t)`.  The coefficient of `u^4` in
the Jacobian gives

```text
2a_2b_3'-3a_2'b_3=0,
a_2=alpha h^2,       b_3=beta h^3.                    (24)
```

The second statement follows by valuations in the UFD `k[t]`.  It also
descends through the nonnormal ring: since `a_2,b_3` belong to

```text
S=k[t^2-1,t(t^2-1)]={f in k[t]: f(1)=f(-1)},           (25)
```

the two coprime powers in `(24)` force `h(1)=h(-1)`, so `h in S`.  The top
forms `(13)` force `deg_t a_2=70`, `deg_t b_3=105`, hence `deg h=35`.
At a root `c` of this nonconstant `h`, the restriction to `t=c` has degree at
most `(1,2)`.  It is immersed because the full Jacobian is a unit, so `(4b)`
makes it injective.  Gwozdziewicz then makes the whole map an automorphism,
contradicting `(4)`.  Therefore

```text
(deg_u A,deg_u B)>=(4,6) in the scaled width ladder.    (26)
```

This is an all-degree obstruction, not a bounded coefficient search.  At the
first remaining width `(4,6)`, the same calculation gives

```text
a_4=alpha r^2,       b_6=beta r^3,       r in S.       (27)
```

Target cap `38` leaves degree at most `32` for the coefficient of `u^6`, so
`deg_t b_6<=96` and

```text
deg_t r<=32.                                             (28)
```

The infinity face uses a coefficient base of `t`-degree `35`, whereas its
first possible `u`-width correction uses a different base of degree at most
`32`.  Thus cap `38` really does require a second Newton scale in *both*
components; quotient reduction repaired `(16)` but did not remove this
every-line width debt.

### 3.2 An exact cap-38 two-scale scaffold

The bound `(28)` is attained by a particularly small target-ring object.  Put

```text
V=lambda uY+mu u^2,
D=Y^12+XY^10 V.                                        (29)
```

Its two source scales have degrees `36` and `34`; the coefficient bases have
`t`-degrees `35` and `32`.  The square has target cap `26`.  The raw cube has
cap `39`, but the same nodal relation gives the exact cap-`38` representative

```text
D^3 = Y^36+3XY^34V+3X^2Y^32V^2
      +(Y^32-X^2Y^30)V^3             in R.             (30)
```

Consequently `(D^2,D^3)` has source degrees `(72,108)`, widths `(4,6)`,
top forms `(K^2,K^3)`, and identically zero mutual bracket.  Its leading
width factor is

```text
r=XY^10=t^10(t^2-1)^11,                               (31)
```

whose distinct roots are `t=-1,0,1`.  Let `eta=3X+1` and define the
**FINITE-EXACT** scaffold

```text
A_*=D^2+u^2-1+(eta/2)Y,
B_*=D^3+u^3-u+(3eta/4)uY.                             (32)
```

On all three degree-drop lines, `D=D_t=0` and

```text
(A_*,B_*)=(u^2-1,u^3-u),
Jac(A_*,B_*)=(3X+2)(3X+1)/2=1.                       (33)
```

Thus every root of the second-scale leading coefficient lands on the minimal
noninjective immersed boundary rather than an injective low-degree curve.
The four conductor points still form one fibre.  This simultaneously passes
the cap, Newton-width, relation-reduction, leading-bracket, root-line,
four-point, and first-jet gates.  Its full Jacobian is nonconstant.

### 3.3 The conductor period and the first-layer no-go

There is an exact global condition between pointwise Jacobian matching and
the endpoint jets.  For any `A,B in R`, set

```text
Phi_(A,B)(u)=integral_(-1)^1 A(u,t) B_t(u,t) dt.       (34)
```

Because the values of `A,B,B_u` agree at `t=+/-1`, integration of

```text
Jac(A,B)=partial_u(AB_t)-partial_t(AB_u)               (35)
```

shows that a normalized Keller pair must obey

```text
Phi_(A,B)=2u+constant.                                 (36)
```

This is the de Rham period, or symplectic action, around the loop created by
normalizing the node.  It is the precise analogue of retaining a cycle
current rather than only Kirchhoff data at each endpoint.

The period is completely explicit.  Write

```text
A=a(u,X)+Yc(u,X),          B=b(u,X)+Yd(u,X),
a=sum a_iX^i,              c=sum c_iX^i,               (37)
```

and similarly for `b,d`.  If

```text
beta_n=(-1)^n 2^(2n+1)(n!)^2/(2n+1)!,
K_ij=-2i beta_(i+j)/(2(i+j)+3),                        (38)
```

then parity and one integration by parts give the **PROVED** formula

```text
Phi_(A,B)=sum_(i,j) K_ij(a_i d_j-b_i c_j).             (39)
```

The gate has independent power but is not sufficient.  The hostile

```text
A=u^2-1+Y/2,
B=u^3-u+(3/4)uY+(105/16)uX^2                          (40)
```

has the nodal four-point fibre, unit endpoint jets, and exactly `Phi=2u`,
yet its Jacobian is nonconstant.

The normalization also turns the full correction problem into two polynomial
PDEs.  Put

```text
L(f)=(3X+2)f+2X(X+1)f_X.                              (41)
```

Then `Jac=E+tO`, where

```text
E=a_uL(d)-L(c)b_u+2X(X+1)(c_u b_X-a_Xd_u),
O=2(a_u b_X-a_Xb_u)+X(c_uL(d)-L(c)d_u).               (42)
```

The Keller equation is exactly `E=1,O=0`.  This already closes the entire
first conductor layer.  With boundary `(17)`, suppose there are no
`(X,Y)^2` terms:

```text
A=a_0+Xa_1+Yc_0,          B=b_0+Xb_1+Yd_0.            (43)
```

The endpoint equations parameterize every first jet by `h,k in k[u]`:

```text
a_1=2uh,                    b_1=(3u^2-1)h,
c_0=1/2+2uk,                d_0=3u/4+(3u^2-1)k.       (44)
```

Equation `(39)` gives `Phi=4h/15`, so `(36)` forces `h'=15/2`.  The
coefficients of `X,X^2` in `E=1` independently give

```text
h[16(3u^2+1)k+12u]=18.                                (45)
```

Now `h=(15/2)u+c`.  If `k` is nonzero, the left side of `(45)` has degree
`deg k+3`; if `k=0`, it has degree two.  It cannot equal `18`.  Hence every
Keller projection through the nodal cylinder needs a genuine second
conductor layer `X^2`, `XY`, or `Y^2`.  This is again an all-degree no-go.

Finally, the period defect of `(32)` can itself be repaired without spending
any scarce degree.  At `lambda=mu=1`, its odd first component is

```text
f_0+u f_1+u^3 f_3,
f_0=eta Y/2,       f_1=2XY^23,       f_3=2X^2Y^21.    (46)
```

For `P_j=X^(j+2)(X+1)`, `j=0,1,2`, the exact moment matrix
`M_rj=integral f_r dP_j`, `r=0,1,3`, is invertible.  Its first two dual
combinations `G_0,G_1` satisfy

```text
integral A_* dG_0=1,          integral A_* dG_1=u.     (47)
```

Writing the degree-six period residual as
`R=sum_(i=1)^6 r_i u^i`, the correction

```text
Delta B=(sum_(i=1)^5 r_i u^i)G_0+r_6u^5G_1            (48)
```

makes `(36)` exact.  It has target cap at most `10`, width at most five, and
contains `X^2(X+1)`, so it preserves `(4c)`, `(26)`, and all three root-line
jets.  Exact evaluation still gives `Jac(A_*,B_*+Delta B)!=(1)` at `(0,2)`.
Thus the period-repaired object is a stronger positive counterexample cell,
not a Keller pair.

### 3.4 The node carries a unipotent jet holonomy

The two normalization branches retain one exact piece of coherent data that
is invisible in the common boundary values.  Put

```text
gamma'(u)=(2u,3u^2-1),
v_+(u)=(A_t,B_t)(u,1),       v_-(u)=(A_t,B_t)(u,-1).  (49)
```

For a normalized Keller jet, both branch frames

```text
M_+=[gamma',v_+],            M_-=[gamma',v_-]          (50)
```

have determinant one.  Their second-column difference is therefore parallel
to their common first column:

```text
v_+-v_-=delta(u)gamma',
M_+=M_- [1 delta;0 1].                                (51)
```

In the normal form `(37)`, the endpoint equations give
`(a_X,b_X)=h(u)gamma'` at `X=0`, while
`v_+-v_-=4(a_X,b_X)`.  Hence

```text
delta=4h.                                              (52)
```

This is an intrinsic polynomial unipotent holonomy of the node.  Higher
conductor layers vanish at `X=0`, so they can change the equations that govern
`h` but cannot erase `(51)--(52)`.  Strong dephasing keeps the squared sizes
of the two local matching channels but discards precisely this branch-relative
shear.

Every unit transverse vector along the minimal boundary has the form

```text
v=(1,3u/2)+k(u)(2u,3u^2-1).                           (53)
```

This frame has an exact elementary-factor ledger.  For
`E_+(z)=[1 z;0 1]` and `E_-(z)=[1 0;z 1]`, direct multiplication gives

```text
[gamma',v]
 =E_-(3u/2-1) E_+(1) E_-(2u-1) E_+(k).               (53a)
```

Thus the canonical boundary is already a three-shear `SL_2` path, and its
arbitrary transverse choice is one final shear.  The relative factor between
the two normalization branches is exactly the last unipotent factor in
`(51)`.  This is a positive factorization control, not an assertion that the
frame extends to a globally integrable polynomial Jacobian matrix.

Consequently, with `T=A_tB_u`, the boundary matching products are

```text
T=(1+2uk)(3u^2-1),
T+1=2u(3u/2+(3u^2-1)k).                               (54)
```

At `u=0` the determinant is purely off diagonal (`T=-1,T+1=0`), whereas at
`u=+/-1/sqrt(3)` it is purely diagonal (`T=0,T+1=1`).  These three channel
switches are independent of `k` and of every higher correction.  They are the
boundary version of the two-path `K_(2,2)` plaquette: the unit determinant is
transported between the two matchings while their coherent difference stays
fixed.

The cheapest normal thickening already exposes why polynomialization is hard.
For the canonical choice

```text
v_0=(1,3u/2),             F_lin(u,r)=gamma(u)+r v_0(u), (55)
```

one has

```text
Jac_(u,r)(F_lin)=1-3r/2.                               (56)
```

A marked reparameterization `r=r(s)`, `r(0)=0`, has unit Jacobian exactly
when

```text
(1-3r/2)r'=1,
s=r-3r^2/4,
r=(2/3)(1-sqrt(1-3s)).                                 (57)
```

Thus the linear normal model classicalizes into a Catalan/Kummer square-root
series rather than a polynomial.  Allowing the general transverse shear in
`(53)` does not make the linear model polynomial-flat: direct differentiation
gives

```text
det(v',v)=k'-3/2-6uk-2(3u^2+1)k^2.                    (58)
```

No polynomial `k` makes `(58)` zero.  If `deg k=n>=1`, the quadratic term has
the unique top degree `2n+2`; a nonzero constant `k` leaves a nonzero quadratic
term, and `k=0` leaves `-3/2`.  This closes only polynomial thickenings that
are linear in the normal parameter.  Quadratic and higher normal terms are
exactly the open resource used by `(29)--(48)`; no all-degree obstruction is
being inferred from `(57)`.

## 4. The open correction problem

The first serious nodal-cylinder search cell is now much narrower:

1. start at target cap `38` from the two-scale base `(29)`--`(32)`, not the
   dead width-`(2,3)` scaffold `(18)`;
2. retain scaled widths at least `(4,6)`, the degree-`32` side base or a legal
   replacement, the nodal root-line specializations, and the four-point fibre;
3. use genuine second-conductor-layer terms, since `(43)`--`(45)` exclude the
   entire first layer;
4. impose the period equation `(36)` from the start--the finite-exact repair
   `(46)`--`(48)` shows that this costs no high-degree budget; and
5. solve `E=1,O=0` coefficient by coefficient, while rejecting any survivor
   injective on even one affine source line and quotienting target shears and
   the image relation `(2)`.

The cell is attractive because the packet obstruction has disappeared and
the first admissible cell now passes high-face resonance, Newton similarity,
all roots of its lower leading coefficient, the conductor period, and the
first-jet gates simultaneously.  What remains is pointwise integrability:
cancel the nonconstant coefficients of `(42)` without losing the every-line
collision curve.  The cheapest next computation is a sparse `X`-adic solve
beginning at conductor order two, followed by a saturated equal-horizontal-
fibre test.  The closed form `(8)` remains the non-potential positive control.

## Reproduction and scope

Run:

```bash
python3 04-computation/nodal_cylinder_keller_projection_scaffold.py
python3 -O 04-computation/nodal_cylinder_keller_projection_scaffold.py
python3 04-computation/nodal_cylinder_width_period_gates.py
python3 -O 04-computation/nodal_cylinder_width_period_gates.py
```

The first companion verifies `(1)`--`(21)` with exact rational symbolic
arithmetic (the existence of a complex root of `alpha^2+1` is certified by the
displayed nonconstant polynomial and a coprimality gate).  The independent
`nodal_cylinder_width_period_gates.py` companion checks the finite-exact parts
of `(22)`--`(48)`: the width ledger, relation-reduced two-scale cube, all three
root-line jets, period pairing signs, PDE split, first-layer equations, exact
moment-dual period repair, and a nonconstant-Jacobian witness after repair.
Equations `(49)`--`(58)` are direct two-by-two determinant and degree
calculations; they do not enlarge the finite-exact claim of either companion.
The degree-`38` floor additionally uses the cited reduced-height theorem and
the current proved proper-power/infinity-root gates; `(22)` uses Lang's cited
Newton-polygon similarity theorem and `(4a)` uses Gwozdziewicz's cited
injectivity-on-one-line theorem.  No claim is made that the open correction
equations are solvable, or that every useful collision surface is a nodal
cylinder.
