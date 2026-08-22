---
id: THM-3586
title: "Nodal-cylinder cap-38 width, period, and second-conductor Keller gates"
status: >
  PROVED structural identities and all-degree consequences + FINITE-EXACT
  cap-38 controls + CITED global inputs.  For the normalization
  (u,t)->(u,t^2-1,t(t^2-1)), every polynomial Keller projection would be
  noninjective on every affine source line, have an immersed noninjective
  conductor restriction and at least a four-point fibre.  In the cited
  reduced (72,108) cell its first possible target cap is 38; cited Newton
  similarity and the every-line gate exclude width (2,3), so the first
  unexcluded width is (4,6) and needs two coefficient scales.  An exact
  conductor period and normalization PDE hold generally; with the minimal
  nodal boundary, unipotent node-jet holonomy and an all-degree no-I^2
  obstruction force a genuine second conductor layer.
  The displayed (4,6) scaffold and its exact moment-dual period repair both
  have nonconstant Jacobian.  The residual second-layer Keller PDE is OPEN;
  no planar Jacobian counterexample is constructed.
source: boxeph / nodal-cylinder projection session, 2026-08-20--21
audit: >
  PASS.  An independent hostile audit rederived the cap-38 relation
  reduction, all-degree width-(2,3) exclusion, degree-32 second scale,
  period and PDE formulas, first-conductor no-go, and moment-dual repair.
  Ordinary, optimized, and stored companion transcripts are byte-identical.
depends_on: []
related:
  - THM-3555-catalan-thickening-universal-cubic-root-cover
  - THM-3556-cusp-square-packet-marked-root-kummer-owner
  - THM-3563-nongraph-conductor-node-cycle-keller-obstruction
  - THM-3567-separated-rational-keller-maximal-observable-nodal-completion-no-go
  - THM-3572-squarefree-danielewski-affine-modification-and-two-bracket-collapse
external:
  - "Gwozdziewicz, Injectivity on one line, Theorem 1.1, arXiv:alg-geom/9305008."
  - "Lang, Newton polygons of Jacobian pairs, DOI:10.1016/0022-4049(91)90128-O."
  - "Shastri, one-point-at-infinity criterion, DOI:10.18910/6794."
  - "Guccione--Guccione--Horruitiner--Valqui, arXiv:2204.14178v1, reduced degree pairs below 125."
scripts:
  - 04-computation/nodal_cylinder_keller_projection_scaffold.py
  - 04-computation/nodal_cylinder_width_period_gates.py
outputs:
  - 05-knowledge/results/nodal_cylinder_keller_projection_scaffold.out
  - 05-knowledge/results/nodal_cylinder_width_period_gates.out
script_sha256:
  - de61c326e049e056f204d21f04fc2265e43742411d26bd6703b504f9ed0f40df
  - a2483f1da3237df5896c9fd508b2d137bca5fb82a67ffec4043accbe088fe050
output_sha256:
  - 0df404c494ef6b579d681203d521da2f8896fdf25cd0b96f0e179b01e64ce897
  - de3514b83157b63558deefee1e8137cd6ba5e99b13cfa813d142c6b4c8fc9f16
hash_basis: raw LF bytes
---

# THM-3586 -- nodal-cylinder cap-38, width, period, and conductor gates

**PROVED structural consequences + FINITE-EXACT controls + CITED global
inputs; OPEN residual Keller PDE.**  This theorem does not construct a
Keller pair.  In particular, neither displayed scaffold in Section 6 has
constant Jacobian.

All algebraic identities are over a characteristic-zero field.  The four
global inputs in Section 3 are invoked over `C`, where their cited forms
apply.  A nonzero constant Jacobian is normalized to one by rescaling the
second component.

## 1. The normalization packet and its collision

Put

```text
X=t^2-1,                 Y=t(t^2-1),
nu(u,t)=(u,X,Y),
H=Y^2-X^2(X+1).                                           (1)
```

Then `nu` is the normalization of

```text
S=Spec C[u,X,Y]/(H).                                      (2)
```

Indeed, `t` is integral, `t=Y/X` off `X=0`, and `C[u,t]` is normal.  More
explicitly,

```text
C[t^2-1,t(t^2-1)]={f in C[t]:f(1)=f(-1)}.               (3)
```

For the reverse inclusion in `(3)`, subtract the common endpoint value,
divide by `t^2-1`, and split the quotient into its even and odd parts.  The
conductor ideal is `I=(X,Y)`, and its inverse image is the pair of lines
`t=+1,-1`.

The tangent minors are

```text
M_uX=2t,                 M_uY=3t^2-1=3X+2,
M_XY=0.                                                    (4)
```

They never vanish together.  Thus `nu` has differential rank two
everywhere, but

```text
nu(u,1)=nu(u,-1)=(u,0,0).                               (5)
```

Consequently any pair `A,B in C[u,X,Y]/(H)` pulls back to a planar map
retaining the collision `(5)`.

## 2. The descended normal gate passes, but one closed certificate is not a wedge

The image normal and tangent normal satisfy

```text
(H_u,H_X,H_Y)(nu)=X(0,-M_uY,M_uX).                      (6)
```

Unlike the cusp-square packet of THM-3556, the multiplier `X` belongs to the
image Jacobian ideal.  Modulo `H`,

```text
X=-((3X+1)/2)H_X-(9/4)YH_Y,                            (7)
```

and division of `(7)` after pullback gives the unit tangent-minor identity

```text
1=((3X+1)/2)M_uY-(9/4)YM_uX.                           (8)
```

Even closedness can be imposed.  The polynomial two-form

```text
omega=-(9/4)Y du wedge dX
      +(3X+1)/2 du wedge dY
      +(15/4)u dX wedge dY                              (9)
```

is closed and pulls back to `du wedge dt`.

This particular certificate is nevertheless not `dA wedge dB` for any
polynomials `A,B in C[u,X,Y]`.  With `V=3X+1`, the vector field dual to four
times `(9)` is

```text
D=15u partial_u-6V partial_V-9Y partial_Y.              (10)
```

If `(9)=dA wedge dB`, then `D(A)=D(B)=0`.  Every invariant monomial obeys

```text
5a=2b+3c             for u^a V^b Y^c.                  (11)
```

At `u=0`, equation `(11)` forces `b=c=0`.  Hence every polynomial invariant
restricts to a constant on that plane; the two gradients are parallel there,
while `D` is generically nonzero.  This contradiction proves the claim.
Other closed certificates are not excluded.

## 3. CITED global inputs and their exact consequences

This section separates imported theorems from the consequences proved here.

### 3.1 Every-line and four-point gates

Gwozdziewicz's **CITED** theorem says that a planar Keller map injective on
one affine source line is an automorphism.  A Keller projection through
`nu` cannot be an automorphism because of `(5)`.  It is therefore
noninjective on every affine source line.

In particular, its common conductor restriction

```text
gamma(u)=F(u,1)=F(u,-1)                                (12)
```

must be noninjective.  It is immersed, since `gamma'(u)=0` would make the
full Jacobian vanish at `(u,1)`.  Thus a collision
`gamma(u_0)=gamma(u_1)`, `u_0!=u_1`, produces four distinct source points

```text
(u_0,+/-1),(u_1,+/-1)                                  (13)
```

in one fibre.

This already has a sharp degree consequence.  If a polynomial plane curve
of degree at most two has `gamma(r)=gamma(s)` with `r!=s`, writing it as
`c_2u^2+c_1u+c_0` gives

```text
(r+s)c_2+c_1=0,       gamma'((r+s)/2)=0.               (14)
```

Hence an immersed noninjective conductor curve has degree at least three.
The bound is attained by the nodal cubic

```text
gamma_0(u)=(u^2-1,u^3-u),
gamma_0(1)=gamma_0(-1)=(0,0),
gamma_0'(u)=(2u,3u^2-1)!=0.                            (15)
```

### 3.2 The first reduced target cap is 38

Give `(u,X,Y)` their pulled-back source degrees `(1,2,3)`.  Define the
target cap of a reduced projection pair using chosen representatives in
`C[u,X,Y]` modulo `H`.  A representative of total degree at most `E` pulls
back to source degree at most `3E`.

The Guccione--Guccione--Horruitiner--Valqui **CITED** sub-`125`
classification says that a hypothetical reduced counterexample in this
range has degree pair `(72,108)`, up to order.  Its top forms therefore have
the standard common-power shape

```text
A_72=alpha K^2,              B_108=beta K^3,
deg K=36.                                                   (16)
```

The common-power step in `(16)` is proved locally, not imported.  For
homogeneous binary forms of degrees `m,n`, write on the chart `t!=0`
`P=t^m p(u/t)` and `Q=t^n q(u/t)`.  Their top Jacobian is

```text
t^(m+n-2)(n p' q-m p q').
```

Its vanishing makes `p^n/q^m` constant.  Unique factorization and
`gcd(72,108)=36` then give exactly the square/cube form `(16)`, including
the appropriate power of `t`.  Thus no unlisted proper-power theorem is a
dependency here.

The homogeneous degree-`108` monomials available at caps `35,36,37,38`
have respective `u`-exponent sets

```text
empty,              {0},              {0,1},
{0,1,2,3}.                                                (17)
```

At cap `36`, `B_108` is a scalar multiple of `t^108`.  At cap `37`, it has
the form `t^107(at+bu)`: if `b!=0`, its root multiplicities are incompatible
with a cube, while `b=0` again has one projective root.  Shastri's **CITED**
one-place-at-infinity criterion excludes that one-root case.  Thus

```text
E>=38.                                                    (18)
```

At equality the first permitted two-root base is, after scaling,

```text
K=t^35(t+lambda u),             lambda!=0.              (19)
```

The statement `(18)` is specifically about a GGHV-reduced pair in the
sub-`125` degree cell and its reduced target representatives.  It is not a
bound on arbitrary unreduced presentations or on higher degree cells.

### 3.3 Newton similarity excludes width `(2,3)`

Lang's **CITED** theorem scales the origin-augmented Newton polygons of a
Keller pair of degrees greater than one.  In the `(72,108)` cell,

```text
Newt(B)=(3/2)Newt(A).                                   (20)
```

Taking maximal `u`-coordinates gives the width ladder

```text
(deg_u A,deg_u B)=(2r,3r),             r>=1.           (21)
```

The first rung is impossible.  At width `(2,3)`, let `a_2(t),b_3(t)` be the
leading `u`-coefficients.  The coefficient of `u^4` in the Jacobian is

```text
2a_2b_3'-3a_2'b_3=0.                                   (22)
```

Valuations in the UFD `C[t]` give

```text
a_2=alpha h^2,             b_3=beta h^3.               (23)
```

Both coefficients lie in the ring `(3)`.  Their coprime powers force
`h(1)=h(-1)`, so `h` also lies in that ring.  The degree-`(72,108)` leading
forms give `deg h=35`; hence `h` has a root `c`.  On the line `t=c`, the
restriction has coordinate degrees at most `(1,2)`.  It is immersed because
the full Jacobian is one, so `(14)` makes it injective.  Gwozdziewicz's
theorem would then make the whole map an automorphism, contradicting `(5)`.
Therefore

```text
(deg_u A,deg_u B)>=(4,6) in the scaled ladder.          (24)
```

At width `(4,6)`, equation `(22)` recurs for the leading coefficients:

```text
a_4=alpha r^2,             b_6=beta r^3,
r in C[t^2-1,t(t^2-1)].                                 (25)
```

For cap `38`, a target monomial carrying `u^6` leaves at most `32` other
target factors, so `deg_t b_6<=96` and

```text
deg_t r<=32.                                             (26)
```

The infinity base `(19)` has coefficient degree `35`.  Thus the first
unexcluded cap-`38` rung needs a second coefficient scale in both
components; relation reduction does not remove this width debt.

## 4. PROVED boundary holonomy, channel switch, period, and PDE

### 4.1 Unipotent node-jet holonomy

Assume the minimal boundary `(15)`.  Put

```text
v_+ = (A_t,B_t)(u,1),       v_-=(A_t,B_t)(u,-1),
M_+=[gamma_0',v_+],         M_-=[gamma_0',v_-].         (27)
```

For a normalized Keller pair, `det M_+=det M_-=1`.  Two unit-area frames
with the same first column differ by a unique unipotent shear:

```text
M_+=M_- [[1,delta],[0,1]].                              (28)
```

Write the unique normal form

```text
A=a(u,X)+Yc(u,X),              B=b(u,X)+Yd(u,X).       (29)
```

At `X=0`,

```text
v_+-v_-=4(a_X,b_X).                                    (30)
```

Equality of the two determinants makes `(a_X,b_X)` parallel to
`gamma_0'`.  Since its components are coprime, there is `ell in C[u]` with

```text
(a_X,b_X)=ell gamma_0',             delta=4ell.         (31)
```

Higher powers of `I` do not change `(30)--(31)`.  This is an intrinsic
unipotent branch-matching law, not a construction of a global Keller pair.

Every unit transverse jet on either branch has the form

```text
v=(1,3u/2)+k(u)(2u,3u^2-1).                            (32)
```

The two summands in the determinant are exactly

```text
A_tB_u=(1+2uk)(3u^2-1),
A_uB_t=2u(3u/2+(3u^2-1)k)=A_tB_u+1.                   (33)
```

At `u=0` the determinant is carried entirely by `-A_tB_u`; after adjoining
`sqrt(3)` (in particular over `C`), at `u=+/-1/sqrt(3)` it is carried
entirely by `A_uB_t`, independently of `k`.  This pure-channel switch is
only a boundary identity.

### 4.2 Linear-normal thickenings force a Catalan/Kummer branch

The simplest thickening of `(15)` already fails.  With `v` as in `(32)`,
put `F_k(u,r)=gamma_0(u)+r v(u)`.  Then

```text
Jac_(u,r)(F_k)=1+r C_k(u),
C_k=k'-3/2-6uk-2(3u^2+1)k^2.                           (34)
```

No polynomial `k` makes `C_k=0`: for nonzero `k` of degree `d`, the unique
top term has degree `2d+2` and coefficient `-6 lc(k)^2`; `k=0` leaves
`-3/2`.  Hence no linear-in-normal polynomial thickening in this family is
Keller.

For the canonical `k=0`, even a scalar reparameterization `r=r(s)` would
have to satisfy

```text
(1-3r/2)r'=1,
s=r-3r^2/4,
r=(2/3)(1-sqrt(1-3s))                                  (35)
```

on the branch marked by `r(0)=0`.  The repair is algebraic but not
polynomial.  Equation `(35)` closes only this linear-normal ansatz; quadratic
and higher normal terms remain available.

### 4.3 The conductor period

For arbitrary `A,B` in the nodal ring, define

```text
Phi_(A,B)(u)=integral_(-1)^1 A(u,t)B_t(u,t)dt.          (36)
```

Endpoint equality from `(3)` and integration by parts give the all-degree
identity

```text
Phi_(A,B)'=integral_(-1)^1 Jac(A,B)dt.                  (37)
```

Thus a normalized Keller pair necessarily satisfies

```text
Phi_(A,B)=2u+constant.                                  (38)
```

The period has a closed coefficient formula.  In `(29)`, write
`a=sum a_iX^i`, `c=sum c_iX^i`, and similarly for `b,d`.  Set

```text
beta_n=(-1)^n 2^(2n+1)(n!)^2/(2n+1)!,
K_ij=-2i beta_(i+j)/(2(i+j)+3).                         (39)
```

Parity and one beta-integral recurrence give

```text
Phi_(A,B)=sum_(i,j) K_ij(a_i d_j-b_i c_j).             (40)
```

This is a necessary global period condition, not a sufficient Keller test.

### 4.4 The exact normalization PDE

Define

```text
L(f)=(3X+2)f+2X(X+1)f_X.                               (41)
```

Direct differentiation of `(29)` splits the Jacobian uniquely as

```text
Jac(A,B)=E+tO,

E=a_uL(d)-L(c)b_u
  +2X(X+1)(c_u b_X-a_Xd_u),

O=2(a_u b_X-a_Xb_u)
  +X(c_uL(d)-L(c)d_u).                                 (42)
```

The normalized Keller equation is exactly

```text
E=1,                       O=0.                         (43)
```

## 5. PROVED all-degree first-conductor-layer no-go

Retain the boundary `(15)` and suppose **exactly** that no `I^2` term is
present:

```text
A=u^2-1+Xa_1(u)+Yc_0(u),
B=u^3-u+Xb_1(u)+Yd_0(u).                               (44)
```

The unit Jacobian on both branches parameterizes every possible first jet
by polynomials `h,k in C[u]`:

```text
a_1=2uh,                    b_1=(3u^2-1)h,
c_0=1/2+2uk,                d_0=3u/4+(3u^2-1)k.       (45)
```

Formula `(40)` becomes

```text
Phi_(A,B)=4h/15.                                        (46)
```

The Keller period `(38)` forces `h'=15/2`.  Independently, the coefficients
of `X` and `X^2` in `E=1` force

```text
h[16(3u^2+1)k+12u]=18.                                 (47)
```

Now `h=(15/2)u+c`.  If `k!=0`, the left side of `(47)` has degree
`deg(k)+3`; if `k=0`, it has degree two.  It cannot be the constant `18`.
Therefore every Keller projection with boundary `(15)` requires a genuine
second conductor term from `I^2`.  This is all-degree in `u`; it does not
exclude those second and higher layers.

## 6. FINITE-EXACT cap-38 width-`(4,6)` scaffold and period repair

The algebraic cap/width/root-line subcell is nonempty; no claim is made that
the full every-affine-line or Keller system is nonempty.  For
`lambda*mu!=0`, put

```text
V=lambda uY+mu u^2,
D=Y^12+XY^10V.                                         (48)
```

The square has target cap `26`.  The raw cube has cap `39`, but modulo `H`
it has the cap-`38` representative

```text
Q_38=Y^36+3XY^34V+3X^2Y^32V^2
     +(Y^32-X^2Y^30)V^3,

Q_38-D^3=Y^30V^3H.                                    (49)
```

After pullback, `(D^2,Q_38)` is literally `(D^2,D^3)`.  It has source
degrees `(72,108)`, widths `(4,6)`, top forms

```text
(t^35(t+lambda u))^2,
(t^35(t+lambda u))^3,                                  (50)
```

and second coefficient base `XY^10` of `t`-degree `32`, attaining `(26)`.

Let `eta=3X+1` and define

```text
A_*=D^2+u^2-1+(eta/2)Y,
B_*=Q_38+u^3-u+(3eta/4)uY.                             (51)
```

On each of `t=-1,0,1`, both `D` and `D_t` vanish, and

```text
(A_*,B_*)=(u^2-1,u^3-u),          Jac(A_*,B_*)=1.     (52)
```

The four points `(u,t)=(+/-1,+/-1)` map to `(0,0)`.  The companion checks
all identities in `(48)--(52)` symbolically and checks at an exact rational
point that the full Jacobian is nonconstant.  Thus `(51)` is a
**FINITE-EXACT scaffold, not a Keller pair**.

The period defect can also be repaired without changing the high cell.
Specialize `lambda=mu=1`.  The period of `(51)` has degree six in `u`.  Put

```text
f_0=eta Y/2,          f_1=2XY^23,        f_3=2X^2Y^21,
P_j=X^(j+2)(X+1),                         0<=j<=2.      (53)
```

The exact `3 x 3` moment matrix

```text
M_rj=integral_(-1)^1 f_r dP_j,             r=0,1,3,   (54)
```

is invertible.  Its first two dual combinations `G_0,G_1` obey

```text
integral A_* dG_0=1,             integral A_* dG_1=u. (55)
```

If the degree-six residual in `(38)` is `sum_(i=1)^6 r_i u^i`, then

```text
Delta B=(sum_(i=1)^5 r_i u^i)G_0+r_6u^5G_1            (56)
```

makes the period exact.  It has target cap at most ten and width at most
five, and its factor `X^2(X+1)` preserves all three root-line first jets.
The corrected pair still has nonconstant Jacobian: the exact witness
`Jac(0,2)-1` has nonzero numerator and denominator reducing modulo `1000003`
to `(759686,317751)`.  This is a stronger **FINITE-EXACT non-Keller
control**, not evidence that `(43)` is solved.

## 7. Failure boundaries, adjacent canon, and the open cell

The status ledger is:

```text
PROVED:       normalization, collision and rank;
              descended normal-ideal and closed-form gates;
              nonpotentiality of the displayed closed form;
              boundary holonomy and pure-channel switch;
              linear-normal no-go and algebraic Catalan repair;
              period formula, normalization PDE, no-I^2 theorem.

CITED INPUTS: Gwozdziewicz every-line theorem;
              GGHV sub-125 degree classification;
              Shastri one-point-at-infinity criterion;
              Lang origin-augmented Newton similarity.

PROVED FROM
THE INPUTS:   every-line noninjectivity and four-point fibre;
              reduced cap >=38 in the sub-125 cell;
              Newton width >=(4,6) and degree-32 second-scale debt.

FINITE-EXACT: cap-38 width-(4,6) scaffold and moment-dual period repair;
              both have nonconstant Jacobian.

OPEN:         E=1,O=0 with genuine I^2 terms, every affine source line
              noninjective, and all cap/width/root-line/fibre gates retained.
                                                                  (57)
```

The neighboring theorems have different scopes:

1. THM-3556 fails the descended normal-multiplier ideal before integrability;
   `(7)--(9)` show that this nodal packet passes that gate.
2. THM-3563 closes one different normalization by a finite three-node cycle;
   the present two-branch conductor instead leaves the period and PDE debts.
3. THM-3567 computes a separated full-field nodal completion which is
   singular and not globally etale; `nu` is everywhere immersive, but no
   plane projection has been shown Keller.
4. THM-3572 makes a smooth non-`A2` affine modification symplectically exact
   and leaves one-bracket compression.  Here `(9)` is one closed unit
   certificate but is itself not a two-potential wedge.  The analogy supplies
   no implication in either direction.

The honest next calculation is an `I`-adic coefficient solve of `(43)`
starting at order two, with `(38)` imposed from the outset and an exact
equal-horizontal-fibre test after every surviving truncation.  Nothing here
settles `JC(2)`.

## 8. Exact verification contract

Reproduce with

```bash
python3 04-computation/nodal_cylinder_keller_projection_scaffold.py
python3 -O 04-computation/nodal_cylinder_keller_projection_scaffold.py
python3 04-computation/nodal_cylinder_width_period_gates.py
python3 -O 04-computation/nodal_cylinder_width_period_gates.py
```

Both ordinary/optimized pairs are byte-identical to their stored outputs.
The companions use explicit failure raises rather than `assert`.  They check
the exact structural identities, cap and width ledgers, relation-reduced
cube, three root-line jets, four-point fibre, beta-period signs, PDE split,
first-layer equations, moment-dual correction, and exact non-Keller
witnesses.  Finite hostile boxes support the displayed all-degree proofs;
they are not extrapolated universes.

**QED for the statements marked PROVED and FINITE-EXACT; the residual Keller
PDE in `(57)` remains OPEN.**
