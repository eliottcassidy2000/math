---
id: THM-3888
title: "F-zero equianharmonic Jacobian and two-section integrality"
status: >
  PROVED + INDEPENDENTLY HOSTILE-AUDITED BASE; FINAL EXACT-SHELL
  STRENGTHENING AWAITING FOCUSED AUDIT.  The
  THM-3881 f=0 square residual is a binary quartic with I=0 and an explicit
  equianharmonic Jacobian.  Over the generic x-line its y-fibration is a
  rational elliptic surface with fibres II,II,II,II,IV and geometric
  Mordell--Weil group of rank six and trivial torsion.  Polynomial residual
  squares give sections integral away from two explicitly marked boundary
  sections whose sum is the second T=0 section.  The converse is asserted
  only for sections factoring through the original affine quartic chart.
  An explicit nonzero k(x)[y] point proves that x-integral descent and the
  origin address are load-bearing.  Within the integral Weierstrass-coordinate
  shell, there are six constant-u sections; exactly two give polynomial
  points of the chosen affine quartic chart, namely the second T=0 section
  and one explicit x-polar hostile section.  Every other polynomial quartic
  point has a genuine Weierstrass denominator.
  Effective two-section S-integral enumeration, x-integral descent, a Keller
  atlas, and JC(2) remain OPEN.
source: jc_zero_debt_lift / post-THM-3885 elliptic reframe, 2026-08-23
audit: >
  PROVED BASE WITH PROVISIONAL EXACT STRENGTHENING.  Independent hostile
  audits passed the elliptic-surface packet at ce633e56c5 and the first
  integral-height delta at 58b095278c.  The companion verifies the binary
  quartic invariants and discriminant, the birational map and inverse, all
  four marked sections and their chord law, the T-divisor factorization, the
  generic squarefree y-discriminant, short-Weierstrass discriminant, infinity
  minimalization, IV component addresses, the normalized cubic-factor chart,
  the explicit descent-hostile point, the 3-division polynomial, and the
  fifteen-case constant-T rational-root exhaustion and all four
  alternate/boundary constant-u sections in 68 active gates.
  Normal and optimized runs
  must byte-match the frozen
  output.  Independent audit must recheck extension to smooth projective
  models, the exact one-way integrality scope at bad fibres, rationality and
  Shioda--Tate inputs, torsion injection at IV, and the local intersection of
  the two boundary sections.  The final low-height rational-root exhaustion
  awaits a focused delta audit.
depends_on:
  - THM-3881-cusp-ideal-residual-transport-rank-two-matrix-factorization
  - THM-3885-cusp-residual-f-zero-arm-dichotomy-and-quadratic-closure
related:
  - THM-3874-three-cusp-quadratic-k3-affine-class-group
  - THM-3884-cusp-residual-total-degree-leading-gauge-filtration
  - THM-3886-cusp-residual-equality-seam-second-layer-trichotomy
script: 04-computation/jc2_f_zero_equianharmonic_jacobian_integrality_thm3888.py
output: 05-knowledge/results/jc2_f_zero_equianharmonic_jacobian_integrality_thm3888.out
script_sha256: 329598af4a09c2ff3aa37db65ef176975f4df218d48c2c4f6ecd74e30ec04e47
output_sha256: b612c4f412c2610197aa70c57aee52879ec3fd021610e60982caa91c07913203
semantic_sha256: 9a50da8fdffac4a46ca3d9074d86f2e1c85943a5150792ec486c41a2c062fbb6
hash_basis: raw LF bytes
---

# THM-3888 -- the f-zero lane is a two-boundary-section problem

**PROVED + INDEPENDENTLY HOSTILE-AUDITED BASE; FINAL EXACT-SHELL
STRENGTHENING AWAITING FOCUSED AUDIT.**
Work over an algebraically closed field `k` of characteristic zero.  Retain
the THM-3881 notation

```text
a=x+1,                         L=9x+4,
F=15x^2+15x+4,                 K=y^2-F,
Delta=a^3L^2-K^2.                                           (1)
```

For the `f=0` lane, a residual square is a pair `T,G in k[x,y]` satisfying

```text
G^2=q(T):=L^4-6aL^2T^2-8KT^3-3a^2T^4.                    (2)
```

The curve `(2)` is not an unstructured quartic.  It is an equianharmonic
elliptic curve with two marked points removed.  This reframe proves the exact
generic elliptic-surface anatomy and identifies the remaining problem as an
effective two-section integrality problem plus descent in `x`.

## 1. The binary quartic has `I=0`

For a binary quartic with affine coefficients

```text
A_4T^4+B_3T^3+C_2T^2+D_1T+E_0,
```

use the classical invariants

```text
I=12A_4E_0-3B_3D_1+C_2^2,
J=72A_4C_2E_0+9B_3C_2D_1-27A_4D_1^2
    -27B_3^2E_0-2C_2^3.                                  (3)
```

For `(2)`, direct substitution gives

```text
I=0,
J=1728L^4Delta,
disc_T(q)=-110592L^8Delta^2.                              (4)
```

Thus the generic smooth projective completion has `j=0`.  The square of
`Delta` in `(4)` is the exact reason the quartic residual and the repo's
recurring `C_3` grammar meet here.

## 2. Exact birational map and inverse

On the dense chart `T!=0`, put

```text
X_0=(G+L^2)/T^2,
Y_0=2(X_0^2+3a^2)T+8K.                                   (5)
```

Reduction modulo `(2)` gives the exact cubic

```text
Y_0^2=8L^2(X_0-a)^3-64Delta.                              (6)
```

After the scaling

```text
X=2L^2(X_0-a),                    Y=L^2Y_0,                (7)
```

this is the short equianharmonic model

```text
E: Y^2=X^3-64L^4Delta.                                    (8)
```

Conversely, on the dense chart where the displayed denominator is nonzero,

```text
D_E=X^2+4aL^2X+16a^2L^4,
T=2L^2(Y-8KL^2)/D_E,
G=(a+X/(2L^2))T^2-L^2.                                   (9)
```

Equations `(5)-(9)` are inverse rational maps.  Since the generic quartic and
`E` are smooth projective genus-one curves and `(T,G)=(0,L^2)` is rational,
the map extends uniquely to an isomorphism of their smooth projective models,
with this point chosen as the origin `O`.

### 2.1. The normalized cubic factor--cofactor chart

There is a more revealing scaling than `(7)`.  Put

```text
u=X/(4L^2),                         v=Y/(8L^2).             (9a)
```

Then `(8)` and its inverse become

```text
v^2=K^2+L^2(u^3-a^3),
T=(v-K)/(u^2+au+a^2),
G=(a+2u)T^2-L^2.                                          (9b)
```

Equivalently,

```text
(v-K)(v+K)=L^2(u-a)(u^2+au+a^2).                         (9c)
```

Thus the two deleted sections are the two nontrivial cube-root addresses of
`u^3=a^3`, while `P_0` is the address `u=a`.  On any Weierstrass section for
which `u,v in k(x)[y]`, the first inverse in `(9b)` is polynomial exactly when

```text
u^2+au+a^2 divides v-K in k(x)[y].                        (9d)
```

This literal divisibility statement is scoped to integral Weierstrass
coordinates.  For a general Mordell--Weil section, `u,v` may have
denominators and `(9d)` must be replaced by its valuation form.  This is why
the first integral height shell is the cheapest exact computation, but not
the whole problem.

### 2.2. The integral Weierstrass shell has six points, only two affine

The first height-shell computation admits an all-degree proof.  Suppose

```text
u,v in k(x)[y],              T=(v-K)/(u^2+au+a^2) in k(x)[y].   (9e)
```

Put `m=deg_y(u)`.  If `m>=2`, equation `(9b)` has leading degree `3m`, since
`deg(K^2)=4`.  Therefore an odd `m` is impossible, while for even `m`

```text
deg_y(v)=3m/2 < 2m=deg_y(u^2+au+a^2).                    (9f)
```

The same strict inequality holds at `m=2`.  Hence the divisibility in `(9e)`
forces `v-K=0`.  Equation `(9b)` then gives `u^3=a^3`, contradicting `m>=2`.
Thus

```text
u,v,T polynomial in y  ==>  deg_y(u)<=1.                 (9g)
```

The two remaining degrees are finite grammars, and both can be completed
exactly.  If `m=1`, then `deg_y(v)=2`, so `(9e)` forces `T=t in k(x)`.
Writing `(2)` as a polynomial in `y` gives

```text
q(t)=-8t^3y^2+C(t),
C(t)=L^4-6aL^2t^2+8Ft^3-3a^2t^4.                        (9h)
```

There is no `y` term.  If `q(t)` is a square and `t!=0`, its square root is
linear in `y`; comparison of the missing linear term forces its constant
term to vanish, and therefore `C(t)=0`.  The polynomial `C` is primitive in
`k[x][t]`.  Its constant and leading coefficients are `L^4` and `-3a^2`,
respectively.  The UFD rational-root theorem therefore leaves only

```text
t=c L^i/a^j,       c in k^*,       0<=i<=4, 0<=j<=2.     (9i)
```

For each of these fifteen pairs, clear the power of `a` in `C(t)`, collect
coefficients in `x`, and take their gcd in `Q[c]`.  The exact companion
obtains the unit ideal for every pair.  Since these identities persist after
extension to the algebraically closed constant field, none of `(9i)` is a
root.  The case `t=0` would force `v=K` and then `u^3=a^3`, contradicting
`m=1`.  Thus the `m=1` shell is empty.

If `m=0`, the right side of `(9c)` is constant in `y`.  If the product
`(v-K)(v+K)` were a nonzero constant, both factors would be units in
`k(x)[y]`, contradicting their nonconstant difference `2K`.  Hence the
product is zero and

```text
u^3=a^3,                         v=+K or -K.               (9j)
```

Write `u=zeta*a`, where `zeta^3=1`.  The choice `zeta=1` gives the two
polynomial quartic points

```text
(u,v)=(a,K):       (T,G)=(0,-L^2),
(u,v)=(a,-K):      (T,G)=(T_*,G_*).                      (9k)
```

There are four further integral Weierstrass sections, and they must not be
discarded by dividing the inverse denominator.  For either nontrivial cube
root,

```text
u^2+au+a^2=a^2(zeta^2+zeta+1)=0,                         (9l)
```

For `v=-K` these are exactly the deleted sections `Q_+,Q_-`.  For `v=+K`
both numerator and denominator of `(9b)` vanish.  Put

```text
s=1+2zeta,                         s^2=-3.
```

Then `X_0=as`, so an affine quartic preimage would have
`G=asT^2-L^2`.  Exact substitution in `(2)` gives

```text
G^2-q(T)=-2T^2(-4KT+aL^2(s-3)).                          (9m)
```

The nonzero preimage is therefore

```text
T=aL^2(s-3)/(4K),                                        (9n)
```

which is not in `k(x)[y]`.  Thus the complete integral-Weierstrass shell has
six sections: two chosen-chart polynomial residual points, two deleted
boundary sections, and two alternate-chart rational points.  Consequently
every **other polynomial quartic point** has a denominator in `u` or `v`;
cancellation through `(9b)` is now the genuinely remaining height problem.

## 3. The four marked sections and the divisor of `T`

Choose `s in k` with `s^2=-3`.  Besides the origin, the second point over
`T=0` is

```text
P_0=(4aL^2,8KL^2) on E.                                  (10)
```

Homogenize `(2)` in weighted coordinates `[T:U:G]` of weights `(1,1,2)`.
The two points at `T=infinity`, namely

```text
[1:0:+sa],                    [1:0:-sa],                  (11)
```

map to the two sections

```text
Q_+=(2aL^2(s-1), -8KL^2),
Q_-=(2aL^2(-s-1),-8KL^2).                                (12)
```

Indeed, in the local coordinate `1/T`, the two square-root branches begin

```text
G=epsilon*s*a*T^2-4K/(epsilon*s*a) T+O(1),
epsilon in {+1,-1};
```

substitution in `(5)` gives `(12)`.  Thus these are the actual two ends of
the quartic chart, not merely two points which happen to lie on `(8)`.

Their `Y`-coordinates agree, so their chord is horizontal.  The short
Weierstrass addition law gives exactly

```text
Q_+ + Q_- = P_0.                                          (13)
```

On the generic elliptic curve the coordinate `T` has two simple zeroes and
two simple poles.  Therefore

```text
div_E(T)=O+P_0-Q_+-Q_-.                                   (14)
```

The elementary factor shadow of `(14)` is

```text
(G-L^2)(G+L^2)
 =-T^2(3a^2T^2+8KT+6aL^2).                               (15)
```

Over `k(x)[y]`, the two left factors are coprime because their difference is
the unit `2L^2`.  Thus every root of `T` is polarized between the two `T=0`
sections, while poles are polarized between `Q_+` and `Q_-`.  This is the
global elliptic version of THM-3885's finite root-partition grammar.

## 4. What polynomiality means -- and what is not yet proved

Regard `(8)` as an elliptic fibration over the `y`-line with coefficient
field `k(x)`.  A polynomial solution

```text
T,G in k(x)[y]                                             (16)
```

defines a rational section of its smooth projective model.  Over every finite
point of the `y`-line, the original coordinates remain in the affine chart
`U=1`; hence the lifted section has no finite contact with either deleted
boundary section `Q_+` or `Q_-`.  Its boundary contact is supported at
`y=infinity`.  In other words, every polynomial residual square supplies a
two-section `S`-integral point for

```text
S={y=infinity},                    boundary=Q_+ union Q_-.  (17)
```

There is a safe converse with the integrality condition stated rather than
inferred: any section which, over `A^1_y`, factors through the original affine
quartic chart has regular pullbacks of `T,G`, hence gives elements of
`k(x)[y]` satisfying `(2)`.  Merely knowing on a resolved model that a section
does not meet the strict transforms of `Q_+` and `Q_-` is **not asserted** to
be sufficient at singular fibres; vertical exceptional data may be needed.

The original JC lane is narrower still.  It requires the coefficients in
`(16)` to lie in `k[x]`, with no poles in the `x` direction, and it retains
the address

```text
T(0,0)=0                                                   (18)
```

and THM-3885's `a=0` arm dichotomy.  These descent and integrality conditions
are invisible to the geometric Mordell--Weil rank computed below.

They are also genuinely necessary.  The normalized integral section

```text
u=a,                         v=-K                           (18a)
```

gives the nonzero exact solution over `k(x)[y]`

```text
T_*=-2K/(3a^2),
G_*=4K^2/(3a^3)-L^2.                                     (18b)
```

Direct substitution gives `G_*^2=q(T_*)`.  But

```text
T_*(0,0)=8/3,                    ord_(a=0)(T_*)=-2,        (18c)
```

and `G_*` has an `a=0` pole of order three.  Thus the generic two-section
integral problem is **not empty**.  This point is excluded precisely by the
global `k[x,y]` descent and origin address, so neither condition may be
dropped from a future closure argument.

## 5. Generic fibre packet over the `y`-line

Pass to `C=overline{k(x)}`.  As a polynomial in `y`,

```text
Delta(y)=a^3L^2-(y^2-F)^2.                                (19)
```

It has four simple roots.  One exact certificate is

```text
disc_y(Delta)=256a^6L^4(1-a)^3(9a-4)^2 !=0               (20)
```

in `k(x)`, using

```text
F^2-a^3L^2=(1-a)^3(9a-4)^2.                              (21)
```

For `(8)`, the short-Weierstrass discriminant is

```text
Delta_E=-2^16 3^3 L^8Delta(y)^2.                          (22)
```

Each simple root of `(19)` has `ord(a_6)=1`, hence is a Kodaira `II` fibre.
At infinity put `u=1/y` and make the degree-one elliptic scaling

```text
X_infinity=u^2X,                    Y_infinity=u^3Y.       (23)
```

The transformed constant coefficient is

```text
a_(6,infinity)
 =64L^4u^2-128FL^4u^4-64L^4(a^3L^2-F^2)u^6.              (24)
```

Its order is exactly two, so the infinity fibre is `IV`.  There are no other
singular fibres.  The complete packet is therefore

```text
II, II, II, II, IV,                  Euler sum=4*2+4=12.   (25)
```

The fundamental line bundle has degree one; the minimal surface is a rational
elliptic surface.  Thus its geometric Neron--Severi rank is ten.  The sole
reducible fibre contributes the root lattice `A_2`, so Shioda--Tate gives

```text
rank E(C(y))=10-2-rank(A_2)=6.                            (26)
```

The torsion subgroup is zero.  At the additive `IV` fibre, characteristic
zero makes the formal group and additive identity component torsion-free, so
torsion injects into the component group `Z/3`.  It remains only to exclude
three-torsion.  For `Y^2=X^3+B`, where `B=-64L^4Delta`,

```text
psi_3=3X(X^3+4B)=3X(X^3-256L^4Delta).                    (27)
```

The branch `X=0` would require `Delta` to be a square in `C(y)`; the other
branch would require `Delta` to be a cube.  Both are impossible because the
four valuations in `(20)` are one.  Hence

```text
E(C(y)) ~= Z^6.                                           (28)
```

## 6. The infinity component addresses

The two boundary sections are not generic anonymous Mordell--Weil points.
In the `IV` chart, blow up with

```text
X_infinity=uX_1,                    Y_infinity=uY_1.       (29)
```

At `u=0`, the two nonidentity components have addresses

```text
Y_1=+8L^2,                          Y_1=-8L^2.             (30)
```

The marked sections specialize as follows:

```text
P_0:          Y_1=+8L^2,
Q_+,Q_-:      Y_1=-8L^2.                                  (31)
```

Moreover, in this smooth blow-up chart

```text
X_1(Q_+)-X_1(Q_-)=4saL^2u.                               (32)
```

Thus `Q_+` and `Q_-` meet once at infinity and separate to first order; they
have no finite meeting forced by the model.  This is the precise local
sidecar that an effective two-section integral-point enumeration must retain.

## 7. Exact surviving problem and cheapest decisive computation

The rank-six result `(28)` explains why Shioda--Tate alone cannot close the
lane: the relevant subset is not the full Mordell--Weil group.  The remaining
problem is

```text
rank-six MW sections
 intersect two-section S-integrality
 intersect k[x]-integral descent
 intersect the origin and a=0 addresses.                  (33)
```

The cheapest decisive next computation is therefore not another unrestricted
quartic coefficient Groebner basis.  It is:

1. compute an explicit Mordell--Weil basis for the `II^4+IV` rational
   elliptic surface and express `Q_+,Q_-,P_0` in that basis;
2. use `(9k)-(9n)` as the completed integral-coordinate base case and enumerate
   the first rational-coordinate denominator shell, subject to no finite
   intersection with `Q_+ union Q_-` and the `IV` address `(31)-(32)`;
3. apply the inverse denominator test

   ```text
   D_E divides 2L^2(Y-8KL^2)                              (34)
   ```

   and then demand that the quotient and reconstructed `G` lie in `k[x,y]`
   and satisfy `(18)`.

A finite shell search is evidence, not a proof, until accompanied by an
effective height bound or a classification of the two-section integral
points.  That height bound, the full all-degree `f=0` closure, the `f!=0`
lane, a Keller atlas, and JC(2) remain **OPEN**.

Reproduce the exact algebraic packet with

```bash
python3 04-computation/jc2_f_zero_equianharmonic_jacobian_integrality_thm3888.py
python3 -O 04-computation/jc2_f_zero_equianharmonic_jacobian_integrality_thm3888.py
```

Both streams must byte-match
`05-knowledge/results/jc2_f_zero_equianharmonic_jacobian_integrality_thm3888.out`.
