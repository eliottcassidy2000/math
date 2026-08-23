---
id: THM-3888
title: "F-zero equianharmonic Jacobian and two-section integrality"
status: >
  PROVED + VERIFIED-EXACT CANDIDATE AWAITING INDEPENDENT HOSTILE AUDIT.  The
  THM-3881 f=0 square residual is a binary quartic with I=0 and an explicit
  equianharmonic Jacobian.  Over the generic x-line its y-fibration is a
  rational elliptic surface with fibres II,II,II,II,IV and geometric
  Mordell--Weil group of rank six and trivial torsion.  Polynomial residual
  squares give sections integral away from two explicitly marked boundary
  sections whose sum is the second T=0 section.  The converse is asserted
  only for sections factoring through the original affine quartic chart.
  Effective two-section S-integral enumeration, x-integral descent, a Keller
  atlas, and JC(2) remain OPEN.
source: jc_zero_debt_lift / post-THM-3885 elliptic reframe, 2026-08-23
audit: >
  PROVISIONAL EXACT PROOF CANDIDATE.  The companion verifies the binary
  quartic invariants and discriminant, the birational map and inverse, all
  four marked sections and their chord law, the T-divisor factorization, the
  generic squarefree y-discriminant, short-Weierstrass discriminant, infinity
  minimalization, IV component addresses, and the 3-division polynomial in
  29 active gates.  Normal and optimized runs must byte-match the frozen
  output.  Independent audit must recheck extension to smooth projective
  models, the exact one-way integrality scope at bad fibres, rationality and
  Shioda--Tate inputs, torsion injection at IV, and the local intersection of
  the two boundary sections.
depends_on:
  - THM-3881-cusp-ideal-residual-transport-rank-two-matrix-factorization
  - THM-3885-cusp-residual-f-zero-arm-dichotomy-and-quadratic-closure
related:
  - THM-3874-three-cusp-quadratic-k3-affine-class-group
  - THM-3884-cusp-residual-total-degree-leading-gauge-filtration
  - THM-3886-cusp-residual-equality-seam-second-layer-trichotomy
script: 04-computation/jc2_f_zero_equianharmonic_jacobian_integrality_thm3888.py
output: 05-knowledge/results/jc2_f_zero_equianharmonic_jacobian_integrality_thm3888.out
script_sha256: c0fa999e83e8ed3062d842ff39a0c4b20c80a9f7f286b1186b6fd1dc36cf8676
output_sha256: 4b26d1194b8cad77d53acd7d42a08ddc35f907785c2592607eb24c5cfbb5f7a2
semantic_sha256: b04048dea4c389b1f96995634ef2b930383984e02267a407ff79d1e0ff52b4ed
hash_basis: raw LF bytes
---

# THM-3888 -- the f-zero lane is a two-boundary-section problem

**PROVED + VERIFIED-EXACT CANDIDATE AWAITING INDEPENDENT HOSTILE AUDIT.**
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
2. enumerate the first height shells subject to no finite intersection with
   `Q_+ union Q_-` and the `IV` address `(31)-(32)`;
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
