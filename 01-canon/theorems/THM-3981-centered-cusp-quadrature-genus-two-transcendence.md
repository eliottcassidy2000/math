---
id: THM-3981
title: "Centered cusp quadrature has a genus-two transcendence obstruction"
status: >
  PROVED + VERIFIED-EXACT / AWAITING INDEPENDENT HOSTILE AUDIT. For every
  centering slope lambda, the generic-Y curve underlying the centered cusp
  quadrature carries a differential that remains nonexact on every finite
  algebraic cover. For lambda!=0 it is an integral genus-two curve and the
  differential is holomorphic with divisor equal to the finite zero-color
  point plus infinity. For lambda=0 the curve is rational and the
  differential has three simple poles with nonzero residues. Consequently
  the canonical centered formal X-coordinate is transcendental over the
  generic-slice function field and cannot terminate algebraically in that
  gauge. This is not an obstruction to other formal gauges or to JC(2).
source: jc-extra-debt-local / post-THM-3979 generic-slice audit, 2026-08-24
depends_on:
  - THM-3979-two-color-formal-cusp-darboux-lifting
related:
  - THM-3977-simultaneous-cusp-arm-family-critical-resultant
  - THM-3980-all-height-canonical-split-formal-cusp-nonalgebraization
script: 04-computation/jc2_centered_cusp_quadrature_genus2_thm3981.py
output: 05-knowledge/results/jc2_centered_cusp_quadrature_genus2_thm3981.out
script_sha256: be8cb9201d8b48ba23cc1db72fc7cb7e7a5722d0a6ff2545163b4c2bbae0dbbc
output_sha256: 5e879daec0413a9177baac4a9a626ab1a834da1adf6aa656bccdb7d988013427
semantic_sha256: d99e424d8aeb80a290bd802dd3944298035a17409a3f0a45f69307ee1b047474
hash_basis: raw LF bytes
---

# THM-3981 -- the centered cusp quadrature is algebraically nonintegrable

**PROVED + VERIFIED-EXACT / AWAITING INDEPENDENT HOSTILE AUDIT.** Let `k` be
an algebraically closed field of characteristic zero, fix `lambda in k`,
and put `K=k(Y)`. Consider the affine `K`-curve

```text
C_0: z(z-1)^2=x^3(Y+(lambda/2)x).                         (1)
```

Let `C` be its smooth projective normalization. If `lambda!=0`, then `C`
is geometrically integral of genus two. The differential

```text
omega=2x dx/((z-1)(3z-1))                                (2)
```

extends to a nonzero holomorphic differential on `C`, with exact divisor

```text
div_C(omega)=P_0+P_infinity,                              (3)
```

where `P_0` is the normalization point `(x,z)=(0,0)` and `P_infinity` is
the unique point above `x=infinity`. If `lambda=0`, the curve is rational
and `(2)` instead has three simple poles with nonzero residues.

More strongly, after base change to an algebraic closure `Kbar`, if
`pi:C' -> C_Kbar` is any finite cover by a smooth projective curve, then
`pi^*omega` is not exact in `Kbar(C')`. Hence the centered formal
quadrature

```text
X^2=2 integral_0^x s/((z-1)(3z-1)) ds                  (4)
```

on the generic `Y`-slice is transcendental over `K(C)`: its formal germ is
not algebraic over the generic-slice function field. The result holds for
every centering slope `lambda`; `lambda=1` is the normal form used in
THM-3979, while general `lambda` matches the linear cusp jet
`2L-lambda*y` occurring across the family of THM-3977.

## 1. Integrality and the degree-three cover

First suppose `lambda!=0`.

Write

```text
phi(z)=z(z-1)^2,
h(x)=x^3(Y+(lambda/2)x).                                  (5)
```

The polynomial `phi(z)-h(x)` is a monic cubic in `z`. If it were reducible
over `Kbar(x)`, cubicity would force a linear factor and hence a root
`z=r(x) in Kbar(x)`. This would factor the degree-four rational map `h` as

```text
h=phi composed with r.
```

Degrees of nonconstant rational maps multiply, giving
`4=3 deg(r)`, an impossibility. Thus the cubic is irreducible even over
`Kbar(x)`; Gauss's lemma proves that `(1)` is geometrically integral. In
particular

```text
x:C -> P1_K
```

has degree three. Notice that the rational-root test is sufficient here
*because the polynomial has degree three*: no unexamined quadratic-factor
case remains.

## 2. Complete ramification ledger

The critical points and values of `phi` are

```text
phi'(z)=(z-1)(3z-1),
phi(1)=0,                    phi(1/3)=4/27.               (6)
```

The zero divisor of `h` consists of `x=0` with multiplicity three and the
simple point `x=-2Y/lambda`. The equation `h=4/27` is

```text
G(x)=27x^3(2Y+lambda*x)-8=0.                              (7)
```

It has four distinct roots over an algebraic closure of `K`. Indeed the
only possible nonzero common root with `G'` is
`x=-3Y/(2lambda)`, and

```text
G(-3Y/(2lambda))=-(729Y^4+128lambda^3)/(16lambda^3)!=0   (8)
```

in `k(Y)`. The normalized fibres and ramification indices are therefore:

At `x=0,z=0`, the `z`-derivative is a unit and `z=x^3` times a unit. At
`x=0,z=1`, writing `s=z-1` gives `s^2` times a unit equal to `x^3` times a
unit; coprimality of `2` and `3` gives one normalization branch with
orders `(ord(x),ord(s))=(2,3)`. At `x=-2Y/lambda,z=1`, the zero of `h` is
simple, so `x+2Y/lambda` is a unit times `s^2`. Finally, at a root `alpha`
of `(7)`, `h'(alpha)!=0` and
`phi(z)-4/27=-(z-1/3)^2+O((z-1/3)^3)`. These local equations give

```text
x=0:
  P_0,       z=0,       ord(x)=1, ord(z)=3,             e=1;
  P_1,       z=1,       ord(x)=2, ord(z-1)=3,           e=2;

x=-2Y/lambda:
  Q_0,       z=0,                                           e=1;
  Q_1,       z=1,       ord(x+2Y/lambda)=2, ord(z-1)=1,  e=2;

G(alpha)=0:
  R_alpha,   z=1/3,     ord(x-alpha)=2, ord(z-1/3)=1,    e=2;
  R'_alpha,  z=4/3,                                         e=1. (9)
```

There are four rows of the third type. At `x=infinity`, let `q=1/x` and
extend its valuation to a normalization place with ramification index `e`.
The right side of `(1)` has value `-4e`; hence `z` has negative value and

```text
3 ord(z)=-4e.                                             (10)
```

Thus `3` divides `e`. Since the whole cover has degree three, there is one
point `P_infinity`, with

```text
e(P_infinity)=3,          ord(x)=-3,        ord(z)=-4.    (11)
```

There is no omitted ramification: the discriminant of the cubic in `z` is

```text
Disc_z(phi(z)-h)=h(4-27h),                               (12)
```

and `(9)--(11)` account for both factors and infinity. The total tame
ramification contribution is

```text
1 + 1 + 4*1 + 2 = 8.                                    (13)
```

Riemann--Hurwitz for the degree-three map now gives

```text
2g(C)-2=3*(-2)+8=2,                 so g(C)=2.           (14)
```

## 3. The quadrature differential is holomorphic

The only possible finite poles of `(2)` occur at `z=1` or `z=1/3`; the
ramification parameters in `(9)` cancel them exactly. The complete order
calculation is

```text
place             ord(x)   ord(dx)   denominator order   ord(omega)
P_0                   1        0              0                1
P_1                   2        1              3                0
Q_1                   0        1              1                0
R_alpha               0        1              1                0
P_infinity           -3       -4             -8                1. (15)
```

For example, at `P_1` one may take a parameter `s` with
`x~s^2,z-1~s^3`; at `Q_1` and `R_alpha`, `x-x_0~s^2`; and at infinity
`x~s^-3,z~s^-4`. All other finite points have regular denominator, and
`dx` can vanish only at the ramification points already listed. Thus
`omega` has no pole. The two displayed simple zeros already have total
degree two, equal to `2g(C)-2`; hence there are no further zeros and `(3)`
is exact.

## 4. No finite algebraic cover creates a primitive when `lambda!=0`

After base change to `Kbar`, let `pi:C' -> C_Kbar` be finite, with `C'`
smooth and projective. Characteristic zero makes the extension separable.
If `Q` lies above `P` with ramification index `e_Q`, then

```text
ord_Q(pi^*omega)=e_Q ord_P(omega)+(e_Q-1)>=0.            (16)
```

The pullback is therefore again a nonzero holomorphic differential. If it
were `df` for some `f in Kbar(C')`, a pole of `f` of order `m>0` would give
a pole of `df` of order `m+1`; the leading coefficient `-m` cannot vanish
in characteristic zero. Hence `f` has no poles, so projectivity makes it a
constant. This contradicts `pi^*omega!=0`.

On the formal branch at `P_0`, fix `Y=y-(lambda/2)x`. The boundary equation
of THM-3979 becomes `(1)`, its volume unit is

```text
w=1/((z-1)(3z-1)),
```

and the defining quadrature satisfies `d(X^2)=omega` at fixed `Y`. If its
formal germ `X` were algebraic over `K(C)`, it would lie on some finite
normal cover, which becomes a smooth projective curve `C'` after completing
and base-changing constants. The difference
`d(X^2)-pi^*omega` is an algebraic differential. Its image in the chosen
completed local differential field is zero, and that map is injective;
hence the difference vanishes globally. This contradicts the preceding
paragraph and proves the transcendence assertion.

## 5. The rational endpoint has a residue obstruction

At `lambda=0`, equation `(1)` is rational. Put

```text
r=x/(z-1).
```

Direct substitution gives the birational parametrization

```text
x=r/(Yr^3-1),                 z=Yr^3/(Yr^3-1),           (17)
```

and `(2)` becomes

```text
omega=-2r dr/(Yr^3-1).                                  (18)
```

Over `Kbar`, the denominator has three distinct roots `rho` satisfying
`Y rho^3=1`. Each is a simple pole, with

```text
Res_(r=rho)(omega)=-2/(3Y rho)=-2rho^2/3!=0.             (19)
```

If `pi:C' -> P1_r` is any finite cover and a point above `rho` has
ramification index `e`, the residue of the pullback is `e` times `(19)`,
still nonzero in characteristic zero. Exact differentials have zero residue
at every place. Thus `(18)` is not exact on any finite algebraic cover,
and the same formal-germ argument as in Section 4 proves transcendence at
`lambda=0`.

This endpoint explains why the nonzero-slope genus-two proof cannot simply
be specialized: the holomorphic obstruction degenerates to a logarithmic
one, but nonintegrability survives.

## 6. Scope

This theorem isolates why the **canonical centered quadrature gauge** in
THM-3979 cannot algebraize or terminate: for nonzero slope it asks for a
primitive of a genus-two holomorphic differential, and at zero slope it
asks for a primitive with nonzero residues. It does not prove that every
formal cusp Darboux gauge is conjugate to this one, does not exclude a
different algebraic or polynomial gauge, and does not produce an
unrestricted Jacobian-conjecture obstruction.

**QED candidate.**

## Reproduction

```bash
python3 04-computation/jc2_centered_cusp_quadrature_genus2_thm3981.py
python3 -O 04-computation/jc2_centered_cusp_quadrature_genus2_thm3981.py
```
