---
id: THM-3981
title: "Centered cusp quadrature has a genus-two transcendence obstruction"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. For every
  centering slope lambda, the generic-Y curve underlying the centered cusp
  quadrature carries a differential that remains nonexact on every finite
  algebraic cover. For lambda!=0 it is an integral genus-two curve and the
  differential is holomorphic with divisor equal to the finite zero-color
  point plus infinity. For lambda=0 the curve is rational and the
  differential has three simple poles with nonzero residues. Consequently
  the canonical centered formal X-coordinate is transcendental over the
  generic-slice function field and cannot terminate algebraically in that
  gauge. A separate lambda=1 supplement proves that this canonical
  obstruction survives every scalar Y-specialization, changing at the
  discriminant from holomorphic to logarithmic. This is not an obstruction
  to other formal gauges or to JC(2).
  A second, all-height supplement proves that the same lambda=1 centered
  quadrature is transcendental over the generic-Y function field for every
  determinantal height n>=2. Odd heights have nonzero seam residues; even
  heights at least four have a nonzero residue pairing with a global
  holomorphic differential. Exceptional scalar fibres are not included in
  that all-height statement.
source: jc-extra-debt-local / post-THM-3979 generic-slice audit, 2026-08-24
audit: >
  PASS (root / jc-cohn3709 and jc-zero-debt-lift, 2026-08-24). Independent
  audits checked geometric integrality, every finite and infinite branch of
  the degree-three cover, the Riemann--Hurwitz genus-two count, the exact
  divisor of omega, persistence of nonexactness under finite pullback, and
  the formal-germ injection. They also rederived the lambda=0 parametrization
  and its three nonzero residues. Normal, optimized, and frozen outputs
  byte-match at CHECKS=56; all hashes agree. The result remains deliberately
  scoped to the canonical centered gauge. A separate root/Huygens audit
  checked the lambda=1 hyperelliptic model and both scalar degeneration
  residues; its 25-gate companion matches in normal and optimized mode.
  A further independent hostile audit checked the all-height lambda=1
  extension: geometric irreducibility (including gcd(3,n+2)=3), the complete
  branch/genus and infinity ledgers, odd seam residues, all terms in the
  even-height residue pairing, and trace descent from the D-place formal
  germ. Its 248-gate companion matches the frozen transcript after LF
  normalization in normal and optimized mode; an independent n<=500 sweep
  passed 2,494 additional gates.
depends_on:
  - THM-3979-two-color-formal-cusp-darboux-lifting
related:
  - THM-3973-exact-volume-simple-cubic-determinantal-affine-plane-completion
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

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.** Let `k` be
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

## 7. Independent lambda-one hyperelliptic and scalar-slice supplement

The preceding proof varies the centering slope while keeping `Y` generic.
There is an orthogonal fixed-slope view which controls every scalar
specialization. Set `lambda=1` and put

```text
u=(z-1)/x,                    v=x+Y-u^3.                 (20)
```

Then the same function field has the hyperelliptic presentation

```text
v^2=(u^3-Y)^2+2u^2,
Disc_u=-512Y^2(729Y^4+128),                              (21)

omega=2du/v.                                             (22)
```

The two hyperelliptic infinities in `(21)` are the points `P_0` and
`P_infinity` in `(3)`, so the two divisor calculations agree. For every
scalar `Y_0` off the discriminant, `(21)` is a smooth genus-two curve and
the holomorphic argument still applies.

At `Y_0=0`, write `v=uV`. The normalization is

```text
V^2=u^4+2,
du/v=du/(uV).                                            (23)
```

The two residues over `u=0` have squares `1/2`. On the other discriminant
locus the unique repeated root is `a=4/(9Y_0)` and

```text
R''(a)/2=-4.                                             (24)
```

Writing `v=(u-a)V` gives two residues whose squares are `-1/4`. Hence the
form remains nonexact on both genus-one normalizations, and trace descent
again forbids a primitive after finite algebraic extension. Thus the
canonical `lambda=1` quadrature is nonalgebraic on **every** scalar
`Y`-slice. This scalar degeneration is not the `lambda=0` rational endpoint
of Section 5; the two axes must not be conflated.

The independent companion
`04-computation/jc2_centered_cusp_quadrature_genus_two_thm3981.py` has 25
optimization-safe gates. Its script/output/semantic hashes are
`9ba76e156c092410e755d2680e5c7413721d729957802bd03da1dc7ef20ff3b8`,
`3f29a9da90581085fae01be3fe8887bdf294f5a1100f247b0278c179ac765b32`,
and `f88ec164728e05717827b54c573afed10668a4c4c448aae3f27f0bd8bcbb15b7`.

## 8. All-height centered-gauge supplement

The preceding theorem concerns height two. There is also a uniform
generic-`Y` statement in the all-height determinantal tower of THM-3973.
For `n>=2`, put

```text
z=1+x^n t,       p=zt,       y=x^(n-1)zt^2,
B_n=k[x,z,p,y],                 Y=y-x/2.                  (25)
```

On the completed color `D=(x,z,p)`, let `z=zeta_n(x,y)` be the branch
with `zeta_n(0,y)=0`. Its volume unit and centered quadrature are

```text
w_n=1/((zeta_n-1)(3zeta_n-1)),
X_n^2=2 integral_0^x s w_n(s,Y+s/2) ds.                 (26)
```

Thus `X_n=x+O(x^2)` and, at fixed `Y`,

```text
d(X_n^2)=2 omega_n,
omega_n=x dx/((z-1)(3z-1)).                             (27)
```

For `6L^2=1`, this supplies the centered cusp pair

```text
A_D=Y^2+2LX_n,                 C_D=Y^3+3LX_nY
```

on the D color. Together with the retained-arm pair
`(A_L,C_L)=(t+2Lx,-x)`, the Hensel idempotents glue it to an exact formal
Darboux pair. Thus `(26)` is an actual all-height formal gauge, not merely
an auxiliary Abelian integral.

The ambient surface function field is the function field over `K=k(Y)`
of the normalization `C_n` of

```text
z(z-1)^2=x^(n+1)(Y+x/2).                                (28)
```

This curve is geometrically integral. Indeed, if the monic cubic in `z`
were reducible over an algebraic closure of `K(x)`, it would have a root
`r in Kbar(x)`. Integrality over `Kbar[x]` makes `r` a polynomial. Degree
at infinity already rules this out unless `3` divides `n+2`. In the
remaining case, the simple zero of the right side at `x=-2Y` forces
`r(-2Y)=0` with multiplicity one: `r=1` would give even multiplicity.
Every zero of `r` or `r-1` lies in `{0,-2Y}`. Since the nonconstant
polynomial `r-1` has a zero, it must vanish at zero; hence `r` has only the
single simple zero `-2Y` and has degree one. Then `3 deg(r)=n+2` forces
`n=1`, contrary to `n>=2`.

Put

```text
d=gcd(3,n+2),       epsilon_n=1 if n is even, 0 otherwise.
```

The degree-three map `C_n -> P1_x` has ramification contribution
`epsilon_n` above `(x,z)=(0,1)`, contribution one above
`(-2Y,1)`, contribution `n+2` at the simple roots of
`x^(n+1)(Y+x/2)=4/27`, and contribution `3-d` at infinity. Therefore

```text
g(C_n)=(n+epsilon_n-d+2)/2,                              (29)
ord_infinity(omega_n)=(2n-2)/d-1>=0.                    (30)
```

Some scalar specializations make this branch ledger degenerate; the
statement here remains over the generic field `k(Y)`. The singular seam in
the affine plane model is
present for every `Y`; smoothness always refers to the normalization.

For odd `n=2q-1`, `q>=2`, extend constants by `a^2=Y`. At the two seam
points `P_epsilon` above `(0,1)`, set `s=z-1` and

```text
s sqrt(1+s)=epsilon a x^q(1+x/(2Y))^(1/2).
```

Since `2 sqrt(1+s)/(2+3s)=1+O(s)` and `s=O(x^q)`, every omitted term in
`omega_n` is regular. Its exact residue is

```text
Res_(P_epsilon)(omega_n)
 =1/(2 epsilon a) binom(-1/2,q-2)(2Y)^(-(q-2)) !=0.     (31)
```

Thus no rational primitive exists at any odd height.

Now take even `n=2q`. For `q=1`, Sections 1--4 already show that
`omega_2` is a nonzero holomorphic genus-two differential. If `q>=2`, the
only pole of `omega_(2q)` is a residue-free second-kind pole at the seam.
The global differential

```text
eta_q=x^q dx/((z-1)(3z-1))                              (32)
```

is holomorphic; its infinity order is `(q+1)/d-1>=0`. With `x=tau^2`
and `a^2=Y`, write `h` for the local primitive of `omega_(2q)`. Direct
expansion gives

```text
omega_(2q)=a^-1 tau^(2-2q)(1+tau^2/(2Y))^(-1/2)d tau
             + terms starting at tau^3 d tau,
eta_q      =a^-1(1+tau^2/(2Y))^(-1/2)d tau
             + terms starting at tau^(2q+1)d tau.       (33)
```

The first omitted contribution to `h eta_q` starts at `tau^4 d tau`, so
the following binomial convolution is the complete residue, not merely a
leading approximation:

```text
Res_P(h eta_q)
 =(-1)^(q-1)(q-2)!/((2q-3)!! Y^(q-1)) !=0.              (34)
```

If `omega_(2q)=df`, then `f eta_q` can have a pole only at `P`, while
`(34)` says its residue there is nonzero, contradicting the global residue
theorem. Hence `omega_n` is nonexact for every `n>=2`.

Finally use the D-place embedding `F=K(C_n) -> K((x))`. If the formal germ
`h=X_n^2/2` were algebraic over `F`, then `E=F(h)` would be a finite
subextension of the completed algebraic closure. In characteristic zero,
injectivity of differentials upgrades the local identity `dh=omega_n` to
an identity in `Omega_(E/K)`. Taking field trace gives

```text
d Tr_(E/F)(h)=[E:F] omega_n,                             (35)
```

contradicting nonexactness. Thus both `X_n^2` and `X_n` are transcendental
over `Frac(B_n)` for every height. This proves nontermination only for the
specific centered gauge `(26)`. It is distinct from the rational split
`K times K` gauge of THM-3980, does not cover exceptional scalar fibres,
and leaves alternative formal gauges, global Keller pairs, and `JC(2)`
open.

The independent companion
`04-computation/jc2_all_height_centered_cusp_quadrature_thm3981_independent.py`
has 248 optimization-safe gates. Its script and frozen-output hashes are
`d790684bbf4a1a00346f3e97bacd33bc5f42bc07252dabd3526d55324d9c8d3f`
and `60e053e4bd8d823af0133a9e80b30309c22329ed1c5bd46bc8e66d4c0a68738e`.
The frozen LF transcript matches normal and optimized runs after line-ending
normalization. Its first output line retains the pre-promotion scratch label
`THM-3980`; the canonical routing and filename are this THM-3981 supplement,
which is a different gauge from actual THM-3980.

**QED.**

## Reproduction

```bash
python3 04-computation/jc2_centered_cusp_quadrature_genus2_thm3981.py
python3 -O 04-computation/jc2_centered_cusp_quadrature_genus2_thm3981.py
python3 04-computation/jc2_all_height_centered_cusp_quadrature_thm3981_independent.py
python3 -O 04-computation/jc2_all_height_centered_cusp_quadrature_thm3981_independent.py
```
