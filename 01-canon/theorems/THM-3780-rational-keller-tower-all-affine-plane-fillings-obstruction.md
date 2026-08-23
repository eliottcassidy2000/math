---
id: THM-3780
title: "Rational Keller tower opposite-boundary unit and all-affine-plane-fillings obstruction"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  The
  THM-3774 tower factors through its normalization plane as the product of
  two reciprocal nonconstant Jacobians x and g=1/x.  Its common etale locus
  is the cylinder G_m x A^1.  For any normal affine model R of the same
  function field on which the same functions U,P are regular and define an
  etale map, the monic cover equation forces its root t into R and the
  cotangent identity forces g to be a unit.  Therefore no affine-plane
  model can make the same pair a polynomial Keller pair.  This is a
  same-pair source-filling obstruction, not JC(2).
source: jc_zero_debt_lift / all-source-fillings audit, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT by root.  The birational inverse and reciprocal
  Jacobian factors, finite normalization equality, branch parametrization,
  shared-cylinder equality, normality/integrality step, etale cotangent-basis
  argument, and nonconstant-unit obstruction were rederived independently.
  Normal, optimized, and frozen exact transcripts agree byte for byte.
depends_on:
  - THM-3774-three-component-rational-keller-cover-tower
related:
  - THM-3600-danielewski-arm-plane-atlas-singular-shear-and-no-filling
  - THM-3773-quadratic-rational-keller-cover-degree-and-target-word-nonpolynomialization
  - THM-3779-three-component-tower-maximal-danielewski-observable
script: 04-computation/jc2_rational_keller_tower_all_affine_plane_fillings_thm3780.py
output: 05-knowledge/results/jc2_rational_keller_tower_all_affine_plane_fillings_thm3780.out
script_sha256: ce55cbd60a58a2e20ab741218a2ae6619653fb530603319b5a7015a2f6d48193
output_sha256: 01b7da80ef343bbf73294c88ed1966e875004a92c191cb4f4e227d1bd2751235
semantic_sha256: 08c4b6f0f8d8874c6b20f28ca8a66684f8f59be73ef6cf4efdcf7852740dfd4a
hash_basis: raw LF bytes
---

# THM-3780 -- the reciprocal boundary unit forbids every affine-plane filling

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  The
all-degree rational Keller tower has a simpler structural description than
its expanded formulas suggest.  It is the composition of a birational chart
with Jacobian `x` and a polynomial fold with Jacobian `g=1/x`.  Their product
is one.  The pole and ramification defects are therefore the two ends of one
unit on a shared cylinder.

The description also closes a natural escape from THM-3774.  Changing the
affine-plane model of the source function field cannot make the same rational
pair polynomial Keller: every normal etale model must keep that nonconstant
unit invertible, whereas an affine plane has only constant units.

Let `k` be an algebraically closed field of characteristic zero, let `m>=1`,
and retain THM-3774's functions

```text
A=1+xy,
B=1+x^(2m+1)A^m,
U=xAB,
P=((m+1)B-1)/(mx)                                    (1)
```

inside `K=k(x,y)`.  Define

```text
t=x^2A,
g=P-[(m+1)/m]t^m.                                    (2)
```

Then

```text
g=1/x,
U=tP-t^(m+1)/m,                                      (3)
x=1/g,
y=tg^3-g.                                             (4)
```

Consequently

```text
K=k(t,P)=k(t,g),                                      (5)

k[x,x^(-1),y]
 =k[t,g,g^(-1)]
 =k[t,P,g^(-1)].                                      (6)
```

The common ring in `(6)` is a Laurent cylinder, the coordinate ring of
`G_m x A^1`.

Here is the all-fillings assertion.  Let `R` be a normal finitely generated
`k`-domain such that

```text
Frac(R)=K,                       U,P in R.             (7)
```

If the morphism

```text
(U,P):Spec R -> A^2                                  (8)
```

is etale, then

```text
t in R,                         g in R^*,             (9)
k[t,P,g^(-1)] subset R.                               (10)
```

In particular, no such `R` is isomorphic to a polynomial ring `k[X,Y]`.
Equivalently, no birational change of affine-plane source model can make the
same two rational functions `(U,P)` into a polynomial constant-Jacobian
pair.

## 1. The reciprocal-Jacobian factorization

The definition of `t` turns the two high-degree expressions in `(1)` into

```text
B=1+xt^m,
P=1/x+[(m+1)/m]t^m,
U=t/x+t^(m+1).                                        (11)
```

Equations `(3)` follow immediately.  Conversely, substituting `(4)` gives

```text
1+xy=tg^2,                       x^2(1+xy)=t,          (12)
P=g+[(m+1)/m]t^m,                                    (13)
```

which proves the field equality `(5)` and both directions of the localized
ring equality `(6)`.

There is a conceptual two-step factorization of the rational target map:

```text
A^2_(x,y)  -- beta -->  A^2_(t,P)  -- nu -->  A^2_(U,P),               (14)

beta(x,y)=(x^2(1+xy), 1/x+[(m+1)/m]t^m),
nu(t,P)=(tP-t^(m+1)/m,P).                             (15)
```

Both arrows are defined and mutually compatible on the cylinder in `(6)`.
Their Jacobians are

```text
J_(x,y)(t,P)=x,
J_(t,P)(U,P)=P-[(m+1)/m]t^m=g.                        (16)
```

Since `beta^*(g)=1/x`, the chain rule gives

```text
J_(x,y)(U,P)=xg=1.                                    (17)
```

Thus the all-degree cancellation of THM-3774 is not an accidental expansion:
it is a reciprocal factorization of one into two nonconstant Jacobians.  The
rationality of `beta` is exactly what permits the reciprocal `1/x`.

## 2. The alternative plane is the normalization and exposes the cusp fold

The relation in `(3)` is equivalent to the monic equation

```text
f_m(t)=t^(m+1)-mPt+mU=0.                              (18)
```

Because `(17)` makes `U,P` algebraically independent, `k[U,P]` is a
polynomial ring.  Equation `(18)` makes `k[t,P]` finite and integral over
it.  By `(5)` the fraction field is `K`, and `k[t,P]` is integrally closed.
It is therefore exactly the normalization of `k[U,P]` in `K`.

On this normalization plane the same pair is polynomial, but its Jacobian is
not a unit:

```text
J_(t,P)(U,P)=g=P-[(m+1)/m]t^m.                        (19)
```

Hence the etale locus of the finite normalization map `nu` is precisely
`D(g)`, the cylinder in `(6)`.  Its ramification curve `g=0` is smooth and
maps by

```text
P=[(m+1)/m]t^m,                 U=t^(m+1).             (20)
```

This is the normalization parametrization of THM-3774's branch binomial

```text
(m+1)^(m+1)U^m-m^(m+1)P^(m+1)=0.                     (21)
```

In particular, the `m=2` cusp is not merely an eliminant: it is the image of
the Jacobian-zero boundary that the original `(x,y)` filling placed at the
opposite end of `G_m`.

## 3. Normality forces the cover root and etaleness forces the unit

Now let `R` satisfy `(7)` and suppose `(8)` is etale.  Since `(18)` is monic
with coefficients in `k[U,P] subset R`, the element `t in K` is integral
over `R`.  Normality means that `R` is integrally closed in `K`, so

```text
t in R.                                                (22)
```

Equation `(2)` then gives `g in R`.

Differentiate the polynomial identity in `(3)` inside `Omega^1_(R/k)`:

```text
dU=t dP+g dt,
dU wedge dP=g dt wedge dP.                             (23)
```

Etaleness of `(8)` identifies `Omega^1_(R/k)` with the free module having
basis `dU,dP`.  In particular, `dU wedge dP` remains nonzero in the
two-dimensional cotangent space at every maximal ideal.  If `g` lay in a
maximal ideal `M`, reducing the second identity in `(23)` modulo `M` would
make its right side zero and its left side a basis wedge, a contradiction.
Thus `g` lies in no maximal ideal, so

```text
g in R^*.                                              (24)
```

Equations `(22),(24)` and `P=g+[(m+1)/m]t^m` prove `(10)`.  The unit is
genuinely nonconstant: `(5)` says that `t,g` are a transcendence basis, and
also `g=1/x` with nonconstant `x`.  But

```text
k[X,Y]^*=k^*.                                         (25)
```

Therefore `R` cannot be an affine-plane coordinate ring.  If one begins
directly with polynomials `U,P in k[X,Y]` having nonzero constant Jacobian,
the Jacobian criterion makes `(8)` etale and `k[X,Y]` is already normal, so
all hypotheses used above are automatic.  This proves the announced
same-pair source-polynomialization obstruction.

## 4. The two completions and the exact failure coordinate

The shared cylinder has two immediate affine-plane completions, and each
pays one side of the reciprocal invoice:

| model | added divisor | what extends | first failure |
|---|---|---|---|
| `k[x,y]` | `x=0` | the source volume `dx wedge dy` | `P` has pole `1/x` |
| `k[t,P]=k[t,g]` | `g=0` | both target functions `U,P` | `J(U,P)=g` vanishes |
| `k[t,P,g^(-1)]` | neither end | target functions and etaleness | nonconstant unit `g` remains |

The first row is THM-3774's rational Keller near miss.  The second is the
polynomial cusp fold in Section 2.  The third proves sharpness of the unit
conclusion: normal affine etale models with the same pair do exist, but they
must retain a nonconstant unit and therefore cannot be affine planes.

The boundary variables obey

```text
xg=1.                                                  (26)
```

Thus adding `x=0` sends `g` to a pole, while adding `g=0` sends `x` to a
pole.  The constant Jacobian in `(17)` exists only on the overlap where both
are units.  This is the exact geometric reason that merely choosing a more
favorable affine source chart cannot repair the tower.

## 5. Assumption boundaries and scope

1. **Normality is load-bearing only at `(22)`.**  It promotes the integral
   root `t` into the model.  For a proposed Keller filling it is not an
   optional restriction: an etale finite-type algebra over the regular
   target plane is regular and hence normal.
2. **Etaleness is load-bearing only at `(24)`.**  Dropping it gives the
   normalization plane `k[t,P]`, where the same functions are polynomial and
   `g=0` is exactly the critical divisor.
3. **Affine-plane units are the terminal obstruction.**  Dropping the
   constant-unit condition gives the sharp cylinder model in the last row of
   the table.
4. **The functions are fixed.**  The theorem excludes alternative source
   fillings for this same pair `(U,P)`.  It does not exclude replacing one or
   both target functions, leaving the field `k(U,P)`, deforming the tower, or
   constructing an unrelated polynomial Keller pair.
5. **No Jacobian-conjecture conclusion follows.**  The original source plane
   does not contain `P`, the normalization plane is ramified, and the common
   etale cylinder is not an affine plane.

The exact companion checks the symbolic all-`m` factorization, inverse,
normalization, differential coefficient, and branch parametrization, along
with eight direct source and opposite-completion controls.  The theorem's
new content is the normal-model quantifier in Section 3: the pole cannot be
hidden by a different affine-plane filling because its reciprocal is forced
to survive as a global unit.  **QED.**
