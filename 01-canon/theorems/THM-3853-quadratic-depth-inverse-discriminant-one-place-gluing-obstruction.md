---
id: THM-3853
title: "Quadratic-depth inverse-discriminant one-place gluing obstruction"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Two explicit
  irreducible discriminants have
  affine normalization A1, glue the four simple rays of the rational
  THM-3808 binary cubic at one affine point, and have one place at infinity.
  Nevertheless neither discriminant is realized by any homogeneous
  quadratic perturbation of all four binary-cubic coefficients with the
  linear packet fixed.  The exact saturated coefficient ideals are [1].
  This is a bounded inverse-discriminant obstruction, not a cubic-cover or
  planar-Jacobian theorem.
source: jc_quartic_c3_construct / inverse binary-cubic discriminant lane, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (root, 2026-08-23).  The audit checked the
  finite birational parametrizations, localization proof of irreducibility,
  four finite addresses and unique infinity place, and the GL2 endpoint
  criterion for unit representation.  It also inspected the full
  12-parameter coefficient universe, scalar elimination, nonzero
  saturation, and both hostile controls.  Normal and optimized exact
  Groebner replays byte-match the frozen 96-gate transcript and both hashes.
depends_on:
  - THM-3808-homogeneous-linear-binary-cubic-veronese-unit-trap
related:
  - THM-3801-cubic-etale-normalization-nonmonogenic-and-companion-sheet-gate
  - THM-3811-ramification-class-unit-criterion-and-nonlinear-cubic-packet
  - THM-3844-two-cusp-polynomial-branch-quadratic-resolvent-design-gate
  - THM-3847-two-place-cubic-deformation-monogenic-unit-debt
  - THM-3850-nonconstant-cubic-profile-irreducible-branch-puncture-formula
script: 04-computation/jc2_inverse_discriminant_quadratic_depth_thm3853.py
output: 05-knowledge/results/jc2_inverse_discriminant_quadratic_depth_thm3853.out
script_sha256: 2e5b4a240a4f81c793adc5cda1a3edc5da24056d5e49324e4617a1745be84bed
output_sha256: f74d5652d24fa1b881d6da4ca0ac3cdcc5da144b42d5011e8c14b27f20c9c0a7
semantic_sha256: 69a0700ffd7dde6653d1ebdc37d3fdefd70a60c176d41c1352b3783c8d566845
hash_basis: raw LF bytes
---

# THM-3853 -- the first one-place inverse discriminants do not terminate at quadratic depth

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  Work over an
algebraically closed field `k` of characteristic zero, and put `R=k[A,C]`.
Start from the rational
homogeneous-linear Delone--Faddeev packet of THM-3808,

```text
f_0(X,Y)=A X^3+C X^2Y+7A XY^2-3A Y^3,                         (1)

Delta_0=A(C+5A)(4C+19A)(3C-17A).                              (2)
```

For either

```text
L=C                         or                         L=A+C,  (3)
```

and any `lambda in k*`, the curve

```text
Gamma_(L,lambda)=V(delta_(L,lambda)),
delta_(L,lambda)=Delta_0+lambda L^5                             (4)
```

is irreducible.  Its affine normalization is `A1`; four distinct finite
normalization points map to the origin, and its projective normalization has
exactly one point over infinity.

Now perturb all four coefficients of `(1)` by arbitrary homogeneous
quadratics:

```text
a=A+q_a,       b=C+q_b,       c=7A+q_c,       d=-3A+q_d,
q_a,q_b,q_c,q_d in k[A,C]_2.                                  (5)
```

Let `Delta(a,b,c,d)` be the binary-cubic discriminant.  Then, for each of
the two lines in `(3)`,

```text
Delta(a,b,c,d) != Delta_0+lambda L^5
for every lambda in k*.                                        (6)
```

Thus the most direct nonmonogenic repair of the homogeneous cone does not
terminate at coefficient depth two.  Any polynomial realization of either
named inverse discriminant with the same linear packet must use a coefficient
term of degree at least three.  This does **not** exclude higher-degree or
formal lifts, other linear packets, or other one-place discriminants.

## 1. The exact monogenicity coordinate

For a finite free cubic order with trace-zero basis `(omega,theta)`, its
intrinsic binary index form is

```text
I(x,y)=det(1,x omega+y theta,(x omega+y theta)^2).               (7)
```

Over any commutative ring, if `I(x,y)` is a unit, then `(x,y)=R`, because
`I(x,y)` belongs to the ideal `(x,y)`.  Choose `r,s` with

```text
xs-yr=1.                                                        (8)
```

The corresponding `SL_2(R)` change sends the first endpoint coefficient of
the transformed index form to `I(x,y)`.  Conversely, a unit endpoint is
already a represented unit.  Hence

```text
I represents a unit
iff some GL_2(R) coordinate has a unit endpoint coefficient.    (9)
```

This identifies the exact old mechanism.  THM-3844 is globally a power
basis, so `(9)` has endpoint `1`; THM-3847 has

```text
I_beta=A X^3+3X^2Y-2CXY^2+2 beta Y^3,
I_beta(0,1)=2 beta in k*.                                       (10)
```

Replacing the scalar endpoint by a polynomial is therefore the cheapest
visible repair, but THM-3850 shows why that particular profile family pays
extra branch punctures.

In `(5)`, all four coefficients belong to the maximal ideal `(A,C)`.  Thus
for every `x,y in R`,

```text
I_(a,b,c,d)(x,y) in (A,C),                                     (11)
```

and the index form cannot represent a scalar unit.  Assertion `(11)` is
only a secure **index-form nonrepresentation** statement.  By itself it
does not prove that a corresponding cubic algebra is a domain, normal, has
constant units on an etale open, or supplies a plane atlas.

Here the Delone--Faddeev index form is the negative of the displayed binary
cubic; the sign does not change `(11)` or the unit-representation question.

## 2. Two exact one-place targets

Use complementary linear coordinates `(M,L)` and write

```text
Delta_0(M,L)=L^4 D_L(M/L).                                     (12)
```

For `L=C,M=A`, and for `L=A+C,M=A`, respectively, one obtains

```text
D_C(t)     =-t(5t+1)(17t-3)(19t+4),
D_(A+C)(t)=-t(4t+1)(15t+4)(20t-3).                             (13)
```

Each polynomial has four distinct roots.  On `(4)`, away from `L=0`, put

```text
t=M/L.                                                         (14)
```

Equation `(4)` is then equivalent to

```text
L=-D_L(t)/lambda,                   M=tL.                       (15)
```

This is a polynomial parametrization by every finite `t`.  Conversely,
`t` is integral over `k[L]` by the degree-four equation
`D_L(t)+lambda L=0`, and `(14)` identifies the function fields.  Therefore
`k[t]` is the finite normalization of the affine curve.

The localization of `delta_(L,lambda)` at `L` is, up to the unit `L^4`, the
linear polynomial `D_L(t)+lambda L`; it is irreducible.  Any additional
factor would be supported on `L=0`, but `L` does not divide `Delta_0` in
either orientation.  Hence `(4)` is irreducible.

At a root of `D_L`, formulas `(15)` give `L=M=0`; the four distinct roots in
`(13)` are the four normalization addresses glued at the origin.  At
`t=infinity`, `L` has degree four and `M` degree five.  Thus the projective
normalization has one and only one point over the line at infinity.  The
affine normalization is exactly `A1`, not a punctured rational curve.

The excluded coordinate `L=A` is a useful typing warning: `A` already
divides `Delta_0`, so `Delta_0+lambda A^5` remains reducible.  It is not a
third one-place orientation.

## 3. The saturated inverse-discriminant computation

Write

```text
q_a=u_0 A^2+u_1 AC+u_2 C^2,
q_b=u_3 A^2+u_4 AC+u_5 C^2,
q_c=u_6 A^2+u_7 AC+u_8 C^2,
q_d=u_9 A^2+u_10 AC+u_11 C^2.                                  (16)
```

The exact discriminant law is

```text
Delta=b^2c^2-4ac^3-4b^3d-27a^2d^2+18abcd.                      (17)
```

Its degree-four part is automatically `(2)`, and it has no terms outside
degrees four through eight.  For a fixed `L` in `(3)`, equate:

- every coefficient of degrees six, seven, and eight to zero; and
- the six degree-five coefficients to one common scalar multiple of `L^5`.

After eliminating the common scalar via the `C^5` row, there are exactly

```text
29 equations in u_0,...,u_11.                                  (18)
```

In both orientations the desired scalar is

```text
lambda=-4u_11.                                                  (19)
```

Let `J_L` be the ideal generated by `(18)` over `Q`, and introduce `eta` to
enforce `lambda!=0`.  Exact Groebner reduction gives

```text
Groebner(J_C+(eta lambda-1))     =[1],
Groebner(J_(A+C)+(eta lambda-1)) =[1].                           (20)
```

Because `(20)` is an identity over `Q`, it remains so after extension to
every characteristic-zero field.  This proves `(6)`.

Both hostile boundaries are active.  At `lambda=0`, the zero perturbation
survives and returns the reducible four-line discriminant `(2)`.  The
tempting perturbation `q_d=C^2` already supplies the desired coefficient
`lambda=-4`, but its complete debt is

```text
Delta-Delta_0
=-4C^5+126A^2C^3+162A^3C^2-27A^2C^4.                           (21)
```

The last three terms cannot be discarded; `(20)` says that no simultaneous
homogeneous-quadratic correction of the other coefficients cancels them
while retaining a nonzero target scalar.

## 4. Scope and next design grammar

The result is deliberately bounded.  It fixes the linear packet `(1)`,
allows every homogeneous quadratic correction in `(16)`, and closes only
the two explicit inverse targets `(3)`--`(4)`.  It does not classify all
quadratic packets or all one-place discriminants.  It also does not turn
the index obstruction `(11)` into normality, nonprincipal different, an
affine etale surface, or a Keller map.

The mechanism is nevertheless sharp.  The desired curve already has the
correct one-place geometry and the square-zero point naturally demanded by
secure nonmonogenicity.  What fails is finite termination: the first
quadratic coefficient that creates the fifth-degree gluing also creates
mixed fifth-degree and sixth-degree debt.  A separate formal deformation
analysis may solve that debt recursively; `(20)` says only that the solution
cannot stop after the pure quadratic layer.  The first open polynomial cell
therefore allows interacting coefficient terms of degree at least three.

The assertion-free companion freezes the two normalizations, the full
12-parameter universe, both 29-equation ideals, the two saturated bases,
and the zero and active controls.  Normal and optimized replays must
byte-match its frozen transcript.

Reproduction:

```bash
python3 04-computation/jc2_inverse_discriminant_quadratic_depth_thm3853.py
python3 -O 04-computation/jc2_inverse_discriminant_quadratic_depth_thm3853.py
```
