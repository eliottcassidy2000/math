---
id: THM-3965
title: "Unit-ideal two-parameter cubic deformations retain multiple infinity places"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. Every constant
  two-parameter deformation
  A U^3+(C+gA)U^2V+(AC-1)UV^2+hAV^3 of the THM-3907 unit-ideal binary
  cubic is either generically reducible or has no one-place discriminant.
  For h!=0, the repeated-root incidence is a quadratic in A whose leading
  coefficient 2t^3+gt^2-h has at least two distinct roots; each produces a
  distinct normalization place over [1:0:0]. The incidence has one reducible
  seam, (g,h)=(0,-4/3), and its two target components have respectively two
  and four infinity places. Thus even rational genus-zero degenerations do
  not repair THM-3907's infinity debt. Bivariate coefficient deformations,
  other binary-cubic grammars, Keller realization, and JC(2) remain open.
source: jc-zero-debt-lift / post-THM-3907 positive-deformation lane, 2026-08-24
audit: >
  INDEPENDENT MATHEMATICAL AND HOSTILE AUDIT PASS. The audit rederived the
  domain and scalar-unit gates, repeated-root incidence, unique reducible
  seam, height-one maximality and S3 packet, and exact two/four infinity-place
  ledgers on the seam. The weak symbolic-factor presence test was hardened to
  an exact divisibility identity. Normal and optimized runs match the frozen
  29-gate output after canonical LF normalization on Windows; raw and semantic
  hashes and documentation checks pass.
depends_on:
  - THM-3801-cubic-etale-normalization-nonmonogenic-and-companion-sheet-gate
  - THM-3907-unit-ideal-nonmonogenic-cubic-six-place-boundary
related:
  - THM-3920-affine-plane-boundary-unibranch-depressed-cubic-chart-obstruction
  - THM-3961-arbitrary-q-hidden-repetition-normality-and-conductor-debt
script: 04-computation/jc2_unit_ideal_two_parameter_infinity_places_thm3965.py
output: 05-knowledge/results/jc2_unit_ideal_two_parameter_infinity_places_thm3965.out
script_sha256: 942b28312b093f84efcf8acb3864953089c104303404536d740e8917c7cffb5c
output_sha256: 3421b597d3eb40e702e4687d92da952d2c0555c5d87acdb64bcc574037f5b55e
semantic_sha256: 3ec72321d50110d452b2b0c88a0097ca9cf0397f7591241525b6385a637c6a7a
hash_basis: raw LF bytes
---

# THM-3965 -- the first unit-ideal deformation cannot merge its infinity roots

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.** Work over an
algebraically closed field `k` of characteristic zero. For constants
`g,h in k`, consider the binary cubic

```text
Phi_(g,h)=A U^3+(C+gA)U^2V+(AC-1)UV^2+hA V^3.           (1)
```

This is the smallest constant-coefficient deformation of the unit-ideal
nonmonogenic cubic in THM-3907 that preserves its scalar-unit hostile fibre.
The complete result is:

1. if `h=0`, the generic cubic is reducible;
2. if `h!=0` and `(g,h)!=(0,-4/3)`, the associated Delone--Faddeev order is
   a normal nonmonogenic `S3` cubic order with irreducible reduced
   discriminant, but the normalization of that discriminant has at least two
   distinct places above the projective point `[1:0:0]`; and
3. at `(g,h)=(0,-4/3)`, the discriminant splits into two reduced rational
   components, whose normalizations have respectively two and four infinity
   places.

Consequently no member of `(1)` supplies an irreducible one-place
discriminant. The statement is about this exact two-parameter family; it is
not a no-atlas theorem for arbitrary binary cubics.

## 1. The unit-ideal, domain, and nonmonogenic gates

The coefficient ideal is the whole target ring because

```text
AC-(AC-1)=1.                                             (2)
```

Nevertheless `Phi_(g,h)` represents no scalar unit. If
`Phi_(g,h)(X,Y)=lambda in k*` for `X,Y in k[A,C]`, restriction to `A=0`
would give

```text
X_0 Y_0(CX_0-Y_0)=lambda                    in k[C].     (3)
```

All three factors would be units. Thus `X_0,Y_0` would be nonzero constants,
while `CX_0-Y_0` would be nonconstant, a contradiction. By the index-form
criterion of THM-3801, every domain member is globally nonmonogenic.

If `h=0`, the dehomogenized cubic has the factor `t=U/V`, so this row fails
the cubic-domain gate. Assume henceforth that `h!=0`. Dehomogenizing by
`V=1` gives

```text
psi=A t^3+(C+gA)t^2+(AC-1)t+hA
   =C t(t+A)+(At^3+gAt^2-t+hA).                         (4)
```

Over `k(A)`, a factorization of this polynomial, which is linear in `C`,
would have a `C`-independent factor dividing both rows in `(4)`. The roots
`t=0,-A` of the first row give values

```text
hA,                       A(-A^3+gA^2+1+h),             (5)
```

and neither is identically zero. The rows are coprime, so `(4)` is
irreducible over `k(A,C)`.

## 2. Exact repeated-root incidence

The finite repeated-root equations are `psi=psi_t=0`. Eliminating the linear
color `C` gives

```text
G(A,t)=alpha(t)A^2+beta(t)A+t^2=0,                       (6)
alpha=2t^3+gt^2-h,                    beta=t^4-2ht.
```

Generically the color is recovered from

```text
(A+2t)C+2Agt+3At^2-1=0.                                (7)
```

There is no lost incidence component where the two `C` coefficients vanish.
Indeed `t(t+A)=A+2t=0` forces `(A,t)=(0,0)`, and there
`psi_t=-1`; this is only the familiar isolated raw-resultant artifact.

The quadratic discriminant is

```text
Disc_A(G)=t^2 D(t),
D=t^6-4(h+2)t^3-4gt^2+4h(h+1).                          (8)
```

Because `alpha(0)=-h!=0`, the coefficients of `G` are primitive over
`k[t]`. Thus `G` is reducible exactly when `D` is a square in `k(t)`. Since
`D` is monic, a rational square is a polynomial square. Write

```text
D=(t^3+a t^2+b t+c)^2.                                  (9)
```

The coefficients of `t^5,t^4,t^3,t^2,t^0` successively force

```text
a=b=0,       c=-2(h+2),       g=0,
(h+2)^2=h(h+1).                                         (10)
```

Hence

```text
G reducible  iff  (g,h)=(0,-4/3),                       (11)
D=(t^3-4/3)^2.
```

Away from `(11)`, the genuine incidence is irreducible. Its map to the
target discriminant is generically one-to-one: a nonzero cubic cannot have
two distinct repeated roots. A repeated root at projective root-infinity
would require `A=C=0`, only an isolated target point. Therefore the target
discriminant has irreducible support. Its homogeneous degree-seven part and
its linear row are

```text
Delta_7=-4A^4C^3,                    [A]Delta=4.          (12)
```

The linear term makes the irreducible discriminant reduced. The associated
finite-free cubic order is consequently maximal in every height-one
localization: an order index changes discriminant valuation by an even
integer. Finite freeness gives `S2`, so the order is normal. Its irreducible
cubic and nonsquare discriminant give generic group `S3`.

## 3. Why alpha counts normalization places, not just formal roots

Compactify `(6)` only in the `A` direction, writing `[X:Z]=[A:1]`. Its
equation is

```text
alpha(t)X^2+beta(t)XZ+t^2Z^2=0.                         (13)
```

For every distinct root `rho` of `alpha`, `(13)` contains the point

```text
(t,[X:Z])=(rho,[1:0]).                                  (14)
```

These are distinct points before normalization, so each has at least one
distinct normalization point above it. They are genuine target-infinity
places. Indeed in the local coordinate `z=1/A`, equation `(7)` reads

```text
C=-(2gt+3t^2-z)/(1+2tz),                                (15)
```

which is regular at `z=0`. The induced target point is therefore

```text
[A:C:1]=[1:zC:z] -> [1:0:0].                            (16)
```

Because the incidence and discriminant have the same function field, the
distinct points in `(14)` are distinct places of the projective
discriminant normalization.

It remains only to count the support of `alpha`. A cubic with one distinct
root would have to satisfy

```text
alpha=2(t-r)^3.                                         (17)
```

The missing `t` coefficient gives `r=0`, and then the constant term gives
`h=0`, contrary to the domain hypothesis. More explicitly,

```text
Disc_t(alpha)=4h(g^3-27h).                              (18)
```

If `alpha` is not squarefree, then `h=g^3/27`, necessarily `g!=0`, and

```text
alpha=2(t+g/3)^2(t-g/6),                                (19)
```

which still has two distinct roots. Thus every irreducible row has at least
two distinct normalization places over `[1:0:0]`. This is the exact
one-place obstruction; it remains valid on the repeated-`alpha` boundary.

## 4. The unique reducible seam also fails componentwise

At `(g,h)=(0,-4/3)`, both the incidence and target discriminant factor:

```text
G=(2A+t)(3At^3+2A+3t)/3,                               (20)

Delta=-(12A^3-3AC-1)(AC^3+12A+3C^2)/3.                 (21)
```

Both target factors are irreducible by primitive linearity in `C` or `A`,
and both occur to exponent one.

The first component has normalization parameter `A`:

```text
C=4A^2-1/(3A).                                          (22)
```

It has two infinity places, at `A=0` and `A=infinity`. The second has
normalization parameter `t`:

```text
A=-3t/(3t^3+2),                    C=2/t.                (23)
```

It has one infinity place at `t=0` and three more at the distinct roots of
`3t^3+2`, hence four in total. Therefore splitting the discriminant does not
produce a hidden one-place component.

## 5. Hostile degenerations and the design lesson

The genus-zero degeneration

```text
(g,h)=(-3/4,-1)
```

has

```text
D=t^2(t-1)^2(t^2+2t+3).                                 (24)
```

Thus the incidence normalization is rational: after removing the square,
it is a conic. Nevertheless

```text
Disc_t(2t^3-3t^2/4+1)=-1701/16,                         (25)
```

so the leading coefficient still supplies three distinct infinity places.
This is a sharp warning that genus zero does not pay the place invoice.

At `(g,h)=(3,1)`, the repeated-leading-coefficient boundary is realized:

```text
alpha=(t+1)^2(2t-1).                                    (26)
```

It supplies exactly two distinct roots to the mechanism, showing where the
lower bound in Section 3 is attained at the level of leading support. The
theorem does not claim that these are the only infinity places of that curve.

The constant two-parameter deformation `(1)` is therefore exhausted. A
viable continuation must let a coefficient vary genuinely in both target
coordinates, change the linear-color plane, or use another nonmonogenic
cubic grammar. No Keller atlas or counterexample is claimed here, and JC(2)
remains **OPEN**. **QED.**

## Reproduction

```bash
python3 04-computation/jc2_unit_ideal_two_parameter_infinity_places_thm3965.py
python3 -O 04-computation/jc2_unit_ideal_two_parameter_infinity_places_thm3965.py
sha256sum 04-computation/jc2_unit_ideal_two_parameter_infinity_places_thm3965.py \
  05-knowledge/results/jc2_unit_ideal_two_parameter_infinity_places_thm3965.out
python3 agents/check_docs.py
```
