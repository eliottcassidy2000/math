---
id: THM-3906
title: "Degree-six common-zero normal cubic two-place boundary"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  The binary cubic
  (C^2+A)U^3+A U^2V+(A+C)UV^2+C V^3 defines a normal finite-flat,
  nonmonogenic S3 cubic order with geometrically irreducible rational sextic
  discriminant.  Its branch has an ordinary four-address common zero and two
  A2 cusps.  It has one projective infinity point but exactly two
  normalization places, of orders two and four.  The containing minimal
  triple-root sextic grammar is always reducible or has the same two-place
  boundary.  Independently, a rational one-place sextic with an ordinary
  quadruple point has at most two finite cusps; the bound is sharp as plane
  geometry.  THM-3911 subsequently obstructs a normal finite-flat,
  generic-irreducible-S3 realization of that sharp control with generic
  discriminant exponent one; THM-3908 excludes the
  coefficient-depth-at-most-two pure-sixth top-form lane.  THM-3913 first
  escapes at degree ten and
  depth three, but its normalization is elliptic.  Other sextic grammars,
  nonnormal orders, Keller realization, and JC(2) remain open.
source: root / post-THM-3890 degree-six escape and cusp-cap session, 2026-08-23
audit: >
  INDEPENDENTLY HOSTILE-AUDITED on 2026-08-23 by three disjoint audits.  One
  audit rederived the discriminant, generic irreducibility, Delone--Faddeev
  normality, projective normalization, singularities, and two Newton places.
  It caught and repaired an initial incomplete genus ledger: the origin
  contributes six, infinity contributes two, and two finite A2 cusps
  contribute one each.  A second audit proved the six-parameter grammar
  dichotomy.  Two independent proofs exclude the tempting four-A2 one-place
  packet, one by derivative integration and one by the degree-two radial
  Wronskian; the latter sharpens the cap to two cusps and supplies a sharp
  curve control.  The assertion-free companion freezes all identities,
  normalization data, singular support, Newton edges, family coefficients,
  cusp cap, full family Newton support, and sharp birational control in 91
  active gates.  Normal and optimized runs
  byte-match the frozen output.
depends_on:
  - THM-3801-cubic-etale-normalization-nonmonogenic-and-companion-sheet-gate
related:
  - THM-3844-two-cusp-polynomial-branch-quadratic-resolvent-design-gate
  - THM-3874-three-cusp-quadratic-k3-affine-class-group
  - THM-3879-rational-torus-sextic-c3-packet-one-place-tradeoff
  - THM-3887-binary-cubic-common-zero-quintic-one-tangent-obstruction
  - THM-3889-maximally-confluent-quadratic-binary-cubic-two-place-obstruction
  - THM-3890-universal-quintic-common-zero-resolvent-class-group-dichotomy
  - THM-3908-quadratic-depth-common-zero-one-point-sextic-two-place-obstruction
  - THM-3911-sharp-one-place-sextic-resolvent-three-torsion-obstruction
  - THM-3913-moving-triple-root-one-place-decic-normal-nonmonogenic-cubic
script: 04-computation/jc2_degree_six_common_zero_normal_cubic_two_place_thm3906.py
output: 05-knowledge/results/jc2_degree_six_common_zero_normal_cubic_two_place_thm3906.out
script_sha256: 3558670eb81467e10f924c68f2222d0ba82b036fb98a68fcb30a76bd2f92036a
output_sha256: a559118dfc582684f0319411cee54f963f5c2139d267e53c325936ce4959d47a
semantic_sha256: 6c584836bf42fa095c6a8b14da1aca1b70d30712bfb3ac152dbd97b3694a9757
hash_basis: raw LF bytes
---

# THM-3906 -- degree six realizes the cubic order but pays two places

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

Let `k` be an algebraically closed field of characteristic zero and put
`R=k[A,C]`.  Consider

```text
Phi=(C^2+A)U^3+A U^2V+(A+C)UV^2+C V^3.                  (1)
```

The associated Delone--Faddeev cubic algebra is a normal finite-flat domain
of rank three, is not monogenic over `R`, and has generic Galois group `S3`.
Its discriminant is the geometrically irreducible sextic

```text
Delta=-3A^4+4A^3C-20A^2C^2-4AC^3
      -4A^3C^2+6A^2C^3-48AC^4-4C^5-27C^6.              (2)
```

The projective branch curve is rational.  It has an ordinary quadruple
point at the common coefficient zero `(0,0)`, two finite ordinary cusps, and
exactly one projective point at infinity but two normalization places above
that point.  Thus degree six escapes the degree-five class-group obstruction
and realizes the cubic-order/nonmonogenic layer.  It fails precisely at the
one-place infinity layer.

## 1. The cubic order is irreducible, nonmonogenic, and normal

All four coefficients in `(1)` lie in `(A,C)`, and the coefficients `A,C`
themselves occur.  Hence their common ideal is exactly `(A,C)`.  In the
index-form convention of THM-3801, every index determinant is `+/-Phi(x,y)`
and therefore lies in `(A,C)`.  It cannot be a scalar unit, so the cubic
algebra is not monogenic.

Dehomogenize with `V=1` and write the cubic-root variable as `t`.  Then

```text
Phi(t,1)=A t(t^2+t+1)+C(Ct^3+t+1).                       (3)
```

The two displayed polynomials in `t` are coprime in `k[C,t]`; their
resultant is `C^3(C^2+C+1)`, which is nonzero.  Thus `(3)`, viewed as a
primitive polynomial linear in `A`, is irreducible.  Gauss's lemma makes the
generic cubic irreducible over `k(A,C)`.

Section 2 proves that `(2)` is irreducible.  It has valuation one at its
height-one prime, so it is nonsquare in `k(A,C)`.  The generic cubic
therefore has Galois group `S3`.

Let `T` be the Delone--Faddeev algebra.  It is finite free over the regular
surface `R`, hence `S2`.  Away from `(Delta)` it is etale.  At the DVR
`R_(Delta)`, the discriminant exponent is one.  The order/maximal-order
discriminant formula has an even index correction, so the index is zero and
`T` is the maximal DVR order.  Every height-one localization is normal;
therefore `T` is `R1`, hence normal.  Generic irreducibility makes it a
domain.

## 2. The repeated root gives the complete normalization

At a generic discriminant point, let `t` be the unique double root of
`Phi(t,1)`.  Solving `Phi=partial_t Phi=0` gives

```text
Q(t)=2t^3+4t^2+2t+1,
C=Q(t)/(t^3(t+2)),
A=-(2t+3)Q(t)/(t^4(t+2)^2).                              (4)
```

For homogeneous parameters `[T:S]`, put

```text
Qh=2T^3+4T^2S+2TS^2+S^3.
```

Formula `(4)` extends to

```text
[T:S] |-> [ -S^2(2T+3S)Qh : ST(T+2S)Qh : T^4(T+2S)^2 ]. (5)
```

All three coordinates have degree six.  They have no common zero: `T=0`,
`T=-2S`, and `S=0` are checked directly, while a root of `Qh` has nonzero
third coordinate.  Thus `(5)` is a degree-six morphism.  The generic double
root recovers `t`, so it is birational onto its image.  Its image has degree
six and lies on `(2)`; consequently `(2)` is the irreducible equation of the
image and `(5)` is its normalization.

This also supplies a geometric proof of absolute irreducibility.  As a
hostile arithmetic check, `Delta(A,1)` reduces modulo seven to a scalar
times the irreducible polynomial `A^4+A+1`, and the rational curve point
`(0,-4/27)` is smooth.

## 3. The exact singularity and genus ledger

The degree-four tangent cone of `(2)` at the origin is

```text
-A(3A^3-4A^2C+20AC^2+4C^3).                             (6)
```

The cubic factor has discriminant `-109744=-16*19^3`, so `(6)` consists of
four distinct lines.  The origin is an ordinary quadruple point with

```text
delta_origin=binomial(4,2)=6.                             (7)
```

Its normalization addresses are `S=0` and the three roots of `Qh`; the
affine cubic `Q` has discriminant `-76`, so these addresses are distinct.

Differentiating `(4)` gives

```text
dA/dt= 2(t^2+3t+3)(4t^3+9t^2+7t+4)/(t^5(t+2)^3),
dC/dt=-2(t^2+t+1)(t^2+3t+3)/(t^4(t+2)^2).                (8)
```

The two roots of `t^2+3t+3` are distinct and give the only other affine
singularities.  In target coordinates they satisfy

```text
27C^2+27C+19=0,                 21A-24C-19=0.             (9)
```

An exact saturation of `(Delta,partial_A Delta,partial_C Delta)` by `C`
gives `(9)` with square scheme structure; on `C=0` the only singular point
is the origin.  The tangent Hessian has rank one at both points.  After the
infinity contribution below, the genus formula leaves total delta two for
these two singularities, so each has delta one.  The rank-one tangent cone
makes them ordinary `A2` cusps rather than nodes.

At infinity the leading form is `-27C^6`, so the projective support is the
single point `[1:0:0]`.  In its chart

```text
x=C/A,                         z=Z/A,
```

the local equation is

```text
h=-27x^6
  +z(-4x^5-48x^4+6x^3-4x^2)
  +z^2(-4x^3-20x^2+4x-3).                                (10)
```

Its `z`-discriminant is

```text
Disc_z(h)=16x^4(x^2-x+1)^3.                              (11)
```

The two primitive Newton edges give

```text
z=-(4/3)x^2+...,                 z=-(27/4)x^4+....         (12)
```

They are two smooth branches, corresponding in `(5)` to `t=-2,0`.  Their
intersection multiplicity, hence the local delta invariant of their union,
is two.  The complete genus formula is therefore

```text
p_a(sextic)=10=6+1+1+2,                 genus(normalization)=0. (13)
```

This corrects the tempting but false inference that the ordinary quadruple
point alone exhausts the sextic genus.

## 4. The minimal triple-root coefficient grammar is always two-place

Consider the six-parameter family

```text
a=C^2+alpha A+alpha1 C,
b=beta A+beta1 C,
c=gamma A+gamma1 C,
d=C.                                                        (14)
```

At `C=0` its discriminant is

```text
gamma^2(beta^2-4alpha gamma) A^4.                          (15)
```

If the coefficient in `(15)` vanishes, then `C` divides the full
discriminant, so the branch is reducible.  Otherwise, after total
homogenization, the local infinity equation has the initial packet

```text
gamma^2(beta^2-4alpha gamma) z^2
  -4gamma^3 x^2 z-27x^6+later terms.                       (16)
```

Its Newton polygon has the two primitive edges

```text
(0,2)--(2,1)--(6,0),                                      (17)
```

and therefore two places of orders two and four.  Thus every member of
`(14)` is reducible or fails the one-place condition.  This closes the
minimal triple-root sextic grammar, not every sextic leading row.

## 5. The finite-cusp cap for a one-place sextic

There is a separate geometry-only obstruction.  Let `Gamma` be an
irreducible rational plane sextic with an ordinary quadruple point and a
unique normalization place at infinity whose image on `Gamma` is smooth.
Choose the normalization parameter so that the infinity place is
`t=infinity`.  After a target-linear change fixing the quadruple point, the
affine normalization has polynomial coordinates

```text
X(t),Y(t),                 deg X=5,       deg Y=6.          (18)
```

Smoothness at infinity forces the degree five in `(18)`: in the `Y`-chart,
`X/Y` must have valuation one.

Let the four ordinary branches above the quadruple point have distinct
parameters `r_i`.  Their squarefree product `G` is the full common factor of
`X,Y`, so

```text
X=G P_1,                         Y=G Q_2,                  (19)
```

with `deg P_1<=1`, `deg Q_2<=2`.  The map

```text
phi=[P_1:Q_2]:P1 -> P1                                  (20)
```

is nonconstant because the four branch tangents are distinct, and it has
degree at most two.  At any finite cusp parameter away from the quadruple
point, `X'=Y'=0`; equation `(19)` then makes

```text
P_1 Q_2'-P_1'Q_2=0.                                     (21)
```

Thus every cusp ramifies `(20)`.  Riemann--Hurwitz gives total ramification
at most `2*2-2=2`.  Consequently

```text
number of distinct finite cusps <=2.                      (22)
```

In particular an ordinary quadruple point plus four `A2` cusps cannot admit
any one-place line.  Those singularities force genus zero.  If the line were
supported at the ordinary quadruple, its pullback would contain four
normalization points.  If it were supported at an `A2` cusp, its contact
order would be two or three, not the required degree six.  Its support must
therefore be smooth, reducing to the case above.

The cap `(22)` is sharp as plane geometry.  Put

```text
G=t^4+t^3-2t-5,
X=tG,                         Y=(t^2+1)G.                  (23)
```

The four roots of `G` give an ordinary quadruple point, infinity is unique
and smooth, and

```text
gcd(X',Y')=(t-1)(t+1).                                   (24)
```

The two critical parameters are `A2` cusps; the relevant second/third
derivative determinants are `1680` and `432`.  This control is not claimed
to be a binary-cubic discriminant.

More generally, when removing an `r`-address finite fibre from a rational
one-place degree-`d` curve leaves a nonconstant radial map, the same argument
bounds its off-fibre cusp ramification by `2(d-r)-2`.  Distinct tangent
directions at an ordinary `r`-fold point guarantee nonconstancy.

## 6. Boundary and replay

The theorem proves that degree six already supports the normal nonmonogenic
`S3` cubic layer.  It closes only the family `(14)` and the four-cusp
rational one-place target.  A different degree-six coefficient grammar,
different singularity packet, positive-genus branch, or non-common-zero
index form may still work.  THM-3911 proves that the sharp control `(23)` has
scalar units and no class-group three-torsion, so it cannot support a normal
finite-flat cubic algebra with generic irreducible `S3` fibre and this branch
as generic-exponent-one discriminant divisor.  Nonnormal orders and different
  sextics remain outside that obstruction.  THM-3908 independently excludes
  the stated coefficient-depth-at-most-two one-point sextic lane.  THM-3913
  pays the unit and three-class invoice at degree ten/depth three, but its
  elliptic normalization has no polynomial plane atlas.  A future one-place
  sextic still requires its own resolvent unit/class-group audit before
  inverse Delone--Faddeev or Keller realization.  `JC(2)` remains **OPEN**.

Reproduce the exact packet with

```bash
python3 04-computation/jc2_degree_six_common_zero_normal_cubic_two_place_thm3906.py
python3 -O 04-computation/jc2_degree_six_common_zero_normal_cubic_two_place_thm3906.py
```

Both streams must byte-match the frozen output named in the metadata.
**QED.**
