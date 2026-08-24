---
id: THM-3907
title: "Unit-ideal nonmonogenic cubic six-place boundary"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  The binary cubic
  A U^3+C U^2V+(AC-1)UV^2+A V^3 has coefficient ideal one but represents no
  nonzero scalar unit.  Its Delone--Faddeev algebra is a normal finite-flat,
  globally nonmonogenic S3 cubic order with absolutely irreducible degree-seven
  discriminant.  The projective branch has two infinity support points and
  exactly six normalization places, three above each.  Thus common coefficient
  zero is not necessary for nonmonogenicity, but this first unit-ideal escape
  does not solve the one-place boundary.  Unit-ideal one-place forms, Keller
  realization, and JC(2) remain open.
source: root / unit-ideal escape audit, 2026-08-23
audit: >
  INDEPENDENTLY HOSTILE-AUDITED on 2026-08-23.  The audit independently
  rederived the unit-ideal witness, specialization obstruction, degree-three
  rational-map irreducibility, repeated-root incidence and absolute
  irreducibility, normality, S3 group, and both infinity Newton packets.  The
  assertion-free companion freezes the discriminant, incidence resultant,
  squarefree sidecar, and the 3+3 branch count in 26 active gates.  Normal and
  optimized runs byte-match the frozen output.
depends_on:
  - THM-3801-cubic-etale-normalization-nonmonogenic-and-companion-sheet-gate
related:
  - THM-3808-homogeneous-linear-binary-cubic-veronese-unit-trap
  - THM-3890-universal-quintic-common-zero-resolvent-class-group-dichotomy
  - THM-3906-degree-six-common-zero-normal-cubic-two-place-boundary
script: 04-computation/jc2_unit_ideal_nonmonogenic_cubic_six_place_thm3907.py
output: 05-knowledge/results/jc2_unit_ideal_nonmonogenic_cubic_six_place_thm3907.out
script_sha256: a25cddf1d6abf44f97186bf9d5c947a6e4b7c2cc55a4c61a01e11256b4e8282d
output_sha256: acafe47d68ed7037d55991bf2a787b4c3854ae3b5e31196eaa5a2fac82f97ff5
semantic_sha256: 784084eee175020d1cd6f33768b3e535084924ce7cfee612eab3d687748852c9
hash_basis: raw LF bytes
---

# THM-3907 -- unit ideal does not mean globally monogenic

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

Let `k` be an algebraically closed field of characteristic zero and put
`R=k[A,C]`.  Consider

```text
Psi=A U^3+C U^2V+(AC-1)UV^2+A V^3.                       (1)
```

The coefficient ideal of `Psi` is `R`, but `Psi` represents no element of
`k*`.  The associated Delone--Faddeev algebra is a normal finite-flat domain
of rank three, is globally nonmonogenic, and has generic Galois group `S3`.
Its discriminant is absolutely irreducible of degree seven.  Its projective
closure has two points at infinity and three normalization places above each.

## 1. Unit coefficient ideal, no scalar-unit value

The coefficient ideal is the unit ideal because

```text
1=AC-(AC-1).                                               (2)
```

Suppose, however, that `Psi(X,Y)=lambda` for some `X,Y in R` and
`lambda in k*`.  Specialize `A=0`, writing `X0,Y0 in k[C]`.  Equation `(1)`
becomes

```text
X0 Y0(CX0-Y0)=lambda.                                     (3)
```

Every factor on the left must be a unit of `k[C]`.  Thus `X0,Y0` are
nonzero constants, while `CX0-Y0` is nonconstant, a contradiction.  Hence
`Psi` represents no scalar unit.  By the index-form criterion in THM-3801,
the cubic algebra is not monogenic.

This proves that the common-zero gate used in THM-3890 is sufficient, not
necessary, for global nonmonogenicity.

## 2. Generic irreducibility and the discriminant

Dehomogenize with `V=1`.  In the cubic-root variable `t`, equation `(1)` is

```text
C(t^2+At)+(At^3-t+A).                                     (4)
```

Over `k(A)`, the numerator and denominator of

```text
C=-(At^3-t+A)/(t^2+At)                                   (5)
```

are coprime.  Indeed the denominator roots `0,-A` give numerator values
`A` and `A(2-A^3)`, not zero as rational functions of `A`.  The rational map
in `(5)` has degree three.  Therefore `(4)` is irreducible over `k(A,C)`.

Direct calculation gives

```text
Delta=-4A^4C^3-27A^4+30A^3C^2+A^2C^4
      -30A^2C-6AC^3+4A+C^2.                              (6)
```

To prove absolute irreducibility, parameterize the genuine repeated-root
incidence.  Such a finite repeated root has `t!=0`: at `t=0`, the cubic and
its derivative specialize to `A` and `AC-1`, which cannot vanish together.
After localizing by `t`, equations `Psi=partial_t Psi=0` imply

```text
C t^2=A(1-2t^3),
g=(1-2t^3)A^2+(2t-t^4)A-t^2=0.                           (7)
```

As a quadratic in `A`,

```text
Disc_A(g)=t^2(t^6-12t^3+8).                              (8)
```

The second factor is squarefree, hence is not a square in `k(t)`.  Thus
`g=0` is irreducible after localization by `t`; the first equation in `(7)`
then makes the genuine incidence its graph and hence irreducible.  The raw
unsaturated equations also contain the extraneous line `A=t=0` with `C`
free, which is not a repeated-root locus.  Elimination records that artifact
in the exact identity

```text
Res_t(g,Ct^2-A(1-2t^3))=-A^3 Delta.                       (9)
```

Every discriminant point except `(0,0)`, the sole repeated-root-at-infinity
case, lies in this irreducible incidence image.  A plane curve cannot gain an
isolated point component, so `V(Delta)` has irreducible support.  Finally,
the linear term `4A` in `(6)` shows that `Delta` is reduced, hence absolutely
irreducible.

The generic cubic is irreducible and `Delta` is nonsquare, so its Galois
group is `S3`.

## 3. Normality of the cubic order

Let `T` be the Delone--Faddeev algebra.  It is finite free over the regular
surface `R`, hence `S2`, and generic irreducibility makes it a domain.  Away
from `(Delta)` it is etale.  At the unique height-one branch prime,
`v_Delta(disc T)=1`.  The DVR discriminant/index formula forces index zero
in the integral closure.  Thus every height-one localization is normal,
making `T` `R1`; Serre's criterion proves normality.

## 4. Two infinity points carry six places

The leading form of `(6)` is

```text
Delta_7=-4A^4C^3.                                        (10)
```

Thus the infinity support is exactly `[1:0:0]` and `[0:1:0]`.

At the first point, use `x=C/A,z=Z/A`.  The tangent cone is

```text
-4x^3-27z^3,                                             (11)
```

a product of three distinct lines.  Hence this point has three
normalization places.

At the second point, use `y=A/C,z=Z/C`.  Its tangent cone is `y^2z`, and
the two Newton edges give

```text
z=4y^2+...,
y=(3+2sqrt(2))z^2+...,
y=(3-2sqrt(2))z^2+....                                   (12)
```

These are three distinct smooth branches.  Therefore

```text
infinity support points=2,              infinity places=3+3=6. (13)
```

## 5. Boundary and replay

The example removes common coefficient zero from the list of necessary
conditions for a normal nonmonogenic cubic order.  It does not solve the
infinity problem: it has two support points and six places.  A unit-ideal
index form with an irreducible one-place discriminant, an associated Keller
atlas, and `JC(2)` all remain **OPEN**.

Reproduce the exact packet with

```bash
python3 04-computation/jc2_unit_ideal_nonmonogenic_cubic_six_place_thm3907.py
python3 -O 04-computation/jc2_unit_ideal_nonmonogenic_cubic_six_place_thm3907.py
```

Both streams must byte-match the frozen output named in the metadata.
**QED.**
