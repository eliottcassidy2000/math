---
id: THM-3887
title: "Binary-cubic common-zero quintic unibranch obstruction"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  If all four
  coefficients of a binary cubic vanish at a plane point, its discriminant
  has multiplicity at least four there.  When the quartic tangent cone is a
  fourth power, the coefficient pencil has rank one and the quintic jet is
  divisible by the cube of that tangent.  Consequently no reduced
  irreducible discriminant of total degree at most five is unibranch at such
  a common zero.  Higher degree, multibranch quintics, absence of a common
  zero, cubic-order realizability, Keller atlases, and JC(2) remain OPEN.
source: root / binary-cubic index-form common-zero lane, 2026-08-23
audit: >
  INDEPENDENTLY HOSTILE-AUDITED on 2026-08-23.  The two projective
  root-multiplicity cases for a
  singular binary cubic are reduced independently to `X^2Y` and `X^3`; their
  discriminant expansions bound finite pencil contact by two and three.
  The rank-one Taylor jet and the degree-five factor are checked symbolically.
  The audit independently replayed the contact calculation, rank-one
  implication, Taylor divisibility, and all 17 assertion-free gates in normal
  and optimized modes.  It also checked every quartic root partition and the
  genus-six finite-address packet used by THM-3890.  Normal and optimized runs
  must byte-match the frozen output.
depends_on: []
related:
  - THM-3801-cubic-etale-normalization-nonmonogenic-and-companion-sheet-gate
  - THM-3855-formal-inverse-discriminant-lift-and-algebraization-gate
  - THM-3890-universal-quintic-common-zero-resolvent-class-group-dichotomy
script: 04-computation/jc2_binary_cubic_common_zero_quintic_thm3887.py
output: 05-knowledge/results/jc2_binary_cubic_common_zero_quintic_thm3887.out
script_sha256: 5e6f79c63bc8d3cb565ba0b9b5f01dac0cdc0443bfb4959063daa1337844701a
output_sha256: 1cfdf303a09e50c95c88c1496169337289631e2d409168416e4ea0cccd17ab7d
semantic_sha256: 67b4b929172cefd97c23b2be58275980b27f39663b1969ee6a4c08b8369a45c5
hash_basis: raw LF bytes
---

# THM-3887 -- a common-zero index gate cannot start with a unibranch quintic

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

Let `k` be an algebraically closed field of characteristic zero, put

```text
R=k[u,v],                 m=(u,v),
```

and consider the binary cubic

```text
Phi(X,Y)=aX^3+bX^2Y+cXY^2+dY^3,       a,b,c,d in m.       (1)
```

Its discriminant is

```text
D=b^2c^2-4ac^3-4b^3d-27a^2d^2+18abcd.                    (2)
```

This is the binary index-form situation behind the cheapest sufficient
nonmonogenicity gate in THM-3801: the four coefficients have a common zero.

## 1. The fourth-order tariff

Because `(2)` is homogeneous of degree four in `(a,b,c,d)`,

```text
D in m^4.                                                   (3)
```

Write `a_1,b_1,c_1,d_1` for the homogeneous linear jets and `D_4` for the
degree-four homogeneous jet of `D`.  Then

```text
D_4=Disc(a_1X^3+b_1X^2Y+c_1XY^2+d_1Y^3).                  (4)
```

Thus the tangent cone is the discriminant of the pencil cut out by the
linear coefficient map

```text
k^2 -> Sym^3(k^2),       (u,v) |-> uF+vG.                 (5)
```

## 2. A genuine cubic pencil has no finite fourth-order discriminant contact

Suppose `F,G` are independent and the pencil is not contained in the cubic
discriminant hypersurface.  At a singular member there are only two root
partitions.

For a double root and a simple root, change binary variables and the pencil
parameter so that the member is `X^2Y`.  Write the transverse member as

```text
alpha X^3+beta X^2Y+gamma XY^2+delta Y^3.
```

Then the discriminant of `X^2Y+tG` begins

```text
-4delta t + (gamma^2-12beta delta)t^2 + ... .              (6)
```

If `delta!=0`, the contact is one.  If `delta=0` and `gamma!=0`, it is two.
If both vanish, every member is divisible by `X^2`, so the whole pencil is
singular.

For a triple root, use `X^3`.  The discriminant of `X^3+tG` begins

```text
-27delta^2 t^2 + (-54alpha delta^2+18beta gamma delta
                  -4gamma^3)t^3 + ... .                   (7)
```

If `delta!=0`, the contact is two.  If `delta=0` and `gamma!=0`, it is three.
If both vanish, every member is divisible by `X^2`, so the pencil is again
contained in the discriminant.  Hence

```text
finite intersection multiplicity of a line with Disc=0 <= 3.              (8)
```

In particular a rank-two pencil cannot have nonzero discriminant
`lambda ell^4`: that would give one singular member with finite contact four.

## 3. The rank-one jet forces a cubic tangent divisor at order five

Assume now

```text
D_4=lambda ell^4,             lambda in k^*.               (9)
```

By Section 2 the linear coefficient map has rank one.  Therefore, for a
fixed nonsingular binary cubic `H`,

```text
(a_1,b_1,c_1,d_1)=ell H,      Disc(H)=lambda.              (10)
```

Let `Q_2` be the vector of quadratic coefficient jets.  Homogeneity of the
discriminant gives the Taylor expansion

```text
Disc(ell H+Q_2+higher)
  =ell^4 Disc(H) + ell^3 (d Disc)_H(Q_2) + terms of degree >=6.             (11)
```

Consequently

```text
ell^3 divides D_5,                                                (12)
```

where `D_5` is the homogeneous quintic jet.  This divisibility is intrinsic:
it does not depend on the chosen coordinates or normalization of `H`.

## 4. No reduced irreducible quintic can be unibranch at the common zero

Suppose `D` is nonzero, reduced, irreducible, has total degree at most five,
and the germ `D=0` at the origin is unibranch.  Equation `(3)` leaves only
two cases.

If `ord_m(D)=5`, then `D` is homogeneous of degree five.  Over the
algebraically closed field it is a product of linear forms, a contradiction.

If `ord_m(D)=4`, unibranchness forces its tangent cone to be a fourth power:

```text
D_4=lambda ell^4.                                          (13)
```

There are no terms above degree five, so `(11)--(12)` give

```text
D=ell^4 lambda+ell^3 q_2=ell^3(lambda ell+q_2),            (14)
```

again contradicting reduced irreducibility.  Therefore

```text
common zero + reduced irreducible + local unibranch
       ==> total degree of the binary-cubic discriminant >= 6.             (15)
```

Equivalently, a viable degree-five common-zero design must be multibranch at
that point.  The generic degree-four tangent cone of a genuine pencil has
several tangent directions, so this theorem does not eliminate quintics; it
changes their required singularity packet from the tempting `(4,5)` cusp to
a multi-address configuration.

## 5. Exact boundary

The theorem uses only the common-zero sufficient gate.  It says nothing
against an intrinsically nonmonogenic binary cubic whose coefficients generate
the unit ideal but whose index form represents no scalar unit.  It also does
not construct or exclude a multibranch quintic, a higher-degree one-place
discriminant, a normal finite-flat cubic order, an etale plane atlas, or a
Jacobian counterexample.  `JC(2)` remains **OPEN**.

Reproduce the exact contact and Taylor packet with

```bash
python3 04-computation/jc2_binary_cubic_common_zero_quintic_thm3887.py
python3 -O 04-computation/jc2_binary_cubic_common_zero_quintic_thm3887.py
```

Both runs must byte-match
`05-knowledge/results/jc2_binary_cubic_common_zero_quintic_thm3887.out`.
