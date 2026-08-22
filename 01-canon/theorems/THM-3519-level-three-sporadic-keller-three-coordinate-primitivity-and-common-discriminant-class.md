---
id: THM-3519
title: "Level-three sporadic Keller three-coordinate primitivity and common discriminant class"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENT CLEAN-ROOM AUDIT.  For the fixed
  sporadic Keller map at level three, each actual source coordinate x,y,z is
  primitive in the generic degree-27 function-field extension.  Their monic
  minimal-polynomial discriminants, and every saturated degree-27 coordinate
  eliminant obtained by target-field scaling, have the common full square
  class [-2J].  A lawful F_101 fibre at target (77,62,4) consists of nine
  separable rank-three quotient blocks; all 27 local coordinate
  characteristic polynomials are cubic and squarefree, all three global
  degree-27 products are squarefree, and independent 27x27 power matrices
  have ranks (27,27,27).  The parent target (93,28,83) is a sharp hostile:
  its y power rank is 26 and two blocks share T+76.  This is fixed-map,
  level-three scope only.
source: codex/keller-level-three-three-coordinate-independent-audit/2026-08-16
audit: >
  The clean-room companion imports neither the candidate nor its parent.  It
  re-derives the inverse section symbolically, proves DQ-xK^2=-3g and the map
  congruences, checks det J_F=-2, enumerates all of F_101^3 once to recover
  both depth-two split trees directly, computes quotient inverses by modular
  linear solves, computes characteristic polynomials by direct determinant,
  and independently checks three 27x27 power-basis determinants.  Normal,
  optimized, and stored transcripts agree after LF normalization.
depends_on:
  - THM-2473-sporadic-keller-branch-tower-depressed-trisection-anatomy
  - THM-3495-level-three-sporadic-keller-norm-divisor-and-three-component-nonproperness
related:
  - THM-2546-integral-coordinate-dichotomy-and-parity-lens-scope
  - THM-3508-level-two-sporadic-keller-three-coordinate-primitive-discriminant-square-class
  - MISTAKE-413
script: 04-computation/keller_level_three_three_coordinate_primitive_independent_audit_20260816.py
output: 05-knowledge/results/keller_level_three_three_coordinate_primitive_independent_audit_20260816.out
script_sha256: 10520f0314a521ca479c0971916b52794084cc287f446bd05bc4eb0ab55e75e2
output_sha256: 859104c1aee350220b54142d8be55270ddbb3630e047ce492bf479d3df2e1e91
semantic_sha256: 51b5945a5f26c3335b6837b169a916d9562acd8448d2dde4f1fe1449c1e157a3
hash_basis: LF-normalized UTF-8 bytes; exact finite-field and symbolic semantic ledger
---

# THM-3519 -- all three level-three coordinate views carry `[-2J]`

**PROVED + VERIFIED-EXACT + INDEPENDENT CLEAN-ROOM AUDIT.**

## 1. Statement

Let `F` be the fixed sporadic Keller map of THM-2473, let

```text
K=Q(a,b,c),
K_3=Q(x,y,z),                    (a,b,c)=F^3(x,y,z).     (1)
```

THM-3495 proves `[K_3:K]=27`.  Then every actual source coordinate is a
primitive element:

```text
K(x)=K(y)=K(z)=K_3.                                      (2)
```

Let `m_x,m_y,m_z in K[T]` be their monic minimal polynomials, and let `J` be
the primitive irreducible third-image polynomial in THM-3495's exact
normalization.  Then

```text
[Disc(m_x)]=[Disc(m_y)]=[Disc(m_z)]=[-2J]
                                      in K*/K*2.          (3)
```

If `E_x,E_y,E_z` are saturated actual degree-27 coordinate eliminants and
each differs from its corresponding monic minimal polynomial by a nonzero
element of `K`, their discriminants also satisfy `(3)`.

The theorem concerns this one map and the third self-composition.  It does
not assert three-coordinate primitivity for another Keller map, for the
outside counterexample family, or at a higher level.

## 2. The inverse formula is derived from the map

For a target `(a,b,c)`, put

```text
L=27a^2c^2-18abc+16a+b^3c-b^2,
g(X)=LX^3+(4-3bc)X-2c.                                  (4)
```

Write `r=1+xy`.  The third coordinate of `F(x,y,z)=(a,b,c)` gives

```text
z=[x(2-3xy)-c]/x^3.                                    (5)
```

Substitution into the first two displayed coordinates of THM-2473 reduces
them exactly to

```text
ax^3=r[x(1+r)-cr^2],
bx^2=2x(1+2r)-3cr^2.                                   (6)
```

Eliminating `r^2,r^3` from `(6)` gives

```text
3ax^2=r(1+bx-r),
Qr=xK,                                                 (7)

K=2+(9ac-b)x,
Q=3c(1+bx)-4x.
```

Define

```text
D=(12a-b^2)x^2+bx+2.                                   (8)
```

Direct expansion gives the load-bearing identity

```text
DQ-xK^2=-3g(x).                                        (9)
```

On the dense open set where the displayed quantities are units, equations
`(7)`--`(9)` give `rK=D`; substituting this in the first relation of `(7)`
gives

```text
y=b-3ax((9ac-b)x+2)/D,
z=(2x-c-3x^2y)/x^3.                                   (10)
```

Thus `(10)` is derived from the original map rather than guessed from a
finite witness.  The symbolic audit additionally substitutes `(10)` into
all three coordinates of `F`; after clearing denominators, the first two
residual numerators have `x`-degrees `9` and `6` and are exactly divisible
by `g`, while the third residual is identically zero.  The symbolic quotient
ledger has hash

```text
96ff553c33cb980eee181f895f730d7a4b8331c7ce59b4b677623284e4c401b3. (11)
```

The same calculation independently recovers `det J_F=-2`, so reduction
modulo `101` remains etale.

## 3. A lawful rank-27 fibre over `F_101`

The clean-room companion evaluates the original three polynomial coordinate
formulas at every point of `F_101^3` once.  It does not use `(10)` for this
census.  It reconstructs the fibre of

```text
(a,b,c)=(77,62,4)                                      (12)
```

through two completely split levels:

```text
level one:
(13,36,2), (39,84,75), (49,74,79),

level two:
((35,40,46),(83,54,93),(84,61,87)),
(( 2,91,44),(35, 7,50),(64,28,37)),
((24, 2,74),(84,50,58),(94,59, 2)).                   (13)
```

The thirteen `L` values at the final target, three first points, and nine
second points are

```text
(31,1,21,85,62,35,85,70,82,99,16,86,95),             (14)
```

so all are nonzero.  At all twelve split inverse edges, both `x` and the
denominator `D` in `(8)` are nonzero.  Hence both split stages lie in the
lawful inverse chart.

For each of the nine second points `q_i=(a_i,b_i,c_i)`, form

```text
A_i=F_101[X]/(g_i),
g_i=L(q_i)X^3+(4-3b_ic_i)X-2c_i.                      (15)
```

Every `g_i` has degree three and derivative gcd one.  In each quotient,
the multiplication determinants of `D` and `X^3` are respectively

```text
(73,7),(62,90),(42,89),(34,82),(42,56),
(35,49),(86,45),(32,100),(46,82),                     (16)
```

and are all nonzero.  Therefore `(10)` defines exact quotient elements
`y_i,z_i`.  Substitution into the original map gives `q_i` coefficient by
coefficient in every quotient.  Thus the specialized third fibre algebra is
the honest finite etale algebra

```text
A=A_1 x ... x A_9,                   dim_F101 A=27.    (17)
```

No claim that the nine last cubics split over `F_101` is used.

## 4. Independent characteristic polynomials and power bases

For each `xi in {X,y_i,z_i}`, the companion builds the multiplication matrix
on `A_i` in the basis `(1,X,X^2)`.  It computes

```text
det(TI-M_xi)                                             (18)
```

directly from the six determinant permutations, not from the candidate's
trace/trace-square formula, and checks Cayley--Hamilton inside `A_i`.
All `9*3=27` local characteristic polynomials are monic cubic and squarefree.
For `X`, `(18)` also equals the monic normalization of `g_i`.

Multiplying the nine local polynomials coordinate by coordinate gives three
monic degree-27 polynomials.  Their derivative gcds are all one, and their
ordered coefficient hashes are

```text
x: 306dc5831a989a0fe953ae5eb6127d26bd6ff07f8cc4e25d54dc0f65e53e5e19
y: aba63853d168e670461e67b67e42d26380e18c62b363c8931aa187f1a97d4e51
z: b8a09c5bfed7affdabbbb3fcf222002a3a7c11c8d55cb19bcc949494725bb51d. (19)
```

This independently reproduces all three candidate products.  As an
orthogonal conclusion check, the audit concatenates the nine quotient
coordinates and forms, for each source coordinate, the full `27 x 27`
matrix with columns

```text
1,xi,xi^2,...,xi^26.                                   (20)
```

The exact ranks and determinants modulo `101` are

```text
ranks        =(27,27,27),
determinants =(25,64,1).                               (21)
```

Thus every coordinate generates the specialized algebra `(17)`.  Equation
`(21)` proves the desired conclusion directly; the squarefree products in
`(19)` are a separately coded spectral cross-check.

## 5. The parent target is a sharp hostile

The same clean-room paths reconstruct the earlier lawful target

```text
(93,28,83).                                             (22)
```

All its local coordinate cubics remain squarefree, but the derivative-gcd
degrees of the three global products are

```text
(0,1,0),                                                (23)
```

and the three power ranks are

```text
(27,26,27).                                             (24)
```

The defect is localized exactly: flattened quotient blocks `4` and `5`,
corresponding to second points `(72,58,17)` and `(80,32,9)`, share the
linear `y` factor

```text
T+76,                    whose root is 25 in F_101.    (25)
```

Hence one lawful fibre can establish `x,z` primitivity while failing the
`y` test by precisely one dimension.  This is the minimal witness against
transporting coordinate primitivity without an explicit separator or power
determinant.

## 6. Good reduction proves generic primitivity

On the generic open set where the cubic leading coefficients and inverse
denominators are units, the three successive monic cubic constructions give
a free rank-27 tower basis over a localization of the target ring.  For a
source coordinate `xi`, express `(20)` in that basis and call its determinant
`d_xi`.  It is a rational function in `(a,b,c)` with integer coefficients
after a common denominator is cleared.

The fibre `(12)` lies in this open model by `(14)`--`(16)`.  Splitting the
first two cubics identifies its rank-27 specialization with `(17)`, and the
specializations of `d_x,d_y,d_z` are the three nonzero values in `(21)`.
Consequently none of the three rational functions is identically zero in
characteristic zero.  Therefore the powers in `(20)` are generically
linearly independent.  Since `[K_3:K]=27`, this proves `(2)`.

This argument uses one nonzero minor on one good fibre.  It does not assert
that every finite fibre separates every coordinate; `(22)`--`(25)` show the
sharp opposite.

## 7. Trace-form congruence gives the common full class

Fix any `K`-basis `e` of `K_3` and let `G_e` be its trace Gram matrix.  For a
primitive coordinate `xi`, let `C_xi` be the matrix whose columns are the
power basis `(20)` expressed in `e`.  Then

```text
Gram(1,xi,...,xi^26)=C_xi^T G_e C_xi,
Disc(m_xi)=det(C_xi)^2 det(G_e).                        (26)
```

Thus the three monic minimal-polynomial discriminants differ only by
literal squares in `K`.  THM-3495 supplies the exact full `x`-coordinate
class

```text
[Disc(m_x)]=[-2J].                                     (27)
```

Combining `(26)` and `(27)` proves `(3)`.  In particular, the constant unit
`[-2]` is retained; no divisor-parity argument discards it.  This is exactly
the MISTAKE-413 guardrail.

If a saturated degree-27 eliminant is `E_xi=c_xi m_xi`, then

```text
Disc(E_xi)=c_xi^(2*27-2) Disc(m_xi)
           =c_xi^52 Disc(m_xi)
           =(c_xi^26)^2 Disc(m_xi).                    (28)
```

Therefore denominator clearing or target-field scaling preserves `(3)`.
The saturation and degree conditions exclude a repeated polynomial attached
to a nonprimitive coordinate view.

## 8. Meaning of the three cubic `-4` factors

THM-2546 gives all three level-one coordinate discriminants in the form

```text
-4 * (square)^2 * L,                                   (29)
```

so their common square class is `[-L]`.  THM-3508 proves the corresponding
level-two class `[H]`; this theorem proves the level-three class `[-2J]`.
The common mechanism is not equality of the three displayed eliminants.
At each fixed level, primitivity makes their power bases bases of one trace
space, and trace-form congruence retains only one permutation-sign class.

This loses labelled roots, full wreath monodromy, exact discriminant
multiplicities, basis volumes, projection-collision divisors, and the other
Jelonek components.  It is a three-view theorem for one fixed compositional
tower, not a classification of dimension-three Keller counterexamples or of
the outside infinite family.

## 9. Reproduction and scope

Run from the repository root:

```text
python -B 04-computation/keller_level_three_three_coordinate_primitive_independent_audit_20260816.py
python -B -O 04-computation/keller_level_three_three_coordinate_primitive_independent_audit_20260816.py
```

Both executions reproduce the stored transcript exactly after LF
normalization.  The companion imports neither candidate computation nor its
parent, uses no `assert` gates, and has semantic digest

```text
51b5945a5f26c3335b6837b169a916d9562acd8448d2dde4f1fe1449c1e157a3. (30)
```

No exact positive discriminant multiplicities, level-four three-coordinate
primitivity, all-level induction, arbitrary Keller theorem, outside-family
classification, `JC(2)`, or LRC consequence follows.

**QED.**
