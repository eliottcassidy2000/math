---
id: THM-2055
title: "The rank-two determinant gate is a polygonal normal-fan problem"
status: >
  PROVED, elementary convex geometry. The maximum determinant in THM-2053 is
  the support norm of the centrally symmetric column polygon for a fixed
  ordered saturated basis. Only geometric hull vertices can own a residual
  direction, and each owner cone cuts the uncertified set to one explicit
  tangent disk. Basis changes give new valid sufficient gates but need not
  preserve this Euclidean polygon/disk picture. This is a carrier reduction,
  not LRC(14): every remaining lattice point still needs an LRC sidecar.
source: codex-2026-07-21-LRC-normal-fan-sail
depends_on:
  - THM-2053
related:
  - THM-2052
  - THM-2054
  - THM-2056
  - THM-2057
  - HYP-2896
  - HYP-2986
  - HYP-8871
  - MISTAKE-229
---

# THM-2055 -- determinant gate = polygonal normal fan

## 1. Fixed-basis scope and the convex body

Keep the notation of THM-2053, and fix an **ordered** `Z`-basis `u,z` of the
saturated lattice `U intersect Z^n`. In this basis, write

```text
c_i=(u_i,z_i) in Z^2,
d=(a,b),
D_U(d)=max_i |a z_i-b u_i|.                            (1)
```

Let `R(a,b)=(-b,a)` be counterclockwise quarter-turn and put

```text
K_U=conv{+c_i,-c_i:1<=i<=n}.                           (2)
```

Then

```text
D_U(d)=h_(K_U)(R d),                                   (3)
```

where `h_K(y)=max_(x in K) x dot y` is the support function. Indeed,
`c_i dot R d=a z_i-b u_i`.

Every object in (1)--(3), including the standard Euclidean norm, the
quarter-turn `R`, the polygon, its polar, its normal fan, and the tangent
disks below, is relative to this chosen basis. This is a load-bearing scope
condition. Explicitly, if

```text
(u',z')^T=A(u,z)^T,             A in GL_2(Z),
```

then `c_i'=A c_i`, while the coefficient vector representing the same speed
row is `d'=A^(-T)d`. The torus and speed row are merely reparameterized, but
the standard coordinate norm and the quantity `det(d,c_i)` need not be
preserved by a nonorthogonal `A`. Thus a new basis supplies another valid
THM-2053 sufficient gate; it does not canonically preserve the numerical
function `D_U`, polygon, radial roof, or disk carrier written here. A row may
be certified in any basis, but basis-specific residual carriers must not be
identified without transporting the metric and coordinates as well.
The subscript in `D_U` follows THM-2053's notation; it does not assert basis
invariance.

If the columns span `R^2`, then `K_U` is centrally symmetric with the origin
in its interior. Consequently `D_U` is a norm: its unit ball is the inverse
quarter-turn of the polar polygon `K_U^o`. In particular, it is generally a
**polyhedral norm**, not a binary quadratic form.

## 2. Hull-vertex deletion

Let `V(K_U)` be the vertex set of `K_U`. Equation (3) immediately gives

```text
D_U(d)=max_(p in V(K_U)) p dot R d.                    (4)
```

Thus a signed point `sigma c_i` that is not a geometric hull vertex can be
omitted from the maximum in the THM-2053 determinant test without changing
that maximum. Central symmetry makes `c_i` a vertex exactly when `-c_i` is a
vertex. If several runner labels realize the same vertex, one geometric copy
suffices for (4), but the full set of labels remains the owner sidecar.

This deletion applies **only** to the determinant maximum and its geometric
carrier. A nonvertex may lie in the interior or on an edge. Its runner label
must remain in the original speed row, the positivity and collision walls,
the transverse residue deck, pair-sum data, and every phase-height or endpoint
calculation. It can never be the unique geometric owner of gate failure, but
it can still be decisive for LRC.

This is the first exact atlas compression: the determinant sidecar depends on
the signed column polygon, not on all thirteen labelled columns separately.

## 3. Normal-fan owner sectors

For `p in V(K_U)`, let

```text
N_p={y in R^2:p dot y=h_(K_U)(y)}
C_p=R^(-1)N_p.                                         (5)
```

The cones `C_p` form the pulled-back normal fan of `K_U` and cover the
parameter plane. On `C_p`,

```text
D_U(d)=p dot R d.                                      (6)
```

Write `p=(r,s)` and `w_p=(s,-r)`. Since `p dot R d=d dot w_p`, failure of the
THM-2053 gate inside `C_p` is exactly

```text
d in C_p,
||d-(91/2)w_p||^2 < (91^2/4)||p||^2.                  (7)
```

Hence the uncertified primitive parameters are

```text
Z^2_prim intersect union_(p in V(K_U))
  (C_p intersect B((91/2)w_p,(91/2)||p||)).            (8)
```

Away from the origin, adjacent cones overlap on their owner-tie rays. On a
tie ray, all tied vertices give the same support value, so the corresponding
inequalities in (7) agree on that ray. A fixed cyclic half-open convention can
assign each nonzero parameter to one cone, but this is bookkeeping rather
than an intrinsic unique owner; later sidecars should retain the whole tied
owner set. The origin lies in every cone but is excluded by primitivity.

Every disk boundary passes through the origin, and every **uncertified** point
in the `p` sector satisfies the owner-local bound

```text
||d||<91||p||.                                         (9)
```

This is sharper than enumerating the round `91L(U)` disk whenever the active
hull vertex is shorter than the longest column.

## 4. Exact radial roof

For a unit vector `e` in the chosen coordinate Euclidean norm and real
`rho>0`, homogeneity of the support function gives

```text
D_U(rho e)=rho D_U(e).
```

Therefore, on this real ray, the **THM-2053 determinant inequality** is
equivalent to the radial condition

```text
rho>=91 D_U(e).                                        (10)
```

Equality is certified. This is not an if-and-only-if criterion for LRC safety:
points below the roof are merely uncertified by this gate. Nor is it an
enumeration of primitive multiples. An oriented rational lattice ray contains
exactly one primitive lattice point (the opposite ray contains its negative),
while an irrational ray contains no nonzero lattice point. The roof is a
real-geometric boundary whose Cartesian pieces are linear and whose
restriction to the unit circle is piecewise sinusoidal.

Each rational owner cone has two rational boundary rays, so ordinary
continued fractions or a Stern--Brocot/Klein sail can address its primitive
lattice points. THM-2056 makes this finite Farey certificate precise.

## 5. One-tail calibration

For `u=(1,...,13)` and `z=e_12`, all axial columns except `(13,0)` lie inside
the signed hull. Thus

```text
V(K_U)={+-(13,0), +-(12,1)},
D_U(a,b)=max(13|b|,|a-12b|).                           (11)
```

The nominal `26` signed-column disks collapse to four vertex sectors. Along
`d=(1,b)` with integral `b!=0`, equation (11) gives
`D_U(1,b)=13|b|`, recovering the `|b|>=1183` terminal from THM-2053.
THM-2057 now closes the whole scaled one-tail plane by missing-clock and
binding witnesses; the determinant residual was only its geometric address.

## 6. Predicate and loss ledger

The normal-fan quotient preserves exactly:

- the THM-2053 determinant value and safe/uncertified gate;
- the active signed hull owner and its tangent disk;
- primitive parameter coordinates and rational slope.

It forgets the non-hull runner labels, the exact phase-height optimum,
pair-sum rulers, clock divisibility, and endpoint ownership. Therefore an
empty residual proves a plane safe, but a nonempty sector is only
**unresolved**. It must be joined to a valid model-specific sidecar, such as
THM-2057's clock/binding certificate in the one-tail plane, or to a separately
verified application of THM-2054, pair-sum, or owner-labelled Euler machinery.
QED.
