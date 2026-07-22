---
id: THM-2055
title: "The rank-two determinant gate is a polygonal normal-fan problem"
status: >
  PROVED, elementary convex geometry. The maximum determinant in THM-2053 is
  the support norm of the centrally symmetric column polygon. Only hull
  vertices can own a residual direction, and each owner cone cuts the
  uncertified set to one explicit tangent disk. This is a carrier reduction,
  not LRC(14): every remaining lattice point still needs a pair-sum, Fejer, or
  owner-labelled Euler certificate.
source: codex-2026-07-21-LRC-normal-fan-sail
depends_on:
  - THM-2053
related:
  - THM-2052
  - THM-2054
  - HYP-2896
  - HYP-2986
  - HYP-8871
  - MISTAKE-225
---

# THM-2055 -- determinant gate = polygonal normal fan

## 1. The correct convex body

Keep the notation of THM-2053. Thus a saturated basis `u,z` of a rational
two-plane has columns

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

If the columns span `R^2`, then `K_U` is centrally symmetric with the origin
in its interior. Consequently `D_U` is a norm: its unit ball is the inverse
quarter-turn of the polar polygon `K_U^o`. In particular, it is generally a
**polyhedral norm**, not a binary quadratic form.

## 2. Hull-vertex deletion

Let `V(K_U)` be the vertex set of `K_U`. Equation (3) immediately gives

```text
D_U(d)=max_(p in V(K_U)) p dot R d.                    (4)
```

Thus any column `+/-c_i` that is not a hull vertex can be deleted from the
THM-2053 determinant test without changing it. Interior columns still matter
to the original LRC row, but they can never own failure of the geodesic gate.

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

The cones overlap only on owner-tie rays; a fixed cyclic half-open convention
makes (8) disjoint. Every disk boundary passes through the origin, and every
point in the `p` sector satisfies the owner-local bound

```text
||d||<91||p||.                                         (9)
```

This is sharper than enumerating the round `91L(U)` disk whenever the active
hull vertex is shorter than the longest column.

## 4. Exact radial roof

For a unit vector `e` and `rho>0`, homogeneity of the support function gives

```text
D_U(rho e)=rho D_U(e).
```

Therefore the determinant certificate is valid precisely beyond the radial
roof

```text
rho>=91 D_U(e).                                        (10)
```

The roof is piecewise sinusoidal, with pieces indexed by the rational normal
fan. This is the correct place to use ordinary continued fractions: each
rational owner cone has two rational boundary rays, so its primitive lattice
points can be addressed by a Stern--Brocot/Klein sail without invoking a
multidimensional continued-fraction analogy.

## 5. One-tail calibration

For `u=(1,...,13)` and `z=e_12`, all axial columns except `(13,0)` lie inside
the signed hull. Thus

```text
V(K_U)={+-(13,0), +-(12,1)},
D_U(a,b)=max(13|b|,|a-12b|).                           (11)
```

The nominal `26` signed-column disks collapse to four vertex sectors. Along
`d=(1,b)`, equation (11) gives `D_U(1,b)=13|b|` for `b!=0`, recovering the
`|b|>=1183` terminal from THM-2053. HYP-2896 then supplies the arithmetic
resonance fan inside the remaining sector intervals.

## 6. Predicate and loss ledger

The normal-fan quotient preserves exactly:

- the THM-2053 determinant value and safe/uncertified gate;
- the active signed hull owner and its tangent disk;
- primitive parameter coordinates and rational slope.

It forgets the non-hull runners, the exact phase-height optimum, pair-sum
rulers, and endpoint ownership. Therefore an empty residual proves a plane
safe, but a nonempty sector is only **unresolved**. It must be joined to
THM-2054 relative decorrelation, HYP-2896-style resonance fans, or the
THM-2047/HYP-2108 owner-labelled Euler gate. QED.
