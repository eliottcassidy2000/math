---
id: THM-3558
title: "Linear collision-plane to smooth target-hypersurface no-go"
status: >
  PROVED + VERIFIED-EXACT.  For the fixed THM-1300 three-dimensional Keller
  map, no affine source plane containing two points of the displayed triple
  collision has image contained in a target hypersurface smooth at the
  collision value.  Pairs involving the fixed point fail at the tangent
  plane.  The moving pair has exactly two tangent-compatible planes over
  Q(sqrt(13)), but their two image branches have unequal target Hessians.
  Hence no affine coordinate plane can descend this collision through any
  smooth target coordinate hypersurface.  Nonlinear source coordinate
  hypersurfaces remain open.
source: kps-s188
depends_on:
  - THM-1300-jacobian-counterexample-dixmier-A3-explicit
related:
  - THM-3546-invariant-graph-keller-descent-criterion
  - THM-3553-fixed-threedimensional-keller-map-polynomial-graph-section-no-go
  - THM-3554-punctured-kummer-collision-surface-normal-form
companion: 04-computation/jacobian_linear_collision_plane_second_jet_no_go_kps_s188.py
output: 05-knowledge/results/jacobian_linear_collision_plane_second_jet_no_go_kps_s188.out
script_sha256: fd56a535bcb48601a5eca0a22aed066d71d1e8ea9baa120a027e6373a4dc8288
output_sha256: 9764c1507dc148ccc9fc6af63d43186b0e3cb0c6ad0d2d1c516e97fd490f6cee
hash_basis: LF-normalized bytes
---

# THM-3558 -- linear collision-plane to smooth target-hypersurface no-go

**PROVED + VERIFIED-EXACT.**  In the fixed THM-1300 map, changing the target
linear projection does not rescue any affine source plane through a collision
pair.  More strongly, no such plane can map into an arbitrary target
hypersurface that is smooth at the common image.  The exceptional first-order
case exists only over `Q(sqrt(13))` and dies at the second jet.

All varieties below are over `C`.

## 1. Statement and local necessity

Put `u=1+xy` and write the THM-1300 Keller map as

```text
F1=u^3 z+y^2u(4+3xy),
F2=y+3xu^2z+3xy^2(4+3xy),
F3=2x-3x^2y-x^3z.                                      (1)
```

Its Jacobian determinant is `-2`, and

```text
p0=(0,0,-1/4),
p+=(1,-3/2,13/2),
p-=(-1,3/2,13/2)                                       (2)
```

all map to

```text
q=(-1/4,0,0).                                           (3)
```

Let `Pi` be any affine plane containing two distinct points in `(2)`.  Then
there is no algebraic hypersurface `T` in target `A^3`, smooth at `q`, such
that

```text
F(Pi) subset T.                                         (4)
```

Here set-theoretic containment is already enough.  Since `F` is etale, the
restriction of `dF` to the tangent plane of `Pi` has rank two.  If `(4)` held,
each of the two source branches would therefore map locally etale onto the
smooth two-dimensional germ `(T,q)`.  In particular:

1. the two transported tangent planes at `q` would be equal; and
2. after choosing two target coordinates on that tangent plane, the two
   image branches would be the same formal graph and hence have equal graph
   Hessians.

The first condition eliminates four of the six ordered collision branches.
The second eliminates the only surviving pair of source-plane parameters.

## 2. Transported normals

Let `J_p` be the Jacobian matrix of `F` at a source point `p`.  If a source
plane has row normal `n`, its image branch at `q` has row normal

```text
m_p=n J_p^(-1).                                         (5)
```

Indeed, `ker(n)` is the source tangent plane, and `(5)` is the unique target
normal whose pullback by `dF_p` is `n`.  Thus two target tangent planes agree
exactly when their normals in `(5)` are proportional.

### 2.1. Planes through `p0,p+`

Every normal to a plane through these two points has the form

```text
n=((6B-27C)/4,B,C).                                     (6)
```

The cross product of the transported normals at `p0` and `p+` is

```text
(-27/32 (B^2-6BC+3C^2),
  27/16 (4B^2-38BC+77C^2),
  -3/4  (8B^2-53BC+C^2)).                               (7)
```

The resultant in `B` of the first two parenthesized quadratics is

```text
-647 C^4.                                               (8)
```

If `C` is nonzero, `(8)` forbids their simultaneous vanishing.  If `C=0`,
the first component in `(7)` vanishes only when `B=0`, which is not a plane
normal.  Hence the two image tangents never agree.

### 2.2. Planes through `p0,p-`

The reflected family has normal

```text
n=((6B+27C)/4,B,C),                                     (9)
```

and transported-normal cross product

```text
(-27/32 (B^2+6BC+3C^2),
 -27/16 (4B^2+38BC+77C^2),
   3/4  (8B^2+53BC+C^2)).                              (10)
```

The corresponding first-two-component resultant is again `-647 C^4`.
The same `C!=0` / `C=0` split excludes every projective normal.

### 2.3. Planes through `p+,p-`

Every plane through the moving pair has equation

```text
3r x+2r y+2c z-13c=0,                  [r:c] in P^1.   (11)
```

For `n=(3r,2r,2c)`, the two target normals are

```text
m+ =(-52c+12r, -3c/2+r/2,  27c/8-9r/8),
m- =(-52c-12r,  3c/2+r/2, -27c/8-9r/8).               (12)
```

Their cross product is

```text
(0,27(r^2-13c^2),12(r^2-13c^2)).                       (13)
```

Therefore the target tangents agree only when `c!=0` and

```text
a=r/c,                         a^2=13.                 (14)
```

These are the two exceptional source planes

```text
3a x+2a y+2z-13=0.                                     (15)
```

## 3. The exceptional planes fail at second order

On `(15)`, substitute

```text
z=(13-a(3x+2y))/2                                      (16)
```

into `(1)`.  Use `(alpha,beta)=(F1,F2)` as local target coordinates and
write the third target coordinate locally as

```text
gamma=F3=g(alpha,beta).                                (17)
```

The two source-to-`(alpha,beta)` Jacobian determinants at `p+` and `p-` are

```text
9(a-3)/8,                    9(a+3)/8,                 (18)
```

respectively.  They are nonzero because `a^2=13`, so `(17)` is legitimate
on both branches.  Exact implicit differentiation gives the same first jet:

```text
grad g+=grad g-=(-32a/9,4/9).                          (19)
```

The target-coordinate Hessians are

```text
H+ = [[320a/9+10816/27, 16a/9+224/9],
      [16a/9+224/9,             -16/27]],

H- = [[320a/9-10816/27, 224/9-16a/9],
      [224/9-16a/9,               16/27]].             (20)
```

Consequently

```text
H+-H-=[[21632/27,32a/9],
       [32a/9,   -32/27]] != 0.                        (21)
```

If both branch germs lay in one smooth target hypersurface germ, `(17)`
would be the unique graph equation of that germ, so its Hessian would be
unique.  Equation `(21)` is the required contradiction.  This proves `(4)`
for the final collision pair and completes the theorem.

## 4. Consequence and surviving counterexample architecture

A coordinate hypersurface in affine three-space is smooth and isomorphic to
`A^2`.  Hence no affine source coordinate plane containing a collision pair
can participate in THM-3546's planar descent mechanism for this fixed map,
regardless of the target polynomial coordinate system.

This is stronger than THM-3553 in the target direction but narrower in the
source direction:

- THM-3553 treats polynomial graphs of every degree, but only the displayed
  target graph chart;
- the present theorem treats every target hypersurface smooth at `q`, but
  only affine source planes.

Their intersection closes the complete linear descent cell.  It does **not**
exclude a nonlinear source coordinate hypersurface, a hypersurface singular
at `q`, a different ambient Keller map, or a direct two-dimensional
construction.  For coordinate-hypersurface descent, the sharp surviving cell
is therefore a nonlinear source coordinate surface whose two transported
branch germs agree to all orders in one smooth target coordinate surface.

## 5. Exact verification

Reproduce with

```bash
python3 04-computation/jacobian_linear_collision_plane_second_jet_no_go_kps_s188.py
python3 -O 04-computation/jacobian_linear_collision_plane_second_jet_no_go_kps_s188.py
```

The ordinary and optimized transcripts agree byte-for-byte.  The companion
checks the ambient determinant and collision, all three normal families,
the two resultants `(8)`, the exceptional quadratic `(14)`, the invertible
local coordinates `(18)`, and the complete first- and second-jet identities
`(19)--(21)` by exact symbolic arithmetic.

**QED.**
