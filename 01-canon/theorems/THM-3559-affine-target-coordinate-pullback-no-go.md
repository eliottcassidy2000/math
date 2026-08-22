---
id: THM-3559
title: "Affine target-coordinate pullback no-go for the fixed Keller map"
status: >
  PROVED + VERIFIED-EXACT.  Let ell be any affine-linear target coordinate
  vanishing at the fixed THM-1300 triple-collision value.  No irreducible
  coordinate factor of ell(F) contains two collision points.  The proof
  classifies the complete factor locus of ell(F)=Kz+L.  Outside two Kummer
  parameters and the F3 row it is primitive; projection to the (x,y)-plane,
  compactly supported Euler characteristic, and nonconstant units exclude an
  affine-plane fibre.  The exceptional Kummer cylinders are G_m x A1 and
  each meets at most one collision point; the F3 row is xC.  Thus a descent
  of the known threefold collision must make both source and target
  coordinate hypersurfaces nonlinear.
source: kps-s188
depends_on:
  - THM-1300-jacobian-counterexample-dixmier-A3-explicit
related:
  - THM-3546-invariant-graph-keller-descent-criterion
  - THM-3554-punctured-kummer-collision-surface-normal-form
  - THM-3558-linear-collision-plane-smooth-target-hypersurface-no-go
companion: 04-computation/jacobian_affine_target_pullback_factor_atlas_kps_s188.py
output: 05-knowledge/results/jacobian_affine_target_pullback_factor_atlas_kps_s188.out
script_sha256: b51e6a258ce7854d71142b93555498de472082c6fa7d37cecbbc5fc67792b148
output_sha256: 6fdd20ee8944e915f05bc58aad7cfa2c07c7bddb53bdd229cbf2531a36cc2c36
hash_basis: LF-normalized bytes
---

# THM-3559 -- affine target-coordinate pullback no-go

**PROVED + VERIFIED-EXACT.**  The coordinate-hypersurface descent of
THM-3546 cannot use a linear coordinate on either side.  THM-3558 closes a
linear source coordinate against every smooth target coordinate surface.
The present theorem closes an affine-linear target coordinate against every
source-coordinate factor of its pullback.

All varieties are over `C`, and Euler characteristic means compactly
supported topological Euler characteristic.

## 1. Complete affine target pencil

Retain the fixed THM-1300 map

```text
u=1+xy,
F1=u^3z+y^2u(4+3xy),
F2=y+3xu^2z+3xy^2(4+3xy),
F3=2x-3x^2y-x^3z.                                      (1)
```

The triple collision points and their common image are

```text
p0=(0,0,-1/4),
p+=(1,-3/2,13/2),
p-=(-1,3/2,13/2),
q=(-1/4,0,0).                                           (2)
```

Every affine-linear target coordinate vanishing at `q` is, up to a nonzero
scalar,

```text
ell=A(X+1/4)+BY+CZ,                 [A:B:C] in P^2.     (3)
```

Set

```text
H=ell(F)=K(x,y)z+L(x,y).                                (4)
```

Direct expansion gives the structured leading coefficient

```text
K=A u^3+3Bx u^2-Cx^3.                                  (5)
```

This binary cubic in `(u,x)` is the source of the complete factor atlas.

## 2. Common-factor classification

When `A!=0`, every reduced component of `K=0` has equation

```text
d_r:=u-rx=0                                             (6)
```

for a root of

```text
P(r)=Ar^3+3Br^2-C.                                      (7)
```

The parametrization

```text
x=t,                    y=r-t^(-1)                     (8)
```

identifies `V(d_r)` with `G_m`.  Exact restriction of the constant term in
`(4)` gives

```text
t L|_(d_r=0)
 =3rP(r)t^3-5P(r)t^2
  +(Ar^2+A/4+4Br)t+Ar+2B.                              (9)
```

Thus `d_r` divides both `K` and `L` exactly when all four coefficients in
`(9)` vanish.  A nonzero solution must have `A!=0`: if `A=0`, the last
coefficient forces `B=0`, and then `(7)` forces `C=0`.  Normalize `A=1`.
The remaining ideal has reduced Groebner basis

```text
C+r/8,                   B+r/2,             r^2-1/4.  (10)
```

Consequently the only common Kummer factors occur at

```text
r= 1/2, [A:B:C]=[1:-1/4:-1/16],
r=-1/2, [A:B:C]=[1: 1/4: 1/16].                       (11)
```

At either parameter the root `r` is double in `K`, but `d_r` divides `H`
only once.  The exact root polynomials are

```text
(2R-1)^2(4R+1)/16,       (2R+1)^2(4R-1)/16.           (12)
```

There is one additional type of common factor.  Equation `(5)` is divisible
by `x` only if `A=0`, while

```text
L(0,y)=By                 when A=0.                    (13)
```

Hence `x` is common only in the pure `F3` row `A=B=0`, where

```text
H=CF3=-Cx(x^2z+3xy-2).                                 (14)
```

Equations `(10)--(14)` classify the full gcd `G=gcd(K,L)`: it is `1`, one
of the two simple factors `d_(+/-1/2)`, or `x` in the pure `F3` row.

## 3. Primitive affine modifications are not planes

Write

```text
H=G R,               R=K0z+L0,             gcd(K0,L0)=1.   (15)
```

Gauss's lemma makes `R` irreducible.  We now prove that `V(R)` is never
isomorphic to `A^2`.

For any primitive polynomial `R=K0z+L0`, project its zero surface `S` to the
`(x,y)`-plane.  Put

```text
D=V(K0)_red,                  Z=V(K0,L0).               (16)
```

Primitivity makes `Z` finite.  The projection is an isomorphism over
`A^2 minus D`, has fibre `A^1` over every point of `Z`, and has empty fibre
over `D minus Z`.  Additivity and multiplicativity of Euler characteristic
therefore give

```text
chi(S)=1-chi(D)+#Z.                                     (17)
```

There are three cases.

### 3.1. `A!=0`

The reduced curve `D` is a nonempty disjoint union of curves `(6)`, each
isomorphic to `G_m`; distinct roots give disjoint curves because `x=0` and
`u=0` cannot hold together.  Thus `chi(D)=0` and

```text
chi(S)=1+#Z.                                            (18)
```

If `Z` is nonempty, `(18)` differs from `chi(A^2)=1`.  If `Z` is empty,
then `K0` has no zero on `S`, so it is a nonconstant unit in `C[S]`; this is
again impossible for `A^2`, whose only units are constants.  This includes
the primitive residual factors at both exceptional parameters `(11)`.

### 3.2. `A=0,B!=0`

Here

```text
K0=x(3B u^2-Cx^2).                                     (19)
```

The reduced divisor `D` is one `A^1`, namely `x=0`, disjoint from at least
one `G_m` component.  By `(13)`, `Z` meets the `A^1` component in exactly
the point `(0,0)`.  If `N` is the number of its points on the `G_m`
components, then

```text
chi(S)=1-1+(1+N)=1+N.                                  (20)
```

For `N>0`, Euler characteristic excludes a plane.  For `N=0`, the nonconstant
factor `3Bu^2-Cx^2` has no zero on `S` and is a unit, again excluding a
plane.

### 3.3. Pure `F3`

After removing the factor `x` in `(14)`, the residual surface is

```text
2-3xy-x^2z=0.                                           (21)
```

Here `D=V(x)` and `Z` is empty, so `(17)` gives `chi(S)=0`.  Equivalently,
THM-3554 identifies `(21)` explicitly with `G_m x A^1`.  It is not a
coordinate plane.

This proves the claim about every primitive residual `R`.

## 4. No factor can carry a collision pair

It remains only to inspect the extracted factors `G`.

Each exceptional Kummer cylinder has equation `d_r=0` and coordinate ring

```text
C[x,x^(-1),z],                                          (22)
```

so it is `G_m x A^1`, not a coordinate hypersurface.  More sharply, its
values at the three collision points are

```text
d_r(p0)=1,
d_r(p+)=-1/2-r,
d_r(p-)=-1/2+r.                                        (23)
```

It contains at most one collision point.  In the pure `F3` row, the factor
`x` is a genuine coordinate but contains only `p0`; the other factor is the
noncoordinate Kummer surface `(21)`, which contains `p+` and `p-`.

Therefore no irreducible coordinate factor of `ell(F)` contains two points
of `(2)`.  In the divisibility language of THM-3546, there are no source and
target coordinates with the target coordinate affine-linear and

```text
rho_s divides ell(F)
```

such that `V(rho_s)` carries a known collision pair.

Combined with THM-3558, this proves a two-sided linear barrier: a descent of
the fixed three-dimensional collision to `JC(2)` must use a nonlinear source
coordinate and a nonlinear target coordinate simultaneously.

## 5. Exact verification

Reproduce with

```bash
python3 04-computation/jacobian_affine_target_pullback_factor_atlas_kps_s188.py
python3 -O 04-computation/jacobian_affine_target_pullback_factor_atlas_kps_s188.py
```

The ordinary and optimized transcripts agree.  The companion checks `(4)`,
`(5)`, the Laurent restriction `(9)`, the complete Groebner basis `(10)`,
both exceptional factorizations including their multiplicities, the pure
`F3` factorization, and the collision incidences `(23)` by exact arithmetic.

**QED.**
