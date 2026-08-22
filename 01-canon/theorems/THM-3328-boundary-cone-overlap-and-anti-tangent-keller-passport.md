---
id: THM-3328
title: "Boundary-cone overlap and anti-tangent Keller passport"
status: >
  PROVED + INDEPENDENTLY AUDITED.  For a C1 local diffeomorphism on a
  polygonal domain, two labelled boundary preimages whose derivative-image
  inward cones overlap force a nonempty open target set with two distinct
  interior inverse branches.  At two smooth positively oriented boundary
  points, failure of such local overlap forces the image tangents to be
  anti-parallel.  For a degree-d polynomial triangle map this gives an exact
  collision/tangency/sign passport and, under explicit component-freeness
  hypotheses, at most 6d^2-6d+3 complex edge-pair collision parameters.
  Combined with THM-3303, a constant-J HFC(3) candidate has an open
  double-sheet cell, a smooth anti-tangent collision, or a distinct-source
  vertex collision with disjoint transformed inward cones.  HFC(3) and JC(2)
  remain OPEN.
audit: >
  Independent proof audits checked the uniform inverse-chart expansion,
  oriented half-plane sign under both Jacobian orientations, the
  anti-parallel iff boundary, the vertex-cone type, same-edge saturation by
  t-s, the conditional Bezout count, and the slit-annulus hostile control.
  No computation is used.
source: root/creative-passports-ii/2026-08-03
depends_on:
  - THM-3303-keller-simplex-null-moments-force-a-boundary-collision
related:
  - THM-3018-factorial-conjecture-as-a-simplex-moment-problem
  - THM-3319-exceptional-quadratic-two-clutch-etale-persistence
---

# THM-3328 -- boundary-cone overlap and anti-tangent Keller passport

**PROVED + INDEPENDENTLY AUDITED.**

THM-3303 forces a boundary collision in the constant-Jacobian homogeneous
factorial sector but leaves open whether that collision creates overlap or a
one-sheeted cut.  This theorem resolves the local alternative and turns the
exception into an exact edge passport.

## 1. Local cone-overlap lemma

Let `Omega` be a connected polygonal domain in `R^2`, let `F` be `C^1` on a
neighbourhood of `closure(Omega)`, and suppose `det DF` never vanishes there.
For `x in boundary(Omega)`, define the open inward tangent cone

```text
C_x={v: x+epsilon*v lies in Omega for every sufficiently small epsilon>0}.
                                                                    (1)
```

Suppose two labelled boundary points satisfy

```text
x!=y,                         F(x)=F(y)=w.              (2)
```

**Cone-overlap lemma.**  If

```text
DF(x)C_x intersects DF(y)C_y,                            (3)
```

then there is a nonempty open set `U`, with `w` in its closure, such that
every `z in U` has two distinct preimages in `Omega`, one near `x` and one
near `y`.

Choose a target direction in `(3)` and a small closed angular sector `K`
whose directions remain a positive distance inside both transformed cones.
The inverse-function theorem gives disjoint inverse charts at `x` and `y`.
Uniformly for unit `v in K`, they have expansions

```text
F_x^(-1)(w+epsilon*v)
 =x+epsilon*DF(x)^(-1)v+o(epsilon),
F_y^(-1)(w+epsilon*v)
 =y+epsilon*DF(y)^(-1)v+o(epsilon).                       (4)
```

The leading directions lie strictly inside `C_x` and `C_y`.  For all small
positive `epsilon`, both inverse images in `(4)` are therefore interior
points.  The interior of the resulting truncated target sector is the
required `U`.  The disjoint source charts make its two inverse branches
distinct.

This is a local theorem: it uses neither polynomiality, moments, nor global
properness.

## 2. Smooth collisions and the anti-tangent exception

Orient every smooth boundary component so that `Omega` lies to its left.  At
two nonvertex collision points, let `tau_x,tau_y` be the oriented boundary
tangents and put

```text
u_x=DF(x)tau_x,                    u_y=DF(y)tau_y.          (5)
```

Both vectors are nonzero.  Since `det DF` has constant sign `sigma` on the
connected domain, the transformed inward half-planes are

```text
H_x={z:sigma*det(u_x,z)>0},
H_y={z:sigma*det(u_y,z)>0}.                                (6)
```

Two open oriented half-planes through the origin are disjoint exactly when

```text
u_y=-lambda*u_x                    for some lambda>0.       (7)
```

Indeed, nonparallel half-planes meet in an open wedge, positive-parallel
half-planes coincide, and negative-parallel tangents select opposite sides.
Combining `(7)` with the cone-overlap lemma proves:

> A smooth boundary collision which creates no open target set of
> multiplicity at least two must have anti-parallel oriented image tangents.

Condition `(7)` is necessary, not sufficient, for a collision to remain
one-sheeted.  Higher normal jets may still make two anti-tangent branches
overlap.

## 3. Exact triangle passport

Let `F=(p,q)` be a real polynomial map of degree at most `d` with constant
nonzero Jacobian.  Parametrize the three triangle edges counterclockwise by

```text
gamma_i(t)=a_i+t*v_i,                    0<=t<=1.           (8)
```

For two smooth edge parameters define

```text
E_ij(t,s)=F(gamma_i(t))-F(gamma_j(s)),
T_ij(t,s)=det(DF(gamma_i(t))v_i,DF(gamma_j(s))v_j),
A_ij(t,s)=dot(DF(gamma_i(t))v_i,DF(gamma_j(s))v_j).        (9)
```

Every smooth collision compatible with the absence of an open double-sheet
cell satisfies the exact semialgebraic system

```text
E_ij(t,s)=0,                    T_ij(t,s)=0,
A_ij(t,s)<0.                                             (10)
```

The derivative vectors are nonzero, so the last two conditions in `(10)` are
equivalent to anti-parallelism.  For `i=j`, the diagonal `t=s` is not a
collision of distinct source points and must be saturated away.  For
`i!=j`, the parameter pair representing the edges' common endpoint is the same
source point and must be discarded.  A distinct-source collision involving a
vertex is a cone case rather than a smooth two-tangent case.

At a vertex the lawful object is the full inward wedge generated by the two
incident edge rays.  If its derivative image intersects the other derivative-
image inward cone, Section 1 again gives an open double-sheet cell.  Only
disjoint transformed cone interiors survive in the one-sheet passport.

## 4. Conditional finite collision bound

For one unordered pair of distinct edges, the two coordinate equations in
`E_ij=0` have total degree at most `d`.  If they have no common component,
Bezout gives at most `d^2` complex projective solutions, hence at most that
many affine collision parameters.

On one edge, both raw coordinate differences are divisible by `t-s`.  Divide
by that factor and assume the two resulting degree-at-most-`d-1` polynomials
have no common component after saturation.  There are then at most
`(d-1)^2` off-diagonal complex collision parameters on that edge.  Under
these component-freeness hypotheses for all three cross-edge and all three
same-edge systems,

```text
number of complex edge-pair collision parameters
 <=3d^2+3(d-1)^2
 =6d^2-6d+3.                                             (11)
```

This is an upper bound on parameter solutions, with multiplicity; it need not
be sharp and may count shared vertices.  A common collision component is not
silently declared finite.  It is a separate positive-dimensional complex
collision stratum; on its distinct-source real smooth locus the
cone/anti-tangent test applies pointwise.

## 5. Constant-J HFC(3) consequence

Let `g=p+iq` be the triangle restriction of a homogeneous `HFC(3)` candidate
and assume `det D(p,q)` is a nonzero constant.  THM-3303 proves that
`F=(p,q)` has a boundary collision.  Apply the present theorem hierarchically.

1. Either some collision has overlapping transformed inward cones, and the
   multiplicity function `N_F` is at least two on a nonempty open target cell.
2. Or no such open cell exists.  Then a collision between two smooth edge
   interiors obeys `(10)`, while any collision involving a vertex has
   disjoint transformed inward cone interiors.

Thus every candidate in this sector has an open double-sheet cell, a smooth
anti-tangent collision, or a distinct-source vertex collision with disjoint
transformed inward cones.  The smooth alternative is encoded by `(10)`; the
vertex alternative is checked by finitely many full-cone comparisons.  This
refines the residual; it does not exclude it.

## 6. Sharp hostile control

The slit-annulus map realizes the exceptional orientation.  On the positively
oriented rectangle `[0,2*pi]x[0,1]`, put

```text
Phi(t,s)=sqrt(rho^2+s)*(cos(t),sin(t)),                    (12)
```

where `rho>0`.  Its Jacobian is the nonzero constant `-1/2`.  The vertical
sides `t=0` and `t=2*pi` map to the same radial seam, but their positively
oriented image tangents are anti-parallel and their transformed inward sides
are opposite.  The interior maps one-to-one onto a slit annulus, so the seam
does not force a local double-sheet cell.

The control is nonpolynomial and its source is not a triangle.  It proves only
that the anti-parallel exception cannot be deleted from the local theorem.
Nothing here proves `HFC(3)`, `JC(2)`, a mate, or a polynomial inverse.

QED.
