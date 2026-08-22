---
id: THM-3303
title: "Keller simplex null moments force a boundary collision"
status: >
  PROVED + INDEPENDENTLY AUDITED.  Let Delta be a nondegenerate closed real
  triangle and let g=p+iq be a complex polynomial on Delta.  If every positive
  area moment of g vanishes and det D(p,q) is a nonzero constant, then the
  Keller map (p,q) is not injective even on the boundary of Delta.  Therefore
  any HFC(3) counterexample in this constant-Jacobian sector is already a real
  noninjective planar Keller map and hence a counterexample to JC(2).
  Equivalently, JC(2) implies HFC(3) on this sector.  This proves neither
  conjecture: outside the boundary-embedded case the joint boundary and sheet-
  multiplicity passport is load-bearing.
source: root/cross-frontier-35-81/2026-08-03
audit: >
  An independent proof audit checked the boundary-degree step, the location of
  the quadrature node, Bergman-kernel covariance, the polynomial-circle edge
  obstruction, and the HFC/JC quantifiers.  It strengthened the initial
  injectivity hypothesis to boundary injectivity and identified the residual
  boundary/multiplicity passport.  No computation is used.
depends_on:
  - THM-3018-factorial-conjecture-as-a-simplex-moment-problem
related:
  - THM-3300-factorial-gaussian-torus-bridge-and-the-archimedes-no-go
  - THM-3328-boundary-cone-overlap-and-anti-tangent-keller-passport
external:
  - "S. R. Bell, The Bergman kernel and quadrature domains in the plane, arXiv:math/0403371 (quadrature-domain background only; the order-one argument is given here)."
  - "M. Fleeman and E. Lundberg, The Bergman analytic content of planar domains, arXiv:1602.03615 (one-point simply connected quadrature-domain background only)."
---

# THM-3303 -- Keller simplex null moments force a boundary collision

**PROVED + INDEPENDENTLY AUDITED.**

This is a typed one-way bridge between the homogeneous Factorial Conjecture in
three variables and the planar Jacobian Conjecture.  Its point is not that the
two conjectures are equivalent.  It isolates a natural HFC(3) sector in which
every would-be counterexample must already contain a planar Keller collision,
and it locates that collision on the boundary of the simplex.

The disk/edge mechanism was anticipated, without a canon proof, in
`07-reflections/fc3-cap-the-decisive-step-null-quadrature-and-the-edge-obstruction-opus-S4.md`.
The constant-Jacobian and boundary-injectivity hypotheses close precisely that
reflection's uniform one-sheet case; they do not close its overlap residual.

## 1. Statement

Let `Delta` be a closed nondegenerate triangle in `R^2`, and write

```text
g=p+iq in C[x,y],             F=(p,q):R^2 -> R^2,       (1)
```

where `p,q` have real coefficients.  Assume

```text
int_Delta g^m dA=0                  for every m>=1,      (2)
det DF=c in R\{0}.                                        (3)
```

Then `F|partial(Delta)` is not injective.  In particular, there are distinct
real boundary points `a,b` with `F(a)=F(b)`, so the complexification
`(p,q):C^2 -> C^2` is a noninjective polynomial map with constant nonzero
Jacobian.

**HFC(3)--JC(2) corollary.**  For a homogeneous
`f in C[X,Y,Z]`, put

```text
g(x,y)=f(x,y,1-x-y).                                  (4)
```

THM-3018's proved polar identity says that `L(f^m)=0` for every `m>=1` if and
only if `(2)` holds on the standard triangle.  Hence a counterexample to
HFC(3) for which `det D(Re g,Im g)` is a nonzero constant would be a
counterexample to `JC(2)`.  Equivalently,

```text
JC(2)  =>  HFC(3) on the constant-nonzero-Jacobian sector.             (5)
```

## 2. Boundary injection would give one sheet

Suppose for contradiction that `F|partial(Delta)` is injective.  Its image
`Gamma` is a Jordan curve.  For `w notin Gamma`, the Brouwer degree of `F` on
`int(Delta)` is the winding number of `Gamma` about `w`, equal to `+1` or `-1`
on the bounded side and zero on the unbounded side.  Nonzero degree supplies a
preimage on the bounded side.  Every preimage has local degree `sign(c)`, so
the winding sign must equal `sign(c)` and the number of preimages is exactly
one there; it is zero on the unbounded side.

No interior point can map to `Gamma`: the local diffeomorphism `F` is open, so
an interior neighbourhood of such a point would map across both sides of
`Gamma`, producing a preimage on the side whose degree is zero.  It follows
that

```text
F:Delta -> closure(Omega),       Omega=F(int(Delta)),                  (6)
```

is a homeomorphism and `Omega` is a bounded Jordan domain.  Its boundary is
the union of three regular polynomial arcs.

## 3. The image is a one-node quadrature domain

The change-of-variables formula and `(2)` give

```text
int_Omega z^m dA(z)=|c| int_Delta g^m dA=0,       m>=1.                (7)
```

Writing `A=area(Omega)>0`, linearity gives the one-node identity

```text
int_Omega r(z)dA(z)=A r(0)        for every holomorphic polynomial r.  (8)
```

The node `0` really lies in `Omega`; it is not being inserted formally.  For
`w` outside `closure(Omega)`, set

```text
C(w)=int_Omega dA(z)/(w-z).                                      (9)
```

Expanding at infinity and using `(7)` gives `C(w)=A/w` there.  The complement
of a closed Jordan domain is connected, so analytic continuation gives that
identity throughout the exterior where `w!=0`.  If `0` were outside or on the
boundary, take exterior points `w -> 0`.  The elementary estimate

```text
|C(w)| <= int_Omega dA(z)/|w-z| <= 2*pi*R                         (10)
```

holds with one fixed `R` for all such nearby `w`, because `Omega` is bounded;
but `|A/w| -> infinity`.  Thus `0 in Omega`.

Holomorphic polynomials are dense in the Bergman space `A^2(Omega)` for a
bounded Jordan domain (one may exhaust by Riemann-map dilates and apply
Mergelyan approximation).  Both sides of `(8)` are continuous there, so `(8)`
extends to every `h in A^2(Omega)`.  Therefore the Bergman kernel satisfies

```text
K_Omega(z,0)=1/A.                                                    (11)
```

Choose a Riemann map `phi:D -> Omega` with `phi(0)=0`.  Kernel covariance and
`K_D(w,0)=1/pi` give

```text
(1/A) phi'(w) conjugate(phi'(0))=1/pi.                              (12)
```

Hence `phi'` is constant, so `Omega` is a disk centred at zero.

## 4. A polynomial triangle edge cannot trace the circle

The homeomorphism `(6)` sends each triangle edge into the circle
`partial(Omega)`.  Parametrize one edge affinely and write its restriction as
`h(t) in C[t]`.  If the radius is `rho`, then

```text
|h(t)|=rho                    for every real t in [0,1].             (13)
```

Let `h*(t)=conjugate(h(conjugate(t)))`.  Equation `(13)` makes
`h(t)h*(t)-rho^2` vanish on an interval and hence identically.  If
`deg h=d>0`, the product has degree `2d` and leading coefficient equal to the
squared modulus of the leading coefficient of `h`, impossible for a constant.
Thus `h` is constant.  This contradicts boundary injectivity (and also the
invertibility of `DF` on the nonzero edge tangent), proving the theorem.

## 5. The surviving boundary/multiplicity passport

Without boundary injectivity, the area formula does not give a uniform domain.
It gives the integer multiplicity function

```text
N_F(z)=#{a in int(Delta):F(a)=z},
int_C z^m N_F(z)dA(z)=|c| int_Delta g^m dA=0,       m>=1.             (14)
```

`N_F` is locally constant away from the three polynomial boundary arcs.  What
the proof excludes is the **joint** one-sheet Jordan case supplied by a
boundary embedding.  A boundary collision need not force `N_F>=2`: an
interior-one-sheet survivor could instead have a cut, a hole, or a boundary
identification.  Other survivors can have genuine overlap cells.  The honest
sidecar is therefore the boundary-incidence passport together with `N_F`, not
the multiplicity function alone.

THM-3328 subsequently makes this local alternative exact: derivative-image
inward-cone overlap forces an open double-sheet cell, while a smooth
one-sheeted survivor must be anti-tangent and vertex survivors retain the
full transformed cone.

The cheapest next problem is to classify the winding/multiplicity cells cut
out by the three edge images, including shared-cut and hole topologies, and
test `(14)` on that passport.  Forgetting the boundary data and rerunning the
disk argument would silently assume the case just excluded.

A smooth hostile control shows why this distinction is real.  For `rho>0`, on
the rectangle `[0,2*pi]x[0,1]`, the local diffeomorphism

```text
(t,s) |-> sqrt(rho^2+s) (cos(t),sin(t))                         (15)
```

has constant Jacobian `-1/2`, identifies its two vertical boundary edges, and
is one-sheeted on the interior of a slit annulus.  Its multiplicity is one
almost everywhere on the annulus, whose positive holomorphic moments all
vanish by rotation.  This is not polynomial and its source is not a triangle,
so it is not a counterexample to the theorem.  It is the sharp control showing
that the area formula alone cannot replace boundary topology by an overlap
claim.

## 6. Scope and failure boundary

Proved: `(1)`--`(3)` force a collision already on `partial(Delta)`, and the
HFC(3)--JC(2) implication `(5)` holds in the stated sector.

Not proved: `HFC(3)`, full `FC(3)`, `JC(2)`, or existence of any polynomial
satisfying `(2)`--`(3)`.  Candidates with zero or nonconstant real Jacobian lie
outside the theorem.  Candidates with constant nonzero Jacobian remain
possible only by being noninjective, in which case the boundary/multiplicity
passport surrounding `(14)` is mandatory.

The one-node quadrature viewpoint is classical; Bell and Fleeman--Lundberg,
cited in the front matter, are background for that language.  The proof above
supplies the order-one disk conclusion directly from Bergman-kernel covariance
and does not import a classification theorem.

**QED.**
