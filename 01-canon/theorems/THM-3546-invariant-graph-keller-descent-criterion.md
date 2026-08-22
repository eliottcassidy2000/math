---
id: THM-3546
title: "Invariant graph Keller descent criterion and collision transfer"
status: >
  PROVED / INDEPENDENTLY AUDITED WITH FINITE-FIELD SCOPE REPAIR.  A
  constant-Jacobian polynomial map in n+1 variables that carries one
  polynomial graph scheme-theoretically to another restricts to a
  constant-Jacobian map in n variables.  In graph coordinates the normal
  component is r times a polynomial; block-triangular differentiation on r=0
  factors the ambient unit Jacobian into tangential and normal factors,
  forcing both to be constants.  Any geometric collision on the source graph
  therefore descends.  The same holds for arbitrary coordinate hypersurfaces
  after polynomial left-right straightening.
source: kps-s184
depends_on: []
related:
  - THM-1350-equivariant-fixed-locus-JC
  - THM-3543-torus-quotient-ramification-square-no-go
---

# THM-3546 -- a collision on a polynomial graph descends one dimension

**PROVED / INDEPENDENTLY AUDITED WITH FINITE-FIELD SCOPE REPAIR.**  This gives
a precise way in which a higher-dimensional Keller counterexample could produce a lower-
dimensional one.  Unlike categorical invariant quotients, restriction to a
polynomial graph loses no tangent direction.  The ambient Jacobian unit then
forces both the tangent and normal factors to be units.

## 1. Graph descent theorem

Let `k` be a field and let

```text
F:A^(n+1)->A^(n+1),              det JF=c in k*.        (1)
```

Write source coordinates `(x,z)`, with `x=(x_1,...,x_n)`, and target
coordinates `(y,w)`.  Let `F=(F_y,F_w)`, `h in k[x]`, and `g in k[y]`.
Suppose the graph containment is scheme-theoretic, namely the identity

```text
F_w(x,h(x))=g(F_y(x,h(x)))              in k[x].       (2)
```

Equivalently, `F` carries the graph subscheme of `h` into that of `g`.
Containment on all geometric points also suffices; over an infinite field,
containment on all `k`-points implies `(2)`.

Define the induced map on the graph by

```text
f(x)=pr_y F(x,h(x)).                                   (3)
```

Then

```text
                         det Jf in k*.                  (4)
```

In particular, if two distinct points of `Graph(h)` over one common field
extension `K/k` have the same image under `F`, then `f_K` is a noninjective
Keller map in `n` variables.

### Proof

Straighten the two graphs with the triangular polynomial automorphisms

```text
Phi_h(x,r)=(x,h(x)+r),
Psi_g(y,w)=(y,w-g(y)).                                  (5)
```

Both have Jacobian determinant one.  Set

```text
F_tilde=Psi_g o F o Phi_h=(P(x,r),R(x,r)).              (6)
```

Condition `(2)` says `R(x,0)=0`.  Since `(r)` is a principal prime ideal in
`k[x,r]`, there is a polynomial `A(x,r)` such that

```text
R(x,r)=r A(x,r).                                       (7)
```

On `r=0`, the Jacobian matrix is block upper triangular:

```text
J F_tilde(x,0) = [ J_x P(x,0)    partial_r P(x,0) ]
                  [      0             A(x,0)      ].   (8)
```

Taking determinants in `(8)` and using `(1)` gives the polynomial identity

```text
c=det J_xP(x,0) * A(x,0)       in k[x].                (9)
```

The only units in `k[x]` are the nonzero constants.  Hence both factors in
`(9)` lie in `k*`.  Since `P(x,0)=f(x)`, this proves `(4)`.  If two points on
the source graph collide, their distinct `x`-coordinates collide under `f`,
proving the final assertion.  QED.

The proof simultaneously shows that the normal multiplier `A(x,0)` is a
nonzero constant; this does not need to be imposed as an extra hypothesis.

## 2. Coordinate-hypersurface form

The graph language is only a convenient chart.  Let `rho_s` and `rho_t` be
polynomial coordinates: each occurs as the last component of a polynomial
automorphism of its ambient affine space.  Suppose

```text
rho_t o F = a rho_s                                  (10)
```

for a polynomial `a` in the source ring, so `F` carries `V(rho_s)` into
`V(rho_t)`.  Choose source and target automorphisms
`alpha_s=(u,rho_s)` and `alpha_t=(v,rho_t)`.  The left-right straightening
`alpha_t o F o alpha_s^(-1)` reduces `(10)` to `(7)`; its determinant is the
ambient constant multiplied by the constant ratio of the two automorphism
Jacobians.  Therefore the restriction

```text
F:V(rho_s) ~= A^n  ->  V(rho_t) ~= A^n               (11)
```

is Keller.  Any collision in `V(rho_s)` descends.

For a proposed descent of the fixed three-variable counterexample, this turns
the vague request for a special surface into a four-gate test for each
proposed pair:

1. choose a source coordinate polynomial `rho_s` vanishing at two known
   colliding preimages;
2. choose a target coordinate polynomial `rho_t` vanishing at their common
   image;
3. solve the exact divisibility `rho_s | rho_t(F)`; and
4. verify that both polynomials are coordinates, rather than merely smooth or
   irreducible.

If these four gates pass, planar noninjectivity and constant Jacobian follow
from the theorem; they do not require a second determinant search.

## 3. Relation to the current fixed map

THM-3543 shows that the categorical torus quotient of the fixed
three-dimensional map contracts an entire divisor and acquires the Jacobian
square `2(2-3v-t)^2`.  That operation is a quotient, not a graph restriction,
so `(9)` does not apply.  The present theorem identifies the missing
coordinate: a successful descent must retain a transverse normal coordinate
until it has produced an actual coordinate hypersurface.

THM-1350 is the linear fixed-locus analogue.  There the equivariance makes the
Jacobian block diagonal along a fixed affine subspace; the same unit-product
argument forces its tangential determinant to be constant.  Here no group
action is required, and the invariant subspace is replaced by a possibly
nonlinear polynomial graph.

For the explicit collision in THM-1300, a graph candidate must at least obey

```text
h(0,0)=-1/4,
h(1,-3/2)=h(-1,3/2)=13/2.                              (12)
```

Interpolation makes `(12)` cheap; the load-bearing condition is the global
graph identity `(2)`, or equivalently the divisibility `(10)`.  Low-degree
exact searches can therefore falsify a sharply defined lane without saying
anything about arbitrary coordinate hypersurfaces.

## 4. Scope

This theorem is a descent criterion, not a planar counterexample.  It does
not assert that the fixed map has any invariant graph, that a smooth
hypersurface is a polynomial coordinate, or that a rational/local graph can
be polynomially straightened.  Failure through a bounded degree does not
exclude higher-degree graphs, nongraph coordinate hypersurfaces, or other
three-dimensional Keller maps.

The finite-field qualifier is essential if containment is stated set-
theoretically.  Over `F_q`, the Keller map

```text
(x,z) -> (z,x^q-x)
```

sends every `F_q`-point of the graph `z=0` to that graph, but the induced map
on rational points is constant; the polynomial identity `(2)` fails.  The
independent audit supplied this hostile and verified the repaired proof, the
coordinate-hypersurface form, and the absence of any characteristic
restriction under the scheme-theoretic hypothesis.
