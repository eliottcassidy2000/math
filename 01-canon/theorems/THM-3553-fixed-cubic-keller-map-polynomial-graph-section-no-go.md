---
id: THM-3553
title: "Fixed cubic Keller map polynomial graph-section no-go"
status: >
  PROVED / INDEPENDENTLY HOSTILE-AUDITED.  For the displayed coordinates of
  the fixed THM-1300 cubic Keller map, restriction to every polynomial source
  graph z=h(x,y) has a nonconstant tangential Jacobian.  The nonzero top row
  is explicit in every degree, including constant h.  Hence no such graph is
  carried scheme-theoretically into any polynomial target graph in the
  displayed target coordinates.  This does not exclude graphs after
  nonlinear ambient coordinate changes or nongraph coordinate hypersurfaces.
source: root-2026-08-18-planar-jacobian-counterexample-hostiles
depends_on:
  - THM-1300-jacobian-counterexample-dixmier-A3-explicit
  - THM-3546-invariant-graph-keller-descent-criterion
related:
  - THM-3543-torus-quotient-ramification-square-no-go
  - THM-3554-punctured-kummer-collision-surface-normal-form
---

# THM-3553 -- fixed cubic Keller map polynomial graph-section no-go

**PROVED / INDEPENDENTLY HOSTILE-AUDITED.**  Let `F=(F_1,F_2,F_3)` be the
fixed THM-1300 map.  In its displayed coordinates, put `u=1+xy`, so

```text
F_1=u^3 z+y^2u(4+3xy),
F_2=y+3xu^2z+3xy^2(4+3xy),
F_3=2x-3x^2y-x^3z.                                    (1)
```

For every `h in C[x,y]`, define the restriction to the source graph

```text
p(x,y)=F_1(x,y,h(x,y)),
q(x,y)=F_2(x,y,h(x,y)).                                (2)
```

Then `Jac(p,q)` is nonconstant.  Consequently, there are no polynomials
`h,g in C[x,y]` satisfying

```text
F_3(x,y,h(x,y))=g(p(x,y),q(x,y)).                      (3)
```

Equivalently, no polynomial source graph `z=h(x,y)` is carried
scheme-theoretically into a polynomial target graph
`Y_3=g(Y_1,Y_2)` in these displayed coordinates.

## 1. Nonconstant source graph

Let `D=deg h>=1` and let `H=h_D` be its nonzero top homogeneous form.  The
terms involving `h` dominate the fixed terms in `(1)`, so

```text
p_(D+6)=x^3y^3H,
q_(D+5)=3x^3y^2H.                                     (4)
```

Put `K=x^3y^2H`.  Thus the two forms in `(4)` are `yK` and `3K`, and direct
differentiation gives

```text
Jac(yK,3K)=-3K K_x
             =-3x^5y^4 H(3H+xH_x).                   (5)
```

The right side is nonzero.  Indeed, on the degree-`D` binary forms the
operator

```text
3+x partial_x                                             (6)
```

sends `x^i y^(D-i)` to `(i+3)x^i y^(D-i)`.  It is injective in
characteristic zero, so `3H+xH_x` cannot vanish.  Since `C[x,y]` is a
domain, `(5)` is nonzero.

The bracket in `(5)` is the unique possible top homogeneous row of
`Jac(p,q)`, of degree

```text
(D+6)+(D+5)-2=2D+9.                                    (7)
```

All brackets involving a lower homogeneous row have smaller degree.
Therefore lower terms cannot cancel `(5)`, and `Jac(p,q)` is nonconstant.

## 2. Constant source graph, including zero

Let `h=c in C`.  The degree-six and degree-five top forms now include the
fixed terms of `(1)`:

```text
p_6=x^2y^3(cx+3y),
q_5=3x^2y^2(cx+3y).                                   (8)
```

With `K=x^2y^2(cx+3y)`, these are again `yK` and `3K`.  Hence

```text
Jac(p_6,q_5)
  =-3K K_x
  =-9x^3y^4(cx+3y)(cx+2y),                            (9)
```

which is nonzero for every `c`.  At the potentially hostile value `c=0`, it
is `-54x^3y^6`.  Thus `Jac(p,q)` has a nonzero degree-nine top row and is
again nonconstant.

## 3. Graph descent contradiction

If `(3)` held, the graph-containment identity in THM-3546 would apply to the
ambient Keller map `F`.  That theorem forces the tangential determinant of
the induced map `(p,q)` to lie in `C*`.  Sections 1--2 show instead that it
has a nonzero positive-degree top row.  This contradiction proves the graph
no-go.

## 4. Scope and hostile boundary

The source and target graph orientations in the theorem are load-bearing.
It does not rule out a graph after nonlinear polynomial source or target
coordinate changes, a nongraph coordinate hypersurface, a rational graph,
or a graph for a different three-dimensional Keller map.

The plane `x=0` is the canonical hostile against a broader claim.  It is not
a graph `z=h(x,y)` over the displayed `(x,y)`-plane, and `(1)` maps it into
`F_3=0`; the induced map

```text
(y,z) -> (z+4y^2,y)                                    (10)
```

is a polynomial automorphism.  Thus arbitrary coordinate-hypersurface
descent remains a live architecture even though every displayed polynomial
graph is now excluded.
