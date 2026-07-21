---
id: THM-2042
title: "The smooth non-coordinate R=x(2-3xq) has an explicit rank-two Hamiltonian mate"
status: >
  PROVED CORE / completion open. In the canonical rank-two Poisson algebra
  C[x,q,p,z] with {p,x}={z,q}=1, the displayed polynomials R and D below
  satisfy {D,R}=1. The plane polynomial R has no critical point but is not a
  coordinate. An explicit tangent momentum L and two joint invariants I,J are
  computed exactly; their nonconstant bracket identifies why the naive
  symplectic completion fails. The existence and quantization of a second
  canonical pair remain HYP-8801.
source: codex-2026-07-21-DC2-JC2
related:
  - THM-1300
  - THM-1345
  - HYP-8801
---

# THM-2042 -- the first rank-two canonical-pair scaffold

Work in

```text
P_2=C[x,q,p,z],             {p,x}={z,q}=1,
```

with all other brackets between generators zero. Put

```text
R=x(2-3xq),
D=((1+3xq)/2)p-3q^2z.                              (1)
```

Then

```text
{D,R}=1.                                           (2)
```

Indeed

```text
R_x=2-6xq,       R_q=-3x^2,
((1+3xq)/2)R_x+(-3q^2)R_q=1.                      (3)
```

Thus (2) is the Bezout identity for the gradient row of `R`.

The polynomial `R` has no critical point by (3), but it is not a coordinate
of `C[x,q]`: its zero fibre

```text
R^{-1}(0)={x=0} union {2-3xq=0}
```

is reducible, whereas every fibre of a polynomial coordinate is isomorphic to
the irreducible affine line. Hence stabilization by the momentum variables has
created a canonical pair from a genuinely non-coordinate plane polynomial.

Define the tangent momentum and the auxiliary base factor

```text
L=3x^2p+(2-6xq)z,       A=3xq-1.                   (4)
```

Direct calculation gives

```text
{R,L}=0,             {D,L}=(9/2)qL,
{D,q}=-3q^2,         {D,A}=-(3/2)qA.               (5)
```

Consequently the two polynomials

```text
I=q^3L^2,                 J=A^3L                   (6)
```

commute with both `R` and `D`. They do not yet form a canonical pair:

```text
{J,I}=-6q^2A^2L^2.                                 (7)
```

Equations (5)-(7) are the exact first centralizer scaffold. They show that the
obstruction is not internal cancellation in (2): it is the nonconstant
symplectic density left on the joint centralizer. A successful rank-two
Poisson counterexample must replace or correct `(I,J)` so that this density is
one while retaining polynomiality and non-surjectivity.

This theorem does not assert that the full four-polynomial construction in the
owner-supplied abstract has been reconstructed; only the displayed `R`-slice
and its forced elementary consequences are used.
