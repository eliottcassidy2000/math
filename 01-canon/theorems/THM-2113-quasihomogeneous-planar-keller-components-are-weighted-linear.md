---
id: THM-2113
title: "Quasi-homogeneous planar Keller components are weighted-linear coordinates"
status: >
  PROVED. Let w_1,w_2 be positive weights and let f in C[x,y] be a
  nonconstant w-quasi-homogeneous polynomial. If there is g in C[x,y] with
  Jac(f,g)=1, then, after possibly interchanging x and y,
  f=a*x+b*y^q with a nonzero. Conversely every such polynomial is a
  coordinate and has an explicit Jacobian mate. Thus the planar Jacobian
  conjecture holds for a component supported on a single positive-weight
  face. The proof is a direct weighted-degree argument; it does not assume
  that a mate makes the generic fiber C, does not use flow completeness, and
  does not prove the nonhomogeneous one-principal-face descent.
source: codex-2026-07-22-JC2-quasihomogeneous-repair
depends_on: []
related:
  - THM-2045
  - THM-2063
  - THM-2071
  - HYP-8950
  - HYP-8955
---

# THM-2113 -- the quasi-homogeneous planar Keller base case

Let `w=(w_1,w_2)` with `w_1,w_2` positive integers. Give `x` weight
`w_1` and `y` weight `w_2`. Suppose that the nonconstant polynomial

```text
f in C[x,y]
```

is weighted homogeneous of degree `delta`, and that

```text
Jac(f,g)=f_x g_y-f_y g_x=1                         (1)
```

for some polynomial `g`.

We prove that, after possibly swapping the variables,

```text
f=a*x+b*y^q,          a!=0,                         (2)
```

where `q>=1` when `b!=0`. In particular `f` is a coordinate: the triangular
map `(f,y)` has constant Jacobian `a` and polynomial inverse.

## 1. The bracket selects one weighted piece

Write the finite weighted-homogeneous decomposition

```text
g=sum_e g_e.
```

For every `e`, either `Jac(f,g_e)=0` or it is weighted homogeneous of degree

```text
delta+e-w_1-w_2.                                    (3)
```

Distinct weighted degrees cannot cancel. Since the right side of (1) has
degree zero, the component

```text
e_0=w_1+w_2-delta                                   (4)
```

must exist and satisfy

```text
Jac(f,g_(e_0))=1.                                   (5)
```

In particular `e_0>=0`. Equality is impossible: with positive weights every
weight-zero polynomial is constant, whose bracket with `f` is zero. Hence

```text
delta<w_1+w_2.                                      (6)
```

This is the exact weight obstruction. No assertion about the topology of a
generic fiber has entered.

## 2. A degree below the weight sum has only axial monomials

Every mixed monomial `x^i y^j` with `i,j>=1` has weight at least
`w_1+w_2`, so (6) excludes it from `f`. At a fixed positive weighted degree
there is at most one pure power on each axis. Thus

```text
f=A*x^p+B*y^q,                                      (7)
```

where a missing term is allowed and every displayed coefficient is nonzero.

If both terms occur and `p,q>=2`, then

```text
delta=p*w_1=q*w_2<w_1+w_2.
```

Consequently

```text
(p-1)w_1<w_2,       (q-1)w_2<w_1.
```

Multiplication gives `(p-1)(q-1)<1`, impossible for integers `p,q>=2`.
Therefore `p=1` or `q=1`, and (7) has the form (2) after a swap.

If only `A*x^p` occurs, (5) reads

```text
A*p*x^(p-1) * partial_y g_(e_0)=1.
```

The polynomial factor `x^(p-1)` must be a unit, so `p=1`. The one-term
`B*y^q` case is identical. This also covers the axial degenerations without
an appeal to primitivity or places at infinity.

## 3. Converse and exact scope

For `f=a*x+b*y^q` with `a!=0`, take `g=y/a`; then `Jac(f,g)=1`. For
`f=a*y+b*x^p` with `a!=0`, take `g=-x/a`. Thus (2) is both necessary and
sufficient.

The argument proves the zero-descent-step, single-positive-weight-face case.
For a polynomial

```text
f=Phi_w + lower w-degree terms,
```

equation (3) couples different pieces of `f` and no longer isolates the
principal face by itself. Showing that a nonlinear principal face cannot be
cancelled by lower faces, or that a weighted-linear principal face can always
be removed by a terminating triangular descent, is the remaining
single-principal-face problem. That is a genuine planar-Jacobian obstruction;
THM-2113 does not supply the missing properness or descent theorem.

QED.
