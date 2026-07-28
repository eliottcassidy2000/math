---
id: THM-2702
title: "Polynomial-graph coordinate-projection Keller classification"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  On an arbitrary
  polynomial graph s_3=f(s_1,s_2), the three
  coordinate-pair projections of the THM-2696 S3<S4 quotient have a complete
  constant-Jacobian classification.  The (A,d) pair gives exactly the
  triangular family f=(kappa/2)x+H(y-x^2/2); the (B,d) pair has no nonzero-
  constant-Jacobian graph; and the (A,B) pair gives exactly
  f=xy/2-x^3/8-kappa/4.  Both surviving families are explicit polynomial
  automorphisms of A2.  Arbitrary polynomial surfaces, arbitrary linear or
  nonlinear target projections, general JC(2), and DC(2) remain open.
source: root/reflection-quotient-polynomial-graph-2026-07-28
depends_on:
  - THM-2696-reflection-completed-s4-relative-different-and-coordinate-invariant-jacobian-gate
  - THM-2699-affine-plane-linear-projection-keller-slice-classification
related:
  - THM-2655-quartic-keller-resolvent-v4-quasietale-torsor-and-kummer-class-group-gate
script: 04-computation/jacobian_s4_polynomial_graph_coordinate_pairs_thm2702.py
output: 05-knowledge/results/jacobian_s4_polynomial_graph_coordinate_pairs_thm2702.out
script_sha256: c7a4cf1ecdd58bbfa0946cd0064bfeb55ca5fbe106b1bc5fe30c2484b6bef040
output_sha256: 7fefddc476c54712439e6c6dd0e6389cc1a3862e765f506bf7f5d491e0768e3a
hash_basis: LF-normalized bytes
---

# THM-2702 -- all three coordinate-pair graph slices are triangular or empty

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2699 closes every affine source plane and every affine-linear target
projection of the reflection quotient.  Polynomial source graphs are already
an infinite nonlinear test family.  For the three natural coordinate-pair
projections, that family can still be classified completely without assuming
the planar Jacobian conjecture.

## 1. The three graph equations

Use THM-2696's quotient coordinates

```text
A=x^2-2y,             B=y^2-2xz,             d=z.        (1)
```

Let `f in C[x,y]` and restrict `(1)` to the polynomial graph

```text
Gamma_f={z=f(x,y)} ~= A2_(x,y).                          (2)
```

Direct differentiation gives

```text
Jac(A,d)= 2(f_x+x f_y),                                  (3)
Jac(B,d)=-2(y f_x+f f_y),                                (4)
Jac(A,B)=-4(x f_x+x^2 f_y+f-xy).                         (5)
```

Fix a desired constant `kappa in C^*`.

## 2. The complete `(A,d)` family

Set

```text
t=y-x^2/2,
D=partial_x+x partial_y.                                 (6)
```

The change `(x,y)<->(x,t)` is a polynomial automorphism, and in these
coordinates

```text
D=partial_x with t fixed.                                (7)
```

Equation `(3)=kappa` is therefore `D(f)=kappa/2`.  Its complete polynomial
solution is

```text
f=(kappa/2)x+H(t),             H in C[t] arbitrary.       (8)
```

On `(8)` the target map is

```text
(A,d)=(-2t,(kappa/2)x+H(t)).                              (9)
```

It is triangular, with inverse

```text
t=-A/2,
x=2(d-H(t))/kappa,
y=t+x^2/2.                                                (10)
```

Thus this infinite nonlinear escape consists entirely of polynomial
automorphisms.  THM-2699 Family I is its affine subfamily.

## 3. The `(B,d)` family is empty

Suppose `(4)=kappa` and put

```text
E=y partial_x+f partial_y,              c=-kappa/2.       (11)
```

Then `E(f)=c!=0`, and

```text
E(x)=y,       E^2(x)=f,       E^3(x)=c,       E^4(x)=0,
E(y)=f,       E^2(y)=c,       E^3(y)=0.                  (12)
```

Since `x,y` generate `C[x,y]`, `(12)` makes `E` a locally nilpotent
derivation.  Hence

```text
exp(T E):C[T,x,y] -> C[T,x,y]                            (13)
```

is a polynomial automorphism with inverse `exp(-T E)`.  Its Jacobian is a
unit of `C[T,x,y]`; it equals one at `T=0`, so it is identically one.  The
coefficient of `T` in that Jacobian is the divergence

```text
div(E)=f_y.                                               (14)
```

Therefore `f_y=0`.  But `(11)` then reads `y f_x=c`, impossible for a
nonzero constant `c`.  Consequently

```text
there is no polynomial graph with Jac(B,d) in C^*.        (15)
```

This is a structural locally-nilpotent-derivation obstruction, not a bounded
degree calculation.

## 4. The unique `(A,B)` family

Let

```text
L=x partial_x+x^2 partial_y.                              (16)
```

In the `(x,t)` coordinates of `(6)`, `L=x partial_x`.  Equation
`(5)=kappa` becomes

```text
(L+1)f=xy-kappa/4
       =xt+x^3/2-kappa/4.                                (17)
```

On every monomial `x^i t^j`, the operator `L+1` acts by the nonzero scalar
`i+1`.  It is therefore injective on `C[x,t]`, and `(17)` has the unique
polynomial solution

```text
f=xt/2+x^3/8-kappa/4
 =xy/2-x^3/8-kappa/4.                                   (18)
```

For `(18)`, the target map simplifies to

```text
(A,B)=(-2t,t^2+(kappa/2)x).                              (19)
```

Again it is triangular.  Its inverse is

```text
t=-A/2,
x=2(B-t^2)/kappa,
y=t+x^2/2.                                                (20)
```

Thus the cubic forced by the `i+1=4` eigenspace is not a counterexample; it
is the unique automorphic graph in this coordinate pair.

## 5. Boundary and consequence

Combining `(8)`, `(15)`, and `(18)` proves:

```text
every nonzero-constant-Jacobian coordinate-pair map on a graph z=f(x,y)
is one of the two displayed triangular automorphisms.                    (21)
```

The THM-2696 constant-different graph `z=xy-c` is a useful hostile: its three
coordinate-pair Jacobians have respective total degrees `2,3,3`, so it lies
in none of the positive families.  This does not contradict its constant
three-dimensional relative different.

The graph coordinate `z`, and the choice of one of the three literal target
coordinate pairs, are load-bearing.  A source surface not globally of the
form `(2)`, an arbitrary linear combination of target coordinates, or a
polynomial target projection has a different PDE.  No arbitrary `S4`
resolvent, `JC(2)`, or `DC(2)` conclusion follows.

## 6. Exact companion

Run

```text
python 04-computation/jacobian_s4_polynomial_graph_coordinate_pairs_thm2702.py
python -O 04-computation/jacobian_s4_polynomial_graph_coordinate_pairs_thm2702.py
```

Both modes must byte-match

```text
05-knowledge/results/jacobian_s4_polynomial_graph_coordinate_pairs_thm2702.out.
```

The companion uses explicit `require` checks.  It derives `(3)`--`(5)`,
checks nonlinear members of `(8)` and both polynomial inverses, diagonalizes
`L+1` through degree six, verifies the LND flow's divergence coefficient,
and checks the constant-different hostile.  The degree-six coefficient audit
has seven free `H(t)` coefficients, twenty forced non-slice coefficients, and
Euler spectrum `1,...,7`.  This finite control supports but does not replace
the all-degree arguments in Sections 2--4.

An independent hostile audit rederived `(3)`--`(5)`, checked the
`D=partial_x` normal form in `(x,t)`, the `L+1` spectrum, the LND generator
nilpotence and `exp(TE)` unit-Jacobian/divergence step, both inverses, graph
and target scope, normal/`-O`/stored replay, and both exact hashes.  No defect
was found.

QED.
