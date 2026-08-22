---
id: THM-3685
title: "Weighted-quartic native polynomial-graph descent no-go"
status: >
  PROVED + VERIFIED-EXACT.  Let G=(A,B,C) be the explicit degree-four
  determinant-minus-six weighted Keller map of THM-3438.  For every
  polynomial h(x,y), the restriction (A(x,y,h),B(x,y,h)) has nonconstant
  Jacobian.  If h has positive degree d and top form H, its Jacobian has
  nonzero top form -13122 W partial_x W, where
  W=x^6 y^3 H^2, of degree 4d+17.  For constant h=c the same formula holds
  with W=x^4 y^3(4y+cx)^2 and degree seventeen.  Consequently no source
  graph z=h can be carried by G into any target graph C=g(A,B): triangular
  graph descent would factor the ambient unit Jacobian through the
  nonunit tangential Jacobian.  This includes the graph h=1-x through
  G's known collision pair.  The result closes the native (A,B)-graph
  descent of this explicit threefold counterexample, not nonlinear source
  coordinate hypersurfaces, other target coordinate pairs, arbitrary
  quartic C3 boundary data, or JC(2).
source: codex-jc-quartic-c3-construction-2026-08-22
depends_on:
  - THM-3438-weighted-lift-keller-degree-spectrum
  - THM-3546-invariant-graph-keller-descent-criterion
related:
  - THM-3059-quartic-twojet-even-jelonek-c3-escape-counterexample
  - THM-3276-derivative-normalized-inverse-different-line-and-all-finite-jet-hostile
  - THM-3441-weighted-quartic-jelonek-components-and-inertia-parity
  - THM-3553-fixed-threedimensional-keller-map-polynomial-graph-section-no-go
  - THM-3559-affine-target-coordinate-pullback-no-go
script: 04-computation/jacobian_weighted_quartic_native_graph_no_go_thm3685.py
output: 05-knowledge/results/jacobian_weighted_quartic_native_graph_no_go_thm3685.out
script_sha256: b2b92fdaa66993eefd955d41c1d608e3775248558eb2d1aef3fd49092ec07d79
output_sha256: 2045b388adad86ae5cffc4af9ee0fbd41452bfbe5272ace7365d7b27d5fe317b
hash_basis: LF-normalized bytes
---

# THM-3685 -- the weighted quartic collision cannot descend on a native polynomial graph

**PROVED + VERIFIED-EXACT.**

This theorem attacks a concrete counterexample construction rather than an
abstract quartic invariant.  The degree-four map in THM-3438 is already a
noninjective polynomial Keller map in three variables.  Its collision points
lie on an extremely simple polynomial graph.  The question is whether some
polynomial graph is invariant enough to restrict that map to two variables.

The answer is no for the native first-two-coordinate projection, in every
polynomial degree.  The obstruction is not the quartic discriminant or a
bounded cofactor jet.  It is the highest homogeneous part of the actual
tangential Jacobian.

All rings below have characteristic zero.

## 1. The explicit weighted quartic

Use source coordinates `(x,y,z)` and put

```text
u=1+3xy,
gamma=1-4xy-x^2z.                                      (1)
```

The THM-3438 quartic map is

```text
A=(2u+u^2-3u^4 gamma^2)/x^2,
B=(1+u-2u^3 gamma^2)/x,
C=x gamma.                                              (2)
```

The displayed quotients cancel in `C[x,y,z]`, and

```text
Jac_(x,y,z)(A,B,C)=-6.                                  (3)
```

It has generic field degree four and the exact collision

```text
G(1,0,0)=G(-1,0,2)=(0,0,1).                            (4)
```

For `h in C[x,y]`, write

```text
A_h=A(x,y,h(x,y)),
B_h=B(x,y,h(x,y)),
J_h=Jac_(x,y)(A_h,B_h).                                 (5)
```

We prove that `J_h` is never constant.

## 2. Every positive-degree graph has a nonzero top Jacobian

Suppose `d=deg h>=1` and let `H` be the degree-`d` homogeneous part of `h`.
The leading homogeneous forms of `(1)` are

```text
u_top=3xy,
gamma_top=-x^2 H.                                      (6)
```

Only the terms containing `u^4 gamma^2` and `u^3 gamma^2` can contribute to
the top forms of `A_h` and `B_h`.  Define

```text
W=x^6 y^3 H^2.                                         (7)
```

Direct substitution into `(2)` gives

```text
(A_h)_top=-243 x^6 y^4 H^2=-243yW,
(B_h)_top=-54 x^6 y^3 H^2=-54W,                        (8)

(A_h)_top=(9/2)y(B_h)_top.                              (9)
```

The apparent proportionality in `(9)` does not kill the Jacobian, because
the proportionality factor is the variable `y`.  In fact

```text
Jac((9/2)yB_top,B_top)
  =-(9/2)B_top partial_x(B_top)
  =-13122 W partial_x W.                               (10)
```

The polynomial in `(10)` is nonzero.  Indeed, if `partial_x W=0`, then in
characteristic zero `W` would belong to `C[y]`; but the nonzero polynomial
`W=x^6y^3H^2` is divisible by `x^6`.  Therefore `(10)` is the nonzero top
homogeneous part of `J_h`.  Its degree is

```text
deg J_h=(2d+9)+(2d+8)=4d+17.                            (11)
```

In particular, no cancellation with lower graph terms can make `J_h` a
unit.

## 3. Constant graphs have the same obstruction

Let `h=c` be constant, including `c=0`.  Now

```text
gamma_top=-x(4y+cx).                                   (12)
```

Put

```text
W_c=x^4 y^3(4y+cx)^2.                                  (13)
```

The same calculation gives

```text
(A_c)_top=-243yW_c,
(B_c)_top=-54W_c,
(J_c)_top=-13122 W_c partial_x W_c.                    (14)
```

Here `(14)` is nonzero of degree seventeen.  Thus `(5)` is nonconstant for
every polynomial `h`, including the zero polynomial.

Combining Sections 2 and 3 proves the exact statement

```text
for every h in C[x,y],       Jac(A_h,B_h) notin C.      (15)
```

## 4. Block triangularity forbids every native graph descent

Suppose, toward a contradiction, that there are polynomials

```text
h in C[x,y],              g in C[A,B]                  (16)
```

such that `G` carries the source graph `z=h(x,y)` into the target graph
`C=g(A,B)`.  Equivalently,

```text
C(x,y,h)=g(A_h,B_h).                                   (17)
```

Use the source normal coordinate

```text
r=z-h(x,y)                                              (18)
```

and the target normal coordinate

```text
s=C-g(A,B).                                             (19)
```

Both are triangular coordinate changes with Jacobian one.  Equation `(17)`
says that `s(x,y,r)` is divisible by `r`; write

```text
s=r m(x,y,r).                                           (20)
```

On `r=0`, the derivative of `(A,B,s)` with respect to `(x,y,r)` is block
upper triangular.  Equations `(3)` and `(20)` therefore give

```text
-6=J_h m(x,y,0)                  in C[x,y].             (21)
```

The only units of `C[x,y]` are the nonzero constants.  Hence `(21)` would
force `J_h` to be constant, contradicting `(15)`.  Therefore

```text
there are no h,g satisfying (17).                       (22)
```

The same conclusion survives any polynomial automorphism of the native
`(A,B)` target plane, since its constant Jacobian only rescales `J_h`.

## 5. The known collision graph is already hostile

The two points in `(4)` lie on the graph

```text
h=1-x,                                                  (23)
```

because its values at `(1,0)` and `(-1,0)` are `0` and `2`.  Its top form is
`H=-x`, so `(7)` gives `W=x^8y^3`.  Equation `(10)` becomes

```text
(J_(1-x))_top=-104976 x^15 y^6.                         (24)
```

Thus the most immediate graph through the collision does not merely fail an
ad hoc target equation: its tangential Jacobian carries a degree-21 debt.
The all-degree proof shows that changing the polynomial graph cannot remove
this native debt.

There is a compatible direct target-plane control.  For constants `a,b,c`,
the pullback

```text
C-aA-bB-c                                               (25)
```

is quadratic in `z` unless `a=b=0`.  After putting `t=xy`, its discriminant
has `x`-coefficients

```text
[x^0] 36a^2(t+1)(3t+1)^5,
[x^1] 12ab(3t+1)^4(5t+4),
[x^2] 4(3t+1)^3(9act+3ac+6b^2t+4b^2),
[x^3] 8bc(3t+1)^3,
[x^4]=[x^5]=0,
[x^6]=1.                                                (26)
```

If this discriminant were a square in `C[x,y]`, its image under the injective
substitution `y=t/x` would be a square in `C[t][x,x^{-1}]`.  No negative
`x`-power can occur in a square root, because its least Laurent term would
square uniquely.  If `a!=0`, the least coefficient is the `x^0` term

```text
36a^2(t+1)(3t+1)^5,
```

which is not a square in `C[t]`: both displayed irreducible factors have odd
valuation.  If `a=0` but `b!=0`, the least coefficient is the `x^2` term

```text
8b^2(3t+1)^3(3t+2),
```

again not a square in `C[t]`.  Hence a square is possible only when `a=b=0`,
when it is exactly `x^6`.  The residual `C-c`
has no polynomial graph factor: for `c!=0` it is primitive linear with
nonunit `z`-coefficient `-x^3`, while for `c=0` its factors are `x` and
`1-4xy-x^2z`, neither a graph `z=h(x,y)`.  This independently checks the
affine target-plane edge of `(22)`.

## 6. Exact frontier and construction reframe

The connection ledger is now

```text
source object:      the explicit degree-four weighted Keller map G;
candidate descent:  source graph z=h -> target graph C=g(A,B);
preserved target:   ambient constant Jacobian and the known collision;
destroyed by graph: the source normal direction;
required sidecar:   unit normal multiplier in (21);
obstruction:        nonzero top tangential Jacobian (10);
minimal hostile:    h=1-x through both collision points.                (27)
```

This closes one construction lane completely: no amount of increasing the
degree of a native polynomial graph repairs the restriction.  A planar
counterexample descended from the weighted quartic must instead use at least
one of the following genuinely different objects:

1. a nonlinear source coordinate hypersurface not presented as `z=h(x,y)`;
2. a target coordinate pair not obtained from the native `(A,B)` plane; or
3. a non-restriction construction that algebraizes the exact inverse-
   different/cofactor incidence directly.

The theorem does **not** exclude those lanes.  In particular it does not
turn THM-3276's local inverse-different line into a global no-go, exclude an
arbitrary quartic `C3` Jelonek component, or prove `JC(2)`.

MISTAKE-416 is not a dependency: its current corrected entry concerns an
unstored historical Keller family, not quartic `C3` graph descent.  The
present proof uses only the literal THM-3438 formula and the graph
block-triangular mechanism.

## 7. Exact reproduction

Run

```bash
python3 -B 04-computation/jacobian_weighted_quartic_native_graph_no_go_thm3685.py
python3 -B -O 04-computation/jacobian_weighted_quartic_native_graph_no_go_thm3685.py
```

Both modes byte-match the stored transcript.  The companion checks
polynomial cancellation, `(3)`, coefficient-generic top forms through degree
five, six full lower-order graph controls, the collision and `(24)`, and the
complete affine-plane discriminant coefficients `(26)`.  Every truth gate
uses an explicit exception and remains active under optimized execution.

**QED.**
