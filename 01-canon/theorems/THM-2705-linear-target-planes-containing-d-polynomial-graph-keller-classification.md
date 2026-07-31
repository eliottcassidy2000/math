---
id: THM-2705
title: "Linear target planes containing d on polynomial graph Keller slices"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For every linear
  target plane containing d in the THM-2696 S3<S4
  quotient, all constant-nonzero-Jacobian polynomial graph slices are
  classified.  The pure A direction is the scaled THM-2702 triangular family;
  the pure B direction is empty; every genuinely mixed direction has exactly
  two affine graph solutions, exchanged by the square-root gauge, and both
  target maps are explicit triangular automorphisms.  Target planes not
  containing d, nonlinear target projections, arbitrary polynomial source
  surfaces, JC(2), and DC(2) remain open.
source: root/reflection-quotient-mixed-target-graph-2026-07-28; root-long-frontiers independent audit 2026-07-28
depends_on:
  - THM-2696-reflection-completed-s4-relative-different-and-coordinate-invariant-jacobian-gate
  - THM-2699-affine-plane-linear-projection-keller-slice-classification
  - THM-2702-polynomial-graph-coordinate-projection-keller-classification
related:
  - THM-2655-quartic-keller-resolvent-v4-quasietale-torsor-and-kummer-class-group-gate
script: 04-computation/jacobian_s4_polynomial_graph_mixed_d_planes_thm2705.py
output: 05-knowledge/results/jacobian_s4_polynomial_graph_mixed_d_planes_thm2705.out
script_sha256: cfaf4135a52d96fecb7644e5ffcbdb6c7f60e446289e2a8c4a8ca281d5f19fe2
output_sha256: 2a7ae1bc65af8303e4bcc661f7aebeb0c4cbb94ec86e846fcaaafeebc5966ebc
hash_basis: LF-normalized bytes
---

# THM-2705 -- every linear target plane containing `d` is classified on polynomial graphs

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2699 classifies affine source planes for arbitrary affine-linear target
projections.  THM-2702 permits arbitrary polynomial graphs but treats only the
three literal coordinate pairs.  The intersection of those two directions is
larger than either theorem states: every linear target plane containing the
third quotient coordinate `d` can be classified on every polynomial graph.

The mechanism is unexpectedly rigid.  A genuinely mixed target produces a
nonlinear first-order PDE, but its top degree in one characteristic coordinate
is an uncancellable square.  The sole quadratic boundary is a polynomial
Riccati equation, and the remaining linear equation has only three possible
slopes, one of which is impossible and two of which give affine automorphisms.

## 1. Statement and row reduction

Work over `C`.  On the graph

```text
z=f(x,y),            f in C[x,y],                         (1)
```

use THM-2696's quotient coordinates

```text
A=x^2-2y,            B=y^2-2xf,            d=f.          (2)
```

Every linear target two-plane containing `d` has, after an invertible target
row operation, the form

```text
(C,d),               C=alpha A+beta B,                   (3)
```

with `(alpha,beta)!=(0,0)`.  Adding a scalar multiple of `d` to `C` does not
change `Jac(C,d)`.  Fix `kappa in C^*`.

Then the complete classification of

```text
Jac(C,d)=kappa                                        (4)
```

is as follows.

1. If `beta=0`, necessarily `alpha!=0`, and

   ```text
   f=(kappa/(2 alpha))x+H(y-x^2/2),       H in C[t].     (5)
   ```

2. If `alpha=0` and `beta!=0`, there is no polynomial solution.

3. If `alpha beta!=0`, choose `q in C^*` and `h in C^*` by

   ```text
   q^2=-alpha/beta,          h=kappa/(2 beta q).          (6)
   ```

   There are exactly two solutions:

   ```text
   f_+=-q^2 x+q y+q^3-h,
   f_-=-q^2 x-q y+h-q^3.                                (7)
   ```

   Replacing `q` by `-q` also replaces `h` by `-h` and swaps the two
   displayed solutions, so `(7)` is independent of the square-root choice.

Every solution in `(5)` and `(7)` makes `(C,d)` a polynomial automorphism of
`A^2`.  Thus this entire graph/target-plane class contains no planar Jacobian
counterexample.

## 2. The master equation

The determinant is

```text
Jac(C,d)
 =2[(alpha-beta y)f_x+(alpha x-beta f)f_y].              (8)
```

When `beta=0`, put

```text
t=y-x^2/2.
```

Then `partial_x+x partial_y` becomes `partial_x` on `C[x,t]`, and `(8)` is

```text
2 alpha partial_x f=kappa.
```

This gives exactly `(5)`.  Moreover

```text
C=-2 alpha t,       d=(kappa/(2 alpha))x+H(t),           (9)
```

which is triangular and has the explicit inverse

```text
t=-C/(2 alpha),
x=(2 alpha/kappa)(d-H(t)),
y=t+x^2/2.                                             (10)
```

The factor `1/alpha` in `(5)` is load-bearing; omitting it would silently
change the prescribed Jacobian.

## 3. The pure `B` direction is empty

Suppose `alpha=0`, `beta!=0`.  With

```text
E=y partial_x+f partial_y,
```

equation `(8)` says

```text
E(f)=-kappa/(2 beta)=c!=0.                              (11)
```

Hence

```text
E(x)=y,       E^2(x)=f,       E^3(x)=c,       E^4(x)=0,
E(y)=f,       E^2(y)=c,       E^3(y)=0.                 (12)
```

Thus `E` is locally nilpotent on the polynomial generators and therefore on
`C[x,y]`.  Its exponential is the polynomial automorphism

```text
X=x+Ty+(T^2/2)f+(T^3/6)c,
Y=y+Tf+(T^2/2)c                                      (13)
```

of `C[T,x,y]`.  Its Jacobian is a unit in `C[T]`, hence is identically one
after evaluation at `T=0`.  The coefficient of `T` in that Jacobian is
`div(E)=f_y`, so `f_y=0`.  But `(11)` then becomes

```text
y f_x=c,
```

which is impossible after setting `y=0`.  This proves case 2.

## 4. The mixed PDE and its degree collapse

Now assume `alpha beta!=0` and choose `(q,h)` as in `(6)`.  Put

```text
t=y-qx,
f=qt+q^3-h+g(x,t).                                    (14)
```

Remembering that `partial_x` at fixed `y` is
`partial_x-q partial_t`, substitution in `(8)` gives

```text
Jac(C,d)/2-beta qh
 = beta[(h-g)g_t-(q^2+qx+t)g_x-qg].                    (15)
```

Since `beta qh=kappa/2`, equation `(4)` is equivalent to

```text
(h-g)g_t-(q^2+qx+t)g_x-qg=0.                           (16)
```

Let `r=deg_t(g)` when `g!=0`, with leading coefficient `a_r(x)`.

- If `r>=3`, the term `-g g_t` has the unique top degree `2r-1`, with
  coefficient `-r a_r^2`.  Every other term has degree at most `r+1`.
  Therefore `a_r=0`, a contradiction.
- If `r=2`, the coefficient of `t^3` is

  ```text
  -(a_2'+2a_2^2).                                      (17)
  ```

  A nonzero polynomial cannot satisfy `a_2'+2a_2^2=0`: for positive degree,
  the two terms have different degrees, and for degree zero the square term
  is nonzero.  Thus `a_2=0`.

Consequently `deg_t(g)<=1`.  Write

```text
g=a(x)t+b(x).                                          (18)
```

The coefficients of `t^2,t,1` in `(16)` give

```text
a'=0,
b'=-a(a+q),
b'(a+2q)=0,
ah-(a+q)b_0-b'q^2=0,                                  (19)
```

where `b=b'x+b_0` after the first three equations.  Hence

```text
a in {0,-q,-2q}.                                       (20)
```

For `a=0`, equation `(19)` gives `b=0`, hence `g=0` and the first solution in
`(7)`.  For `a=-q`, its last equation reads `-qh=0`, impossible because
`q,h!=0`.  For `a=-2q`,

```text
b=-2q^2x+2h-2q^3,                                     (21)
```

which gives the second solution in `(7)`.  No nonlinear polynomial solution
has survived.

## 5. The two mixed maps are triangular

For `f_+`, retain `t=y-qx`.  Direct substitution gives

```text
d=qt+q^3-h,
C=beta(t^2+2q^2t+2hx).                                (22)
```

Thus

```text
t=(d-q^3+h)/q,
x=(C/beta-t^2-2q^2t)/(2h),
y=t+qx.                                                (23)
```

For `f_-`, put `u=y+qx`.  Then

```text
d=-qu+h-q^3,
C=beta(u^2+2q^2u-2hx),                                (24)
```

with inverse

```text
u=(h-q^3-d)/q,
x=(u^2+2q^2u-C/beta)/(2h),
y=u-qx.                                                (25)
```

These are precisely the two square-root gauges of THM-2699's affine Family
III.  The PDE proof adds the substantive statement that no nonlinear graph
can sit above them.

## 6. Boundary and non-implications

This theorem uses both restrictions in its title.

- The source is a global graph `z=f(x,y)` in the fixed quotient coordinates.
  An arbitrary polynomial surface is not covered.
- The target plane must contain `d`.  A plane not containing `d` cannot in
  general be row-reduced to `(alpha A+beta B,d)`.
- Only affine-linear target coordinates are classified.  Nonlinear target
  functions have different equations.
- The theorem classifies a large, natural slice of the THM-2696 `S4`
  resolvent geometry.  It neither proves nor disproves `JC(2)` or `DC(2)`.

In particular, the next honest target is a plane not containing `d`, or a
polynomial source surface not representable as `(1)`, not another pass over
the mixed PDE `(16)`.

## 7. Exact companion

Run

```text
python 04-computation/jacobian_s4_polynomial_graph_mixed_d_planes_thm2705.py
python -O 04-computation/jacobian_s4_polynomial_graph_mixed_d_planes_thm2705.py
```

Both modes must byte-match

```text
05-knowledge/results/jacobian_s4_polynomial_graph_mixed_d_planes_thm2705.out.
```

The companion uses explicit `require` checks.  It verifies `(8)`, the scaled
pure-`A` family and inverse, the exact mixed reduction `(15)`, every
coefficient in `(19)`, the all-degree inequalities behind the top-degree and
Riccati gates, both mixed Jacobians and inverses, and the square-root gauge
swap.  These exact controls support but do not replace the all-degree proof in
Sections 2--4.

An independent hostile audit rederived the target-plane row reduction, the
scaled pure-`A` characteristic calculation, and the pure-`B` locally
nilpotent derivation/divergence contradiction.  It separately checked the
mixed PDE, including the `r=2` polynomial Riccati boundary, every equation in
`(19)`, the square-root gauge, both inverses, and the stated scope.  Normal and
optimized execution byte-matched the recorded output SHA-256
`2a7ae1bc65af8303e4bcc661f7aebeb0c4cbb94ec86e846fcaaafeebc5966ebc`;
the script SHA-256 is
`cfaf4135a52d96fecb7644e5ffcbdb6c7f60e446289e2a8c4a8ca281d5f19fe2`.
The stored transcript has nine lines; the script has zero Python `assert`
nodes and eighteen explicit `require` calls.

QED.
