---
id: THM-2715
title: "Nonlinear d-target shear polynomial-graph Keller classification"
status: >
  PROVED + VERIFIED-EXACT + TWO INDEPENDENT HOSTILE AUDITS.  For every H in
  C[d], all nonzero-constant-Jacobian polynomial
  graph solutions for the nonlinear target pair (A,B+H(d)) are classified.
  Linear H gives the unique shifted cubic of THM-2709.  Quadratic H gives
  exactly two explicit quadratic sections on one sharp coefficient wall;
  both induce the same triangular automorphism.  Degree at least three gives
  no polynomial graph.  Arbitrary nonlinear target pairs, nongraph source
  surfaces, JC(2), and DC(2) remain open.
source: root-nonlinear-graph-target-scout-2026-07-28
depends_on:
  - THM-2696-reflection-completed-s4-relative-different-and-coordinate-invariant-jacobian-gate
  - THM-2709-complete-linear-target-plane-polynomial-graph-keller-classification
related:
  - THM-2702-polynomial-graph-coordinate-projection-keller-classification
  - THM-2705-linear-target-planes-containing-d-polynomial-graph-keller-classification
script: 04-computation/jacobian_s4_nonlinear_d_target_shear_thm2715.py
output: 05-knowledge/results/jacobian_s4_nonlinear_d_target_shear_thm2715.out
script_sha256: 96501cb96639e1081bcc0056233b3137a8385a69bf94a2d50270a4a6ab568be5
output_sha256: f6872ba42e9bc085cfe8cc87146abc0ec9b0cdc956b5f8be72fbe9118f7133a0
hash_basis: LF-normalized bytes
---

# THM-2715 -- nonlinear `d`-target shears have only a quadratic wall

**PROVED + VERIFIED-EXACT + TWO INDEPENDENT HOSTILE AUDITS.**

THM-2709 classifies every *linear* target plane on a polynomial graph in the
THM-2696 quotient.  The first nonlinear family is not harder by degree-by-
degree elimination.  Its PDE has an exact primitive, and that primitive
leaves only one quadratic wall.

## 1. Statement

Work over `C`, put

```text
z=f(x,y),
A=x^2-2y,              B=y^2-2xz,              d=z,       (1)
```

and fix `H in C[d]`.  Suppose

```text
Jac(A,B+H(d))=kappa!=0.                                (2)
```

Write

```text
t=y-x^2/2,                    c=-kappa/4!=0.            (3)
```

Then exactly the following cases occur.

1. If `H(d)=beta d+gamma` (including `beta=0`), there is one graph:

   ```text
   f=(x+beta/2)t/2
     +(x+beta/2)(x^2+beta^2/4)/8+c.                    (4)
   ```

2. If

   ```text
   H(d)=alpha d^2+beta d+gamma,              alpha!=0, (5)
   ```

   a graph exists if and only if

   ```text
   beta+2 alpha c=0,              equivalently beta=alpha kappa/2. (6)
   ```

   For either of the two roots

   ```text
   p^2=-1/(4 alpha),                                      (7)
   ```

   the corresponding graph is

   ```text
   f_p=p x^2-4p^2x+2pt+8p^3+c.                          (8)
   ```

   These are the only two graphs.

3. If `deg H>=3`, no polynomial graph satisfies `(2)`.

Every graph in cases 1--2 makes `(A,B+H(d))` a polynomial automorphism of
the graph plane.  Thus this entire nonlinear one-coordinate target family
contains no counterexample to the planar Jacobian conjecture.

## 2. The exact primitive

For a polynomial `F(x,t)=f(x,t+x^2/2)`, differentiation at fixed `t` gives

```text
F_x=f_x+x f_y.                                           (9)
```

The coordinate determinants from THM-2709 are

```text
Jac(A,B)=-4(xf_x+x^2f_y+f-xy),
Jac(A,d)= 2(f_x+x f_y).                                  (10)
```

Consequently `(2)` is

```text
xF_x+F-xt-x^3/2-H'(F)F_x/2=c.                          (11)
```

The left side is an exact derivative.  Integrating in `x` over `C[t]`
gives

```text
xF-H(F)/2=x^2t/2+x^4/8+cx+K(t)                         (12)
```

for a unique `K in C[t]` once `F` is fixed.  Conversely differentiating
`(12)` recovers `(11)`, so no solution is lost.

There is also a target-level identity independent of the degree of `H`.
Using `(12)` and `y=t+x^2/2`,

```text
B+H(d)=t^2-2cx-2K(t).                                  (13)
```

Hence every solution of `(12)` is automatically triangular whenever
`c!=0`; classification of the graphs is the only remaining task.

## 3. Degree collapse

Let `n=deg H>=2` and `m=deg_x F`.  The right side of `(12)` has `x`-degree
four.  The two potentially leading terms on the left have degrees

```text
m+1,                         nm.                       (14)
```

The case `m=0` is impossible.  If the degrees in `(14)` agree, then
`(n-1)m=1`, hence `(n,m)=(2,1)`; cancellation can only lower the degree
below four.  Otherwise the larger term cannot cancel.  Therefore

```text
(n,m)=(2,2) or (4,1).                                  (15)
```

This already excludes every degree except two and four.

For degree four, write `F=P(t)x+R(t)` and let `a_4,a_3` be the top two
coefficients of `H`.  The `x^4` row of `(12)` says

```text
-a_4 P(t)^4/2=1/8,                                    (16)
```

so `P` is a nonzero constant.  The `x^3` row then says

```text
4a_4R(t)+a_3=0,                                       (17)
```

so `R` is constant too.  The left `x^2` coefficient is now constant, while
the right coefficient is `t/2`, a contradiction.  This removes degree four
and proves case 3.

## 4. The quadratic wall and the two sections

Let `H` be `(5)` and initially write

```text
F=P(t)x^2+Q(t)x+R(t).                                 (18)
```

The `x^4,x^3,x^2,x` rows of `(12)` successively give

```text
P^2=-1/(4alpha),
Q=1/alpha,
R(t)=(1/alpha-beta P-t)/(2alpha P),
beta+2alpha c=0.                                      (19)
```

Indeed the first row makes the polynomial `P(t)` a nonzero constant and
the second then makes `Q(t)=1/alpha` constant, so the displayed comparison
does not silently assume constant leading coefficients.

The constant row merely defines

```text
K(t)=-H(R(t))/2.                                      (20)
```

Thus `(6)` is necessary and sufficient.  Renaming `P=p` and using `(6)`--
`(7)` turns `(18)`--`(19)` into `(8)`.  The two square roots in `(7)` give
two distinct polynomials and exhaust the solutions.

Equation `(13)` simplifies for either root to the same target coordinate

```text
U=A=-2t,
V=B+H(d)
  =-2cx+2t/alpha+gamma-alpha c^2-1/alpha^2.            (21)
```

Since `c!=0`, the inverse is polynomial:

```text
t=-U/2,
x=(2t/alpha+gamma-alpha c^2-1/alpha^2-V)/(2c),
y=t+x^2/2.                                             (22)
```

Its Jacobian is `-4c=kappa`.

For a concrete sharp control, take

```text
alpha=-1,        c=1,        beta=2.                   (23)
```

The two graphs are

```text
f=-x-y,                     f=y-x+2,                   (24)
```

and both give

```text
(U,V)=(-2t,-2x-2t+gamma).                              (25)
```

Keeping `alpha=-1,c=1` but changing `beta` away from `2` destroys every
polynomial graph.  At `c=0`, the quadratic wall can still carry the two
formal sections, but `(21)` loses its `x` coordinate and its Jacobian is
zero.  These are the sharp coefficient and nonzero-response boundaries.

## 5. The linear branch and triangularity

For `H(d)=beta d+gamma`, `(12)` is divisible by `x-beta/2`.  Evaluation at
that root determines `K(t)` uniquely, and division gives `(4)`.  Equivalently
this is THM-2709 with `b=-beta/2`.  Direct substitution gives

```text
U=-2t,
V=(t+beta^2/8)^2-2c(x-beta/2)+gamma,                  (26)
```

so this branch is triangular as well.

## 6. Scope and next target

The classification concerns the infinite-dimensional but special target
family

```text
(A,B+H(d))
```

on one global polynomial graph in the fixed THM-2696 quotient.  It does not
classify a pair with both target coordinates nonlinear, a polynomial such as
`H(A,d)` with genuine mixed dependence, a source surface not graphical in
this chart, an arbitrary planar Keller map, or a Weyl-algebra endomorphism.
It proves neither `JC(2)` nor `DC(2)`.

The mechanism nevertheless identifies a sharper next probe than “nonlinear
targets”: allow the second coordinate to depend jointly on `A` and `d`.
The one-variable shear is exhausted, and every survivor here is triangular.

## 7. Exact companion

Run

```text
python 04-computation/jacobian_s4_nonlinear_d_target_shear_thm2715.py
python -O 04-computation/jacobian_s4_nonlinear_d_target_shear_thm2715.py
```

and compare both transcripts with

```text
05-knowledge/results/jacobian_s4_nonlinear_d_target_shear_thm2715.out.
```

The companion uses exact symbolic arithmetic and explicit failure raises.
It checks the coordinate/PDE primitive, the complete quadratic coefficient
system, both sections, the common triangular target and inverse, the degree-
four terminal contradiction, and the linear positive control.  The general
degree argument is Section 3, not a finite computation.

QED.
