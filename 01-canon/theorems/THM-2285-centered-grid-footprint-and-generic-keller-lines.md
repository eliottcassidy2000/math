---
id: THM-2285
title: "Centered grid footprint and generic Keller lines"
status: >
  PROVED + INDEPENDENTLY AUDITED. A nonzero polynomial with separate
  coordinate degrees A_i has at least the product of s_i nonzero values
  on every product grid whose ith side has A_i+s_i points; the bound and
  centered integer radii are sharp for arbitrary polynomials. Consequently
  every monic planar Keller map has a bounded vertical integer target line
  with at most deg_y(P) nongeneric fibres and a bounded horizontal line
  with at most deg_y(Q) nongeneric fibres. Enlarged centered rectangles
  contain an explicit quadratic footprint of generic integer fibres, so
  these targets have density one. This does not prove generic degree one,
  properness, JC(2), or DC(2).
source: codex-2026-07-25-centered-grid-footprint
depends_on:
  - THM-2241-monic-transverse-response-depth-and-resultant-nonproper-quotient
  - THM-2280-centered-polynomial-grid-avoidance-and-bounded-generic-keller-fibre
related:
  - THM-2270-simultaneous-balanced-cut-relation-and-six-uniform-orientation
---

# THM-2285 -- centered grid footprint and generic Keller lines

**PROVED + INDEPENDENTLY AUDITED.**

THM-2280 finds one bounded grid point outside a finite polynomial bank.
The same induction retains much more information: every extra sample in
every coordinate multiplies the guaranteed nonvanishing footprint.

For a monic planar Keller map, the polynomial being sampled is the top
resultant coefficient. The footprint theorem therefore produces infinitely
many integer targets with the exact generic fibre count on two bounded
coordinate lines, and quadratically many such targets in centered boxes.

## 1. Separate-degree footprint theorem

Let `K` be a field. Take a nonzero polynomial

```text
F in K[T_1,...,T_r]                                   (1)
```

and nonnegative integers `A_i` such that

```text
deg_(T_i)F<=A_i.                                      (2)
```

For each coordinate, let

```text
E_i subset K,
|E_i|=A_i+s_i,
s_i>=1.                                               (3)
```

Then

```text
#{t in E_1 x ... x E_r:F(t)!=0}
 >=product_(i=1)^r s_i.                               (4)
```

### Proof

Induct on `r`. For `r=1`, a nonzero polynomial of degree at most `A_1`
has at most `A_1` roots among the `A_1+s_1` sample points, leaving at
least `s_1` nonzeros.

For the induction step, expand in the last variable:

```text
F=sum_(j=0)^(A_r) F_j(T_1,...,T_(r-1))T_r^j.         (5)
```

Choose one nonzero coefficient `F_j`. It still satisfies

```text
deg_(T_i)F_j<=A_i,                    i<r.            (6)
```

By induction, there are at least

```text
product_(i<r)s_i                                     (7)
```

base points at which `F_j` is nonzero. At each such base point, the slice
of `F` in `T_r` is a nonzero polynomial of degree at most `A_r`, so it is
nonzero at at least `s_r` points of `E_r`. Multiplying (7) by `s_r`
proves (4). QED.

This proof uses only the one-variable root bound. No total-degree estimate,
irreducibility, or generic-position hypothesis is hidden in it.

## 2. Exact sharpness

The lower bound in (4) is sharp for every choice of finite grids and
separate degree bounds.

Choose subsets

```text
Z_i subset E_i,
|Z_i|=A_i,                                             (8)
```

and put

```text
F_0(T_1,...,T_r)
 =product_i product_(z in Z_i)(T_i-z).                (9)
```

Empty products are allowed when `A_i=0`. The separate degree in coordinate
`i` is exactly `A_i`, and

```text
F_0(t)!=0
 iff t_i in E_i minus Z_i for every i.                (10)
```

There are exactly

```text
product_i(|E_i|-|Z_i|)=product_i s_i                 (11)
```

such points.

Thus neither the product in (4) nor any side length `A_i+s_i` can be
improved in a theorem using only the separate degrees.

## 3. Centered integer grids and polynomial banks

Use the centered consecutive set from THM-2280:

```text
S(n)={-floor(n/2),...,ceil(n/2)}.                    (12)
```

It has `n+1` points and radius `ceil(n/2)`.

Assume now that `K` has characteristic zero. Apply (4) to

```text
E_i=S(A_i+s_i-1).                                    (13)
```

The integer images in `K` are distinct. Therefore `F` is nonzero at at
least `product_i s_i` integer points satisfying

```text
|t_i|<=ceil((A_i+s_i-1)/2).                          (14)
```

This radius is optimal for a centered integer carrier of the required
size. If `A_i+s_i` distinct integers all have absolute value at most `R`,
then

```text
A_i+s_i<=2R+1,
R>=ceil((A_i+s_i-1)/2).                             (14a)
```

and the set in (13) attains equality.

For a finite bank

```text
f_1,...,f_m in K[T_1,...,T_r] minus {0},             (15)
```

put

```text
F=product_j f_j,
A_i=sum_j deg_(T_i)f_j.                              (16)
```

The polynomial ring is an integral domain, so `F!=0`, and its separate
degrees are exactly the sums in (16). Equations (4) and (14) give at least

```text
product_i s_i                                        (17)
```

centered integer points at which every member of the bank is nonzero.

For two variables, if

```text
A=sum_j deg_U f_j,
B=sum_j deg_V f_j,                                   (18)
```

then the rectangle

```text
S(A+s-1) x S(B+t-1)                                  (19)
```

contains at least `st` simultaneous nonzeros.

### A bounded generic vertical and horizontal line

There is also a one-dimensional footprint that is infinite rather than
finite. Expand the nonzero product `F(U,V)` in powers of `V` and choose a
nonzero coefficient `g(U)`. Since

```text
deg g<=A,                                             (20)
```

the univariate case on `S(A)` gives an integer

```text
u_0 in S(A),
|u_0|<=ceil(A/2),
g(u_0)!=0.                                            (21)
```

Hence `F(u_0,V)` is a nonzero polynomial of degree at most `B`. It vanishes
at no more than `B` scalar values. Consequently

```text
f_j(u_0,v)!=0 for every j
```

for all but at most `B` integers `v`.

Interchanging the variables gives an integer

```text
v_0 in S(B),
|v_0|<=ceil(B/2),                                    (22)
```

such that all but at most `A` integers `u` avoid the complete bank on the
horizontal line `V=v_0`.

## 4. Quantitative density in square grids

Let

```text
I_R={-R,-R+1,...,R},
m=2R+1.                                               (23)
```

If

```text
m>A,
m>B,                                                  (24)
```

apply (4) with

```text
s=m-A,
t=m-B.                                                (25)
```

The bank is simultaneously nonzero at at least

```text
(m-A)(m-B)                                            (26)
```

points of `I_R^2`. Thus the exceptional count is at most

```text
m^2-(m-A)(m-B)
 =(A+B)m-AB.                                          (27)
```

For fixed degree data, the exceptional proportion is therefore

```text
O(1/R).                                               (28)
```

The simultaneous nonvanishing points have density one among integer target
points. Equation (27) is a universal separate-degree bound; it does not
claim that every particular algebraic curve contains that many lattice
points.

## 5. Application to monic planar Keller maps

Let

```text
P,Q in C[x,y],
Jac(P,Q)=kappa in C*,

d=deg_y P>=1,
e=deg_y Q>=0,                                         (29)
```

and assume that `P` is monic in `y`. Let

```text
Phi=(P,Q):A^2_C -> A^2_C.                             (30)
```

THM-2241 and THM-2280 define the nonzero top resultant coefficient

```text
c(U,V) in C[U,V] minus {0}                           (31)
```

and prove

```text
deg_U c<=e,
deg_V c<=d,                                           (32)

c(u,v)!=0
 iff #Phi^(-1)(u,v)=N,                               (33)
```

where `N` is the generic geometric fibre cardinality. The Keller condition
makes every local multiplicity one.

Write the actual separate degrees as

```text
a=deg_U c<=e,
b=deg_V c<=d.                                         (34)
```

### 5.1 Bounded generic target lines

Apply (21) to `c`. There is an integer

```text
u_0,
|u_0|<=ceil(a/2)<=ceil(e/2),                          (35)
```

such that

```text
#Phi^(-1)(u_0,v)=N                                   (36)
```

for all but at most

```text
b<=d                                                  (37)
```

integers `v`.

Likewise, there is an integer

```text
v_0,
|v_0|<=ceil(b/2)<=ceil(d/2),                          (38)
```

such that

```text
#Phi^(-1)(u,v_0)=N                                   (39)
```

for all but at most

```text
a<=e                                                  (40)
```

integers `u`.

Thus one bounded vertical line and one bounded horizontal line each carry
infinitely many exact generic fibres.

### 5.2 Enlarged centered target rectangles

For every `s,t>=1`, the centered rectangle

```text
S(a+s-1) x S(b+t-1)                                  (41)
```

contains at least `st` integer targets whose fibres have exactly `N`
points. The coarser rectangle

```text
S(e+s-1) x S(d+t-1)                                  (42)
```

also contains at least `st` such targets: apply the footprint theorem
directly with the upper bounds in (32).

### 5.3 Density-one generic integer fibres

If

```text
m=2R+1>a,b,                                           (43)
```

then at least

```text
(m-a)(m-b)                                            (44)
```

targets in `I_R^2` have exact fibre cardinality `N`, and at most

```text
(a+b)m-ab                                             (45)
```

do not. In particular,

```text
#{(u,v) in I_R^2:#Phi^(-1)(u,v)=N}/m^2
 ->1.                                                 (46)
```

This is a lattice-target density statement, not a Zariski-density
replacement: the generic locus was already Zariski open by construction.

## 6. Boundaries and non-implications

The grid theorem is universally sharp by Section 2, but that sharpness is
for arbitrary polynomials with the same separate degrees. No claim is made
that a Keller resultant realizes every sharp grid pattern.

The vertical-line conclusion also cannot force properness. A nonconstant
polynomial in one variable has only finitely many zeros, so a punctured
target line can still lie in the generic locus while the global exceptional
curve `V(c)` is nonempty. Finite monodromy over that punctured line can be
nontrivial, and nothing in the footprint forces

```text
N=1.                                                   (47)
```

The load-bearing scope remains exactly that of THM-2280:

```text
base field C for the Keller application;
P monic in y;
the selected target orientation;
c(U,V) the nonzero top-X resultant coefficient;
generic geometric fibre cardinality, not a quotient continuation state.
                                                                  (48)
```

The centered bounds are not invariant under arbitrary complex coordinate
changes. The theorem does not make `c` constant, prove properness, decide
the singular trigonal boundary of THM-2262, continue the grade-six DC gauge
ladder of THM-2240, or prove `JC(2)` or `DC(2)`.

## 7. Connection and loss ledger

```text
source:
  THM-2280's nonzero polynomial bank and the top resultant coefficient
  of a monic planar Keller map;

target:
  a sharp quantitative grid footprint, bounded generic target lines,
  and density-one exact generic integer fibres;

map:
  induct on variables using one nonzero coefficient and the univariate
  root bound, then apply the result to c(U,V);

preserved:
  simultaneous nonvanishing, separate coordinate degrees, centered
  integer radii, and exact generic fibre cardinality;

destroyed:
  the geometry and monodromy of V(c), coordinate covariance, sheet
  labels, nonproperness data away from the sampled targets, and N itself;

cheapest hostile probes:
  products of grid-linear factors, the boundary e=0, a nonconstant
  punctured target line, and a bank whose pullback vanishes identically;

needed sidecar:
  an independent argument making c constant or forcing N=1.             (49)
```

No computation is needed: every assertion follows from the root bound,
separate-degree induction, and the proved resultant equivalence (33). QED.
