---
id: THM-3544
title: "Planar Keller target-pencil total-degree-six floor"
status: >
  PROVED.  If a complex planar Keller pair is not a polynomial automorphism,
  then every nonzero member of its target pencil has total degree at least
  six.  Degrees at most four are reduced to the complete linear-through-cubic
  source-fibre closures by placing a projective root of the top binary form
  on a fibre axis.  In degree five, the power-free-face theorem forces the
  ordinary top form to be a fifth power of one linear form; the all-direction
  fibre floor then exposes a weighted top face lambda*x^5+c*y^4, which is
  squarefree and contradicts that same theorem.  Tame shears show the
  nonautomorphic hypothesis is essential, and a sextic square shows that the
  two inherited gates alone genuinely stop at six.  This is a pencil-degree
  floor, not a construction, degree-pair classification, or proof of JC(2).
source: boxeph-2026-08-18-planar-jacobian-dephasing
depends_on:
  - THM-2063-one-fiber-linear-planar-keller-pairs
  - THM-2071-quadratic-fiber-square-parity-gate
  - THM-2102-power-free-weight-face-and-first-defect-descent
  - THM-2118-all-degree-cubic-faber-boundary-flux-coprimality
  - THM-2740-polynomial-coordinate-first-target-triangularity-and-mixed-resolvent-shear-closure
related:
  - THM-2063-one-fiber-linear-planar-keller-pairs
  - THM-2071-quadratic-fiber-square-parity-gate
  - THM-2102-power-free-weight-face-and-first-defect-descent
  - THM-2118-all-degree-cubic-faber-boundary-flux-coprimality
  - THM-2740-polynomial-coordinate-first-target-triangularity-and-mixed-resolvent-shear-closure
---

# THM-3544 -- planar Keller target-pencil total-degree-six floor

**PROVED.**

## 1. Statement

Let

```text
F=(P,Q):C^2 -> C^2,             Jac(P,Q)=kappa in C*.  (1)
```

If `F` is not a polynomial automorphism, then for every `(s,t)!=(0,0)`,

```text
R=sP+tQ                 is nonconstant and deg(R)>=6.  (2)
```

The quantifier is over the complete projective target pencil, not only the
two displayed coordinates.  Total degree and the conclusion are invariant
under invertible affine source and target changes.

## 2. The inherited all-direction fibre floor

Fix `(s,t)!=(0,0)`.  Choose `(u,v)` so that

```text
Delta=sv-tu !=0,
S=uP+vQ.                                                  (3)
```

Then

```text
Jac(R,S)=Delta kappa in C*.                              (4)
```

In particular `R` is nonconstant.  It is also not a polynomial coordinate:
if it were, THM-2740 would make `(R,S)`, and hence `(P,Q)`, a polynomial
automorphism.

Let `(ell,m)` be any invertible affine linear source coordinate system.
The current canon closes each of the first three fibre degrees:

```text
deg_m R<=1   -> tame by THM-2063,
deg_m R=2    -> tame by THM-2071,
deg_m R=3    -> automorphism by THM-2118.                (5)
```

Therefore a nonautomorphic pair satisfies the simultaneous curvature floor

```text
deg_m R>=4                                                    (6)
```

for every nonzero target-pencil member and every linear source-fibre
direction.  THM-2084 is a precursor to the last line of `(5)`, but its
finite reduced-degree statement alone is not used as a closure.

## 3. Degrees at most four collapse to the fibre closures

Write `d=deg R` and let `R_d` be its nonzero top homogeneous binary form.
Assume `1<=d<=4`.  Over `C`, `R_d` has a projective root `[b] in P^1(C)`.
Choose a second vector `a` and make the invertible linear source change

```text
(x,y)=a ell+b m.                                         (7)
```

The coefficient of `m^d` in the transformed polynomial is exactly `R_d(b)`,
so

```text
deg_m R<=d-1<=3.                                         (8)
```

This contradicts `(6)`.  Multiplicity can only lower the contribution of
`R_d`; lower homogeneous layers may restore fibre degree `d-1`, but the
universal bound `(8)` is all that is used.  A root at infinity in the
original chart is still an ordinary projective vector and is made the
`m`-axis by `(7)`.  No translation or unproved genericity is being used.

## 4. The quintic is pinched by two different weights

It remains to exclude `d=5`.  Apply THM-2102 to ordinary total-degree weight
`(1,1)`.  Since `R` is not a coordinate, its top form must be a proper power.
A homogeneous quintic has prime degree, so

```text
R_5=lambda L^5,              lambda in C*,             (9)
```

for one linear form `L`.  Make a linear source change with `L=x`.  Then

```text
R=lambda x^5+(terms of total degree <=4).               (10)
```

By `(6)` in the `y` direction, `deg_y R>=4`; equation `(10)` gives the
reverse inequality.  Thus `deg_y R=4`.  Its `y^4` coefficient is a nonzero
constant `c`: a term `x^i y^4` with `i>0` would have total degree at least
five, but `(9)` says the only degree-five term is `lambda x^5`.  Hence

```text
[y^4]R=c in C*.                                          (11)
```

Now use positive weights

```text
wt(x)=4,                     wt(y)=5.                   (12)
```

Both terms in `(10)--(11)` have weight twenty.  Every other term has total
degree at most four, so

```text
4i+5j=4(i+j)+j<=20,                                    (13)
```

with equality only for `(i,j)=(0,4)`; among degree-five terms, only
`(i,j)=(5,0)` occurs.  Therefore the complete weighted top form is

```text
in_(4,5)(R)=lambda x^5+c y^4.                           (14)
```

This binomial is not a proper power.  Indeed

```text
gcd(lambda x^5+c y^4, 5lambda x^4, 4c y^3)=1,          (15)
```

whereas in characteristic zero a nonconstant power `H^e`, `e>=2`, shares
the factor `H` with both partial derivatives.  Equation `(14)` contradicts
THM-2102's conclusion that every positive-weight top face of a hard Keller
component is a proper power.  The quintic case is impossible, completing
the proof of `(2)`.

## 5. Sharpness and failure boundary

The nonautomorphic hypothesis is indispensable.  For every `d>=1`,

```text
(P,Q)=(x^d+y,x),                Jac(P,Q)=-1             (16)
```

is the tame automorphism with inverse `(p,q) |-> (q,p-q^d)` and has a pencil
member of degree `d`.  At `d=5`, it escapes exactly through
`deg_y(x^5+y)=1`, the first line of `(5)`.

The argument also genuinely stops at six.  The homogeneous sextic

```text
H=[xy(x-y)]^2                                             (17)
```

has fibre degree six generically and degree four along each of its three
root directions; every positive-weight exposed face is a square.  It is not
a Keller component—its two partial derivatives share `xy(x-y)`—but it is a
sharp hostile for the combination of the two mechanisms used above.  Any
improvement beyond six needs an additional Keller equation, response,
nonproperness, or coefficient sidecar rather than another repetition of the
projective-root plus proper-power-face argument.

The theorem gives no upper bound, no list of degree pairs, and no existence
claim at degree six.  It proves neither a planar counterexample nor `JC(2)`.

**QED.**
