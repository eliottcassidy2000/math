---
id: THM-3630
title: "Russell-cylinder non-even formal survivor and no finite jet bound"
status: >
  RESERVED / PROVISIONAL PROOF CANDIDATE / PENDING INDEPENDENT AUDIT.
  On the non-even rho=0 tangent stratum, an explicit three-branch formal
  packet admits an all-order decomposable formal target pair (F,w) with
  common Jacobian 12. Hermite and target-jet interpolation therefore give
  polynomial folds with regular decomposable target-pair jets surviving every
  prescribed finite collision cutoff. This constructs neither one all-order
  polynomial nor a global constant-Jacobian pair or Keller map.
source: root / audit_thm3627 all-order non-even continuation, 2026-08-21
audit: PENDING -- provisional theorem and exact companion require hostile audit
depends_on:
  - THM-3624-russell-cylinder-noneven-fold-weighted-cokernel-boundary
  - THM-3619-russell-cylinder-even-fold-higher-jet-staircase
related:
  - THM-3605-russell-cylinder-graph-slice-puncture-no-filling
  - THM-3627-russell-cylinder-noneven-hostile-degree-six-closure
script: 04-computation/jc2_russell_cylinder_noneven_formal_survivor_no_finite_jet_thm3630.py
output: 05-knowledge/results/jc2_russell_cylinder_noneven_formal_survivor_no_finite_jet_thm3630.out
script_sha256: e057bef58af3f6f6a979becff2f7e9e63f738de8934e396642addd9a3ec285af
output_sha256: adc573378b8de821d14d83cd4b0781949bbf06515189ff1414cb069edf5fc320
hash_basis: raw LF bytes
---

# THM-3630 -- Russell-cylinder non-even formal survivor and no finite jet bound

**RESERVED / PROVISIONAL PROOF CANDIDATE / PENDING INDEPENDENT AUDIT.**
This theorem identifies the sharp limitation of the local collision-jet
route on the principal non-even stratum left by THM-3624 and THM-3627. The
single hostile polynomial of THM-3627 closes at source degree six, but the
whole polynomial family cannot be closed at any uniform finite jet order.

All rings, formal germs, and differential forms are over `C`.

## 0. Setup and strict scope

Use the regular local compiler chart from THM-3619,

```text
D=1+x^2q,
a=q/D^2,                  c=xD(D+2),
Jac_(x,q)(a,c)=-3,                                      (1)
```

and the quadratic stable-coordinate fold

```text
q=Q(x)+w^2.                                             (2)
```

At the retained collision, the three source points are `x=-1,0,1` and the
common target is `(a,c,w)=(-3/4,0,0)`. Fix

```text
u in C \ {0,9/4,-9/4}.                                  (3)
```

The exclusion `u!=0` merely makes the packet genuinely non-even; the formal
construction itself also works at `u=0`. The exclusions `u!=+-9/4` are
load-bearing: they are exactly the exceptional tangent values already closed
at degree two in THM-3624.

This theorem uses a **three-point formal multi-germ** of `Q`, namely an
element of

```text
C[[x+1]] direct_product C[[x]] direct_product C[[x-1]]. (4)
```

It is not one polynomial or one globally continued analytic function.

## 1. The exact formal survivor packet

Put

```text
Q_infinity(x)=-3/4-9/(4x^2).                            (5)
```

This rational function is regular at `x=+-1`. Prescribe

```text
Q_-(x)=Q_infinity(x)       in C[[x+1]],
Q_0(x)=-3/4+ux+O(x^2)     in C[[x]],
Q_+(x)=Q_infinity(x)       in C[[x-1]].                 (6)
```

The higher coefficients of `Q_0` are completely arbitrary. The value and
first-derivative packet is

```text
x=-1: (-3,-9/2),       x=0: (-3/4,u),
x=+1: (-3, 9/2),                                      (7)
```

so `(6)` lies on the nonzero `rho=0` stratum of THM-3624.

Let

```text
s=w^2,              g=(1-4s/3)^(-1/2).                 (8)
```

The three moving points

```text
x_-(w)=-g,                 x_0(w)=0,
x_+(w)= g                                                   (9)
```

lie on the respective fold germs `(2),(6)` and have the exact common target
curve

```text
Gamma: (a,c,w)=(-3/4+w^2,0,w).                         (10)
```

Indeed `Q_infinity(+-g)+s=-3/g^2`, so the two side points have `D=-2`,
while the middle point has `D=1`.

Differentiate at fixed `w` along each fold surface. On `Gamma`,

```text
                 a_x                         c_x
x=-g         -9/(4g^3)                       3
x= 0              u                          3
x=+g          9/(4g^3)                       3.          (11)
```

Thus `c,w` are formal coordinates on every branch surface. If

```text
z=a+3/4-w^2,                                             (12)
```

the three surface images have unique graph equations

```text
S_i: z=f_i(c,w),
f_i(c,w)=lambda_i(w)c+O(c^2),                           (13)
```

where

```text
lambda_-=-3/(4g^3),       lambda_0=u/3,
lambda_+= 3/(4g^3).                                     (14)
```

Their pairwise differences are units times `c`: the constants
`-3/4,u/3,3/4` are distinct by `(3)`.

## 2. Exact formal decomposable-pair gluing

On `S_i`, write the desired source area form in the surface coordinates
`c,w`:

```text
12 dx wedge dw=J_i(c,w) dc wedge dw,
J_i=12/c_(x,i).                                         (15)
```

Equation `(11)` gives the decisive common boundary value

```text
J_-(0,w)=J_0(0,w)=J_+(0,w)=4.                          (16)
```

Take the formal primitive in the transverse coordinate

```text
A_i(c,w)=integral_(0)^c J_i(r,w) dr
        =4c+O(c^2).                                     (17)
```

Write `phi_i=f_i` for the three graph series `(13)`, in any fixed branch
order `1,2,3`. Their differences are

```text
phi_j-phi_i=c times a unit.                             (18)
```

The direct Newton divided differences are

```text
D_21=(A_2-A_1)/(phi_2-phi_1),
D_31=(A_3-A_1)/(phi_3-phi_1),                           (19)

D_321=(D_31-D_21)/(phi_3-phi_2).                       (20)
```

All three are regular. Indeed `(17)` gives

```text
A_j-A_1 in c^2 C[[c,w]],                               (21)
```

while every denominator in `(19)` is `c` times a unit. Thus
`D_21,D_31` lie in `c C[[c,w]]`; their difference is again divisible by `c`,
and the last denominator in `(20)` is `c` times a unit.

Now define the single ambient formal target function

```text
F(z,c,w)=A_1+D_21(z-phi_1)
               +D_321(z-phi_1)(z-phi_2).               (22)
```

Direct substitution gives the exact Newton identities

```text
F(phi_i(c,w),c,w)=A_i(c,w),               i=1,2,3.     (23)
```

Take the second formal target function to be `G=w`. On `S_i`, equations
`(15),(17),(23)` give

```text
Phi_(Q_i)^*(dF wedge dw)
  =dA_i wedge dw
  =J_i dc wedge dw
  =12 dx wedge dw.                                     (24)
```

Thus `(F,w)` is an all-order **decomposable formal target pair** with common
Jacobian `12` on the three source germs. The earlier normal-degree
Vandermonde view gives the same extension order by order; `(19)--(23)` are
its direct Newton interpolation.

## 3. No uniform finite collision-jet bound for polynomials

Fix any source cutoff `N`. The pullback through total source degree `N`
depends on only finitely many derivatives of `Q` at `-1,0,1` (matching
through order `N+1` is more than sufficient). By ordinary three-point
Hermite interpolation there is a polynomial `Q_N` whose jets through that
order equal the three formal germs `(6)`. It may be chosen of degree at most

```text
3(N+2)-1=3N+5.                                         (25)
```

Take a sufficiently long target jet of the formal `F` in `(22)`. The global
coordinate ring of the affine target surjects onto every finite local Artin
quotient, so that jet has a regular target-function representative `F_N`.
The second function `w` is already regular; after the Russell isomorphism it
is the polynomial target function

```text
w=Y(B+2)/8-CS.                                         (26)
```

Therefore `(F_N,w)` is a decomposable regular target-pair jet whose source
Jacobian equals `12` through degree `N` at all three collision preimages.
Since `u!=0`, every such `Q_N` is non-even.

Consequently:

```text
for every N, some non-even polynomial rho=0 fold has a
regular decomposable target-pair jet surviving through N.          (27)
```

There is therefore no **uniform finite collision-jet cutoff**, even in the
closed decomposable subsystem, that closes the general non-even `rho=0`
polynomial family. This does not say that one fixed polynomial survives
every order.

## 4. Exact finite controls

For compact controls, match `(6)` only through side/central derivative order
`n` and take the unique Hermite polynomial of degree at most `3n+2`, with
`u=1` and all central derivatives of orders `2,...,n` set to zero. Exact
arithmetic gives

| matched order `n` | polynomial degree | source cutoff `N=2n-2` | `rank P_N` | `rank[P_N|tau_N]` |
|---:|---:|---:|---:|---:|
| 2 | 8 | 2 | 15 | 15 |
| 3 | 11 | 4 | 40 | 40 |
| 4 | 14 | 6 | 77 | 77 |
| 5 | 17 | 8 | 126 | 126 |

Here `P_N` is the enlarged arbitrary-two-form map in the polynomial local
coordinates `(c,e+3,w)` used in THM-3627, and `tau_N` is the common constant
column normalized to `12`. These ranks are independent controls in the
larger system; the decomposable finite-jet conclusion comes from the Newton
construction `(17)--(26)`.

There is also a hostile repair chain starting from the exact polynomial
`Q_h` of THM-3627. Put

```text
H_4=x^5(x+1)^5(x-1)^4,
lambda_4=-58212503/2249728,
Q_[4]=Q_h+lambda_4 H_4.                                 (28)
```

The perturbation preserves every derivative through order four at `x=-1,0`,
every derivative below order four at `x=1`, and changes only the new right
fourth derivative by

```text
768 lambda_4=-174637509/8788.                           (29)
```

It repairs the THM-3627 degree-six failure:

```text
N=6: rank=augmented rank=77,
N=7: rank=augmented rank=100.                           (30)
```

Next put

```text
H_5=x^6(x+1)^6(x-1)^5,
lambda_5=2428928805/58492928,
Q_[5]=Q_[4]+lambda_5 H_5.                               (31)
```

This preserves all previous jets and changes only the new right fifth
derivative by

```text
7680 lambda_5=36433932075/114244.                       (32)
```

Then

```text
N=8: rank=augmented rank=126.                           (33)
```

These repairs are finite controls of the one-new-relation/two-side-jet
freedom. They are not used in the formal proof of Sections 1--3.

## 5. Honest remaining frontier

The theorem proves an all-order survivor only in the product of the three
completed source local rings, and proves arbitrarily long survival only by a
sequence of polynomials depending on the cutoff. It does **not** provide

- one polynomial `Q` surviving all orders;
- a convergent or algebraic globalization of the formal multi-germ;
- one global regular target pair having constant Jacobian on one polynomial
  fold;
- a Keller map or a Jacobian-conjecture counterexample.

Whether every fixed non-even polynomial is eventually closed remains
**OPEN**. A terminal argument cannot simply retain the `Q_infinity` right-
hand side: after parity is removed, lower independent side errors feed later
invoices. Already for the THM-3624 packet at `u=1`, the pure moving-tangent
determinant predicts the degree-four constant `11178`, whereas exact
elimination gives `10449`.

## 6. Exact reproduction

The deterministic companion verifies the exact common curve, all transverse
derivatives and exceptional slope boundaries, the Vandermonde and direct
Newton interpolation identities, every divided-difference regularity gate,
the decomposable pullback, the four Hermite controls, both hostile repairs,
optimized replay, and an AST gate excluding inactive Python `assert`
statements.

Reproduce with

```bash
python3 04-computation/jc2_russell_cylinder_noneven_formal_survivor_no_finite_jet_thm3630.py
python3 -O 04-computation/jc2_russell_cylinder_noneven_formal_survivor_no_finite_jet_thm3630.py
```

Both streams must be byte-identical to
`05-knowledge/results/jc2_russell_cylinder_noneven_formal_survivor_no_finite_jet_thm3630.out`.
