---
id: THM-2316
title: "Degree-eighteen C--D ratio-bank closure"
status: >
  PROVED + VERIFIED-EXACT. In the genuine nonsplit polynomial
  exact-square-prefix degree-eighteen branch of THM-2262/2297, all four
  C--D ratio points in THM-2311's two-sparse bank are empty. The rational
  linear point has an absolutely irreducible genus-four trigonal spectrum:
  ten simple branches plus one smooth totally ramified cubic fibre. The
  other three points form one irreducible cubic orbit. Uniformly over its
  ratio field, their spectra are absolutely irreducible of genus three:
  ten simple branches plus one ordinary node whose two normalization
  branches are unramified. Infinity is unramified throughout. Every
  rational Keller trajectory is therefore constant and gives the inherited
  nonsplit-deck contradiction. The C--D bank is empty; this does not prove
  JC(2).
source: codex-2026-07-25-degree18-cd-linear-closure
depends_on:
  - THM-2262-degree-eighteen-trigonal-spectral-discriminant-reduction
  - THM-2297-degree-eighteen-target-translation-normal-form
  - THM-2311-degree-eighteen-two-sparse-weighted-ratio-bank
related:
  - THM-2314-degree-eighteen-bd-linear-ratio-closure
script: 04-computation/jc2_degree18_cd_linear_ratio_closure_thm2316.py
output: 05-knowledge/results/jc2_degree18_cd_linear_ratio_closure_thm2316.out
script_sha256: f61104412da79fd49c7fd90b21b404799a0c0eb70091c583204ff42940ebbb3a
output_sha256: f76fb74cb3b4e535de03011a48e73892b42875ac219047c3a57c6370a2f61af1
hash_basis: working-tree bytes (LF)
---

# THM-2316 -- the full C--D ratio bank has positive genus

**PROVED + VERIFIED-EXACT.**

THM-2311 reduces the exactly two-sparse degree-eighteen branch to a finite
weighted-projective bank. On the `C`--`D` line, its factors are

```text
22143375+6397664t,

387420489-8964338040t+54880100352t^2
 +16544432128t^3,

t=D^3/C^4.                                           (1)
```

At the linear root, no number-field normalization is needed. The identity

```text
C=D=1/t
```

is a rational representative. Its spectral normalization is an absolutely
irreducible genus-four curve. The cubic factor is one irreducible ratio
orbit, and its three normalizations all have genus three.

## 1. A rational weighted representative

Use THM-2297's target-translation normal form

```text
(B,C,D,W),                         weights (2,3,4,5),

G_0(u,y;B,C,D,W)=0,               wt(u,y)=(2,1).    (2)
```

Put

```text
r=1/t=-6397664/22143375,

(B,C,D,W)=(0,r,r,0).                               (3)
```

Then

```text
D^3/C^4=1/r=t.                                      (4)
```

THM-2311 proves that this ratio classifies the off-axis orbit on the
weighted `C`--`D` line, so (3) represents the required point. Substitution
in `G_0`, followed by multiplication by the nonzero scalar `91125/49`,
gives the integral equation

```text
H(u,y)
 =-48427561125u^3
  +2989355625u^2y^2
  -258339375uy^4
  +28427563008u
  +2095875y^6
  +233971712y^3
  -1052872704y^2
 =0.                                                 (5)
```

The constant rescaling changes neither the curve nor its normalization.

## 2. Absolute irreducibility

The leading `u` coefficient in (5) is a nonzero constant. If `H` were
reducible over `C[y,u]`, its cubic `u`-degree would give a root in `C(y)`.
After dividing by the leading coefficient, that root is integral over
`C[y]`; since `C[y]` is integrally closed, it belongs to `C[y]`.

Let its `y`-degree be `d`. If `d>2`, substitution in (5) gives `y`-degree
`3d` from the cubic term, strictly larger than the possible degrees

```text
2d+2,             d+4,             6                (6)
```

from the other terms. The leading coefficient cannot cancel. Hence every
possible polynomial root has the form

```text
u=ay^2+by+c.                                        (7)
```

Substitute (7) into (5) and equate the seven coefficients of
`1,y,...,y^6`. Exact Buchberger reduction over `Q` gives the reduced
Groebner basis

```text
{1}.                                                 (8)
```

The coefficient ideal remains the unit ideal after extending constants to
`C`. Thus (7) is impossible and `H` is absolutely irreducible. Its
projective normalization is a connected degree-three cover of the
`y`-line.

## 3. Ten simple branch values

For the unscaled normalized curve `G_0` at (3), the exact `u`-discriminant
is

```text
Delta(y)
 =-464767295827968/1953125
   *(15y-52)^2 h_10(y),                             (9)
```

where

```text
h_10(y)
 =2864799140625y^10
  +19862607375000y^9
  +103285558350000y^8
  +367474086720000y^7
  +527020367328000y^6
  -762214033152000y^5
  -6205500108288000y^4
  -55515871936512000y^3
  -144341267034931200y^2
  -333588706036285440y
  -578220423796228096.                              (10)
```

Exact Euclidean reduction gives

```text
gcd(h_10,h_10')=1,
gcd(h_10,15y-52)=1.                                 (11)
```

Thus (10) has ten distinct roots, all simple roots of the discriminant.
At every such value, the normalization has exactly one simple ramification
point of index two. These fibres contribute

```text
10*(2-1)=10                                         (12)
```

to the ramification divisor.

## 4. The squared factor is smooth total ramification

At the remaining finite branch value

```text
y_0=52/15
```

the whole cubic is

```text
G_0(u,y_0)
 =-49/2460375*(10935u-2704)^3.                     (13)
```

Put

```text
u_0=2704/10935.
```

Direct differentiation gives

```text
partial_y G_0(u_0,y_0)
 =-2384639688704/2278125 !=0.                       (14)
```

The curve is therefore smooth at `(u_0,y_0)`. Since the first two
`u`-derivatives vanish there and the third does not, the implicit local
solution has

```text
y-y_0 = unit*(u-u_0)^3+higher terms.                (15)
```

There is one normalization point above `y_0`, with ramification index
three. It contributes `3-1=2`, exactly explaining the squared factor in
(9).

## 5. Infinity is unramified

Use the weighted infinity chart

```text
z=1/y,                         v=u/y^2.              (16)
```

At `z=0`,

```text
z^6G_0(v/z^2,1/z)
 =L_infinity(v)+O(z),

L_infinity(v)
 =1127-138915v+1607445v^2-26040609v^3.              (17)
```

Its discriminant is

```text
Disc_v L_infinity
 =-153384762202971019112448 !=0.                    (18)
```

Hence the normalization has three distinct smooth points above infinity,
and `z` is a local parameter at each. All three are unramified for the
projection to `P^1_y`.

Equations (9)--(11) list every finite branch value. The total ramification
is therefore

```text
10+2=12.                                            (19)
```

Riemann--Hurwitz for the connected degree-three cover gives

```text
2g-2=3*(-2)+12=6,

g=4.                                                 (20)
```

## 6. The cubic ratio orbit has genus three

The other three ratios are the roots of

```text
P(t)
 =16544432128t^3+54880100352t^2
  -8964338040t+387420489.                           (21)
```

Modulo five, its monic reduction is

```text
t^3-t^2-2,
```

which has no root in `F_5`. Thus (21) is irreducible over `Q`. Let

```text
K=Q(alpha),                         P(alpha)=0.
```

The rational-function-free representative

```text
(B,C,D,W)=(0,alpha^2,alpha^3,0)                     (22)
```

has

```text
D^3/C^4=alpha.
```

The irreducibility argument of Section 2 works over `K`: every reducible
cubic would have a root `ay^2+by+c`. Exact Buchberger reduction over
`K` again gives

```text
{1}.                                                 (23)
```

Hence the spectral curve is absolutely irreducible for every embedding of
`K` into `C`.

Write

```text
Y
 =(2130734520alpha-175906971+557240320alpha^2)
   /71030268.                                        (24)
```

The exact discriminant factorization in `K[y]` is

```text
Delta_alpha(y)
 =-153384762202971019112448
   *(y-Y)^2 h_(10,alpha)(y),                         (25)
```

where `h_(10,alpha)` is monic of degree ten and exact Euclidean reduction
gives

```text
gcd(h_(10,alpha),h'_(10,alpha))=1,
gcd(h_(10,alpha),y-Y)=1.                            (26)
```

Thus (25) again has ten simple branch values. The squared factor now has a
different geometry from Section 4. Put

```text
R
 =(11851365288alpha-1192022163+5827039232alpha^2)
   /12375051136,

S
 =-(251919849384alpha-20715038739+34173669376alpha^2)
   /111375460224.                                    (27)
```

Then

```text
G_alpha(u,Y)
 =-26040609(u-R)^2(u-S),                            (28)
```

and at `(R,Y)`,

```text
G_alpha=partial_u G_alpha=partial_y G_alpha=0.
```

For local coordinates `U=u-R`, `V=y-Y`, write the tangent cone as

```text
A U^2+B UV+C V^2.
```

Its `U^2` coefficient and discriminant are the following nonzero
power-basis elements of `K`:

```text
A
 =-59049(179291068488alpha-15721619103
          +43308511232alpha^2)/126276032 !=0,

B^2-4AC
 =31381059609
  *(23902870780036903800843
    -559844552503271371280616alpha
    +3545334969099076026159104alpha^2)
  /2466893777661329408 !=0.                         (29)
```

They are nonzero because their displayed degree-at-most-two
representatives are nonzero in the cubic field. Therefore `(R,Y)` is an
ordinary node. Since `A!=0`, neither tangent is the vertical fibre
`V=0`; both normalization branches use `V` as a local parameter and are
unramified over the `y`-line. The third branch at `u=S` is also unramified.

The ten roots in (26) are consequently the only finite ramification,
each with index two. The infinity cubic is still (17), so infinity is
unramified. Hence

```text
sum_P(e_P-1)=10,

2g-2=3*(-2)+10=4,

g=3.                                                 (30)
```

All statements are identities in `K`; applying any of its three complex
embeddings proves the same genus-three result for every root of (21).

## 7. Every C--D Keller trajectory is impossible

A putative survivor in the retained branch supplies

```text
(u,y) in C(x)^2,                  G_0(u,y)=0.         (31)
```

If (31) is nonconstant, it extends to a nonconstant morphism from `P^1`
to the relevant projective normalization. Riemann--Hurwitz forbids a
nonconstant map from genus zero to genus four or genus three. Thus `u,y`
are constant.

The wall `y=0` is already empty by THM-2262. Its inherited first-flux
identity

```text
Z=T^2=-2N_2/(5103y)                                 (32)
```

then makes `Z`, nonzero `T`, and `q` constant. The genuine nonsplit deck
fixes the algebraically closed constant field but sends `q` to `-q`,
contradicting `q!=0`. Hence all four ratios in (1) are empty.

## 8. Exact gain and stopping boundary

The entire `C`--`D` bank of THM-2311 is empty:

```text
4 -> 0.                                             (33)
```

Together with THM-2314's full `B`--`D` closure, the exactly two-sparse bank
shrinks from `31` to

```text
BC:9,             BW:8,             DW:4,

total:21.                                           (34)
```

The connection ledger is

```text
source:
  all four ratios in THM-2311's C--D bank;

map:
  use C=D=1/t at the rational point and
  (C,D)=(alpha^2,alpha^3) over the irreducible cubic orbit;

preserved:
  labelled weighted orbit, spectral isomorphism class, normalization,
  branch indices, and genus;

temporarily forgotten:
  the scaled third flux, Keller one-form, and whole-Faber sidecar;

why restoration is unnecessary:
  genera four and three already make every rational spectral trajectory
  constant;

hostile control:
  C=D=1 has squarefree branch discriminant, detecting an accidentally
  generic or identically repeated calculation.                        (35)
```

This theorem closes only the `C`--`D` bank inside the genuine nonsplit
polynomial exact-square-prefix degree-eighteen branch. The `BC`, `BW`,
and `DW` banks, every three-/four-sparse singular stratum, split/even
descent, other Newton edges, `JC(2)`, and `DC(2)` remain open.

## 9. Exact reproduction

Run

```bash
python3 04-computation/jc2_degree18_cd_linear_ratio_closure_thm2316.py
python3 -O 04-computation/jc2_degree18_cd_linear_ratio_closure_thm2316.py
```

Both executions are byte-identical to the stored output. The companion
verifies the rational representative, integral curve (5), complete
discriminant (9)--(10), squarefreeness and coprimality (11), exceptional
fibre (13)--(14), cubic-field irreducibility, factorization (25),
ordinary node (27)--(29), both unit Groebner bases, separable infinity,
ramification arithmetic, and squarefree hostile control. The
Riemann--Hurwitz and deck steps are the mathematical proof above, not
computer assumptions. QED.
