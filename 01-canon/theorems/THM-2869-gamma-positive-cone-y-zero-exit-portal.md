---
id: THM-2869
title: "Gamma positive-cone y-zero exit portal"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT.  On the boundary
  U=d_1+x d_3, V=d_2, the Gamma quadratic/cubic divisibility equations
  have exactly one solution with alpha>0 and x>0.  Its Gamma shape is the
  unique positive root of an explicit degree-23 polynomial, isolated in
  (23.2446,23.2447), and its x-coordinate is the unique positive linear
  subresultant root, isolated in (0.71257,0.71259).  This classifies a
  boundary portal; it does not prove that the THM-2865 local branch
  reaches it.  Independent hostile audit is pending.
source: root/audit-2809-gamma-y-zero-exit-2026-07-28
depends_on:
  - THM-2824-arbitrary-three-slot-factorial-moment-divisibility-and-atomic-orientation-boundary
  - THM-2865-gamma-transverse-null-holotopy-and-uniform-fourth-exit
related:
  - THM-2830-disjoint-positive-adjacent-cone-factorial-moment-three-detection
  - THM-2846-arbitrary-positive-cone-moment-three-transverse-boundary
  - THM-2853-gamma-adjacent-tensor-cycle-weighted-positivity
script: 04-computation/gmc_gamma_y_zero_exit_portal_thm2869.py
output: 05-knowledge/results/gmc_gamma_y_zero_exit_portal_thm2869.out
script_sha256: 044cb2bd0821fba508d46381786b52de714e12b845dc5bdae9fa64cac9fe8d4c
output_sha256: f09c39a908ba64b2cf14ea5dfd875ee3347b0d50401ab7bffd384401c0f71dc6
hash_basis: LF-normalized bytes
---

# THM-2869 -- Gamma positive-cone y-zero exit portal

**PROVED CANDIDATE + VERIFIED-EXACT.**

THM-2865 gives a local positive branch

```text
U=d_1+x d_3,                  V=d_2+y d_3
```

near the factorial shape `alpha=1`.  Numerical continuation suggests that
its `y`-coordinate eventually reaches zero.  This theorem proves the exact
algebraic statement which is available without trusting that continuation:
the full boundary face `y=0` contains exactly one positive Gamma shape and
exactly one positive `x` at which the quadratic moment divides the cubic.

The word “portal” records this necessary boundary location.  It does not
assert that any particular interior branch reaches it.

## 1. Boundary equations

For `alpha>0`, put

```text
L_alpha(s^m)=(alpha)_m,
f_n=s^n/(alpha)_n,
d_n=f_(n+1)-f_n,                                      (1)

U=d_1+x d_3,                     V=d_2.                (2)
```

Use the exact Gamma-adjacent tensor formula of THM-2865 to form

```text
g11=L_alpha(U^2),    g12=L_alpha(UV),    g22=L_alpha(V^2),

t111=L_alpha(U^3),  t112=L_alpha(U^2V),
t122=L_alpha(UV^2), t222=L_alpha(V^3).                 (3)
```

The division-free quadratic/cubic invariants are

```text
I1=3t112 g11 g22-t222 g11^2-2t111 g12 g22,
I2=3t122 g11 g22-2t222 g12 g11-t111 g22^2.            (4)
```

Before setting `y=0`, both have the common positive denominator

```text
alpha^4(alpha+1)^4(alpha+2)^4(alpha+3)^3.              (5)
```

Let `F(alpha,x)` and `H(alpha,x)` be the resulting integer numerators
after restriction to `y=0`.  Their `x`-degrees are `4` and `3`, and their
leading coefficients never vanish for `alpha>0`.

## 2. Exact resultant and the unique positive shape

Exact elimination gives

```text
Res_x(F,H)
 =-7077888
   (alpha+2)^6 (alpha+3)^22 (alpha+4)^4
   (alpha+5)^2 (alpha+8)^2
   P9(alpha) P23(alpha),                              (6)
```

where the descending coefficient list of `P9` is

```text
(8,
 3364,
 145397,
 2870386,
 32586944,
 229937674,
 1025612391,
 2810119176,
 4307331852,
 2819687976).                                          (7)
```

Every coefficient in `(7)` is positive, so `P9(alpha)>0` for
`alpha>0`.

The descending coefficient list of `P23` is

```text
(405,
 67862,
 4458573,
 152480084,
 2640768818,
 739093820,
-1227760066990,
-36903041464160,
-665697843962479,
-8704796296181594,
-88073949221238223,
-711703149831960260,
-4673300030645893432,
-25164153182568166184,
-111485361916416386032,
-405718990297087577216,
-1204903266566658216448,
-2884528384378398958208,
-5458050914397962088448,
-7916087911348550367232,
-8378757135131489026048,
-5946662811236683546624,
-2388466752014210826240,
-339470562155929010176).                              (8)
```

There is exactly one sign variation in `(8)`.  Descartes' rule gives at
most one positive root, counted with multiplicity.  Exact evaluation and
Sturm counting give

```text
P23(116223/5000)<0<P23(232447/10000),

# roots in
(116223/5000,232447/10000)=1.                          (9)
```

The polynomial is squarefree.  Hence it has one simple positive root

```text
alpha_* in (23.2446,23.2447).                         (10)
```

All other factors in `(6)` are nonzero for `alpha>0`.  Therefore every
positive boundary solution must have `alpha=alpha_*`.

## 3. Specialization-safe linear subresultant

Remove the `x`-contents of `F,H`.  Their exact subresultant degree profile
is

```text
4,3,2,1,0.                                             (11)
```

The degree-one member has the form

```text
12(alpha+2)^3(alpha+3)^2(alpha+5)^2(alpha+8)
  [D(alpha)x-N(alpha)],                                (12)
```

where `D,N` both have degree `26`, and all `27` coefficients of each are
strictly positive.  Moreover,

```text
gcd(D,P23)=1.                                          (13)
```

At `alpha_*`, the resultant is zero while the linear member `(12)` is
nonzero.  Subresultant specialization therefore gives a gcd of exact
degree one and forces

```text
x_*=N(alpha_*)/D(alpha_*)>0.                           (14)
```

For a direct converse check, substitute `(14)` into both primitive
boundary polynomials, clear denominators, and reduce in `alpha`.
The two remainders modulo `P23` are exactly zero.  Thus `(alpha_*,x_*)`
is a genuine common zero, not merely a necessary resultant root.

Exact Bernstein conversion on the isolating interval `(9)` gives

```text
71257/100000 < x_* < 71259/100000.                    (15)
```

Equations `(10),(14)--(15)` prove existence and uniqueness of the positive
boundary pair.

## 4. Moment consequence and scope

The Gamma integral is strictly positive on squares of nonzero real
polynomials.  Since `U` and `V` in `(2)` are linearly independent,
their Gram quadratic is positive definite.  Let `z` be either of its
nonreal projective roots and put

```text
K=U+zV.
```

The common vanishing of `(4)` is exactly divisibility of the cubic by the
quadratic, while every adjacent difference has zero first moment.
Consequently

```text
L_(alpha_*)(K)
=L_(alpha_*)(K^2)
=L_(alpha_*)(K^3)
=0.                                                    (16)
```

This is an exact positive-cone boundary hostile for a Gamma functional.
It is not a GMC2 counterexample: `alpha_*` is not the one-complex-Gaussian
shape `alpha=1`.  No assertion is made about the fourth or higher moments
at the portal.

Most importantly, THM-2869 does not prove the global continuation

```text
THM-2865 local branch  --->  (alpha_*,x_*,0).          (17)
```

Numerics strongly suggest `(17)`, but a proof needs a connected
continuation corridor with no singularity or branch switch.  The exact
result here is the uniqueness of the only positive `y=0` destination.

## 5. Exact verification

Run

```text
python 04-computation/gmc_gamma_y_zero_exit_portal_thm2869.py
python -O 04-computation/gmc_gamma_y_zero_exit_portal_thm2869.py
```

The companion pins the promoted THM-2865 tensor builder.  It reconstructs
`(4)--(6)`, checks the coefficient and Descartes data in `(7)--(8)`,
isolates the root by exact Sturm counting, verifies `(11)--(15)`, and
checks that the restored common-root numerators vanish modulo `P23`.

Normal, optimized, and stored-output replay must be byte-identical after
LF normalization.  LF-normalized SHA-256:

```text
script  044cb2bd0821fba508d46381786b52de714e12b845dc5bdae9fa64cac9fe8d4c
output  f09c39a908ba64b2cf14ea5dfd875ee3347b0d50401ab7bffd384401c0f71dc6
```

Independent hostile audit is pending.

## 6. Connection contract

```text
source:
  the Gamma positive-cone plane U=d1+x d3, V=d2+y d3;

target:
  the codimension-one face y=0;

map:
  exact resultant projection to alpha followed by a specialization-safe
  linear subresultant recovering x;

preserved:
  positivity of alpha and x, quadratic/cubic divisibility, and the
  nonreal Gram-root moment-three null;

destroyed / not proved:
  connected ancestry from the THM-2865 local branch and all moments
  beyond order three;

needed sidecar:
  an exact continuation corridor from alpha=11/10 to alpha_* with
  nonvanishing Jacobian and y>0 before the terminal face;

cheapest decisive test:
  cover the numerical branch by overlapping rational moving boxes and
  certify face orientation plus Jacobian sign in each box.
```

No Gaussian Moment Conjecture case is closed by this theorem.
