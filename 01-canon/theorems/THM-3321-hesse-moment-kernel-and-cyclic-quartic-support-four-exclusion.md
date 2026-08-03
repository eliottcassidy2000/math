---
id: THM-3321
title: "Hesse moment kernel and cyclic-quartic support-four exclusion"
status: >
  PROVED (closed mixed-moment formula and normalization boundary) +
  FINITE-EXACT (guarded homogeneous Macaulay certificates), with independent
  immutable audit pending.  In the degree-four cyclic eigenspace of THM-3310,
  every mixed simplex moment is a coefficient of
  (1-X^3-Y^3-3XY)^(-1), giving an exact one-sum formula and lattice recurrence.
  For each of the five coefficient hyperplanes, moments M3 through M21 span
  the entire degree-21 homogeneous piece over F_101.  Five displayed maximal
  minors are nonzero modulo 101; the denominator guard 101>86 lifts those
  minors to characteristic zero and excludes every projective coefficient
  support <=4.  The formal torus of the Hesse surface does not preserve the
  simplex moment functional, so this torus supplies no second coefficient
  normalization beyond the always-lawful projective scaling.  No
  classification of all continuous coefficient symmetries is claimed.  Full
  support five and FC(3) remain OPEN.
audit: >
  The companion proves the Hesse denominator identity in exact Q(omega),
  compares the closed formula with an independent six-index Dirichlet
  expansion for all 66 pairs a+b<=10, and checks the recurrence on 7225
  lattice points.  It rebuilds all five degree-18 through degree-22 Macaulay
  maps at p=101, uses known-empty and known-nonempty projective controls, and
  reconstructs explicit 2024-by-2024 maximal minors with determinants
  27,67,42,50,82.  Normal and optimized transcripts byte-match the frozen
  output; no floating point, randomness, affine-only chart, or Python assert
  is used.
source: root/creative-synthesis-next/2026-08-03
depends_on:
  - THM-3310-degree-four-cyclic-eigenspace-on-the-triangle
related:
  - THM-3018-factorial-conjecture-as-a-simplex-moment-problem
  - THM-3300-factorial-gaussian-torus-bridge-and-the-archimedes-no-go
  - THM-3303-keller-simplex-null-moments-force-a-boundary-collision
  - THM-3323-cyclic-quartic-support-five-exact-degree-21-rank
script: 04-computation/degree_four_support4_macaulay_scout_20260803.py
output: 05-knowledge/results/degree_four_support4_macaulay_scout_20260803.out
script_sha256: 74bb6fd1c0019f5fe934fa27744bfab68bfbfdc6bb746857500b91c8507792e1
output_sha256: fec8bc92300373612a1b18c9a1a26818232cdc771ed95b646d685dbd64b3357e
hash_basis: LF-normalized bytes
---

# THM-3321 -- Hesse moment kernel and cyclic-quartic support-four exclusion

**PROVED + FINITE-EXACT; INDEPENDENT IMMUTABLE AUDIT PENDING.**

Retain THM-3310's coordinates

```text
z=s1+omega*s2+omega^2*s3,       zbar=s1+omega^2*s2+omega*s3
```

and degree-four `omega`-eigenvector

```text
g=A*zbar+B*z^2+C*z*zbar^2+D*z^3*zbar+E*zbar^4.             (1)
```

## 1. Closed mixed moments

For nonnegative `a,b`, put `N=a+b` and let angle brackets denote the uniform
probability integral on the triangle.  Polarizing the Dirichlet simplex
formula gives

```text
<(Xz+Yzbar)^N>
 =2*N!/(N+2)! * h_N(X*lambda+Y*lambda^(-1): lambda^3=1).    (2)
```

The complete-homogeneous generating function has denominator

```text
product_(lambda^3=1)(1-X*lambda-Y*lambda^(-1))
 =1-X^3-Y^3-3XY.                                           (3)
```

Define

```text
C(a,b)=[X^aY^b](1-X^3-Y^3-3XY)^(-1).                      (4)
```

Comparing the coefficient of `X^aY^b` in `(2)` proves

```text
mu(a,b)=<z^a zbar^b>
       =2*a!*b!*C(a,b)/(a+b+2)!,                           (5)

C(a,b)=sum multinomial(i+j+k;i,j,k)*3^k,                   (6)
```

where the sum is over `i,j,k>=0` with

```text
3i+k=a,                  3j+k=b.                           (7)
```

Consequently `mu(a,b)=0` exactly when `a-b` is nonzero modulo three.  The
coefficients also satisfy the constant-work lattice recurrence

```text
C(a,b)=C(a-3,b)+C(a,b-3)+3C(a-1,b-1)+[a=b=0],             (8)
```

with negative indices zero.  Equations `(6)` and `(8)` are respectively a
one-dimensional closed sum and an efficient exact sequence compiler.

## 2. The Hesse torus is not a moment symmetry

Put

```text
r=z*zbar,       u=z^3,       v=zbar^3,       uv=r^3,
z*g=A*r+B*u+C*r^2+D*r*u+E*r*v.                             (9)
```

The formal torus `z->t z`, `zbar->t^(-1)zbar` has basis weights
`(-1,2,-1,2,-4)`, but the five pure terms of `M_3=<g^3>` are all nonzero and
have three distinct weights `(-3,6,-3,6,-12)`.  Hence this Hesse torus
preserves the surface equation in `(9)` but not the moment functional.  It
supplies no second coefficient normalization.  We use only the always-lawful
global projective scaling; setting two coefficients to one would require an
additional symmetry and is not justified.

## 3. Homogeneous support-four certificate

Let `M_m(A,B,C,D,E)=<g^m>`.  Cyclic character makes the other moment orders
automatic, so retain `m=3,6,...,21`.  For each coefficient `J`, set

```text
I_J=(M_3,M_6,M_9,M_12,M_15,M_18,M_21)|_(J=0)              (10)
```

in the four remaining coefficient variables.  At homogeneous degree `21`,
the Macaulay multiplication map has `2926` rows and

```text
dim R_21=binomial(24,3)=2024                               (11)
```

columns.  Exact arithmetic modulo `101` gives rank `2024` for every deletion
`J=A,B,C,D,E`.  With rows ordered by `M_3,...,M_21` and lexicographic weak
compositions, the same selected zero-based row ranges

```text
0-1690, 2146-2310, 2315-2320, 2601-2716, 2821-2866        (12)
```

give square `2024`-minors whose determinants modulo `101` are

```text
deleted J:       A    B    C    D    E
determinant:    27   67   42   50   82.                    (13)
```

All moment entries have simplex degree at most `4*21`; therefore every
denominator divides `(4*21+2)!`.  The guard

```text
101>4*21+2=86                                              (14)
```

makes reduction valid.  Each nonzero determinant in `(13)` is the reduction
of a nonzero rational determinant.  Thus over `Q`

```text
R_21 subset I_J                                            (15)
```

for all five `J`.  In particular, every pure twenty-first power lies in
`I_J`.  A common projective zero on `J=0` would force all four remaining
coordinates to vanish, which is impossible.  This includes all points at
infinity and every lower-support boundary.

Every coefficient vector of support at most four lies on at least one of the
five hyperplanes.  Therefore no nonzero `g` of the form `(1)` with coefficient
support at most four can satisfy all factorial moment conditions.  In fact,
the seven conditions through `M_21` already exclude it.

## 4. Boundary

Full support five is the only open coefficient chart in this eigenspace.  We
use the lawful projective normalization; the named Hesse torus supplies no
further reduction.  The degree-21 five-variable Macaulay map has

```text
13972 rows,       12650 columns.                            (16)
```

It **cannot** have full column rank.  For the first five generator degrees
`3,6,9,12,15`, the complete-intersection Hilbert coefficient in degree `21`
is `1705`, a universal lower bound for the quotient by five forms in five
variables.  Adding `M_18` can remove at most `dim R_3-1=34` dimensions because
the Koszul relation with `M_3` kills one multiplier; `M_21` removes at most
one.  Therefore

```text
dim (R/I)_21 >=1670,              rank Mac_21 <=10980.      (17)
```

[THM-3323](THM-3323-cyclic-quartic-support-five-exact-degree-21-rank.md)
proves that this bound is sharp for the actual moment forms:

```text
rank_Q Mac_21=10980,              dim_Q (R/I)_21=1670.      (18)
```

Its explicit maximal minor modulo `101` lifts to characteristic zero, and a
second guarded prime gives the same rank.  Thus there is no hidden rank defect
at this slice, but the slice cannot prove projective emptiness.  The formal
product-series coefficient is `39` in degree `28` and `-354` in degree `29`,
so degree `29` is only the first full-rank candidate not excluded by this
count; no generic-Hilbert-series or degree-29 rank claim is made.  A proof
needs a guarded exact certificate at a sufficient degree or a different
saturation/affine-chart argument.

Nothing here proves `FC(3)` outside this degree-four cyclic eigenspace, treats
non-eigenvectors, or supplies a Jacobian-conjecture consequence.

QED.
