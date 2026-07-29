---
id: THM-2948
title: "Pure resultant negative seam and squarefree core atlas"
status: >
  PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT
  HOSTILE AUDIT.  After dividing the original degree-seven Macaulay
  determinant by the full selected-chart factor
  q200^6*c300*K, the calibrated bare resultant on every normalized
  four-slot support of widths six through ten has repeated-root
  radical
  prod_(j=1)^a(2n+2j-1) prod_(r=a)^M(n+r).  This is an exact
  110-family statement about the genuine resultant, not the common
  positive Pluecker factor times the resultant.  Every displayed root
  is negative; after its full seam multiplicity is removed, the core
  is squarefree, seam-coprime and coefficientwise positive.  Together
  with THM-2945 this gives an independent first-window SFC(4) closure
  through width ten.  No all-width transversality or seam formula is
  claimed.
source: codex-gmc-uniform-width-extension-2026-07-29
depends_on:
  - THM-2925-general-width-terminal-pole-cancellation-and-macaulay-degree-law
  - THM-2942-macaulay-extraneous-flag-factor-and-pluecker-mutation
  - THM-2943-width-seven-eight-two-chart-macaulay-resultant-closure
  - THM-2945-nonnegative-complete-intersection-norm-and-repeated-divisor-gate
related:
  - THM-2927-general-width-flagged-macaulay-leading-coefficient
  - THM-2944-width-nine-ten-two-chart-macaulay-resultant-closure
  - THM-2946-full-macaulay-maximal-minor-gcd-and-chart-free-resultant
script: 04-computation/gmc_pure_resultant_negative_seam_atlas_thm2948.py
output: 05-knowledge/results/gmc_pure_resultant_negative_seam_atlas_thm2948.out
script_sha256: e691f020d8447382356b4f96d667b14532d1b71363092c3ac1baf96538f006af
output_sha256: bb7c92cecc716a05465079631c62b812cea783ece0b11744b653d2b69858786d
calibrated_constructor_sha256: d2f8afeba7dd6c7950405a4845d7bf112b6c9872dd8161146446be8bbdaae0ba
hash_basis: LF-normalized bytes
---

# THM-2948 -- pure resultant negative seam and squarefree core atlas

**PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT
HOSTILE AUDIT.**

## 1. Statement

Let

```text
L:C[s] -> C,                         L(s^j)=j!,          (1)
```

and fix

```text
6<=M<=10,                 1<=a<b<M.                    (2)
```

On the translated support

```text
(n,n+a,n+b,n+M),                                      (3)
```

normalize by `f_j=s^j/j!`, eliminate the mean against the top
coordinate, and clear the order-two, order-three and order-four
mean-difference forms by the positive-denominator polynomials of
PROVED THM-2925.  Call the resulting ternary forms

```text
Q_(M,a,b)(n),          C_(M,a,b)(n),          F_(M,a,b)(n). (4)
```

Write

```text
R_(M,a,b)(n)=Res(Q_(M,a,b),C_(M,a,b),F_(M,a,b))       (5)
```

for the primitive positive-leading associate in `Z[n]` under the
resultant convention fixed below.  Then

```text
deg R_(M,a,b)=46M-26.                                 (6)
```

Define the negative Gamma seam

```text
S_(M,a)(n)
 =product_(j=1)^a (2n+2j-1)
  product_(r=a)^M (n+r).                              (7)
```

For every one of the

```text
C(5,2)+C(6,2)+C(7,2)+C(8,2)+C(9,2)
 =10+15+21+28+36=110                                  (8)
```

families in `(2)`,

```text
rad gcd(R_(M,a,b),partial_n R_(M,a,b))
 =primitive associate of S_(M,a).                     (9)
```

In particular, the radical is independent of the middle offset `b`
and every repeated root is strictly negative.

Remove from `R_(M,a,b)` the full multiplicity of every linear factor
displayed in `(7)`, and call the primitive positive-leading quotient
`U_(M,a,b)`.  Then

```text
gcd(U,partial_n U)=1,              gcd(U,S)=1,         (10)
```

and, in the same finite atlas,

```text
U(0)>0,                 every coefficient of U is >=0. (11)
```

Thus `(10)` is genuinely a squarefree wall-removed resultant core,
not merely a squarefree part of a chart cofactor.

## 2. The load-bearing bare-resultant calibration

Use the original 36-row chart of PROVED THM-2942.  In its coefficient
notation,

```text
Delta_0
 =q200^6*c300*K(Q,C)*Res(Q,C,F),                      (12)
```

where

```text
K
 =c120*q200^2-c210*q110*q200
  -c300*q020*q200+c300*q110^2.                        (13)
```

The companion defines `(5)` by the exact division

```text
R=Delta_0/(q200^6*c300*K).                            (14)
```

It checks zero remainder in every family before taking a derivative
or a gcd.  Equation `(14)` is load-bearing.  In particular, the script
does **not** use

```text
gcd(Delta_0,Delta_1)/(q200^5*c300),                   (15)
```

because `(15)` retains the generally nonconstant positive factor

```text
G_extr=gcd(K,P_alt).                                  (16)
```

The repeated-divisor formula in `(9)` is therefore a statement about
the calibrated genuine resultant.  Negative repeated roots introduced
by `G_extr` are not attributed to the norm.

The source/target map is now explicit:

```text
selected Macaulay determinant
 -> divide the selected flag
 -> genuine complete-intersection resultant
 -> derivative gcd on the depth line.                 (17)
```

The first arrow forgets the chosen Pluecker coordinate while preserving
the common-projective-zero predicate.  The last arrow retains exactly
the ramified/repeated depth divisor.  PROVED THM-2945 is the sidecar
which turns this derivative datum into a nonnegative-ray zero test.

## 3. Why these factors are Gamma seams

For offsets `e_i` the normalized order-`m` tensor coefficient before
mean-difference inclusion--exclusion is

```text
R_e(n)
 =(mn+1)_(sum e_i)/product_i(n+1)_(e_i).              (18)
```

Thus the coefficient family has two inherited kinds of analytic
boundary:

```text
n=-r                                      denominator seams,
2n=-(2j-1)                               order-two seams. (19)
```

PROVED THM-2925 shows how the full mean-difference tensor cancels the
forbidden terminal principal part and supplies polynomial forms after
clearing.  The exact resultant atlas determines which boundaries
remain repeated after that cancellation:

```text
j=1,...,a,                    r=a,...,M.              (20)
```

The smallest positive support offset `a`, rather than the middle
offset `b`, is therefore the first active boundary in both strings.
This explains the `b`-independence in `(9)` without confusing it with
an all-width proof: the pole law predicts the candidate seam, while
the absence of any additional repeated factor is the finite exact
content of this theorem.

Every factor in `(7)` occurs in `R` with multiplicity at least two.
The companion divides each one repeatedly until it no longer divides,
then checks `(10)` directly.  It does not infer core squarefreeness
from a coefficient-sign test.

## 4. Consequence on the physical ray

PROVED THM-2945 orients the complete-intersection resultant so that,
after the strictly positive THM-2925 clearing,

```text
R_(M,a,b)(n)>=0                         for n>=0.       (21)
```

The exact endpoint check gives

```text
R_(M,a,b)(0)>0.                                      (22)
```

All roots of `(7)` lie in `(-infinity,0)`.  Equations `(9)`,
`(21)` and the repeated-divisor gate of THM-2945 therefore imply

```text
R_(M,a,b)(n)>0                         for n>=0.       (23)
```

Hence `Q,C,F` have no common projective zero at any nonnegative
integer depth.  By the standard mean-elimination and projective
reduction, every nonzero polynomial supported on

```text
(n,n+a,n+b,n+M),             6<=M<=10,               (24)
```

has a nonzero factorial moment among orders one through four.

Widths six through eight already have independent closures in proved
canon.  For widths nine and ten, `(23)` is an independent
bare-resultant route parallel to the two-chart common-content route of
THM-2944.

Coefficient positivity in `(11)` gives a second finite proof that the
wall-removed core has no nonnegative real zero.  It is recorded as an
audited boundary feature, not used to extrapolate `(23)` to larger
width.

## 5. Exact finite proof

For every family in `(8)`, the companion:

1. rebuilds every coefficient of `Q,C,F` from the proved THM-2925
   normalized-tensor constructor;
2. interpolates `Delta_0` through the proved degree bound `58M-36`;
3. performs the exact division `(14)` and checks `(6)`;
4. computes `gcd(R,R')`, its squarefree radical, and checks `(9)`;
5. removes every full seam multiplicity and checks both gcds in
   `(10)` and the sign statement `(11)`; and
6. checks one determinant beyond the interpolation grid for each of
   the `110` families.

The degree ranges of the repeated divisor and wall-removed core are:

| width | repeated degree | core degree |
|---:|---:|---:|
| 6 | 60--122 | 121--183 |
| 7 | 65--144 | 144--223 |
| 8 | 67--169 | 164--266 |
| 9 | 72--194 | 184--306 |
| 10 | 74--218 | 205--349 |

The width-specific records include the full seam multiplicity vector
and separate digests of `R`, `gcd(R,R')`, its radical and the core.
No factorization over a floating field is used.

## 6. Exact remaining all-width obstruction

Formula `(18)` explains why the factors in `(7)` are natural.  It
does not prove that they exhaust the repeated divisor at arbitrary
width.

Away from the seam, let

```text
A_n=C[x,y,z]/(Q_n,C_n),                 dim A_n=6.     (25)
```

On a local splitting cover, the resultant is the product of the six
values of `F_n` on the moving `Q_n=C_n` intersection.  A non-seam
repeated root can arise in either of two ways:

1. two distinct moving points have `F_n=0` at the same depth; or
2. one moving point has a tangential zero of `F_n`.

Thus the honest all-width target is:

```text
after saturation by S_(M,a), every zero of the incidence projection
{Q=C=F=0} -> A^1_n is unique and transverse.          (26)
```

Equivalently, the wall-removed norm must be squarefree.  The cheapest
decisive symbolic test is a Jacobian/Fitting saturation of the
incidence ideal by `(7)`, together with a collision test on two copies
of the length-six complete intersection.  Neither THM-2925's pole law
nor selected-chart coprimality proves `(26)`.

This theorem therefore leaves open:

```text
the all-width seam formula,
uniform wall-removed transversality,
width at least eleven,
shifted moment windows, and SFC(5).                    (27)
```

## 7. Reproduction

Run

```text
python 04-computation/gmc_pure_resultant_negative_seam_atlas_thm2948.py
python -O 04-computation/gmc_pure_resultant_negative_seam_atlas_thm2948.py
```

Both modes must byte-match the stored transcript and the LF-normalized
hashes in the front matter.

No promotion is requested before an independent hostile audit checks
the calibration `(14)`, the finite universe, every gcd direction, the
positive-ray consequence and the all-width stopping boundary.
