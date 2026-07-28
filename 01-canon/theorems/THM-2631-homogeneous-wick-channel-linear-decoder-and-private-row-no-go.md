---
id: THM-2631
title: "Homogeneous Wick-channel linear decoder and private-row no-go"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  A Gaussian moment polynomial is one homogeneous scalar row across all
  balanced Wick channels of its degree.  A universal scalar linear
  combination of raw moment polynomials recovers a prescribed level-tagged
  channel monomial if and only if that moment degree has exactly one balanced
  channel.  Equivalently, the nonnegative channel-incidence bank has a left
  inverse exactly when every decoded column has a private row.  Signed scalar
  coefficients do not improve this: degree grading prevents other moment
  levels from cancelling channels in the selected degree.  For
  P=a Z^6+b W^2+c W^18, channels occur exactly at m=4j and are
  (j+2t,3(j-t),t), 0<=t<=j, so every nonempty positive level has
  j+1>=2 channels.
  This strengthens MISTAKE-211 to an all-level linear-decoder obstruction but
  leaves nonlinear moment-ideal, resultant, cumulant, confluent-Hankel, and
  whole-face Frobenius mechanisms untouched.  NC2 and GMC(2) were already
  proved by THM-2022; this theorem is a method boundary, not a new closure.
source: wild-holotopy-mining-2026-07-28-wick-channel-decoder
depends_on:
  - THM-2022-gmc2-frobenius-lowest-balanced-face
related:
  - THM-1645-gmc2-angular-layer-is-dvdk-the-gap-is-purely-radial
  - THM-1770-the-localisation-lemma-as-first-return-renewal-and-the-pair-only-closure
  - THM-2020-gmc2-finite-place-channel-separation
  - THM-2624-two-clock-root-tomography-and-disjoint-carrier-holotopy-boundary
  - THM-2639-gmc-equal-mass-two-rung-persistent-collision-certificate
  - HYP-8765-gmc2-radial-channel-return-tower
script: 04-computation/gmc2_wick_channel_linear_decoder_thm2631.py
output: 05-knowledge/results/gmc2_wick_channel_linear_decoder_thm2631.out
script_sha256: cba25dfb0a8385172555ffd5c4aa008bb0fb02c2e17eb2a1a1c00f5730cf8392
output_sha256: 403d6196e6a7ef39e0e3c548908a8f2520add266a837de56429a37295496d846
hash_basis: LF-normalized bytes
---

# THM-2631 -- raw moments linearly decode a channel exactly at singleton levels

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

MISTAKE-211 showed that distinct first-return channels are not separate
equations inside one scalar Gaussian moment.  The obstruction is not confined
to the first return.  Homogeneous degree makes it categorical: adding any
number of other raw scalar moments still cannot linearly separate two channels
that occur in the same degree.

This theorem identifies the precise linear boundary.  It also explains why
THM-2022 succeeds by preserving and amplifying a complete balanced face rather
than trying to select one of its channels.

## 1. Exact Wick rows

Let `Z` be a circular complex Gaussian and write `W=Zbar`.  For a polynomial
with exact support

```text
P=sum_(i=1)^k c_i Z^(a_i) W^(b_i),       c_i!=0,
q_i=a_i-b_i,
```

put

```text
B_m={r in N^k: |r|=m and q dot r=0},
A(r)=sum_i a_i r_i=sum_i b_i r_i.                         (1)
```

The direct multinomial and Wick expansion, also equation (1) of THM-2022,
is

```text
M_m(c):=E[P^m]
       =sum_(r in B_m) w_(m,r)c^r,

w_(m,r)=m! A(r)! / prod_i r_i! >0.                        (2)
```

Thus `M_m` is one homogeneous polynomial row of total coefficient degree
`m`.  Its columns are the level-tagged channels `(m,r)`, and every channel in
that row has a strictly positive rational-integer weight.

## 2. The scalar linear-decoder theorem

Fix a finite set `S` of moment degrees and a channel `r0 in B_m`, `m in S`.
A **universal scalar linear decoder** for `(m,r0)` means scalars
`lambda_n`, `n in S`, independent of the coefficient variables, such that

```text
sum_(n in S) lambda_n M_n(c)=c^r0                         (3)
```

as a polynomial identity in `Q[c_1,...,c_k]`.  More generally one may replace
the right side by a nonzero scalar multiple of `c^r0`; this makes no
difference.

> **Theorem.**  A decoder (3) exists if and only if `B_m={r0}`.  In that case
> the unique relevant coefficient is `lambda_m=1/w_(m,r0)`; all nonzero
> moment rows of other degrees have coefficient zero.

Indeed, the polynomial ring is graded by total coefficient degree.  Taking
the degree-`m` component of (3) gives

```text
lambda_m sum_(r in B_m) w_(m,r)c^r=c^r0.                  (4)
```

Distinct multiplicity vectors give distinct coefficient monomials.  Every
weight in (4) is nonzero, so (4) can hold only when the row contains the one
channel `r0`.  Conversely a singleton row is decoded by division by its Wick
weight.  Every other nonempty homogeneous component on the left of (3) must
vanish separately.

The same proof allows a formal degreewise sum of arbitrarily many raw
moments: its degree-`m` component is still (4), and every other nonempty
homogeneous component separately forces its scalar coefficient to vanish.
Coefficients on identically zero moment rows are irrelevant.  Hence no later
or earlier raw moment can repair a collision at level `m`.

## 3. Positive matrices, private rows, and signed cancellation

For a finite degree bank, form its level-tagged incidence matrix

```text
A_(n,(m,r)) = 1_(n=m) w_(m,r).                            (5)
```

For an arbitrary nonnegative real matrix `A`, there exists a nonnegative left
inverse `L>=0`, `LA=I`, if and only if every column `j` has a **private row**:

```text
A_(k_j,j)>0,              A_(k_j,i)=0 for every i!=j.     (6)
```

To prove necessity, inspect the off-diagonal zeros in `LA=I`.  Since every
summand is nonnegative, each row used by the `j`-th decoder must vanish in
every other column; the diagonal one forces at least one such row to be
positive at `j`.  Conversely, placing `1/A_(k_j,j)` on one private row for
each column constructs a nonnegative left inverse.

In (5), the sole row at level `m` is positive on every column `(m,r)`.  It is
private exactly when `|B_m|=1`.  Therefore a bank decoding all its channels
has a nonnegative left inverse exactly when every nonempty selected moment
level is singleton.

Signed scalar coefficients do not help.  The block at a colliding level is a
single row of rank one with at least two columns, and grading prevents rows
from other levels from entering that block.  This differs sharply from
THM-2624: there multiple independent rows can give full signed rank even
without a positive inverse.  Here there is only one row per homogeneous
degree, so even signed linear tomography fails.

## 4. The MISTAKE-211 support is hostile at every return

Take

```text
P=a Z^6+b W^2+c W^18.                                    (7)
```

Its charge vector is `(6,-2,-18)`.  If
`r=(r1,r2,r3) in B_m`, balance and total degree say

```text
3r1-r2-9r3=0,             r1+r2+r3=m.                    (8)
```

Eliminating `r2` gives

```text
m=4(r1-2r3).                                               (9)
```

Consequently there are no channels unless `m=4j`.  At `m=4j` every solution
is uniquely

```text
(r1,r2,r3)=(j+2t, 3(j-t), t),             0<=t<=j.        (10)
```

Thus

```text
|B_(4j)|=j+1 >=2                                           (11)
```

at every nonempty positive moment level.  Explicitly,

```text
M_(4j)=sum_(t=0)^j
  (4j)! [6(j+2t)]!
  ----------------------------------
  (j+2t)! [3(j-t)]! t!
  a^(j+2t)b^(3(j-t))c^t.                                  (12)
```

There is therefore no universal scalar linear decoder for **any** balanced
channel of (7), even if all raw moments are supplied.

At the first return `j=1`, equation (12) is the exact MISTAKE-211 identity

```text
M_4=4*6!*a*b^3+4*18!*a^3*c.                              (13)
```

The torus coefficient choice

```text
a=b=1,                 c=-6!/18!=-1/8,892,185,702,400     (14)
```

makes `M_4=0` by cross-channel cancellation.  Equation (11) is stronger than
that one cancellation: it proves there is no later singleton return at which
a linear method could restart.

## 5. What the no-go does not exclude

The theorem concerns scalar linear combinations of the raw moment
polynomials.  It does **not** exclude any operation that changes the graded
problem, including:

- polynomial multipliers followed by moment-ideal elimination or saturation;
- resultants, localized cumulants, and products or determinants of moments;
- confluent factorial-Hankel or Vandermonde systems that create multiple
  equations inside one aligned degree;
- finite-place normalization and channel congruences such as THM-2020; or
- THM-2022's whole-face Frobenius amplification.

These are not cosmetic exceptions.  Multiplying equations can align degrees;
resultants and cumulants are nonlinear; Hankel determinants use products of
rows; and THM-2022 deliberately retains the complete tied face, reduces it
modulo a good prime, and obtains its nonzero Frobenius power.  None asks a
single homogeneous scalar row to behave as several equations.

THM-2639 realizes the first option sharply.  Its free equal-mass two-ray face
has `k+1` channels at every nonempty level `k*ell`, so the present linear
decoder fails forever.  Nevertheless a polynomial multiple of `T_ell`
eliminates the mixed and `y^2` terms of `T_(2ell)` and leaves a nonzero torus
monomial.  Thus the linear no-go and nonlinear two-rung certificate meet at
an exact boundary rather than competing conclusions.

In the holotopy language suggested by THM-2624, the balanced channels form a
multi-point fibre over one degree, while `M_m` collapses that entire fibre to
one weighted scalar.  Different degrees are disjoint graded components, not
overlapping charts with transition data.  Stacking them therefore supplies no
descent map inside the fibre.  A viable effective HYP-8765 refinement must
manufacture new aligned equations or preserve the whole face; it cannot rely
on raw-moment linear isolation.

## 6. Exact companion and scope

The companion enumerates the channels of (7) through degree `48`, verifies
the parametrization (10) for `j=1,...,12`, replays (13)--(14), and exhausts
the private-row criterion on every binary matrix of sizes
`2x2,2x3,3x2,3x3`.  The exact decodable censuses are

```text
(16,2), (64,0), (64,18), (512,6),                         (15)
```

where each pair is `(all matrices, matrices with a nonnegative left inverse)`.

An independent hostile audit rederived the Wick weights, the graded decoder
equivalence, both directions of the private-row criterion, the complete
all-positive-level parametrization (10), the exact fourth-moment
cancellation, and every binary census in (15).  It also replayed normal and
optimized execution against the stored transcript and confirmed both
LF-normalized hashes.  The audit's two scope repairs -- excluding the
singleton degree-zero row and treating an infinite bank degreewise -- are
incorporated above.

This theorem proves no new NC2 or GMC(2) conclusion.  Those statements are
already PROVED by THM-2022.  It strengthens the corrected interpretation of
MISTAKE-211 and narrows the effective-method frontier: raw scalar moments do
not furnish a linear channel selector unless the decoded degree was already
singleton.  THM-2639 is the proved nonlinear survivor.

`QED.`
