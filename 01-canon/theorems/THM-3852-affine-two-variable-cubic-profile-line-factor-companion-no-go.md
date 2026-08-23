---
id: THM-3852
title: "Affine two-variable cubic profiles: every line has a nonpolynomial companion"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Affine profiles
  b(A,C) whose depressed-cubic
  branch contains an affine line are classified completely.  After a
  diagonal symplectic normalization, the nonvertical locus is one
  transverse parameter.  Its residual component always has normalization
  missing at least two projective points, including at the unique
  nonreduced line-collision parameter.
source: jc_sparse_direct_search / two-variable cubic-profile escape lane, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (root / jc-cohn3709, 2026-08-23).  The
  audit rederived both complete line coefficient ideals, the symplectic
  scaling, the residual factorization, and the unique repeated-line
  parameter.  It separately checked the localization normalization at
  `tau=0`, irreducibility via exhaustive second-line elimination, distinct
  infinity supports even at the quadratic double root, and smoothness and
  two-place normalization of the collision conic.  The companion verifies the complete
  vertical and nonvertical coefficient ideals, the symplectic rescaling,
  the residual factorization, the unique second-line parameter, the
  infinity packet, and the smooth collision conic.  Normal and optimized
  runs byte-match the frozen 39-gate transcript and both hashes.
related:
  - THM-3847-two-place-cubic-deformation-monogenic-unit-debt
  - THM-3850-nonconstant-cubic-profile-irreducible-branch-puncture-formula
script: 04-computation/jc2_affine_two_variable_profile_line_companion_thm3852.py
output: 05-knowledge/results/jc2_affine_two_variable_profile_line_companion_thm3852.out
script_sha256: e926174d1592121d758b0f54e88ee2d49cdf27149f1d7eb012f1ce86a67558f0
output_sha256: 69002904c2d33d5d0979deefab8640c4dc16b08dddf9eb6a7fb37df30c4036c0
semantic_sha256: 49dd5f040d11c421e395851164f3cf0a5974e460bdd00ce6df824d78ce662a02
hash_basis: raw LF bytes
---

# THM-3852 -- every affine-profile line has a bad companion

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  Work over an
algebraically closed field `k` of
characteristic zero.  For an affine profile

```text
b(A,C)=b0+bA A+bC C                                               (1)
```

put

```text
Delta_b=-27A^2b^2+8AC^3-54ACb+9C^2-54b.                          (2)
```

Suppose the reduced plane curve `V(Delta_b)_red` contains an affine-line
component.  Then exactly one of the following alternatives occurs.

1. **Vertical boundary.**  The line is `C=0` and

   ```text
   b=bC C.                                                        (3)
   ```

2. **Nonvertical family.**  There are unique `n in k*` and `lambda in k`
   such that the line and profile are

   ```text
   L_n=A-(n^2/6)C-n,

   b=lambda L_n-2/n^2-4C/(3n).                                   (4)
   ```

Conversely every profile in `(3)` or `(4)` has the displayed line as a
branch component.

In every case the line is accompanied in `V(Delta_b)_red` by an irreducible
component whose affine normalization omits at least two distinct points of
its smooth projective normalization.  Hence an affine-line component can
occur, but the complete branch packet never consists only of polynomial
curves.

This is a classification inside the affine-profile grammar `(1)`.  It does
not classify arbitrary two-variable profiles and, by itself, does not prove
a statement about every prospective planar Keller map.  Its exact output is
the branch-passport obstruction above.

## 1. Complete line classification

First let the putative line be nonvertical, so it has a unique equation

```text
A=mC+n.                                                          (5)
```

Along this line the affine profile has the form

```text
b|_(5)=e+dC,
e=b0+bA n,                    d=bA m+bC.                          (6)
```

Set `(5)` and `(6)` in `(2)` and equate all coefficients of `C` to zero.
The exact coefficient ideal in `k[e,d,m,n]` has the following Groebner
basis (the displayed squares are retained because this is an ideal identity,
not merely a radical computation):

```text
27d^3n+54d^2+16e,
(3dn+4)^2,
18dm+3dn^2+8n,
dn^3-4m+2n^2,
(6m-n^2)^2.                                                       (7)
```

Because `k` is a field of characteristic zero, `(7)` forces

```text
n!=0,             d=-4/(3n),
m=n^2/6,          e=-2/n^2.                                     (8)
```

Conversely direct substitution of `(8)` kills every coefficient.  Taking
`lambda=bA` and solving `(6)` for `b0,bC` gives exactly `(4)`.  This proves
the nonvertical classification, including uniqueness of `n` and `lambda`
from the line and the `A` coefficient of `b`.

For a vertical line `C=c`, regard `(2)` as a polynomial in `A`.  Its `A^4`
coefficient is `-27bA^2`, so `bA=0`.  With

```text
e0=b0+bC c,
```

the `A^2` coefficient is `-27e0^2`, so `e0=0`; the remaining `A`
coefficient is `8c^3`, so `c=0` and then `b0=0`.  This is precisely `(3)`.
The converse is immediate because every term of `Delta_(bC C)` is divisible
by `C`.

## 2. Symplectic normalization of the nonvertical family

Use the diagonal symplectic change

```text
A=(n/6) A',                    C=(6/n) C',
b'=(n^2/36)b,                  tau=(n^3/216)lambda.              (9)
```

It preserves `AC`, sends `L_n` to a nonzero scalar multiple of

```text
L=A'-6C'-6,                                                       (10)
```

and transforms `(4)` into

```text
b'=tau L-1/18-2C'/9.                                             (11)
```

Moreover

```text
Delta_b(A,C)=(36/n^2) Delta_(b')(A',C'),                         (12)
```

so neither components nor their normalization punctures are changed.  Drop
the primes.  Exact expansion gives

```text
Delta_b=-(1/12)L R_tau,                                          (13)

R_tau=
 324tau^2 A^3
-1944tau^2 A^2C-144tau A^2C
-1944tau^2 A^2-36tau A^2
+16AC^2+648tau AC+8AC+A+18C+648tau+6.                            (14)
```

Thus the displayed line is always present.  The rest of the proof shows
that its companion never has one place at infinity.

## 3. The one-variable seam `tau=0`

At zero transverse parameter,

```text
R_0=A(4C+1)^2+6(3C+1).                                          (15)
```

The two coefficients in `(15)` are coprime, so this degree-one polynomial
in `A` is irreducible.  Its coordinate ring is

```text
k[C,(4C+1)^(-1)],                                                (16)
```

because Bezout for `(4C+1)^2` and `6(3C+1)` expresses the inverse of the
former inside the quotient ring.  Hence its affine normalization is
`P1_C` with the point `C=-1/4` and the point `C=infinity` deleted.  This is
also the direct affine-line slice of the all-one-variable obstruction in
THM-3850.

## 4. The generic residual cubic

Assume `tau!=0`.  Then `R_tau` has total degree three, with homogeneous
part at infinity

```text
(R_tau)_3
=4A(81tau^2 A^2-486tau^2 AC-36tau AC+4C^2).                     (17)
```

A reducible affine cubic over `k` has an affine-line factor.  No vertical
line can divide `R_tau`, because its `A^3` coefficient is `324tau^2`.
For a nonvertical second line `A=MC+N`, equating the coefficients of
`R_tau(MC+N,C)` gives the exact ideal with Groebner basis

```text
9M-233280tau^2+12528tau-206,
9N- 93312tau^2+ 5400tau-122,
(54tau-1)^3.                                                      (18)
```

Therefore a second line exists exactly when

```text
tau=1/54,                    M=N=6,                              (19)
```

and it is then the original line `(10)`.  It follows that `R_tau` is
irreducible whenever

```text
tau!=0,1/54.                                                      (20)
```

Let `P0=[A:C:Z]=[0:1:0]`.  Equation `(17)` contains `P0`.  Its residual
homogeneous quadratic takes the nonzero value `16` at `(A,C)=(0,1)`, and,
because `k` is algebraically closed, has at least one projective zero.
Thus the projective closure of the irreducible curve `R_tau=0` meets the
line at infinity in at least two **distinct projective points**.  Each has
at least one preimage in the projective normalization.  Preimages of
distinct projective points cannot collide, even when an infinity point is
singular or two branches over one support merge.  The affine normalization
therefore omits at least two distinct points.

For reference, the discriminant of the quadratic factor in `(17)` is

```text
139968 C^2 tau^3(27tau+4).                                      (21)
```

At the hostile double-root value `tau=-4/27`, the quadratic infinity packet
coalesces, but it still does not contain `P0`; the two-support conclusion
survives.

## 5. The unique line collision is nonreduced, but its conic is bad

At the exceptional value `(19)`, exact factorization sharpens `(13)` to

```text
Delta_b=-(1/108)L^2 K,
K=A^2-24AC-6A-27.                                                (22)
```

Thus the branch scheme has a doubled line, while its reduced branch is the
union of that line and `K=0`.  The discriminant of `K` as a monic quadratic
in `A` is

```text
144(4C^2+2C+1),                                                  (23)
```

which is not a square in `k[C]`; hence `K` is irreducible.  Its projective
closure

```text
A^2-24AC-6AZ-27Z^2=0                                             (24)
```

is smooth: its homogeneous gradient vanishes only at the irrelevant cone
vertex `(0,0,0)`.  At infinity `(24)` becomes

```text
A(A-24C)=0,                                                      (25)
```

so this smooth conic has two distinct points at infinity.  Its affine
normalization is therefore `P1` minus exactly two points.  The nonreduced
collision does not remove the companion obstruction.

## 6. Vertical boundary and positive controls

For `(3)` with `bC!=0`, THM-3850 gives the exact reduced packet: the line
`C=0 ~= A1` is accompanied by an irreducible rational component whose
normalization is `P1` minus three points.  At the constant endpoint `bC=0`,

```text
Delta_0=C^2(8AC+9).                                              (26)
```

The reduced packet is again a line together with the hyperbola
`8AC+9=0`; the latter is `G_m`, or `P1` minus its two points at infinity.

These vertical examples and the nonvertical line in `(13)` are positive
controls: genuine affine-line branch components do occur.  Equations
`(15)`, `(17)`, `(22)`, and `(26)` identify the information they inevitably
lose--a residual component with at least two missing projective points.

## 7. Exact replay

Run

```bash
python3 04-computation/jc2_affine_two_variable_profile_line_companion_thm3852.py
python3 -O 04-computation/jc2_affine_two_variable_profile_line_companion_thm3852.py
```

Both commands must byte-match
`05-knowledge/results/jc2_affine_two_variable_profile_line_companion_thm3852.out`.
The assertion-free companion performs 39 exact gates, including equality of
both coefficient ideals rather than only testing the displayed solutions.
The raw-LF SHA-256 hashes are recorded in the metadata.
