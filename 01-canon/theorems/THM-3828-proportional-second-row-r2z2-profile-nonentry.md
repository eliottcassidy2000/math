---
id: THM-3828
title: "Proportional second-row R2Z2 profiles do not enter the cubic pseudo-plane Darboux packet"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  In the first
  canonical second-row extension of
  THM-3821, suppose the nonzero r^2z^2 profile X in the first coordinate and
  its mate Y are proportional, Y=lambda X, and the preceding rz^2 profiles
  are aligned in the same target direction, T=lambda S.  Then the pair is not
  Darboux.  The N!=0 branch enters an exact 10,5,4,3 ladder and fails at a
  half-integral root address/first integral; N=0 fails the arm division.  The
  genuinely two-sided L=T-lambda S!=0 branch has a new 10/7 tower and remains
  OPEN, as do one-sided X=0,Y!=0 and further second-row slots.
source: jc_zero_debt_lift / cubic-pseudoplane second-row profile lane, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (root / jc-cohn-boundary, 2026-08-23).
  The audit independently checked the canonical-bucket typing, the
  nonvanishing arm argument, both UFD valuation ladders, the optional H
  rung and integrated S profile, the fixed-base terminal factor, both local
  coefficients at a nonzero root, the origin and first-integral closure,
  the separately recomputed N=0 branch, and the exact OPEN boundary.  Normal
  and optimized runs byte-match the frozen transcript and the raw hashes
  agree.  The deterministic companion has 43 active
  gates checking the Poisson Casimir, monic normal form, top Wronskian,
  seven descending canonical buckets, target-difference typing, the 5/2 and
  combined 10/5/4 valuation families, optional 10/3 rung, integrated S
  profile, full fixed-base terminal equation, nonzero-root address and
  half-integral next coefficient, confluent first integral, separately
  recomputed N=0 arm remainder, and the complementary OPEN 10/7 tower.
depends_on:
  - THM-3821-cubic-pseudoplane-rz2-odd-ladder-terminal-riccati-gate
related:
  - THM-3814-nodal-rz-kummer-profile-degree-gate
  - THM-3812-nodal-arm-coefficient-second-normal-profile-nonentry
script: 04-computation/jc2_cubic_pseudoplane_proportional_second_row_thm3828.py
output: 05-knowledge/results/jc2_cubic_pseudoplane_proportional_second_row_thm3828.out
script_sha256: 4987fa7ba479b256dca27eb8d2cb051f4c54c62af54612783243da879b88f03e
output_sha256: 316bea85e4e71be29c25ab847f91d127e05a49b45404d0bc1edfb232283dd538
semantic_sha256: 170613b16552465adbb4e26734b0f50828b1e77e8c123d59ab2e2c80e38dfeb1
hash_basis: raw LF bytes
---

# THM-3828 -- the proportional second row is empty

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  Let `k` be an
algebraically closed field of characteristic zero and put

```text
B=k[r,z,e]/(r^2e-z^3+r),
{r,z}=3r^2,       {r,e}=9z^2,       {z,e}=3+6re.       (1)
```

For arbitrary profiles

```text
f,g,h,kappa,p,q,S,T,X,Y in k[e]
```

define the first canonical second-row extension of THM-3821:

```text
A=e^2-z/3+r g+z^2 f+rz p+rz^2 S+r^2z^2 X,
C=e^3-e-ez/2+r h+z^2 kappa+rz q+rz^2 T+r^2z^2 Y.      (2)
```

Suppose `X!=0`.  If

```text
Y=lambda X,                  T=lambda S               (3)
```

for some `lambda in k`, then

```text
{A,C}!=1.                                               (4)
```

Consequently, a hypothetical Darboux pair with `X!=0` first has
`Y=lambda X` forced by its top bucket, but it must lie in the complementary
branch

```text
L=T-lambda S !=0.                                      (5)
```

That OPEN branch necessarily begins with the new `10/7` tower in `(29)`
below.  The theorem does not exclude it.

## 1. Canonical normal form and target differences

The monic presentation

```text
B=k[r,e][z]/(z^3-r^2e-r)                              (6)
```

is free over `k[r,e]` with basis `(1,z,z^2)`.  Thus reduction of
`{A,C}-1` has no quotient loss.  Its highest coefficient is

```text
[r^7]=-30e^2(-XY'+YX').                               (7)
```

If a hypothetical pair has `X!=0`, `(7)` says `(Y/X)'=0` in `k(e)`, hence
`Y=lambda X`.  This explains why the first condition in `(3)` is intrinsic,
not a chosen specialization.

Assume the aligned condition `T=lambda S` and put

```text
M=kappa-lambda f,
N=q-lambda p,
H=h-lambda g,
D=3e^2-2lambda e-1.                                  (8)
```

The arm and next four nontrivial canonical buckets are exactly

```text
D f=2eM-1/12,                                         (9)

[r^5]/(6e)=-2eMX'+5eXM'+2MX,                         (10)

[r^4z^2]/(15e)=2XN'-NX',                             (11)

[r^4z]/3=10eXH'-3eHX'-2XH,                           (12)

[r^4]/3=
 -4e^2MS'+7e^2SM'+2eMS
 -6eMX'+18eXM'+4MX.                                  (13)
```

Each is obtained from a distinct coefficient in the unique normal form.
No division by a profile or truncation of lower buckets has occurred.

The polynomial `M` cannot vanish: `(9)` would make the nonconstant quadratic
`D` times a polynomial equal the nonzero scalar `-1/12`.

## 2. The N-nonzero branch enters 10,5,4,3

Divide `(10)` by the nonzero rational function `XM`.  It is equivalent to

```text
(X^2/(e^2M^5))'=0.                                   (14)
```

The characteristic-zero constant field of `k(e)` is `k`.  Unique
factorization gives first

```text
X=alpha e v^5,                  M=beta v^2            (15)
```

with `alpha,beta!=0` and `v in k[e]`.

Suppose now `N!=0`.  Equation `(11)` gives

```text
(N^2/X)'=0.                                           (16)
```

Combining `(14),(16)`, at a prime other than `e` the valuation equation is
`4 ord(N)=5 ord(M)`, while at `e` it is
`4 ord_e(N)=2+5 ord_e(M)`.  Hence, after absorbing nonzero constants,

```text
X=alpha e^6v^10,
M=beta e^2v^4,
N=delta e^3v^5,                alpha beta delta!=0.  (17)
```

The unused equation `(12)` says either `H=0` or, by the same UFD argument,

```text
H=theta e^2v^3.                                      (18)
```

Thus the second row has exposed the ladder `10,5,4,3`, distinct from
THM-3821's `7,5,4,3,1` anatomy.

Substitution of `(17)` in `(13)` leaves

```text
3alpha e^6v^10v'+alpha e^5v^11
 +7eSv'-evS'+4Sv=0.                                 (19)
```

Equivalently,

```text
(S/(e^4v^7)-alpha e v^3)'=0,
S=alpha e^5v^10+gamma e^4v^7                       (20)
```

for `gamma in k`.

## 3. The fixed-base terminal equation

The next coefficient is the first place where one must not silently replace
`C` by `C-lambda A`: that shear changes the fixed restrictions in `(2)`.
Substitution of `(17),(20)` in the actual `r^3z^2` bucket gives

```text
0=
 3delta e^2v^5(3ev'+v)
 -2e(3e-2lambda)v'
 +2lambda v.                                         (21)
```

Let `rho!=0` be a root of `v` of multiplicity `m`, and write
`v=(e-rho)^m u` with `u(rho)!=0`.  The order-`m-1` coefficient in `(21)`
first forces

```text
3rho-2lambda=0.                                      (22)
```

At that apparent address, the coefficient one order farther is exactly

```text
3rho u(rho)(1-2m).                                   (23)
```

The first term in `(21)` begins only at order `6m-1>m`, so no hidden term
can cancel `(23)`.  Its sole formal resonance is `m=1/2`, impossible for a
polynomial root.  Hence `v` has no nonzero root.

If `v` is constant, `(21)` has nonzero `e^2` coefficient
`3delta v^6`.  Otherwise algebraic closure gives `v=c e^d` with `d>=1`;
absorb `c` into the constants.  The coefficient of `e^d` in `(21)` is

```text
2lambda(2d+1),                                       (24)
```

so `lambda=0`.  Divide `(21)` by `3e^2`.  The surviving equation is

```text
delta v^5(3ev'+v)-2v'=0,

((1+delta e v^5)/v^2)'=0.                            (25)
```

Therefore `1+delta e v^5=c_0v^2` for a constant `c_0`.  A nonconstant
polynomial over `k` has a root, where this identity reads `1=0`.  This
closes the entire `N!=0` branch.

## 4. The N-zero specialization is recomputed

Set `N=0`; do not divide `(11)` by `N`.  The actual `r^3z^2` bucket now is

```text
e(3e-2lambda)X'-2(9e-lambda)X=0,

(X/[e(3e-2lambda)^5])'=0.                            (26)
```

Thus `X=chi e(3e-2lambda)^5`.  Comparing with `(15)` gives, after absorbing
constants,

```text
v=3e-2lambda,             M=beta(3e-2lambda)^2.
```

Polynomial division of the numerator in the arm `(9)` gives

```text
rem_D(2beta e(3e-2lambda)^2-1/12)
  =(72beta e-48beta lambda-1)/12.                    (27)
```

Its linear coefficient is `6beta!=0`, so `(9)` has no polynomial solution.
This closes `N=0` and proves `(4)`.

## 5. Exact surviving boundary

Return only to `Y=lambda X`, and now put `L=T-lambda S`.  Before imposing
`L=0`, the `r^6` bucket is

```text
7eLX'-10eXL'-2LX=0,

(X^7/(e^2L^10))'=0.                                 (28)
```

If `L!=0`, UFD valuations give the exact necessary tower

```text
X=alpha e^6w^10,               L=beta e^4w^7.        (29)
```

This `10/7` branch is **OPEN**.  Equations `(28)--(29)` are recorded to make
the theorem's equality boundary constructive: the aligned branch fails, but
the first misaligned branch carries a genuinely new divisor sidecar rather
than merely more coefficients.

The theorem also does not cover `X=0,Y!=0`, added `r^2` or `r^2z` profiles,
higher canonical rows, different arm immersions, arbitrary Darboux pairs on
`B`, or planar Jacobian maps in general.

## 6. Exact controls

The companion named in the metadata recomputes the full bracket in `(2)`,
reduces it in the monic basis `(1,z,z^2)`, freezes the top Wronskian and every
bucket `(9)--(13),(21),(26),(28)`, verifies all UFD exponent families and
their tower substitutions, checks the integrated profile `(20)`, and derives
the fixed-base terminal without a target-shear loss.  It separately checks
the first and second local root coordinates `(22)--(23)`, the origin address
`(24)`, first integral `(25)`, the exceptional `N=0` remainder `(27)`, and
the OPEN `10/7` tower `(29)`.  It uses no finite-field or bounded-degree
inference.  **QED.**
