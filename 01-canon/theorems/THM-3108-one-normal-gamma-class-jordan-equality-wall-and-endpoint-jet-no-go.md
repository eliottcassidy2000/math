---
id: THM-3108
title: "One-normal Gamma-class Jordan equality wall and endpoint-jet no-go"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  The
  length-two nilpotent attachment has an exact persistent-zero equality
  wall z^2=s even when the repair is a unit, the error is o(z), and every
  scalar term is a proper positive-ratio Gamma term.  For the actual
  one-normal carrier T_C=tau_C^(m+1), the two cancelling norm terms are the
  identical hypergeometric sequence T_C^2.  A cyclic upper quotient
  v^(m+1)=u realizes the wall with the degree-compatible relation
  e^2=tau_C^2 u v^(m-1).  Hence finite Gamma-sum quasianalyticity together
  with the special Jordan type or unshifted Smith bars cannot by itself
  prove that a physical zero-child resultant is nonidentically zero; the
  physical endpoint/characteristic jet remains load-bearing.
source: root/multiscale-newton-flag-2026-08-02
depends_on:
  - THM-3099-finite-gamma-sum-quasianalyticity-and-remote-resultant-zero-dichotomy
  - THM-3105-weighted-jordan-repair-newton-attachment-spectrum
related:
  - THM-2960-local-smith-jet-fitting-barcode-and-negative-depth-chamber-atlas
  - THM-3101-zero-primary-next-moment-unit-suspension-and-conditional-codimension-five-spectrum
script: 04-computation/gmc_one_normal_gamma_jordan_equality_wall_thm3108.py
output: 05-knowledge/results/gmc_one_normal_gamma_jordan_equality_wall_thm3108.out
script_sha256: 4eac76ead23aecf8ea56957ccbb4179e336202f4b6b0106e285518a3cc457873
output_sha256: 28854479072b7ee7c87897e3e699e0b1086d455f98f3a76c0c59f05b80f235fc
hash_basis: LF-normalized bytes
---

# THM-3108 -- one-normal Gamma-class Jordan equality wall

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-3099 reduces a fixed one-parameter remote resultant to a sharp
dichotomy: exact identity or eventual nonvanishing.  THM-3105 identifies the
weighted shifted-determinant spectrum needed to inspect the identity branch.
Those two statements do not, by themselves, rule that branch out.  At a
length-two Jordan block the repair term and a single closing error can be
the same Gamma monomial and cancel exactly.

The obstruction below uses the actual one-normal repair carrier and also
fits the normal degrees after imposing the upper cyclic relation.  It is a
formal local complete-intersection hostile, not an example arising from a
factorial support.  Its purpose is to isolate the extra datum that a positive
physical theorem must compute.

## 1. The exact one-normal Gamma carrier

Fix

```text
m>=3,                 p=m+1,                 n>=0.     (1)
```

For an integer remote offset `C>=0`, put

```text
U_r(C)=(rn+1)_(rC)/(n+1)_C^r,

tau_C=U_m(C)/U_p(C)^(m/p),
T_C=tau_C^p=U_m(C)^p/U_p(C)^m.                        (2)
```

The fractional power in `tau_C` is only a positive normalization.  The
scalar norm carrier `T_C` is the proper Gamma term

```text
T_C=
 [((m(n+C))!/(mn)!)^p]
 [((pn)!/(p(n+C))!)^m].                               (3)
```

Indeed, the powers of `(n+C)!/n!` cancel between numerator and denominator.
Its consecutive quotient is the positive rational function

```text
T_(C+1)/T_C
 =[(m(n+C)+1)_m]^p/[(p(n+C)+1)_p]^m.                 (4)
```

Stirling gives

```text
T_C~kappa_(m,n)^p (m/p)^(mpC) C^(1/2),
kappa_(m,n)>0,                                        (5)
```

so `T_C->0`.  Both `T_C` and `T_C^2` therefore lie exactly in THM-3099's
positive rational-ratio class.

## 2. The minimal shifted-attachment identity

Let

```text
J_2=[0 1],                 E_21=[0 0].                (6)
    [0 0]                       [1 0]
```

For independent variables `s,z`,

```text
det(J_2+s E_21+z I)=z^2-s.                            (7)
```

On the monomial arc

```text
z=T_C,                       s=T_C^2,                 (8)
```

the error satisfies `s=o(z)`, but

```text
det [T_C   1 ]
    [T_C^2 T_C]=T_C^2-T_C^2=0                        (9)
```

for every `C`.  Before repair, `J_2+T_C^2E_21` has special fibre `J_2`, a
unit entry, determinant `-T_C^2`, and hence Smith exponents `(0,2)` over the
`T_C`-adic valuation.  Formula `(7)` is exactly THM-3105's longest-block
attachment face.  The coefficient `1` on its wall is not recorded by the
special Jordan type or by the unshifted Smith multiset.

As a Gamma sum, `(9)` lies in THM-3099's exact-identity alternative: its two
terms are not merely asymptotic equivalents but the identical proper
hypergeometric sequence with coefficients `+1,-1`.  Rational-shift
quasianalyticity correctly preserves, rather than removes, this dependence.

## 3. One-normal cyclic and degree compatibility

The equality wall is compatible with the quotient degrees in the
one-normal zero-child model.  Work over a characteristic-zero field and fix
a unit `u`.  In the upper algebra

```text
B=K[v]/(v^p-u)                                        (10)
```

the coordinate `v` is a unit and

```text
v^(2m)=u v^(m-1),                         p=m+1.       (11)
```

Adjoin a lower nilpotent coordinate by

```text
A_C=B[e]/(e^2-tau_C^2 u v^(m-1)).                     (12)
```

At `tau_C=0`, the element `e` is nonzero nilpotent of Jordan length two;
in particular the scheme-zero condition `eA_0=0` fails.  The upper constant
`u`, the analogue of the next-moment value, is still a unit.  Equations
`(11)--(12)` give the exact factorization

```text
(e-tau_C v^m)(e+tau_C v^m)=0.                         (13)
```

Equivalently, on the rank-two `B`-module with basis `1,e`, multiplication by
`e+tau_Cv^m` has determinant

```text
tau_C^2 v^(2m)-tau_C^2 u v^(m-1)=0.                  (14)
```

Taking the norm from the rank-`p` upper algebra turns each of the two terms
in `(14)` into the same unit multiple of

```text
tau_C^(2p)=T_C^2.                                     (15)
```

Thus the matrix hostile `(9)` is the scalar shadow of a one-normal cyclic
quotient, not an arbitrary all-orders-flat analytic function.

There is also an unscaled carrier invoice.  With the physical upper scaling

```text
lambda_C=U_p(C)^(-1/p),
E_C=U_m(C)^2/U_p(C),                                  (16)
```

one has

```text
U_m lambda_C^m=tau_C,
E_C lambda_C^(m-1)=tau_C^2.                           (17)
```

The second equality uses `p+m-1=2m`.  The coefficient `E_C` is itself a
positive Gamma ratio with rational consecutive quotient, while its normal
monomial has degree only `m-1`; the missing power is supplied lawfully by
the upper relation `(11)`.  A bare degree count therefore does not exclude
the collision.

## 4. Exact consequence and surviving physical question

The following implication is false, even in the degree-compatible
one-normal local model:

```text
proper positive-ratio Gamma coefficients
+ nilpotent special fibre
+ known unshifted Smith bars
+ unit next repair
+ entrywise error o(repair)

=> the repaired determinant is not an exact identity.                 (18)
```

The first failed step is replacing shifted characteristic coefficients by
their valuations.  At a Jordan block of length `lambda`, one error edge can
replace `lambda` repair edges.  For `lambda=2`, `(7)--(9)` show the sharp
wall `s=z^2`; `s=o(z)` is irrelevant there.

This does not prove that an actual factorial moment deformation realizes
the coefficient and sign in `(12)`.  That realization is precisely the
remaining positive target.  On every physical zero-primary component one
must compute a nonzero endpoint minor or characteristic-coefficient jet and
show that its exact Gamma class does not equal the cancellation wall.  A
single exact good append would alternatively select THM-3099's nonidentity
branch, but no such sample is supplied here.

The scheme-zero hypothesis of THM-3101 removes the hostile at its source:
the old multiplication operator is zero on the selected primary algebra,
so no fixed Jordan edge remains to trade one error for two repairs.  The
present theorem neither weakens that hypothesis nor gives a bad factorial
support, an SFC/GMC counterexample, a new good support, NC2, LRC(14), JC(2),
or DC(2).

## 5. Exact companion

The dependency-free companion verifies the factorial cancellation `(3)`,
the rational quotient `(4)`, exact carrier values over a finite bank, the
matrix identity and its two sign walls, derives Smith valuations `(0,2)`
over the exact polynomial ring, checks the shifted determinant is identically
zero there, and verifies the cyclic
identity `(11)--(15)`, and the unscaled normalization `(16)--(17)` for
`3<=m<=10`, `0<=n<=3`, and `1<=C<=5`.

Run

```text
python 04-computation/gmc_one_normal_gamma_jordan_equality_wall_thm3108.py
python -O 04-computation/gmc_one_normal_gamma_jordan_equality_wall_thm3108.py
```

Both executions must byte-match the stored transcript after LF
normalization.

An independent hostile audit rederived the factorial cancellation and
Stirling exponent, the shifted determinant and Smith bars, the cyclic
factorization and both individual norm terms, and the unscaled degree
invoice.  It also confirmed that the statement never applies norm
additively and that its formal-versus-physical boundary is explicit.

**QED.**
