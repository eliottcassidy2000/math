---
id: THM-3798
title: "Cubic pseudo-plane common-AP step-three support peel"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  On every
  THM-3785 cubic pseudo-plane over an algebraically closed characteristic-zero
  field, no exact 2x4 or 4x2 Euler-support pair whose two supports are
  arithmetic progressions of common step three has nonzero scalar bracket.
  The four exhaustive scalar-placement rows all retain the squarefree arm
  divisor; the (-1,-7) row in fact retains its cube.  Common step at least
  four, disjoint-pair and three-chain grammars, arbitrary 2x4 support,
  arbitrary Darboux pairs, and JC(2) remain open.
source: root / post-THM-3796 endpoint-power lane, 2026-08-23
audit: >
  PASS.  An independent hostile audit rederived the complete sign table,
  endpoint UFD law with field units, squarefree divisor boundary, both
  adjacent ODE systems, support-shrink seams, output swap, and the seven
  step-four hostiles.  It added the upper-ODE Delta^3 refinement and ran
  6,105 active gates, including 1,940 UFD valuation cases, with no Python
  assert gate.  The original fraction-field divisibility checks were
  vacuous; MISTAKE-448 records their replacement by explicit polynomial
  quotients and independent Euclidean-remainder paths.
depends_on:
  - THM-3785-linear-higher-pole-russell-pseudoplane-maximal-observable
  - THM-3796-cubic-pseudoplane-first-two-by-four-support-peel
related:
  - THM-3787-cubic-pseudoplane-complete-low-support-darboux-nonentry
script: 04-computation/jc2_cubic_pseudoplane_common_ap_d3_thm3798.py
output: 05-knowledge/results/jc2_cubic_pseudoplane_common_ap_d3_thm3798.out
script_sha256: 6b6a1867ec9311bf108dea2147c39504e66fec78b92c62f0f302043c2db8dad4
output_sha256: 9c65319ccf272948767aa4083435299ad072ed660c1b7b62d592dcf464f0c522
semantic_sha256: 597dc0e9819318b6d6b1ce6beec3a9821f7d0be6cd42840e43bb7acb47615cbc
hash_basis: raw LF bytes
---

# THM-3798 -- common step three is empty

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  Let `k` be
algebraically closed of characteristic zero, let `c in k*`, and retain the
THM-3785 cubic pseudo-plane

```text
B=k[R,Z,E]/(R^2E-Z^3+c^3R)                           (1)
```

with Euler weights `(3,1,-3)` and Poisson-bracket shift `+2`.  Suppose
`F,G in B` have exact active supports

```text
supp(F)={a,a+3},
supp(G)={b,b+3,b+6,b+9}.                              (2)
```

An active weight-zero scalar may be translated away and is not counted,
whereas a nonconstant weight-zero profile remains active.  Then

```text
{F,G} is not a nonzero scalar.                        (3)
```

By skew-symmetry, the same holds for exact `4 x 2` support after swapping
the outputs.

## 1. Four exhaustive scalar rows

Use THM-3785's homogeneous-profile notation

```text
A_u=x^u f(w),                 w=c+x^3y,
[u,f;v,g]=u f g'-v f'g,
Delta=w^3-c^3.                                      (4)
```

Every negative-weight profile is divisible by `Delta`.  Moreover

```text
gcd(Delta,Delta')=1,             disc_w(Delta)=-27c^6 !=0,            (5)
```

so `Delta` is squarefree.

The contribution addresses of (2) are

```text
0,3,6,9,12                  with multiplicities 1,2,2,2,1.           (6)
```

A scalar cannot occupy an endpoint, since that would be a homogeneous
Darboux pair forbidden by THM-3785.  If its relative address is
`q in {3,6,9}`, the scalar weight equation is

```text
a+b=-q-2.                                             (7)
```

The bottom singleton commutes.  For an active commuting singleton, the two
weights have equal strict signs or are both zero: opposite signs would make
a positive product of profiles constant despite the negative profile's
`Delta` factor, and exactly one zero weight would force the zero-weight
profile to be a removable scalar.  Equation (7) therefore forces
`a,b<0`.

The top singleton also commutes.  Its weights are `a+3,b+9`, with sum
`10-q`.  Applying the same sign/zero law on the exhaustive interval
`-q-1<=a<=-1` leaves exactly

```text
q=3:  (a,b)=(-2,-3), (-1,-4),
q=6:  (a,b)=(-2,-6), (-1,-7),
q=9:  no row.                                         (8)
```

There is no bounded-search assumption in (8).

## 2. Endpoint power law

Suppose `u,v` are nonzero integers of the same sign and
`[u,f;v,g]=0`.  Put `U=|u|`, `V=|v|`, and
`delta_0=gcd(U,V)`.  In `k(w)`,

```text
(g^U/f^V)'=0,                    g^U=lambda f^V.       (9)
```

At each irreducible polynomial, unique factorization and coprimality of
`U/delta_0,V/delta_0` give, after absorbing field units using algebraic
closedness,

```text
f=t^(U/delta_0),              g=lambda_0 t^(V/delta_0).              (10)
```

No unit or root choice is lost.  If a negative profile is a positive power
of `t`, then `Delta` divides that power; squarefreeness in (5) implies

```text
Delta | t.                                             (11)
```

Write the two profiles of `F` as `p,q`, bottom to top, and those of
`G` as `g_0,g_1,g_2,g_3`.

## 3. Scalar address three

### Row `(-2,-3)`

Bottom commutation and (10) give

```text
p=t^2,                       g_0=lambda t^3.           (12)
```

Since `p` is negative, (11) gives `Delta|t`.  The scalar bucket factors
before any division:

```text
S=[-2,t^2;0,g_1]+[1,q;-3,lambda t^3]
 =t^2[-2g_1'+3lambda q t'+3lambda q't].               (13)
```

Hence `Delta^2|S), so `S` cannot be a nonzero scalar.

### Row `(-1,-4)`

Now `g_0=lambda p^4`.  Since both `p` and `g_1` have weight `-1`,
write `p=Delta P`, `g_1=Delta N`.  Derivatives of `Delta` cancel:

```text
[-1,p;-1,g_1]=Delta^2(P'N-PN').                      (14)
```

The other scalar summand is

```text
[2,q;-4,lambda p^4]
 =8lambda q p^3p'+4lambda q'p^4,                     (15)
```

which retains at least `Delta^3).  Thus `Delta^2|S`.

## 4. Scalar address six

### Row `(-2,-6)`

The two endpoint equations give

```text
g_0=lambda p^3,                  g_3=mu q^3.           (16)
```

The lower and upper adjacent nonscalar collisions have the complete
polynomial solutions

```text
g_1=3lambda p^2q+H,             [-2,p;-3,H]=0,
g_2=3mu p q^2+nu,               nu in k.              (17)
```

Substitution into the scalar bucket gives

```text
S=6(lambda-mu)p(pq^2)'+[1,q;-3,H].                   (18)
```

If `H=0`, (18) retains the factor `p), hence `Delta`; if
`lambda=mu` it is zero.  If `H!=0`, the UFD law gives

```text
p=t^2,                         H=kappa t^3,            (19)
```

so `Delta|t`, and every term in (18) is divisible by `Delta^2`.
The arbitrary integration constant `nu`, including `nu=0`, creates no
seam.

### Row `(-1,-7)`

Bottom commutation gives `g_0=lambda p^7`.  The lower adjacent equation has
the full polynomial solution

```text
g_1=p^4(7lambda p^2q+nu),          nu in k.            (20)
```

For the scalar bucket,

```text
S=[-1,p;-1,g_2]+[2,q;-4,g_1].                         (21)
```

The first summand has the `Delta^2) Wronskian factor (14), and the second
retains at least `p^3`, hence `Delta^3).  This already proves
`Delta^2|S`.

The independent audit also closes the unused upper seam.  Top commutation
gives `g_3=mu q`.  With `L=g_2-mu p`, the upper adjacent collision is

```text
2qL'+q'L=0,                       (qL^2)'=0.           (22)
```

The exact positive weight-two profile is `q=w^2Q(w^3)`, hence is
nonconstant.  If `qL^2` were a nonzero field constant, `q` would be a
unit; if it is zero, the domain property gives `L=0).  Therefore
`g_2=mu p`.  The first summand in (21) vanishes, sharpening this whole row
to

```text
Delta^3 | S.                                          (23)
```

In all four rows, `S` belongs to `B_0=k[w^3]`.  The nonunit
`Delta=w^3-c^3` divides `S), so `S` is zero or nonconstant, never a
nonzero scalar.  This proves (3).  The identity

```text
[v,g;u,f]=-[u,f;v,g]                                  (24)
```

proves the swapped `4 x 2` statement without a new support seam.

## 5. Sharp failure boundary

At common step four, the same exhaustive sign/zero gate leaves seven rows:

```text
q=4:  (-3,-3), (-2,-4), (-1,-5),
q=8:  (-3,-7), (-2,-8), (-1,-9),
q=12: (-3,-11).                                      (25)
```

These are necessary sign survivors only.  No coefficient system is
constructed or excluded here.  The first disjoint-pair control
`r=4,V={0,1,4,5}` from THM-3796 also remains open.

## 6. Verification and correction lineage

The canonical companion checks the exhaustive tables, symbolic identities,
squarefreeness, endpoint powers, zero integration constants, and explicit
denominator-free divisibility quotients in `24` active gates.  It contains
no Python `assert` gate.  Normal and optimized executions LF-normalize to
the frozen transcript.

The independent audit at

```text
.scratch/thm3798_hostile_audit_20260823/
```

uses both formal first-jet polynomial algebra and actual representatives in
the exact weight modules.  It runs `6,105` active gates, including
`1,940` UFD valuation cases, and checks the stronger upper seam (22).

An earlier companion divided by `Delta^m` in the fraction field and then
multiplied back, a tautology that did not certify polynomial divisibility.
That artifact is **SUPERSEDED**.  MISTAKE-448 records the repair to explicit
polynomial quotients, denominator-one gates, and independent remainder
checks.

Run

```bash
python3 -B 04-computation/jc2_cubic_pseudoplane_common_ap_d3_thm3798.py
python3 -B -O 04-computation/jc2_cubic_pseudoplane_common_ap_d3_thm3798.py
```

Both streams equal
`05-knowledge/results/jc2_cubic_pseudoplane_common_ap_d3_thm3798.out`.

This theorem closes only the exact common-AP step-three cell.  It proves no
common-step-four result, arbitrary `2 x 4` peel, Darboux pair, Keller map,
or Jacobian counterexample.  `JC(2)` remains **OPEN**.
