---
id: THM-1605
title: "Historical outside-family comparison: the fixed m=2 orbit count and the missing higher-family definition"
status: >
  PARTIALLY VERIFIED / HISTORICAL RECORD; SUPERSEDED AS AN EXPLICIT FAMILY
  INPUT.  The fixed sporadic map's degree, determinant, three-point fibre,
  and the abstract 1+2k involution-orbit mechanism are exact.  No literal
  definition or source artifact for the historically reported E_m is
  preserved in this repository, so its m>=3 members, degree list, fibre
  claims, and identification with any current weighted family are not
  verified here.  Use THM-3517 for the first explicit odd weighted m=3 test.
author: opus-2026-07-20-S415; corrected by codex-2026-08-16
depends_on:
  - THM-1300-jacobian-counterexample-dixmier-A3-explicit
  - THM-1350-equivariant-fixed-locus-JC
related:
  - THM-3517-weighted-odd-family-m3-three-coordinate-quintic-and-sign-blind-jelonek-component
  - MISTAKE-416
script: 04-computation/family_comparison_opus_S415.py
output: 05-knowledge/results/family_comparison_opus_S415.out
---

# THM-1605 -- corrected extent-versus-mechanism record

> **CURRENT CORRECTION (2026-08-16).**  Earlier versions promoted a prose
> description of an outside family `E_m` as though its defining formula had
> been imported and checked.  It had not.  The repository contains neither a
> literal `E_m` definition nor a recoverable primary-source artifact for that
> family.  Only the fixed-map side and the abstract orbit mechanism were
> verified.  The historical higher-family claims below are therefore reports,
> not proved canon.  See MISTAKE-416 and the disjoint explicit family in
> THM-3517.

This file also shares the legacy number `THM-1605` with the unrelated toral
nullcone theorem.  Every citation must include this slug.

## 1. What was actually checked

For the fixed sporadic map of THM-1300, exact computation gives

```text
ordinary map degree = 7,
det JF = -2,
F^(-1)(1,0,0)
 = {(0,0,1),(i/2,3i,-26),(-i/2,-3i,-26)}.             (1)
```

More generally, over the target `(a,0,0)`, the selected fibre coordinate
obeys

```text
16a x^3+4x=4x(4ax^2+1),                              (2)
```

so the roots are

```text
x=0,              x=+-i/(2 sqrt(a)).                 (3)
```

For the involution `sigma(x,y,z)=(-x,-y,z)`, (3) is one fixed point plus one
free two-orbit.  This proves the fixed numerical split

```text
3=1+2.                                                (4)
```

The companion script checks only these fixed-map facts and the arithmetic
identity `2m-1=1+2(m-1)`.  It does not instantiate or verify an outside map
for any `m>=3`.

## 2. The proved mechanism that survives

THM-1350 is independent of any proposed family formula.  If a Keller map is
`sigma/tau`-equivariant and `dim Fix(sigma)<=1`, then its restriction from
`Fix(sigma)` to `Fix(tau)` is a constant-Jacobian map in dimension at most
one.  `JC(1)` makes that restriction injective.  Therefore a fibre over a
`tau`-fixed target has at most one `sigma`-fixed point, and all remaining
points occur in free two-orbits.  Whenever exactly one fixed point is present,

```text
|fibre|=1+2k.                                         (5)
```

This is a necessity theorem under its stated equivariance and fixed-locus
hypotheses.  It does not construct a map with every odd fibre cardinality,
and those hypotheses must be checked separately for any candidate family.

## 3. What remains only a historical report

The old session prose attributed the following data to an outside family:

```text
E_m:C^3->C^3,       m>=2,
generic/fixed fibre count 2m-1,
reported ordinary degrees 7,13,26,43,64,... .         (6)
```

It also stated that its `m=2` member was the fixed sporadic map.  In the
current repository, none of the following is available:

- a coordinate formula for `E_m`;
- a source citation or retained source file containing that formula;
- an exact `m=3` determinant, eliminant, or fibre replay; or
- a conjugacy exhibiting literal equality with the fixed map.

Accordingly, (6), the all-odd existence claim, and the map identification are
not proved dependencies.  The numerical coincidence at `m=2` is evidence of
the historical comparison, not a reconstruction of the missing object.

## 4. The explicit replacement test is a distinct family

THM-3438 and THM-3448 contain a literal weighted construction.  THM-3517
reindexes its cyclic subfamily by `ell=2m-3` and calls it `E_m^cyc`.  It has

```text
generic degree 2m-1,
ordinary coordinate degrees (10m-13,10m-14,4),
global monodromy S_(2m-1).                             (7)
```

At `m=3`, THM-3517 computes all three coordinate quintics and the exact
Jelonek set.  But (7) begins `7,17,27,...` in its first coordinate, which
does not match the degree string in (6).  This disagreement is a decisive
typing barrier: `E_m^cyc` is a lawful explicit test family, not a recovered
formula for the historical `E_m`.

## 5. Safe citation rule

This file may be cited for exactly two things:

1. the fixed sporadic calculation (1)--(4); and
2. the conditional involution-orbit grammar (5), through THM-1350.

It may not be cited for an explicit infinite family, a verified `m>=3`
member, a degree classification, a map classification, or a `JC(2)` claim.
Use THM-3517 for the explicit weighted odd family and retain its stated
monodromy and effectivity boundaries.
