---
id: THM-2979
title: "Projected k3 z275 ten-body status and located-torsion closure"
status: >
  PROVED + VERIFIED-EXACT.  The ten projected k=3 atlas rows at z1=275
  contain exactly 2,033 attained denominator states.  Crude fibre capacity
  removes 273 and the common 16-status table removes 1,695 with independently
  replayed exact rational Farkas checks.  The remaining 65 states force
  exactly one high label and are all closed by a uniform order-seven
  clean-cell residue collision.  Hence the projected k=3 cap is z1<=274 and
  the necessary-row ledger is 375,703.  This is not LRC(14).
source: codex-lrc14-k3-z275-ten-body-closure-2026-07-30
audit: >
  Normal and optimized executions are byte-identical.  The referee freezes
  the complete ten-row atlas order, all stage counts, 1,695 exact dual
  rechecks, positive feasible and incompatible hostile status controls, five
  positive two-high gaps, 65 zero-high hostile passes, all literal negative
  low amplitudes, 2,585,952 high-unit recurrence checks across the combined
  package, and a concrete exact order-seven cell pair in every terminal.
depends_on:
  - THM-2941-critical-seven-slot-scalar-wall-and-balanced-boundary
related:
  - MISTAKE-333
script: 04-computation/lrc14_j7_k3_z275_to_z272_septimal_torsion_descent_thm2941.py
output: 05-knowledge/results/lrc14_j7_k3_z275_to_z272_septimal_torsion_descent_thm2941.out
script_sha256: 44d9369fe3151d78f5a281cb8ecc26d859807b635980552d8e5c95c0b5ff826e
output_sha256: d31d9ec1f45a3aa5d087a2b8724fc75480b0d0b4d2b99e544832b0947d408b67
profile_sha256: 1139218bdb4e9d19a3a64b22103030ef682e835803f5caf07c1659b0ce788f5b
semantic_sha256: 0960c65f29707e6ce1b295f35f0118960ce420c1d27df58365b91999c1faca9b
hash_basis: LF-normalized bytes
---

# THM-2979 -- projected k3 z275 ten-body status and located-torsion closure

**PROVED + VERIFIED-EXACT.**

## Statement

In the lossless projected `k=3` scalar atlas inherited from THM-2941, all ten
body rows with first drift `z_1=275` are empty.  Consequently

```text
z_1 <= 274,                                                (1)
```

and removing these ten rows lowers the necessary-row ledger from `375,713`
to `375,703`.  This is a scoped projected-sector theorem, not LRC(14).

## Exact finite reduction

The ten rows contain exactly

```text
2,033 = 273 crude + 1,695 common-status + 65 residual       (2)
```

attained denominator states.  The crude exclusions are exact fibre-capacity
inequalities.  Each common-status exclusion rebuilds the same finite
16-pattern feasibility problem and verifies an exact rational Farkas dual.
Positive feasible and incompatible hostile controls remain active.

The `65` residuals occur on precisely five bodies.  Their inherited
`HIGH-TAIL` gate is essential: the zero-high hostile relaxation passes all
`65`.  On each body, a duplicate-permitting upper bound for two or more high
labels misses the scalar wall by a positive exact gap.  Hence every possible
packet has exactly one high label.  Literal enumeration retains all low
labels, including negative-amplitude ones, and leaves exactly one low pair on
each survivor body:

```text
E=(1,5,7,9,11,13)       low pair (287,351), C=115050, L=630630
E=(1,5,9,11,13,14)      low pair (351,378), C=240580, L=1261260
E=(1,8,10,11,12,14)     low pair (312,364), C=24978,  L=129360
E=(2,5,7,9,11,13)       low pair (287,351), C=232068, L=1261260
E=(2,5,9,11,13,14)      low pair (351,378), C=244176, L=1261260.
```

Here `C` is the number of complete body cells missed by the first and two low
labels.  In every row `7C>L`.

## Uniform septimal closure

Fix one residual case and let `d` be its high denominator.  Exact enumeration
gives `7|d`.  Let `S` be the set of residues modulo `d` represented by those
`C` clean cells.  A residue fibre contains at most `L/d` cells, whence

```text
C <= |S| L/d  and  7C>L  imply  |S|>d/7.                  (3)
```

There are only `d/7` cosets modulo `d/7`; therefore two distinct elements of
`S` lie in one coset.  Their nonzero difference has exact order seven modulo
`d` (seven is prime).  Every high label has the unit-ray form

```text
z=(L/d)u+hL,  gcd(u,d)=1.                                 (4)
```

Multiplication by `u` preserves order seven and the height term cancels.
The two high phases thus have circular separation at least `1/7`.  Two
strict-open danger arcs of radius `1/14` cannot both contain them; at equality
their excluded endpoints merely touch.  At least one of the two complete
cells stays safe for every primitive direction and every ray height.  All
`65` residual cases retain full projected drift-safe mass `1`, strictly above
the three-aligned union cap `36/91`, and are empty.

## Reproducibility and scope boundary

The referee checks all ten `z_1=275` rows and also continues through the next
29 occupied rows.  THM-2979 uses only the first level and therefore proves
the cap `(1)`, not the later cap.  The `z_1=274,273,272` continuation is
recorded separately as a THM-2941 addendum.

The MISTAKE-333 repair boundary is explicit.  Every solver-selected Farkas
basis is checked exactly but excluded from stored row and semantic digests;
only certificate-free problem data and the verified-instance count are
hashed.  Normal and optimized Python executions are byte-identical.
