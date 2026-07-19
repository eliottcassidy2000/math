---
id: THM-1267
title: Centered-spoke survivor protrusion law
status: RESERVED IN PROGRESS -- combining the slowest centered spoke with the full five-comb survivor on its complete safe component.  The exact center-offset identity and the strict 89/42 ratio improvement are independently audited; the stronger endpoint-density and Kakeya-envelope refinements remain under audit before theorem status is claimed
source: codex-2026-07-19 H-drift / centered-component synthesis
depends_on: [THM-1198, THM-1240, THM-1241, THM-1252]
related: [THM-1199, THM-1236, THM-1244]
---

# THM-1267 -- centered-spoke survivor protrusion law

Let `d=d_1` be the slowest fast speed in a hypothetical six-comb cover of a
complete `c`-safe gap.  Apply the centered-spoke construction to `d`, and let
`S` be the complete `d`-safe component through that spoke.  On `S` the five
remaining, genuinely faster combs have the THM-1198 survivor; inside
`S intersect G` those five combs must cover.  The survivor is consequently
forced into the one endpoint tail `S minus G`.

Writing `rho` for the nearest-integer error of the centered spoke, the exact
geometry is

```text
center(S)-center(G)=+-rho/d,
|S minus G|=max(0,(rho+3/7)/d-3/(7c)).
```

The coarse survivor-length consequence already gives the independently
audited strict inequality `42d<89c`.  This theorem will retain the stronger
six-bin endpoint-mass profile and the exact one-comb envelope before freezing
its final constant.  It is a functional density/position conversion, not an
additional overlap invoice.

