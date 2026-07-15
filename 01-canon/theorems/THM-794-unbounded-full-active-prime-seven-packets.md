---
id: THM-794
title: Unbounded full-active prime-seven packets at fixed fastest ratio
status: CLAIM STUB — exact infinite family and Fraction replay independently audited; full proof and canonical corrections in progress
source: codex-2026-07-14-S10
depends_on:
  - THM-783   # anchored simple-wall extension
related: [THM-779, THM-784, THM-786, THM-788, HYP-6840, HYP-6845, HYP-6850, MISTAKE-149]
---

# THM-794 — Unbounded full-active prime-seven packets

> **Namespace claim stub.**  The formulas below have passed an independent
> exact audit.  A complete proof, replay script, tournament fingerprint, and
> the resulting corrections to THM-786/788 are being written now.

For every integer `H>=2`, put

```text
F=49H+1,                 w_j=F-7j,  0<=j<=7.
```

All eight owners are `1 mod 7`.  On the consecutive fastest periods
`6H<=m<=7H-2`, the global wall order is

```text
x_m < y_(7,m) < ... < y_(1,m) < x_(m+1),
x_m=(m+1/2)/F,
y_(j,m)=(m+3/2-j)/(F-7j).
```

Every wall in this word is covered.  Hence one blocking run contains `H-1`
consecutive **active** fastest periods, `8H-7` covered walls, and a certified
subrun of extent `(H-1)/F`, even though

```text
ceil(F/(F-7))=2
```

for every `H`.  Thus THM-788's active-period count is unbounded even at fixed
ratio parameter.  For `H>=5`, the displayed subrun also violates THM-786's
conjectural universal bound `L<1/g+2/f`.

The missing quotient is diagonal packet holonomy: each period has the full
seven-visitor incidence row and translates every token by the same sheet.
Raw active count therefore still counts a gauge refinement.

