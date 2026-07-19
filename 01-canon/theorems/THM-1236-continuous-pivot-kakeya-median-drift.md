---
id: THM-1236
title: THE CONTINUOUS-PIVOT KAKEYA MEDIAN LAW — every auxiliary full-lap scale pays drift, with the fourth speed as the exact constrained optimizer
status: RESERVED/PROOF IN PROGRESS. The moving-arc argument of THM-1235 works for every real auxiliary pivot y>=7c/6, not only runner speeds. The geometric proof, constrained L1 optimization, exact referee, and Lean arithmetic consumers are being audited.
source: codex-2026-07-19-S78 continuation
depends_on: [THM-1235]
related: [THM-1198, THM-1219]
---

# THM-1236 — continuous-pivot Kakeya median drift

For a six-comb cover of a complete `c`-slow gap by ordered speeds
`d1<...<d6`, the claimed strengthening is

```text
sum_i |di-y| > y/14                    for every real y>=7c/6.
```

The proposed proof reparametrizes the initial subneedle by the auxiliary lap
`L*=6y/(7c)`, freezes the resulting six moving arcs, and pays relative drift
`|di/y-1|`.  The remaining work is to record every large-radius/open-cover
branch without strengthening a strict endpoint, and to prove that the exact
constrained optimizer of `sum_i|di-y|-y/14` is `max(7c/6,d4)`.

No LRC(14) conclusion is claimed in this reservation.
