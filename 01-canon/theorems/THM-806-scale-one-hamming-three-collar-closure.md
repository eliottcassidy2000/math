---
id: THM-806
title: Scale-one Hamming-three collar closure
status: CLAIMED (proof and corrected full finite replay in progress)
source: codex-2026-07-15-S10 (continuation of THM-804)
depends_on:
  - LRCUpTo13  # only the 9-, 10-, and 11-speed bounds
related: [THM-770, THM-795, THM-800, THM-804, HYP-6775, HYP-6800, HYP-6820]
---

# THM-806 — Scale-one Hamming-three collar closure

This identifier reserves the scale-one chart isolated by THM-804.  The target
statement is that every proper labelled triple lift

```text
([12]\{r,s,t}) union {r+13i,s+13j,t+13k},   i,j,k>=1,
```

has a strict `1/13` witness.  The proof route now established on scratch work
has two parts: an oriented owner-collar handoff lemma forces one replacement
to be at most `24`; settled lower-runner bounds then reduce the other two,
after ordering, to `v<=381` and `w<=12v`.  An exact component-containment
replay is being completed over that finite box.  This stub claims only the
identifier and records the remaining audit debt; it is not yet a canonical
proof.
