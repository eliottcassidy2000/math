---
id: THM-774
title: Folded-diamond obstruction for the two-sheet equality packet
status: CLAIMED (number reserved; exact proof and sharpness audit in progress)
source: codex-2026-07-14-S3
depends_on:
  - THM-769
related:
  - THM-772
  - HYP-6820
---

# THM-774 — Folded-diamond obstruction

This stub reserves the theorem number for a metric refinement of THM-769's
two-exception equality packet.  Write

```text
A = 2U union {x,y},       |U|=10,       x>y positive odd,
a = (x+y)/2,              b = (x-y)/2.
```

The intended exact statement is that the two eligibility and opposite-colour
conditions of THM-769 are together equivalent to

```text
||a tau|| + ||b tau|| >= 11/13.
```

Thus a tight packet forces the entire strict loose set `G_U` into a closed
folded diamond.  The diamond has Lebesgue measure at most `8/117`, uniformly
over distinct odd `x,y`, with equality for reduced ratio `x:y=9:1` (and its
common odd dilates).  Consequently every tight two-sheet packet must satisfy

```text
measure(G_U) <= 8/117.
```

The proof carrier is the pair of signed errors from the nearest half-integers
to `a tau` and `b tau`: the two odd-runner errors are their sum and difference,
so the maximum error is the `l1` norm.  The sharp measure bound is an exact
half-grid tooth count after reducing by `gcd(a,b)`.  Until the endpoint
conventions, exact measure formula, finite exceptional cases, and independent
replay have been checked, cite this file only as a claimed reduction.
