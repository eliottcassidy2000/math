---
id: HYP-3436
title: LRC14 minimal bad-core cover extractor
status: RESERVED / active exact interval-cover classification; not a proof
source: codex-2026-06-29 continuation of HYP-3435 and the HYP-3422/HYP-3425 two-adic bad-core chain
tangent: T1397
technique: LTI-397
tournament_technique: LTT-297
script: 04-computation/lrc14_minimal_bad_core_cover_extractor_codex_20260629.py
result: 05-knowledge/results/lrc14_minimal_bad_core_cover_extractor_codex_20260629.out
reflection: 07-reflections/lrc14-minimal-bad-core-cover-extractor-codex-20260629.md
related:
  - HYP-3435
  - HYP-3434
  - HYP-3433
  - HYP-3432
  - HYP-3431
  - HYP-3430
  - HYP-3429
  - HYP-3428
  - HYP-3427
  - HYP-3426
  - HYP-3425
  - HYP-3424
  - HYP-3423
  - HYP-3422
  - HYP-3417
  - HYP-3129
  - HYP-2963
  - THM-523
  - OPEN-Q-108
---

# HYP-3436: LRC14 Minimal Bad-Core Cover Extractor

Reserved active packet.  The task is to invert HYP-3435's positive
two-adic branch-cover certificate.

Instead of only reporting a surviving branch witness in

```text
E_safe(1/14) cap (branch0_good union branch1_good),
```

this packet will classify the complementary two-color obstruction

```text
E_safe(1/14) cap B0_odd cap B1_odd.
```

For each `E_safe` component, the intended executable extractor records exact
bad-core components, the minimal branch-0 and branch-1 odd-owner subsets that
cover each bad component, endpoint owners, survivor gaps, and a Tournament
Analysis over proof obligations rather than runners.  The proof-facing hope is
that a full counterexample would have to glue these local minimal covers into a
global bad-core cover, while HYP-3426/HYP-3427/HYP-3429/HYP-3434 already show
that the needed overlap tax and endpoint sidecars are highly constrained.

This stub honestly reserves HYP-3436/T1397/LTI-397/LTT-297 while the exact
script and stored readout are being built.
