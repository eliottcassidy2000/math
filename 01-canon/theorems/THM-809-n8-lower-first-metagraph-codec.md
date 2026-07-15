---
id: THM-809
title: The lower-first B3/S2 metagraph address is injective on all n=8 complement lines
status: FINITE-EXACT (all 1,048,576 n=8 complement lines)
source: codex-2026-07-15-S13
depends_on: [THM-553, THM-796, THM-801]
related: [THM-785, THM-790, HYP-6880]
verification:
  - 04-computation/mobius_cech_n8_frontier_codex_S13.cpp
  - 05-knowledge/results/mobius_cech_n8_frontier_codex_S13.out
---

# THM-809 — lower-face recursion already decides `n=8`

Reserved with the exact completed census.  The full statement and preservation
proof are being written from the checkpointed result:

```text
B3 lower node-pair/colour collision excess       418
+ tau=3                                          252
+ tau=4                                          148
+ tau=5                                           74
+ tau=6                                           52
+ tau=7                                            0
+ fixed tau=8                                      0
```

Thus the lower-only address `Lambda` is injective on all `2^20` lines.  Since
`Omega+B2` adds the ordered upper node pair, its `n=8` injectivity follows.
The theorem stub makes no all-size or continuation-completeness claim.
