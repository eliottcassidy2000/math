---
id: THM-980
title: Scale-sixteen Hamming-six owner-feasibility deficit
status: CLAIMED — a literal C++ reconstruction currently finds that the complete 13,806,600,192-context primitive proper AP-centred common-scale-sixteen Hamming-six bank reduces to 2,540 scalar rows and that no row is owner-locally feasible at more than two of its six owners; independent replay and frozen-build validation are in progress
source: codex-2026-07-17-S66 scale-sixteen exact C++ certificate in progress
depends_on: [THM-765, THM-810, THM-823, THM-859, THM-860]
related: [THM-976, THM-977, THM-978, HYP-6820]
verification:
  - 04-computation/lrc13_scale_sixteen_hamming_six_frontier_scout_codex_c16.cpp
---

# THM-980 — scale sixteen has at least four impossible owners

This namespace is reserved for the exact scale-sixteen continuation of the
primitive proper AP-centred common-scale Hamming-six classification.

For `c=16`, the effective orders are `1,2,4,8,16`, with sixteen literal
`(D,e)` states.  Leave-one-out lcm enumeration currently gives 5,385
hereditary order words and 14,942,208 literal state words per support, hence

```text
924*14,942,208 = 13,806,600,192
```

labelled raw contexts.  Unit-independent scalar owner capacity leaves 2,540
contexts on 404 supports.  Exact owner-local union-mask reachability currently
gives

```text
0 feasible owners: 2,006 contexts,
1 feasible owner :   432 contexts,
2 feasible owners:   102 contexts.
```

Thus the primary implementation finds no context feasible at all six owners;
indeed every row has at least four impossible owners.  Across all 15,240 owner
rows its maximum reachable-union histogram is

```text
 9 sheets:   144,  10 sheets:   468,  11 sheets:   876,
12 sheets: 2,316,  13 sheets: 4,068,  14 sheets: 4,740,
15 sheets: 1,992,  16 sheets:   636.
```

The faithful prospective certificate is the ordered six-owner vector of local
maximum unions (or its thresholded feasibility subset), before any shared-unit
global replay.  This is a **pre-nerve** obstruction: because at least one owner
obligation set is already empty, the nerve of globally covariant unit words
never needs to be formed.  A tournament that ranks owners by
`(feasible, maximum union)` is useful telemetry but loses the absolute
sixteen-sheet threshold unless the full pair observable is retained.

The current C++ program reconstructs CRT masks literally, traverses the full
hereditary state grammar, and checks all `924*5,385=4,975,740` labelled divisor
contexts.  Promotion to `PROVED FINITE-EXACT` requires an independently written
referee, frozen hashes, and optimized/unoptimized/sanitized replay.  Until
those checks land, every census above is explicitly provisional.

This claim concerns only the primitive proper AP-centred common-scale-sixteen
Hamming-six face under THM-860.  It does not close `c>=17`, the H5 bank,
non-AP-centred/deep-sheet continuations, or global sporadic emptiness.
