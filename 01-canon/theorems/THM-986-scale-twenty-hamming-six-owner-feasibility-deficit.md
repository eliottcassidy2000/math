---
id: THM-986
title: Scale-twenty Hamming-six owner-feasibility deficit
status: CLAIMED — a frozen exact C++ certificate exhausts the 52,583,731,200-context primitive proper AP-centred common-scale-twenty Hamming-six bank and finds no scalar row feasible at more than two owners; an independently derived algebraic-CRT referee is in progress
source: codex-2026-07-17-S66 scale-twenty continuation
depends_on: [THM-765, THM-810, THM-823, THM-859, THM-860, THM-982, THM-983]
related: [THM-978, THM-980, THM-981, HYP-6820]
---

# THM-986 — scale twenty has at least four impossible owners

This namespace is reserved for the exact scale-twenty continuation of the
primitive proper AP-centred common-scale Hamming-six classification.

For `c=20`, the effective orders are `1,2,4,5,10,20`, with twenty literal
`(D,e)` states.  Hereditary leave-one-out lcm enumeration gives 26,961 divisor
words and 56,908,800 literal state words per support, hence

```text
924*56,908,800 = 52,583,731,200
```

unquotiented labelled state contexts.  The primary exact certificate finds
12,584 scalar contexts on 830 supports across 65 order-multiplicity patterns.
Owner-local set-union reachability gives the provisional exact census

```text
0 feasible owners: 9,632 contexts,
1 feasible owner : 2,184 contexts,
2 feasible owners:   768 contexts.
```

Thus every scalar row has at least four empty owner projections.  The 75,504
owner rows have maximum-union histogram

```text
12:1320, 13:2904, 14:9720, 15:15924, 16:22800,
17:12876, 18:5700, 19:540, 20:3720.
```

The C++ source/output SHA-256 values are respectively
`6e1fcbafd13e6a9535aaa039c0ba336697c970195da783b059300e6329fc6c98`
and `dbe98f42aa4bef3a13283126120caf51f964420c6d5b6f13c0a697db3d05fbf3`;
optimized, unoptimized, and ASan/UBSan outputs agree byte-for-byte.

The faithful carrier is again the labelled owner feasibility/max-union vector.
The completed owner-summary tournament is transitive in every survivor and
forgets both the absolute twenty-sheet threshold and the four-owner deficit.
Provider, divisor, residue, or isolated-sheet vertices destroy the shared-unit
incidence earlier still.

Promotion to `PROVED STRUCTURAL + FINITE-EXACT` requires the in-progress
independent algebraic-CRT Python referee.  This claim concerns only the
primitive proper AP-centred common-scale-twenty Hamming-six face.  It does not
close scale 21 or higher, the H5 bank, non-AP/deep-sheet continuations, or
global sporadic emptiness.
