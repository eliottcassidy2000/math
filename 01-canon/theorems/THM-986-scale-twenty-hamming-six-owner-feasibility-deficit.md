---
id: THM-986
title: Scale-twenty Hamming-six owner-feasibility deficit
status: PROVED STRUCTURAL + FINITE-EXACT — independent literal-search C++ and algebraic-CRT Python certificates exhaust the 52,583,731,200-context primitive proper AP-centred common-scale-twenty Hamming-six bank and prove that no scalar row is owner-locally feasible at more than two owners
source: codex-2026-07-17-S66 scale-twenty continuation
depends_on: [THM-765, THM-810, THM-823, THM-859, THM-860, THM-982, THM-983]
related: [THM-978, THM-980, THM-981, HYP-6820]
verification:
  - 04-computation/lrc13_scale_twenty_hamming_six_frontier_scout_codex_c20.cpp
  - 05-knowledge/results/lrc13_scale_twenty_hamming_six_frontier_scout_codex_c20.out
  - 04-computation/lrc13_scale_twenty_hamming_six_referee_codex_c20.py
  - 05-knowledge/results/lrc13_scale_twenty_hamming_six_referee_codex_c20.out
---

# THM-986 — scale twenty has at least four impossible owners

The primitive proper AP-centred common-scale-twenty Hamming-six sheet bank is
empty before the global covariant unit-fibre or metric-continuation stages.

For `c=20`, the effective orders are `1,2,4,5,10,20`, with twenty literal
`(D,e)` states.  Hereditary leave-one-out lcm enumeration gives 26,961 divisor
words and 56,908,800 literal state words per support, hence

```text
924*56,908,800 = 52,583,731,200
```

unquotiented labelled state contexts.  The primary exact certificate finds
12,584 scalar contexts on 830 supports across 65 order-multiplicity patterns.
Owner-local set-union reachability gives the exact census

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

## Independent replay

The primary certificate constructs every CRT representative by literal search,
checks scalar cardinalities against a separate period-count formula, traverses
the full literal unit grammar, and computes owner-local unions with a packed
marker table.  The independent Python referee instead solves the CRT
algebraically, enumerates grammar words with Cartesian products, and carries
each owner bank as an immutable set.  It hashes every reachable mask in sorted
order, not only the reported maximum.  The implementations agree on all
censuses, including 65 scalar multiplicity patterns, the all-support context
histogram (with 94 scalar-empty supports), multiplication orbits
`2:1,6:2,12:68`, the minimum-owner histogram, 4,061 owner maximum vectors, and
zero all-six contexts.

| artifact | SHA-256 |
|---|---|
| `04-computation/lrc13_scale_twenty_hamming_six_frontier_scout_codex_c20.cpp` | `6e1fcbafd13e6a9535aaa039c0ba336697c970195da783b059300e6329fc6c98` |
| `05-knowledge/results/lrc13_scale_twenty_hamming_six_frontier_scout_codex_c20.out` | `dbe98f42aa4bef3a13283126120caf51f964420c6d5b6f13c0a697db3d05fbf3` |
| `04-computation/lrc13_scale_twenty_hamming_six_referee_codex_c20.py` | `69cccaaa0c17232fd79c6c12e3e13d2475d4ba606f5a84c8f8b9417bff01424f` |
| `05-knowledge/results/lrc13_scale_twenty_hamming_six_referee_codex_c20.out` | `50566b452a7d0c19f93989e68fb822d4220f2d8ac2b7974d46ecadddc29087cc` |

Optimized, unoptimized, and ASan/UBSan C++ outputs agree byte-for-byte.  Normal
and optimized Python outputs also agree byte-for-byte.

The faithful carrier is again the labelled owner feasibility/max-union vector.
The completed owner-summary tournament is transitive in every survivor and
forgets both the absolute twenty-sheet threshold and the four-owner deficit.
Provider, divisor, residue, or isolated-sheet vertices destroy the shared-unit
incidence earlier still.

For every context a global unit word would induce a feasible choice in all six
owner projections, so the exact maximum of two is terminal.  The theorem
concerns only the primitive proper AP-centred common-scale-twenty Hamming-six
face.  It does not
close scale 21 or higher, the H5 bank, non-AP/deep-sheet continuations, or
global sporadic emptiness.
