---
id: HYP-2943
title: LRC14 regular-solid and Euclidean tiling recursion carrier
status: STUB / namespace claimed for S141 exploration
source: codex-2026-06-24-S141
related:
  - HYP-2942
  - HYP-2941
  - HYP-2940
  - HYP-2938
  - HYP-2894
  - HYP-2892
  - HYP-2908
  - THM-572
results:
  - 04-computation/lrc14_platonic_tiling_recursion_codex_s141.py
  - 05-knowledge/results/lrc14_platonic_tiling_recursion_codex_s141.out
---

# HYP-2943: LRC14 regular-solid and Euclidean tiling recursion carrier

S141 is reserved for the POKE/LRC exploration prompted by Platonic,
Archimedean, and Johnson solids plus the Euclidean square/triangular/hexagonal
tilings.

The working question is whether the prompt's recursion signals

```text
square self-recursion:      4, 9, 16, 25, ...
triangle self-recursion:    4, 16, ...
triangle/hex dual exchange: 6
hex self-recursion:         7, 49, ...
```

can be made into a labelled proof carrier for LRC14, without collapsing
branch-local incidence into scalar numerology.  The intended computation will
compare local vertex configurations, curvature/defect, duality, self-similar
patch counts, and Tournament Analysis fingerprints against the existing
HYP-2942 q=3 unital branch-local rule.
