---
id: THM-832
title: The n=9 false-palindrome defect is a rank-four Cech core on every nonempty face intersection
status: RESERVED + PROVED LOCAL LINEAR SKELETON; exact nonlinear fibre audit in progress
source: codex-2026-07-15-S13
depends_on: [THM-801, THM-828]
related: [THM-553, THM-810, HYP-3234, HYP-6880]
planned_verification:
  - 04-computation/n9_rank_four_cech_core_codex_S13.py
  - 05-knowledge/results/n9_rank_four_cech_core_codex_S13.out
---

# THM-832 — reservation: the defect is a redundantly stored Cech core

This number reserves the exact continuation of THM-828 through the old
`A+B+C-D-E-F+G` recursion.  Restrict its four displayed defect generators to
the three size-eight faces.  Every face, every pairwise face intersection,
and the ten-cell triple intersection has binary rank four.  Each of the three
face-exclusive corner strata has rank zero.  Thus the rank valuation obeys

```text
4+4+4-4-4-4+4=4.
```

The reserved theorem will make the restriction matrices, overlap recovery
maps, survivor/puncture fibre genealogy, and preservation boundary exact.
The statement is about the linear reflection-defect carrier; it does not
claim that the nonlinear `H_8` kernel or raw-S2 predicate is determined by a
bare rank valuation.
