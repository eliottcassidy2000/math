---
id: THM-833
title: Centered-CF transport preserves the n=9 rank-four defect but exposes two dimensions to n=10 raw S2
status: RESERVED + FINITE-EXACT COMPUTATION COMPLETE; canonical proof integration in progress
source: codex-2026-07-15-S13/engine
depends_on: [THM-812, THM-813, THM-828, THM-832]
related: [THM-825, THM-829, HYP-6880]
verification:
  - 04-computation/continued_fraction_n9_defect_transport_codex_S13.py
  - 05-knowledge/results/continued_fraction_n9_defect_transport_codex_S13.out
  - 05-knowledge/results/continued_fraction_n9_defect_transport_codex_S13.json
---

# THM-833 — reservation: centered-CF transport of the defect core

This number reserves the exact recursive coordinate-copy audit.  The derived
map `X_9->X_10` is injective and reflection-equivariant.  It sends THM-828's
four-dimensional defect space injectively to rank four, but its intersection
with the target raw-S2 syndrome has rank two.  Of the 58 literal source pairs,
only ten remain raw-S2-equal after transport; positional moments separate
those ten.  All 58 literal reflection orbits and projected coloured cells
remain distinct in this restricted bank.

The completed theorem will state the derived copy word, layer genealogy,
sector multiplicities, node audit, and preservation boundary.  It will not
claim that an arbitrary continued-fraction word preserves the codec or the
LRC metric predicate.
