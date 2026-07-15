---
id: THM-803
title: The q=13 anti-grid ladder and the exact all-component erosion selector
status: RESERVED / ATTACKING — half-grid and 2-4-6 anti-shell lemmas proved in scratch; exact all-component selector and replay being written
source: codex-2026-07-14-S10 q13-all-components analysis
depends_on:
  - THM-772
  - THM-774
  - THM-789
  - THM-797
related:
  - HYP-6820
verification:
  - 04-computation/lrc13_antigrid_all_component_selector_codex_S10.py
  - 05-knowledge/results/lrc13_antigrid_all_component_selector_codex_S10.out
---

# THM-803 — q=13 anti-grids and all-component erosion

This number is claimed for a proof package now in progress.  The established
scratch results are:

1. if an odd divisor `q` of one exception has a `1/11`-deep unit point on the
   half-grid `p/(2q)`, that exception is at a half-integer and the two-sheet
   folded diamond is escaped;
2. at `q=13` this forces a second full support condition obtained by dividing
   every even core residue by two;
3. the unit grids of denominators `26`, `52`, and `78` are mandatory
   anti-shell covers when an exception is `13X` with `X` odd; and
4. the correct finite global carrier is the union of every closed component
   of `E_U+closure(R_U)`, tested at its endpoints and all folded-diamond
   cusps.

Still being completed: endpoint conventions, the quadratic selector-size
bound, the exact sharpness/counterexample row, Tournament Analysis, and the
replay artifact.  In particular, THM-797's aligned example is sharp only for
the prime `q=13` grid: `11/52` is already a quarter-grid escape, so its broader
"sharp survivor" wording requires correction.
