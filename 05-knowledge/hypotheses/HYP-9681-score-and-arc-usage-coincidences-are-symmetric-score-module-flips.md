---
id: HYP-9681
title: "Every coincidence of the joint observer (score-sequence pair, arc-usage multiset) between distinct merged tournament classes is a converse flip of a strong module whose internal score sequence is self-reversed and whose seam weights the flip preserves"
status: >
  OPEN (CONJECTURE).
  Evidence:
  - THM-5155: every coincidence through order 9 has this form; the only ones are 6 run-flanked core flips at order
    9.
  - Session 320's strong order-24 pair, which reverses a vertex-transitive, hence regular, module, is of this form
    as far as checked (scores trivially self-reversed).
  - The order-12 interior block reversal of the S29 audit (C3 => v => C3 => S4 => w against
    C3 => v => S4 => C3 => w) changes the score sequence, so it is not a coincidence of this observer.
session: kind-pasteur-2026-11-13-S30
source: 05-knowledge/results/kps_s30_observer_lattice_20261113/README.md
related:
  - THM-5155
  - THM-5154
  - HYP-9680
---

# HYP-9681 - score and arc-usage coincidences are symmetric-score module flips

Neither the arc-usage multiset nor the score sequence alone determines a tournament class. Through order 9, the
only way to fool both at once is to reverse a module whose own score sequence reads the same backwards. The
conjecture is that this holds at every order, with "module" in place of "run-flanked core", once the seam weights
are preserved.
