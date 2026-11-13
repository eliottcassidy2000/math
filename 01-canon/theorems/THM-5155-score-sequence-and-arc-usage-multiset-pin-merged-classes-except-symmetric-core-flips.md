---
id: THM-5155
title: "Constraining the merged tournament classes from several observers at once. (A) PROVED: replacing a strong core S, flanked on both sides by nonempty singleton runs, by its converse S^op preserves the arc-usage multiset (THM-5154), and it preserves the global score sequence iff the score multiset of S is self-reversed (s -> |S| - 1 - s). (B) FINITE-EXACT, every class of orders 3-9: the joint observer (score-sequence pair {scores(T), scores(T^op)}, arc-usage multiset) separates all merged classes through order 8. At order 9 it fails on exactly 6 pairs, each TT2 => S => TT1 against TT2 => S^op => TT1 for a rigid, non-self-converse strong 6-tournament S with a self-reversed score sequence. (C) FINITE-EXACT: at order 9 the merged classes number 97144, against 1520 values of H, 184 symmetrised score classes, 50511 slot polynomials, 97121 arc multisets, 29742 of (H, C, |Aut|, c3), 63686 of (scores, Q) and 97138 of (scores, arc multiset)"
status: >
  PROVED: (A), from the component weight formulas of THM-5154 (A) and the score shift by the singletons after a core.
  FINITE-EXACT: (B) and (C) over every class of orders 3-9.
  - Inputs: classgen lists with self-converse flags, |Aut| and converse codes; fprof profiles; kps_s29 arc multisets.
  - The six residual pairs were checked structurally: runs (2, 1); cores converse but not isomorphic; core score
    sequences (1,2,2,3,3,4) four times, (2,2,2,3,3,3) once, (1,1,2,3,4,4) once.
  Exploratory, not pre-registered.
  NOT CLAIMED: completeness beyond order 9; Chan–Pak (OPEN).
session: kind-pasteur-2026-11-13-S30
source: 05-knowledge/results/kps_s30_observer_lattice_20261113/README.md
proof: 05-knowledge/results/kps_s30_observer_lattice_20261113/README.md
related:
  - THM-5154
  - THM-5153
  - HYP-9680
  - HYP-9681
---

# THM-5155 - score sequence and arc-usage multiset pin the merged classes except symmetric-score core flips

The packet README has the observer table and the residual list. Replay: `scripts/observers.py` and
`scripts/score_arcs.py` in a directory with cls_n.txt, prof_n.txt and arcms_n.txt.
