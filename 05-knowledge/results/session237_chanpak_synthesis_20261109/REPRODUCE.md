# Reproduction

Run from repository root with Python3; no third-party packages, downloads,
private data or precomputed tournament-class lists are required.

    python3 05-knowledge/results/session237_chanpak_synthesis_20261109/replay.py

This runs the compact producer and independent programs normally and under
python -O, compares outputs, and tests --fail in both modes. It additionally
checks the preserved full six-vertex symmetry-probe normal/optimized outputs.
For a fresh replay of that slower census too, add --full.

The exact scopes are:

- Cycle producer: all1099 labelled graphs of orders1..5, 16 K1 completions,
  independent path/packing controls. Auditor: all1100 graphs including empty,
  112 randomized completions and a separate graph-deletion polynomial engine.
- Arithmetic producer: root automata M2..40, 393216 recurrence checks,
  256 literal tournaments, 2990 literal pumping controls. Auditor: independent
  Moore refinement, 3570 pumping witnesses, joint-image checks, 48 literal tails.
- Symmetry producer: all202013 labelled rooted tournaments of orders1..6,
  4275 literal short-tail checks. Independent literal auditor: all5405 rooted
  labelled tournaments through order5; actual groups and supplied witnesses.
- Endpoint producer: all157526 positive nonsink rooted cases through order6,
  equality count3294; all1023 apex masks of lengths0..9. Auditor: all3446
  positive nonsink cases through order5, 566 zero-endpoint controls, 127 apex masks.
- Endpoint arithmetic extension: moduli1..50, 524288 evaluations and 256
  literal endpoint DPs. Independent observer minimization: moduli1..40 for
  E, normalized full count (H-1)/2, and actual H.

Universes are labelled unless expressly stated otherwise. The graphs of
orders up to37 in the cycle construction are checked structurally, not all
counted by exponential Hamilton DP. Long words are proved by induction and
pumping; their infinite claims are not inferred from finite data. Normal and
optimized output comparisons do not replace independent mathematical audit.

All executable mathematical checks use explicit exceptions. --fail deliberately
raises before research work; its nonzero exit verifies the failure path in both
Python modes. Final REPLAY.json and SHA256SUMS record the exact artifacts.
