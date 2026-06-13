---
id: HYP-1800
status: EXPLORATORY_ENGINEERING
source: opus-2026-05-30-S5
related:
  - HYP-1797
  - HYP-1798
  - HYP-1776
---

# HYP-1800: OCF Active Ranking Acquisition

## Statement

For active pairwise ranking, the next comparison should be chosen by expected
drop in Hamiltonian-path ambiguity and Omega residue complexity:

```text
score(i,j) = E[Delta H after querying i vs j]
             plus tie-breakers from Omega cycle packets and residue rank.
```

The output should be a ranking plus an ambiguity certificate, not just a sorted
list.

## Evidence

The application probe identified OCF-guided active ranking as the strongest
near-term practical route.  The residue/phase/incidence contrast adds a
diagnostic layer: some ambiguity is localized around near-kill residues, some
is broad phase/fiber ambiguity, and incomplete/weighted data first needs a
completeness-defect warning.

## Predictions

1. Expected `H`-drop acquisition should need fewer pairwise queries than
   uncertainty-only heuristics on small exact benchmarks.
2. Omega cycle packets should provide human-readable explanations for unstable
   rankings.
3. Completeness defect should predict when thresholded weighted comparisons
   make hard tournament invariants unstable.

## Test Plan

1. Build a small `H_topk` prototype for exact `n <= 9`.
2. Compare query counts against random, score-margin, and entropy heuristics.
3. Emit JSON with `H`, top Omega packets, deletion residue, and completeness
   defect.

## Sources

- `07-reflections/applications-probe-2026-05-30.md`
- `07-reflections/applied-residue-phase-incidence-programs.md`
- `04-computation/tournament_tda.py`
