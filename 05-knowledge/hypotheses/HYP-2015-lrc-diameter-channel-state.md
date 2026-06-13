---
id: HYP-2015
status: OPEN
source: codex-2026-06-01-S532
related:
  - HYP-2007
  - HYP-2009
  - HYP-2011
  - HYP-2012
  - HYP-2013
  - HYP-2014
---

# HYP-2015: n=14 diameter-channel state refines the independent-pair metric

**Claim.** HYP-2012 identifies independent pairs as the multi-channel metric.
The n=14 refinement is that one must also track the state and matching shape of
those independent pairs inside a fixed scaffold, especially in the diameter
channel.

**K4 model.** With the scaffold `0->1->2->3` plus `0->3` fixed, the two
independent crossing pair-arcs `(0,2)` and `(1,3)` determine the ordinary
isomorphism class of the four-vertex tournament.  Their four states hit all
four unlabeled K4 classes:

```text
00 -> scores (1,1,2,2), H=5, c3=2, SCC=(4)
01 -> scores (0,2,2,2), H=3, c3=1, SCC=(3,1)
10 -> scores (1,1,1,3), H=3, c3=1, SCC=(3,1)
11 -> scores (0,1,2,3), H=1, c3=0, SCC=(1,1,1,1)
```

Across all K4 choices of two disjoint free arcs and fixed orientations of the
remaining four arcs, 24/48 scaffolds have this complete-coordinate property.

**n=14 evidence.** In the clasp-deleted regular 14-gon, channels `d=1..6` each
have 12 hidden chords, maximum matching size 6, and 7 maximum matchings.  The
diameter channel `d=7` has 6 hidden chords, maximum matching size 6, and exactly
1 maximum matching: all diameter chords are mutually independent.  Thus the
diameter channel is a pure six-bit independent-pair state plus one singleton
runner.

**Prediction.** A wall-only n=14 counterexample must carry a compatible state
vector on the six independent diameter pairs.  If that state is not the
AP-balanced state, then either the outside clasp opens or endpoint/resonance
debt exports to a deeper labelled layer.  Equivalently, "wall-only => AP" should
factor through a finite classifier of independent-pair channel states under a
fixed non-diameter scaffold.

**Next tests.**

1. Add `diameter_pair_state` to n=14 hard-row and finite-check survivor audits.
2. Correlate that state with source measure, endpoint debt, and
   `abs(debt/credit)`.
3. Search for K4 windows inside the n=14 permutohedral handoff graph: fixed
   scaffold plus two independent pair toggles determining a local marked class.

**Files.** `04-computation/lrc_independent_pair_channels_s532.py`;
`05-knowledge/results/lrc_independent_pair_channels_s532.out`;
`07-reflections/lrc-independent-pair-channel-state-s532.md`.
