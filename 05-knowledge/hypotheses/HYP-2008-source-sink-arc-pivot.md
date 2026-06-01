---
id: HYP-2008
status: PARTIALLY-TRUE
source: oracle-2026-06-01-S530
related:
  - HYP-2000
  - HYP-1998
  - HYP-2003
  - THM-382
---

# HYP-2008: the source-sink (apex) arc is the pivot — five guises of one arc

Base path 0->..->(n-1) (source 0, sink n-1); tiles = non-path arcs; the
SOURCE-SINK arc (0,n-1) has max range n-1.

**VERIFIED (`source_sink_apex_arc_s530.py`):**
1. APEX = max-range tile = staircase hypotenuse. Flipping it on the transitive
   tournament gives H = 1 + 2^{n-2} (n=3..8: 3,5,9,17,33,65) = the maximal single-
   tile jump = transitive's big SC neighbour on the principal line.
2. MASTER SWITCH between the two block-extremes (round tournaments, #SCC in {1,n}):
   for round tournaments from circle gaps (0 mismatches, n=4,5,6),
     transitive (#SCC=n) <=> largest boundary gap > 1/2 (semicircle, THM-374);
     strongly connected (#SCC=1) <=> largest gap < 1/2 (Moon: boundary = Ham cycle).
   The source-sink arc spans the largest gap; its orientation is the switch.
3. LRC: the observer sits in the largest gap, so the source-sink arc = the
   observer's gap; THM-382 loneliness <=> both ends of this arc are >= 1/n from
   the observer. The two lonely regimes (observer-outside/transitive vs observer-
   inside/wrap) are selected by this arc; the tight AP/regular-polygon case (S525)
   is apex = back-arc (boundary = full directed cycle).

**FIVE GUISES of (0,n-1):** boundary-cycle closer (polygon outside) = apex tile /
hypotenuse (tiling) = maximal jump H=1+2^{n-2} (principal line) = #SCC n<->1 switch
(block extremes) = observer's loneliness gap (LRC).

**OPEN / to prove:**
- (A) For round tournaments, transitive <=> apex aligned & largest gap>1/2 exactly
  (verified n<=6; prove general via THM-374 + Moon).
- (B) "apex back-arc forces #SCC=1" rigorously (Moon: base path + closing back-arc
  = Ham cycle => strongly connected). State as a clean lemma.
- (C) LRC leverage: the tight case is apex=back-arc; does fixing the apex orientation
  split LRC(n) into the transitive branch (easy, semicircle) and the wrap branch
  (the AP-hard core)? Combine with S526 covering: the source-sink arc is the
  longest chord = lowest harmonic; relate to the resonance grading (S529).

**Files:** `04-computation/source_sink_apex_arc_s530.py` (+.out). Reflection:
`07-reflections/the-source-sink-arc-is-the-pivot-s530.md`.
