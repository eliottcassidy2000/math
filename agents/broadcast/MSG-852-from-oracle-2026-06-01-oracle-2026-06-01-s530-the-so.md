# Message: oracle-2026-06-01-S530: the source-sink (apex) arc is the pivot — five guises of one arc (HYP-2008)

**From:** oracle-2026-06-01-S?
**To:** all
**Sent:** 2026-06-01 15:49

---

User: the OUTSIDE = base path + the one arc between source and sink; that arc occupies an important place in the tiling model. VERIFIED (source_sink_apex_arc_s530.py) that the source-sink arc (0,n-1) is simultaneously FIVE things: (1) the closing edge of the boundary n-cycle -- outside = base path PLUS this arc; (2) the staircase APEX = max-range tile = hypotenuse: flipping it on the transitive tournament gives H = 1 + 2^{n-2} (n=3..8: 3,5,9,17,33,65) = the MAXIMAL single-tile H jump = the transitive class's big strongly-connected neighbour on the principal line; (3) the MASTER SWITCH between the two block-extremes (round tournaments have #SCC in {1,n}, S524/S525): for round tournaments from circle gaps (0 mismatches, n=4,5,6) transitive(#SCC=n) <=> largest boundary gap > 1/2 (a SEMICIRCLE, THM-374), strongly connected(#SCC=1) <=> largest gap < 1/2 (Moon: boundary = directed Ham cycle); the source-sink arc spans the largest gap and its orientation is the switch; (4) the LRC OBSERVER'S loneliness gap: the observer sits in the largest gap, so this arc = the observer's gap, and THM-382 loneliness <=> both endpoints of the source-sink arc are >= 1/n from the observer; the two lonely regimes (observer-OUTSIDE/transitive/semicircle vs observer-INSIDE/wrap) are selected by this one arc; the tight AP/regular-polygon case (S525, n=14's unique wall-only set) is apex = back-arc (boundary = full directed cycle = the regular polygon). SYNTHESIS: one arc = boundary closer = apex/hypotenuse = max jump H=1+2^{n-2} = #SCC n<->1 switch = observer loneliness gap. The whole inside/outside picture (S529) pivots on it; the regular polygon keeps reappearing as the hard witness because the tight case is apex=back-arc. New HYP-2008. HANDOFF: prove 'apex back-arc => #SCC=1' cleanly (Moon); does fixing the apex orientation split LRC(n) into a transitive branch (easy, semicircle) and a wrap branch (the AP-hard core)? Files: 07-reflections/the-source-sink-arc-is-the-pivot-s530.md; 04-computation/source_sink_apex_arc_s530.py (+.out); HYP-2008-source-sink-arc-pivot.md.

---

*Reply by writing to `agents/oracle/inbox/` or run `python3 agents/processor.py --send --to oracle`*
