---
id: THM-743
title: The vertex-term pair-sector sharpening of THM-742 — per-event cost is the PAIR DIFFERENCE delta_v (the arcs' relative speed), buried events are free, and the strand side is sign-favorable for the lower bound; C2's vertex part 4max(J)V_R -> Sum_v^exposed delta_v (5676 -> 2172, shape 1); W0 452 -> 339 / 584 -> 513 (cumulative from crude: 1948 -> 339, 5.7x)
status: PROVED + VERIFIED-EXACT (zero violations over the W-spread, both shapes; exposed-event census exact)
source: opus-2026-07-13-S276 (owner prompt: attack the second-order vertex term with the perspective frame)
depends_on:
  - THM-742 (the bound being sharpened; all other terms unchanged)
related:
  - klein-S294 (HYP-6570 windowed-overlap Farey resonances -- the same pair-sector geometry: event locations are the difference-runner grids (k +- 1/7)/delta)
  - THM-739/740; the owner's (n-1)^2 = 13 + 2 T(12) frame (vertices ARE the double-counted pair sector)
---

# THM-743 — the vertex-term pair-sector sharpening

## Statement

In THM-742's bound |L - Area| <= C1/W + C2/W^2, the vertex part of C2 improves from
4 max(J) V_R to

>  Sum over EXPOSED pair events v of delta_v      (two-sided, stated with margin)
>  (1/2) Sum_v^exposed delta_v                     (lower-bound-only side, derived)

where a PAIR EVENT is a breakpoint u* interior to G_B at which the endpoints of arcs j and j'
coincide (possible only at delta u* == 0 or +-1/7 mod 1, delta = |j - j'| -- the
difference-runner's grid), and an event is EXPOSED iff its meeting point is not strictly
interior to a third arc.  So C2 = 3 Sum_{j in J} j^2 + Sum_v^exposed delta_v (stated).

## Proof (delta from THM-742's step 5)

1. **Per-event cost is the local relative speed.**  The fiber interval (safe gap or bad-run
   overlap) that the event annihilates closes at rate delta_v = |j - j'|, so within the single
   strip containing the zero-crossing the linear model's overshoot is a triangle of area
   <= (delta_v/W)(1/W)/2 per side, i.e. <= delta_v/(2W) on the strip-average (g) side and
   <= delta_v/W on the strand (f) side.  (THM-742 charged 2max(J)/W per breakpoint.)
2. **The strand side is sign-favorable for the lower bound.**  Truth >= linear model always
   (max(0, l) >= l for a collapsing interval; for merging bad runs the independent-run model
   double-counts overlap, understating the safe measure).  Hence for L >= Area - bound (the
   LRC-relevant direction) only the g-side triangle costs: (1/2) Sum delta_v.
3. **Buried events are free.**  If the meeting point lies strictly inside a third arc, the
   union is locally that third arc: no boundary change, no run-structure change, no capping.
4. All other THM-742 terms are unchanged. QED

## Verified numbers (exact census + zero-violation re-verification, W in {10..800})

| shape | Sum delta (all) | Sum delta (EXPOSED) | 742 charged | C2: 742 -> 743 (stated) | lower-side | W0: 742 -> 743 | C1 floor |
|---|---|---|---|---|---|---|---|
| 1 | 4312 | **2172** | 5676 | 8712 -> **3690** | 1845 | 452 -> **339** | 194 |
| 2 | 1518 | **1116** | 2376 | 4086 -> **1971** | 985.5 | 584 -> **513** | 428 |

Cumulative from crude THM-739: shape 1 W0 1948 -> 339 (**5.7x**), shape 2 2676 -> 513 (5.2x).

## The perspective structure (the census histogram)

delta-histogram (all/exposed), shape 1:  1: 22/22, 2: 30/30, 3: 54/38, 4: 72/44, 5: 84/54,
6: 90/50, 7: 90/58, 8: 84/32, 9: 72/26, 10: 54/18, 11: 30/14.

**Adjacent cluster members' events are few and ALWAYS exposed** (their arcs are neighbors in
the union, so their interactions happen on its boundary); **distant pairs' events are numerous
but mostly buried** (they meet inside the union's bulk).  The double-counted pair sector
carries the cost structure: cheap-and-visible near the diagonal, expensive-but-hidden far from
it.  Event locations (k +- 1/7)/delta are the difference set's three-distance grids -- the same
Farey/pair-sector geometry as klein-S294's windowed overlaps.

## What binds now

C1 = 2(#comp + Xi) (THM-742's first-order term): floors W0 at 194 / 428.  Next constants:
(i) Xi's 1/8-caps can use the EXACT segment-end heights ((1/2 - x)^2-type, arrangement-known);
(ii) cross-line joint equidistribution (the Q_s class) remains the deep reserve.

## Files

04-computation/lrc14_vertex_pair_sharpening_thm743_opus_S276.py (+.out): event census with
exposure test, delta-histograms, exact re-verification, new W0 solves.
