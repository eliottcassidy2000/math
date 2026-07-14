---
id: THM-754
title: The 7-clock PARTITION theorem + the clean-slot criterion -- the six unit cells a/7 +- 1/14 and the origin cell TILE the circle exactly (7 x 1/14 = 1/2: k = 7 is the UNIQUE self-dual clock), so the (7,a)-slot family is a PARTITION of LRC(14), not a restriction; tight bodies' witnesses live on the CELL CORNERS; the clean-slot criterion (explicit delta = 1/(14 c_min), thresholds (1,3,5) c_min) is THM-748's mod-7 sibling -- rigorous, fires on 16% of the covering census, silent on the extremals: the full-window k=7 survival IS (A), whose 2x self-similar gap = klein-S306's (H) in slot coordinates
status: PROVED (both parts; the partition is exact arithmetic, the criterion's witness verified exactly) + HONEST (the criterion is a minor tile; the k=7 full window is the problem itself, re-coordinatized -- that is the DISCOVERY, not a defect)
source: opus-2026-07-14-S288 (owner prompt: prove the k=7 slot lemma with the perspective frame)
depends_on: []
related:
  - opus-S287 (HYP-6685: the scarcity collapse to k=7 -- EXPLAINED here)
  - klein THM-748 (the k=2 parity mechanism; the clean-slot criterion is its mod-7 sibling)
  - klein-S295 (the middle-reach localization -- re-derived as the tiling corollary)
  - klein-S306 / THM-753 (the one-step peel's hypothesis (H) = the same residue in peel coordinates)
---

# THM-754 -- the 7-clock partition theorem and the clean-slot criterion

## (i) The partition theorem (why k = 7 is terminal)

For every a: a/7 + 1/14 = (a+1)/7 - 1/14 (exact).  Hence the seven closed cells
[a/7 - 1/14, a/7 + 1/14] (a = 0..6) TILE the circle, overlapping only at corners; the six unit
cells cover [1/14, 13/14].  k = 7 is the UNIQUE clock with k x (1/14) = 1/2: for k >= 8 the
slot windows are genuine restrictions (they can all die -- and do, in the S287 scarcity map);
for k <= 6 they overlap and miss cell middles.  COROLLARIES:
  1. "Some (7,a)-slot survives" <=> the 1/14-safe set meets [1/14, 13/14] <=> (for bodies whose
     smallest speed kills the open origin cell) L-positivity itself.  The (7,a)-family is a
     PARTITION of LRC(14) -- the S287 collapse is explained: the 7-clock's slots are the last
     standing because they cannot all die without the body literally covering.
  2. TIGHT BODIES' WITNESSES LIVE ON CELL CORNERS: at M = 1/14 exactly (AP, GW), the safe set
     is finitely many points and each is a corner t = (2a+1)/14 of the tiling (verified: the
     AP's witness 1/14; both extremals' 5 feasible slots are the corner-adjacent ones).
  3. klein-S295's localization (L({1} u C) > 0 <=> G(C) reaches the middle) is the tiling
     corollary for bodies containing 1.

## (ii) The clean-slot criterion (THM-748's mod-7 sibling; rigorous, explicit)

Let C7 = multiples of 7 in V, c_min = min C7, delta* = 1/(14 c_min), and fix a unit a mod 7.
If (a) every carrier c has ||c delta*|| >= 1/14, and (b) every non-carrier v satisfies
v <= c_min when rho((va) mod 7) = 1/7, v <= 3 c_min when rho = 2/7, v <= 5 c_min when rho = 3/7,
then t* = a/7 + delta* has min clearance >= 1/14 EXACTLY (clearance algebra: carriers see
c delta*; non-carriers see rho - v delta* >= rho - v/(14 c_min) >= 1/14 by the thresholds).
M(V) >= 1/14 with the explicit rational witness t*.  PROOF: three lines per case.  QED.

VERIFIED: fires on 489/3000 (16%) of random covering bodies with exact witnesses; SILENT on the
classic extremals (their second residue-lifts always exceed c_min) -- an honest minor tile.

## (iii) The boundary (what the k = 7 "lemma" really is)

The full-window (7,a)-survival -- all delta in (0, 1/14], the combs of the large non-carriers
reopening windows -- is NOT a lemma below (A): by (i) it IS (A), re-coordinatized into six
delta-windows permuted by the unit group (the pair sector).  The union-bound gap there is the
familiar 2x, SCALE-SELF-SIMILAR (the same 2x at the global scale and inside each window):
this is the named residue, and it equals klein-S306's one-step-peel hypothesis (H)
(disc_v < 6 |G'_P|^2 for moderate non-aligned v) in slot coordinates.  The frame's final
sentence: the problem's last inequality lives on the unique clock whose half-cell is the
origin band -- the incompressible frame measuring itself.

## Files

05-knowledge/results/lrc14_seven_clock_thm754_opus_S288.out (tiling check, criterion witnesses,
census, extremal silence).
