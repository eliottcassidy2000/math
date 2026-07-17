# THM-944 — The below-pair dilate count and the B5 race scoreboard (death-star-2026-07-17-S38)

**Status:** PROVED (Lean, kernel-pure — `TournamentH7/LRCB5Race.lean`; verify the
build report in the session log). Source: HYP-7182. Assembles THM-940/942/943 with
THM-939's trap into the formal race for B5-positivity on the dense core.

## Statement

1. `dilatePairs_card_le` (**the counting lemma**): an injective positive 13-family
   carries at most 12 doubling pairs — each pair is determined by its top, and the
   minimum is never a top.
2. `dilate_top_below_pair` (**the trap confinement**): on a chain-dense family the
   top of every doubling pair sits at element index ≤ j+1 — dilate pairs live
   entirely in the bottom block. (THM-939's A1 instantiated at the explicit (1, −2)
   coefficient vector.) With (1): **a ChainDenseCore family has ≤ 12 dilate pairs,
   all confined below the dense pair** — the directive's count, formal.
3. `deviation_lower` / `deviation_upper_of_mem`: proved two-sided per-subset bounds —
   the floor −(q−1)/7^|T| unconditionally; the ceiling (q−1)/7 − (q−1)/7^|T| through
   THM-942A's exact singleton whenever T contains a gcd-1 element.
4. `B5_race_scoreboard` (**the race, formal**): for q ≥ 14 with all speeds gcd-1 and
   named odd-tail hypotheses tail3/tail5 (the exact sums they bound):
       B5 ≥ (q−1)·2052/16807 − 78(q−1)/49 − 715(q−1)/2401 − tail3 − tail5.
   k = 1 HELPS (singles ≤ 0 enter negatively); k = 2 and k = 4 pay at worst the
   proved trivial floors; the odd tails are exactly what equidistribution owes.
   Sharper per-stratum k = 2 floors (THM-943A's exact dilate mass at 14 | q, with
   the count from (1)+(2)) slot in without touching the frame.

## Honest accounting

The trivial even floors do NOT make the right side positive (78/49 alone exceeds the
equilibrium) — closing the race outright still needs the odd-tail bounds, which is
the equidistribution heart, now isolated as exactly two named rational quantities
per family. What IS closed: the dilate content of the k = 2 term (exact via
THM-943A), its location (below the pair, via the trap), and its count (≤ 12).

## Referee

`b5race_referee_deathstar_S38.py`: counting PASS (30k), confinement PASS (3766
chain-dense families), scoreboard inequality PASS exactly (tails set to their true
sums) on gcd-filtered random (v, q).
