---
id: HYP-4020
title: THE FAR-COUNT-7 DISPATCH -- integrates the LRC(14) endgame into ONE reduction, splitting the remaining content at the union-bound wall j=7. lrc14_of_farcut_split (cite : LRCUpTo13) reduces LRC14Statement to exactly four legs: (0) the LRC(<=13) citation; (1) hwindow = bounded families |v|<=22 (the machine-checked window census); (2) hle6 = COMPRESSED families (no ratio-13 dominant) with <=6 far entries (mac-mini lonely_of_simul_peel / kps block-6 territory -- the union bound holds); (3) hge7 = COMPRESSED families with >=7 far entries (the JointRateCore target -- the last open crux). The DOMINANT case (some entry >= 13x every other) is discharged UNCONDITIONALLY by klein-S114 top_ratio_lonely_13'. j=7 = 1/(2*(1/14)) is exactly the union-bound wall (each danger arc has measure 1/7, so seven tile the circle -- OPEN-Q-108). Sorry-free, foundational-axioms-only (#print axioms lrc14_of_farcut_split = [propext, Classical.choice, Quot.sound]); corpus green. INTEGRATES the fleet's closers (klein ratio-13 peel + kps window census/block engine + mac-mini simulpeel + mac-mini JointRateCore target) into a single theorem that names the two remaining residuals precisely and localizes mac-mini's joint-rate finish to the >=7-far compressed leg
status: VERIFIED (Lean, LRCFarCutSplit.lean, sorry-free, foundational-axioms-only, corpus green). Pure-logic dispatch over lrc14_of_covering_lonely + top_ratio_lonely_13' (S114); the far-count split farCount22 v = #{i : 22 < |v i|} at <=6 vs >=7. HONEST: this INTEGRATES + LOCALIZES; it does NOT close LRC(14). The four legs: citation (owner-sanctioned), window census (machine-checked, done), <=6-far compressed (union bound / simulpeel -- reducible), >=7-far compressed (the joint-rate crux, mac-mini actively wiring JointRateCore's per-cell identity to the concrete drifting arcs). Everything except leg (3) [and the finite middle band inside legs 1-2] is proved or machine-checked.
source: klein-2026-07-02-S115
depends_on:
  - HYP-4019   # S114: top_ratio_lonely_13' (dominant case, ratio 13) -- the DOMINANT leg
  - HYP-3875   # mac-mini-S18: lonely_of_simul_peel (the <=6-far leg)
  - HYP-3874   # mac-mini-S17: JointRateCore (the >=7-far leg target)
related:
  - HYP-3978   # kps-S21: block-6 (near-equal blocks <=6, the <=6 territory)
  - HYP-3974   # kps-S17: window census / peel20
  - OPEN-Q-108 # the r=2 residual after the simultaneous peel; the j=7 union-bound wall
results:
  - 04-computation/lean/TournamentH7/TournamentH7/LRCFarCutSplit.lean
---

# HYP-4020 — the far-count-7 dispatch (integrating the endgame)

## The reduction (LRCFarCutSplit.lean, sorry-free, foundational axioms only)
`lrc14_of_farcut_split (cite : LRCUpTo13) (hwindow) (hle6) (hge7) : LRC14Statement`, where
`farCount22 v := #{i : 22 < |v i|}`:
- **hwindow** — bounded families `|v| ≤ 22` (the machine-checked window census; kps/klein);
- **hle6** — compressed families (no ratio-13 dominant) with `farCount22 ≤ 6`;
- **hge7** — compressed families with `farCount22 ≥ 7`.

The **dominant** case (some `|v i0| ≥ 13·|v i|` for all `i ≠ i0`) is discharged **unconditionally** by
`top_ratio_lonely_13'` (klein-S114). Verified `#print axioms = [propext, Classical.choice, Quot.sound]`.

## Why `j = 7` is the right cut
`7 = 1/(2·(1/14))` is exactly the **union-bound wall** (OPEN-Q-108): each runner's danger arc has measure
`2·(1/14) = 1/7`, so up to six far arcs leave room (the union bound / `lonely_of_simul_peel` / block-6 all
work), but **seven arcs can tile** the circle — the union bound dies, and only the **JointRateCore** Δ-free
telescoping breaks past. So the split at `farCount22 = 7` cleanly separates the *reducible* leg (`hle6`) from
the *genuinely-open* leg (`hge7`).

## What this integrates and localizes
One theorem now assembles the whole fleet's endgame:
| leg | status | owner |
|---|---|---|
| citation | owner-sanctioned | policy |
| dominant (ratio ≥ 13) | **PROVED** unconditionally | klein-S114 |
| window census (`≤22`) | machine-checked | kps/klein |
| `≤ 6` far compressed | union bound / simulpeel / block-6 | mac-mini S18 / kps S21 |
| `≥ 7` far compressed | **the last open crux** | mac-mini JointRateCore (active) |

So LRC(14) is reduced to: the citation + the window census + the `≥ 7`-far compressed families. mac-mini is
wiring `JointRateCore`'s per-cell identity to the concrete drifting arcs — that is exactly leg `hge7`.

## Honest scope
This is **integration and localization**, not a closure: it does not finish LRC(14). Its value is a single
sorry-free theorem that (a) uses the ratio-13 peel to discharge the dominant case outright, (b) pins the two
remaining residuals at the meaningful `j=7` threshold, and (c) names the joint-rate target precisely so the
finish is one clean lemma (`hge7`) away.
