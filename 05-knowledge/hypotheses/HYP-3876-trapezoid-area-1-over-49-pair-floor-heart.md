---
id: HYP-3876
title: THE TRAPEZOID AREA = 1/49, SPEED-INDEPENDENT -- the analytic heart of the drifting pair-floor. klein's two-tooth overlap trap(w1,w2,r) (LRCSpreadPairFloor Stage 1) integrates to EXACTLY 4h^2 = 1/49 over its support, for EVERY pair 0 < w1 <= w2. This is the mean the discrete pair-mass Sum_teeth trap(r_m) concentrates on in the dense-events regime (klein-S118 handed mac-mini the dense case). Lean: trap_integral, [propext, Classical.choice, Quot.sound].
status: CONFIRMED (Lean-formalized, LRCTrapArea.lean, trap_integral = 1/49 SORRY-FREE, [propext, Classical.choice, Quot.sound], builds standalone). The AREA (leading term) is exact; the discrete-sum discretization error (L/49-err) is the next piece. NOTE: self-contained (re-declares trap) because klein's LRCSpreadPairFloor does NOT build against the pinned mathlib v4.30.0 -- a fleet-wide mathlib-drift issue flagged to klein/kps.
source: mac-mini-2026-07-02-S19
related:
  - HYP-4023   # klein-S118 LRCSpreadPairFloor (Stages 1-3) -- this extends it with the aggregate area
  - HYP-4022   # klein-S117 ledger assembly: consumes the pair-floor credits >= (c-1)(L/49) - E
  - HYP-3874   # mac-mini JointRateCore (the joint/pair rate; the pair-floor is its j=2 case)
  - HYP-4021   # klein-S116 path-Hunter (the 7-wall the pair-floor crosses)
results:
  - 04-computation/lean/TournamentH7/TournamentH7/LRCTrapArea.lean
  - 04-computation/trap_area_verification_macmini_20260702.py
---

# HYP-3876 -- the trapezoid area = 1/49, the pair-floor's analytic heart

The LRC(14) proof reduces (klein-S115 far-cut dispatch) to two analytic things; the first is the
**pair-floor** `pairCredits >= (c-1)(L/49) - E` (klein-S117 `hledger_pos_of_bounds`). The commensurate case
is exact `1/49` (LRCCommensuration); the **drifting** case is klein's LRCSpreadPairFloor, and klein-S118
handed the **dense-events regime** (`14D > w1`) to mac-mini. This file lands its analytic heart.

## The theorem (LRCTrapArea.lean, `trap_integral`)
klein's `trap w1 w2 r = max 0 (min (2h/w2) ((h(w1+w2) - |r|)/(w1 w2)))` (`h = 1/14`) is the EXACT overlap
of two danger teeth at lattice residue `r` (her Stage-1 `clipLen_tooth_tooth`). Its integral over the
support is **exactly `1/49`, independent of the speeds**:
```
    ∫_{-S}^{S} trap w1 w2 r dr  =  4 h²  =  1/49 ,     S = (1/14)(w1+w2).
```

## The proof (sorry-free, pure real analysis)
The trapezoid is characterized piecewise (three lemmas, `trap_eq_plateau / _slope_pos / _slope_neg`):
- **plateau** `|r| <= h(w2-w1)`: `trap = 2h/w2` (the narrower tooth width);
- **slopes** `h(w2-w1) <= |r| <= h(w1+w2)`: `trap = (h(w1+w2) - |r|)/(w1 w2)` (linear down to 0).
Split `[-S,S]` at `±S0` (`S0 = h(w2-w1)`); `integral_const` on the plateau, `integral_add + integral_id`
on the two slopes; then `ring`:
- plateau `2 S0 · (2h/w2) = 4h²(w2-w1)/w2`;
- two triangles `(S-S0)²/(w1 w2) = 4h² w1/w2`;
- total `4h²(w2-w1)/w2 + 4h² w1/w2 = 4h² · w2/w2 = 4h² = 1/49`.
The speed-independence is the beautiful part: the plateau shrinks and the triangles grow so their sum is
constant. Verified numerically (5 diverse `(w1,w2)` pairs, ratio exactly `1.0000`).

## Why it matters
`4h² = 1/49` is the DENSITY of the pair overlap. In the dense-events regime the discrete pair-mass
`Sum_{w2-teeth} trap(r_m)` (residues `r_m` walking by `-D`, klein's Stage-3 `walk_one_wrap`) is a Riemann
sum of `trap` and concentrates on `L · (area) = L/49`. So this area is the LEADING TERM of the pair-floor
`pairmass >= L/49 - err`; the remaining `err` is the discretization/discrepancy of the residue walk (the
next piece). This is the exact `1/49` klein cited ("the trapezoid aggregate = 4h² per unit exact") as the
target, now a sorry-free Lean theorem.

## Honest scope
`trap_integral` (the area = 1/49) is exact and sorry-free -- the leading term of the drifting pair-floor.
It does NOT alone close the pair-floor: the discrete-sum-to-integral error bound (the `err`) remains, as
does wiring the aggregate into the Hunter ledger's `credits` slot. But the analytic heart -- WHY the density
is `1/49` regardless of speeds -- is now proved.
