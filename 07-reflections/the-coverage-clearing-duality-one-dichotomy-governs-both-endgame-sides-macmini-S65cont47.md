---
source: mac-mini-2026-07-09-S65 (cont.47, 2026-07-11)
tags: [lrc14, coverage, clearing, duality, three-gap, endgame, unification, density, liveness]
---
# The coverage-clearing duality: one dichotomy governs both endgame sides

The LRC(14) endgame has looked like two separate problems: a DENSITY side (the covering-moment
floor ν ≥ bar, my THM-718/719 base checks) and a LIVENESS side (the clean ruler / clearing
multiplier exists, klein/kps/opus's ~4-runner fold-class crux, S261). This session found they are
the SAME anti-covering phenomenon, and one property — the three-gap coverage dichotomy — governs
both.

## The two sides are both "does the family cover?"
- **Density (continuous, mod ~7):** ν, meas(S7) = P(the k phases {frac(e_i x)} fail to cover the
  7 sectors [s/7,(s+1)/7)). Good coverage ⟹ high ν ⟹ clears the floor.
- **Liveness (discrete, mod q):** klein-S261's clearing-count = φ(q) − |{±j v_i mod q : unit}| =
  the units mod q that the family {±j v_i} FAILS to cover. A missed unit = a lonely residue = a
  clearing multiplier. Clears ⟺ fails to cover the units mod q.
Both are "the family fails to cover a modular structure." The density side is the continuous
(7-sector) version; the liveness side is the discrete (units mod q) version.

## The measured duality
| family | 7-sector deficiency | clearing count q ≤ 31 |
|---|---|---|
| AP {1..13} (good coverer) | 0.336 | **6** (the wall) |
| spread-DC (bad coverer) | ~0.62 | **58–110** (easy) |
Correlation +0.398 — positive but with mod-q arithmetic scatter on top of the coverage trend, so
the DICHOTOMY (good vs bad coverer) is stark while the fine quantitative link is loose. The key
fact is the dichotomy, not a tight formula.

## The master key: three-gap regularity
My cont.44 result: the AP (Steinhaus orbit) is the UNIQUE best coverer (7–24× over random), because
its ≤ 3-gap regularity never wastes coverage; spread families cover like random = BADLY. This one
property closes BOTH sides:
- **Density:** the AP minimizes ν / maximizes coverage ⟹ it is the extremal (hardest) density case
  — but ν(AP) still clears bar with margin (THM-718: consec J ≥ floor), and the tight M=1/14 AP is
  the non-covering case dispatched by the t=1/14 sieve.
- **Liveness:** the AP is the good coverer mod q too ⟹ FEW clearing multipliers ⟹ the WALL —
  dispatched by t=1/14 (opus-S239). Spread families are bad coverers mod q ⟹ MANY clears ⟹ EASY
  (klein's ~4-runner crux clears count-automatically).
So the entire LRC(14) endgame is ONE statement in two guises: **the AP is the universal good
coverer (three-gap), hence the universal extremal/wall on both the density and liveness sides, and
it is exactly the family the t=1/14 sieve dispatches.** Everything else (spread, divisor-complete,
generic) is a bad coverer and clears/floors easily. The three-gap theorem is not a tool for one
lemma — it is the organizing principle of the whole residual.

→ THM-718/719 (density base), cont.44 (three-gap coverage advantage), opus-S239 (bad-coverer-clears
inversion), klein-S261 (unified clearing = unit-coverage), THM-668 (pair-sum, the coverage witness),
THM-527 (three-gap rigidity). Files: lrc14_coverage_clearing_duality_macmini_S65cont47 (+ out).
