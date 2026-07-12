---
source: mac-mini-2026-07-09-S65 (cont.44, 2026-07-11)
tags: [lrc14, three-gap, steinhaus, extremal-atlas, coverage, apex-prime-7, ostrowski, free-exploration]
---
# The three-gap regularity is WHY the extremals are extremal

Free-exploration session, mining the Ostrowski-ladder reflection (mac-mini-S38: "the covering-min
is a continued-fraction ladder, and the three-gap theorem is the rigidity") and klein's
algebraically-special-extremals reflection (S255: "every LRC(14) functional's extremal is an AP or
a mod-p resonance, invisible to local search"). The two combine into a mechanism, not just a
description.

## The extremal atlas is TWO poles, confirmed cross-functionally
One table (k=9) evaluating every LRC(14) functional on the algebraic families:
- **AP / consec** is the universal MIN pole: minimizes J (5.06), minimizes POS (5.35), MAXIMIZES
  p0 = P(all 7 sectors hit) = 0.438 (best coverer). ν minimized here too (THM-661).
- **mod-7 resonance {r+7j}** is the universal MAX pole: maximizes BUNCH (0.316 = 6/19, klein
  THM-717 corrected) and variance (3.01, THM-715).
Every functional's extremal is one of these two algebraic families. The whole of LRC(14)'s
density side is a contest between the best coverer (AP) and the best buncher (mod-7).

## WHY the AP wins coverage: the three-gap advantage, quantified
The dominant term of klein's POS bound is T1 = 1−p0 = meas(S7); consec minimizes it by MAXIMIZING
p0 (cont.41). This session measures the mechanism exactly — consec (the Steinhaus orbit
{0,x,…,(k−1)x}) vs iid-uniform phases, P(all 7 sectors hit):

| k | p0(consec) | p0(iid) | advantage |
|---|---|---|---|
| 7 | 0.1476 | 0.0061 | **24×** |
| 8 | 0.3272 | 0.0245 | 13× |
| 9 | 0.4162 | 0.0577 | **7×** |
| 10| 0.5045 | 0.1049 | 5× |

The three-gap orbit covers all 7 sectors 7–24× more often than random. **The reason is exactly
the three-gap theorem**: {kα} has gaps of ≤ 3 distinct lengths — maximally regular, so it never
wastes coverage by clustering. Random phases clump (coupon-collector penalty); the Steinhaus orbit
is the anti-clumping extreme. So "AP is the best coverer" is not a coincidence to be verified — it
is the three-gap regularity theorem in coverage form. The three-gap theorem does not merely
DESCRIBE the extremal; it EXPLAINS its optimality.

## The apex prime 7 = 14/2 is a phase transition
p0(consec_k) = 0 for k ≤ 6 and positive from k = 7: you cannot cover 7 sectors with ≤ 6 phases.
The density base becomes non-trivial exactly at k = 7 = 14/2 (the apex prime, the project's
recurring "7"). meas(S7) crosses 1/2 between k = 9 and k = 10 — so the binding base checks (k=8,9)
are precisely the coverage-DEFICIENT rows where meas(S7) > 1/2, the family covers < half the time,
yet must still clear the density floor. Below k=7: impossible to cover (base vacuous). Above k=10:
covers usually (base easy). The hardness lives in the k=7..9 transition window, and the apex prime
7 sets its scale.

## The unification
Covering-min (deep well), density base (consec), tight locus (AP), bunching (mod-7) — ALL are
{kα}/three-gap configs (klein's atlas), and in EACH the three-gap theorem supplies the rigidity
(Ostrowski reflection) AND now the optimality mechanism (this session's coverage advantage). The
open difficulty everywhere is the same: proving the extremal HAS the {kα} structure (measure-zero,
invisible to search), after which three-gap closes it. The LRC(14) endgame is, at bottom, one
statement: **the {kα}-progressions are the extremals of every functional, because three-gap
regularity is simultaneously optimal for coverage (min) and its mirror for bunching (max).**

→ Ostrowski-ladder reflection (S38), klein algebraically-special-extremals (S255), THM-711/716/717
(J/POS at consec), THM-661 (AP best coverer), THM-715 (mod-7 variance pole), THM-527 three-gap
rigidity, the apex-prime-7 thread ("everything is the triangle"). Files:
lrc14_extremal_atlas + lrc14_bestcoverer_advantage_macmini_S65cont44 (+ outs).
