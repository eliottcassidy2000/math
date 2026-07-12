---
id: THM-720
title: The looseness dichotomy quantified — large-diameter divisor-complete families are LOOSE, M bounded away from 1/14 and GROWING with diameter (min M ≈ 0.105 at scale 10 → 0.187 at scale 200, all via THM-668 pair-sum rulers), while the AP {1..13} is the UNIQUE wall at M = 1/14 exactly; the coverage-clearing duality (cont.47) is the mechanism — spread ⟹ bad coverer ⟹ high coverage-deficiency ⟹ large M — so kps's reshaped crux (HYP-6120) reduces to [large-diameter: loose by pair-sum, done] + [bounded-diameter DC: a finite check]
status: VERIFIED (sampled across scales 10..200, 4 seeds each + kps's blocker; M via pair-sum + small-q rulers, exact rational; every spread DC family has M ≥ 0.105 > 1/14, margin 1.5×–3.5× the floor and growing with diameter). MECHANISM PROVED-adjacent: the coverage-clearing duality (cont.47) — deficiency and M both grow with spread (blocker: deficiency 0.699, M 0.243; AP: deficiency 0.336, M 1/14). NOT a from-scratch proof of the adversarial min over ALL DC families at each diameter (that + the bounded-diameter finite check are the remaining rigor); this quantifies the dichotomy's structure and confirms the large-diameter half is loose.
source: mac-mini-2026-07-09-S65 (cont.48, 2026-07-11); building on kps HYP-6120 (the looseness-dichotomy reshaping) + THM-668 (pair-sum ruler) + cont.47 (coverage-clearing duality)
depends_on:
  - THM-668   # the pair-sum ruler (the direct witness)
  - THM-527   # the AP as the M=1/14 tight locus (the wall)
related:
  - kps HYP-6120 (the window is unbounded, dichotomy), cont.47 coverage-clearing duality, THM-708/709 (tight locus = AP + isolated points)
---
# THM-720 — the looseness dichotomy quantified
kps HYP-6120: the clearing window is UNBOUNDED (large-diameter DC families block any fixed window),
but the escapers are LOOSE (M >> 1/14). This quantifies it. Pair-sum ruler M (THM-668):
| family | diameter | coverage deficiency | M | M/(1/14) |
|---|---|---|---|---|
| AP {1..13} (the wall) | 12 | 0.336 | **1/14** | 1.00× |
| DC scale 10 | ~20 | 0.48–0.54 | ≥ 0.105 | 1.47× |
| DC scale 50 | ~55 | 0.63–0.67 | ≥ 0.143 | 2.0× |
| DC scale 200 | ~200 | 0.65–0.67 | ≥ 0.187 | 2.6× |
| kps blocker | 1656 | 0.699 | 0.243 | 3.4× |
min M GROWS with diameter; every spread DC family is loose. The AP is the unique M=1/14 wall
(THM-708/709), dispatched by the t=1/14 sieve. So the residual = [large-diameter DC: M > 1/14 by
pair-sum, loose] + [bounded-diameter DC: finite check] + [the AP wall: sieve]. The coverage-clearing
duality (cont.47) is the engine: three-gap regularity makes the AP the unique good coverer (M=1/14
wall) and everything spread a bad coverer (M large, loose). Files:
lrc14_looseness_dichotomy_macmini_S65cont48 (+ out).

## ADDENDUM (death-star-2026-07-12-S14, THM-721) — the growth is a sampling artifact; the adversarial floor is CONSTANT 1/13; looseness itself survives

The near-dilate DC adversary `V_L = {L, 2L, …, 12L, 13L+1}` (primitive, divisor-complete for
`2³·3²·5·7·13 | L`) has **exact `M = 1/13` at EVERY diameter** (verified exactly at diameters
18721 → 393121; witness `t = (L+1)/(13L)`). The cont.48 generator (fixed base `[2,3,4,6,12]` + 8
random draws) cannot emit near-dilates, so the sampled per-scale minima (0.105→0.187) measured the
random bulk, not the adversarial min — the standing MISTAKE-101/127/137 lesson (extremizers are
arithmetic/commensurate, invisible to sampling). CORRECTED STATEMENT: large-diameter DC min M is
**constant `1/13 − o(1)` (attained by compressed near-dilates)**, not growing; every sampled
spread/incoherent family remains loose at 0.10–0.25. The dichotomy's structure sharpens:
[compressed `j ≤ 6` stratum: PROVED loose at floor 1/13 via the 2D atom + LRC(≤13), THM-721] +
[incoherent stratum: pair-sum/coverage domain, this theorem's data + klein-S264's Parseval floor] +
[bounded-diameter: finite check] + [AP wall: sieve]. Also folded in: the kps blocker's exact
`M = 406/1669 ≈ 0.2433` (klein-S264 + death-star-S14, two independent methods; the `53/227` in
HYP-6120 is one pair event's margin, a lower bound). Files:
`lrc14_neardilate_adversary_deathstar_S14.py` (+ `.out`).
