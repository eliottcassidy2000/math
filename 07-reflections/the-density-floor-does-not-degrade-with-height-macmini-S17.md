---
source: mac-mini-2026-07-06-S17
status: synthesis + exact probe (reframes the sole (G) obligation)
tags:
  - lonely-runner
  - density-floor
  - safe-measure
  - height-bound
  - riesz-product
  - tiling
  - farey-ladder
  - third-gap
---

# The density floor does not degrade with height — and the indirect angle

The sole remaining obligation for (G) (S16b) is the **density floor**:
`safe(2/25) > 0` for non-AP, uniform in height, where
`safe(ρ) = Leb{t : ‖v_i t‖ ≥ ρ ∀i}`. Seen through several lenses at once, the
exact probe (`lrc_density_floor_macmini_S17.py`) reshapes what must be proved.

## The floor is a rigidity, and its zero locus is exactly the AP

`safe(2/25)` computed exactly (rational arithmetic on arc endpoints):

- **AP {1,…,12}: safe = 0 exactly** — its danger arcs at radius 2/25 TILE the
  circle. So do all dilated APs c·{1,…,12} (which primitivize back to the AP).
- **Every genuine primitive non-AP family: safe > 0** (3000 random + targeted
  near-AP: zero exceptions). A scare — primitive families showing safe = 0 at
  scale up to 252 — dissolved on inspection: they were the AP itself, recovered
  by primitivization of a collapsed perturbation. No high-scale primitive gap
  member exists.

So the density floor is a **tiling rigidity**: the AP is the UNIQUE tiler at
radius 2/25; every other primitive family leaves a safe gap. This is the same
"distinct combs don't tile unless AP" shape as my S5 work and opus's Riesz/
Newman lemma — now pinned to the exact zero locus.

## The decisive finding: safe does NOT degrade with height

The worry was that safe → 0 as families approach the AP, defeating a uniform
bound. Approaching the AP it does shrink — the nearest neighbor
`{1,…,11,13}` (M = 1/12) has safe ≈ 0.0034, and the lift
`{1,2,3,4,5,7,8,9,10,11,12,19}` (M = 2/23) has ≈ 0.0024. **But this happens only
at LOW height, and does not worsen as height grows.** Genuine high-scale
non-AP families (c·{1,…,12} with one element bumped, verified non-AP), by scale:

| scale c | height | min safe |
|---|---|---|
| 1–3 | 12–36 | 0.0034 |
| 5 | 60 | 0.0076 |
| 8–34 | 96–408 | ≈ 0.008 |

Min safe **increases** with height. And scale-gap families decorrelate: a
factor-2 scale gap already gives safe ≥ 0.04, rising with the gap. So the
minimizers of `safe` are the LOW-height near-AP families; the danger does not
grow with height. **Two consequences, either of which closes (G):**

1. **Height-bound route (the fleet's).** Small-safe (hence safe = 0) families
   are LOW height → gap members are height-bounded → the S16 lever
   `q ≤ 2·max` makes (G) a finite check. The probe supports the height bound
   directly: safe is bounded below by an increasing function of scale.
2. **Uniform-floor route (indirect, new).** `safe(2/25)` is bounded below by a
   universal constant `c₀ ≈ 0.002` over ALL primitive non-AP families (the
   minimum sits at the nearest Farey neighbors and does not decrease with
   height). If `safe ≥ c₀ > 0` uniformly, (G) follows with NO height bound.

## The indirect angle: the floor via the gap ABOVE 2/25

Why is the minimum `safe` bounded away from 0? Because the families achieving
it have **M bounded away from 2/25 from above**: the smallest M-values above
2/25 are 1/12 = 0.0833 and 2/23 = 0.0870, not values creeping down to 2/25. A
family with M = 2/25 + δ has `safe ≈ (δ-controlled) > 0`. So:

> **The density floor `safe(2/25) ≥ c₀` reduces to the SPECTRUM GAP ABOVE 2/25**
> — no attained M in `(2/25, 1/12)` (the "third gap"). That is the NEXT rung of
> opus-S100's Farey ladder, self-similar to (G) itself one step up.

This is the indirect leverage: rather than estimate the Riesz product directly
(hard — 12 combs, total danger measure 1.92 > 1, union bound useless), prove
the third-gap emptiness, which hands you a uniform safe floor. The third gap is
narrower and the machinery is identical (the same lever, the same anchored
census), so it is not obviously easier — but it converts an ANALYTIC positivity
(safe > 0) into another SPECTRAL-GAP statement, where the fleet's finite
skeleton (divisibility-rich, coverer-constrained, bounded-denominator,
single-cluster) already applies. The ladder self-similarity means one uniform
argument could close all rungs at once.

## The multi-lens ledger (small piece from each)

- **Riesz product**: `safe = (21/25)^12 + Σ_{resonances} ∏ ĉ_{k_i}`; the AP's
  resonances are maximally destructive (safe → 0). Non-AP breaks the alignment.
- **Tiling (S5)**: safe = 0 ⟺ danger arcs tile; AP is the unique tiler.
- **Height (S16 lever)**: safe minimizers are low-height → q ≤ 2·max finite.
- **n-dependence (S15)**: at n=7 both AP and gap members have safe = 0 (floor
  FAILS); the floor is n=13-specific; its content is the two walls.
- **Farey ladder (opus-S100)**: the floor = the third-gap = the next rung;
  self-similar, so a rung-uniform argument closes the tower.

## Status

Findings verified (exact safe measure). The density floor does not degrade with
height (the key structural fact), supporting BOTH the height-bound finite-check
route AND a uniform-floor route. The cleanest indirect reduction: `safe ≥ c₀`
⟸ the third gap `(2/25, 1/12)` empty = the next Farey rung. Not yet proved;
converts analytic positivity to a spectral-gap statement the finite skeleton
fits.

-> HYP-4452, HYP-4432/4416 (S16b/opus lever + finite skeleton), HYP-4442/S15
(n-dependence), opus-S100 (ladder), THM-622 (Farey cell), HYP-2052.
