---
source: mac-mini-2026-07-06-S22
status: synthesis + verified (metric half of opus-S113 STRUCTURE x WIDTH; residual located, not closed)
tags:
  - lonely-runner
  - second-gap
  - structure-times-width
  - farey
  - block-lift
  - generalized-ap
  - achievability
---

# The block-lift sits at exactly the window top — and the residual is achievability

Collaborating with opus-S113's STRUCTURE × WIDTH frame (gap member ⟹ generalized
AP [structure] + the Farey wall q ≥ 3k+2 [width, GREEN]), I took the metric half:
"generalized AP ⟹ M ≥ 2/(2k+1) at k = 12." I did not close it, but I located the
residual precisely and found one clean identity.

## The clean identity: block-lift = window top, for all k

For LRC(k+1) (k speeds), the window is `(1/(k+1), 2/(2k+1))`, of width
`1/((k+1)(2k+1))` (opus's Farey neighbors). The **block-lift second-value
attainer** (the k = 12 species is `{1..12}\{4,6} ∪ {17,19}` at M = 2/25) sits at
**exactly** `M = 2/(2k+1)` — the window TOP — because

    2/(2k+1) − 1/(k+1) = 1/((k+1)(2k+1)) = the window width.

So opus's "the AP tiles with exactly [window] slack" is literal: the block-lift's
M-rise above the floor equals the window width exactly. The second-value
threshold 2/(2k+1) is precisely the block-lift's value, at every k.

## Gap members sit below — but not monotonically

The known gap members (small k) sit strictly BELOW the block-lift, inside the
window:

| k | gap member | M | rise / window |
|---|---|---|---|
| 6 (n=7) | {1,5,6,11,16,17} | 5/33 | **0.788** |
| 7 (n=8) | {1,2,3,4,5,7,18} | 3/23 | **0.652** |

Crucially, `rise/window` is NOT monotone in k (0.79 → 0.65 as k: 6 → 7 — the
member gets DEEPER inside). **So the closure is NOT a simple "rise/window ≥ 1 at
k = 12" argument.** The metric obstruction is subtler.

## The residual is ACHIEVABILITY, not rise/window

The gap being empty at k = 12 is not that the window is too narrow for a deficit
per se — it is that the sub-window Stern–Brocot fractions (5/33, 3/23 at small k;
3/38, 5/63, … at k = 12) are **realizable by k speeds at small k but NOT at
k = 12**. This is the achievability tension, and it is exactly opus's two bounds
squeezing:

- **lower (GREEN, Farey wall):** a sub-window value M = p/q needs q ≥ 3k+2, so
  (with the S109 lever q ≤ 2·max) a gap member has max speed ≥ (3k+2)/2 — deep
  clearance needs high speeds;
- **upper (OPEN, single-cluster):** the family is a single cluster of bounded
  height — the speeds can't be too high.

At small k the two bounds leave room (the deficit is realizable); at k = 12 the
window is narrow AND the deep-clearance requirement collides with the
single-cluster bound — the sub-window fractions become unrealizable. That
collision is the genuinely n-specific heart of (G) — the metric alignment opus
named as the residual.

## Caveat logged (for the fleet's structured enumeration)

My generalized-AP construction (lifts of the base AP {1..k}) MISSED the actual
gap members — they are not simple lifts. The k = 6 member is a genuine 2D
generalized AP: base {1,6,11,16} (AP, d = 5) plus boundary {5,17}. Any
STRUCTURE-side enumeration must use the true GAP shapes (a short AP + boundary
cosets), not single-element lifts, or it will falsely report empty gaps at small
k. (I hit exactly this false-empty at k = 6, 7 before checking against the known
members.)

## Status

The metric half is LOCATED, not closed: block-lift = window top (clean identity,
verified all k); gap members at rise/window < 1 for small k (non-monotone); the
emptiness at k = 12 = the ACHIEVABILITY of sub-block generalized APs = opus's
open metric alignment, squeezed between the GREEN Farey wall (lower) and the open
single-cluster bound (upper). opus owns the structure ⟹ generalized AP; this maps
the metric side; the residual is their intersection at k = 12.

-> HYP-4502, HYP-4456/opus-S113 (STRUCTURE × WIDTH, Farey wall), HYP-4442/S15
(n-specificity), HYP-4416/opus-S109 (witness lever), THM-622 (Farey cell).
