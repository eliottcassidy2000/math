---
source: opus-2026-07-08-S160
status: the L=9 correlated remainder CLOSED by kps-S89's EXACT conditional method extended to 2
  points (no grid aliasing): binding min D3 = 0.467183 (block_9 scale 4 + {3,40}), EXACT, margin
  +0.136 >= bar; large scale -> D3_inf^{(9)} = 0.524 (box bound). The conditional D3 matches exact
  Farey to 1e-5. This is the aliasing-immune replacement for S158's fixed-grid search.
tags:
  - lrc14
  - covering-floor
  - D3
  - conditional-exact
  - longest-AP
  - L9
---

# The L=9 exact conditional closure (no aliasing)

**opus-2026-07-08-S160.** Owner: work the natural next step (the L=9 correlated remainder). The clean
tool turned out to be kps-S89's EXACT conditional method (which closed L=10), extended from one
outlier point to two. This is aliasing-immune -- it replaces S158's fixed-`NG` search, which (as S157
found) aliases past prim-diam ~1500.

## The exact conditional structure (kps-S89 extended to 2 points)

For the L=9 family `E = block_9(scale d) u {p,q}`, condition on `a = frac(dx)`: the AP_9 phases are
`{frac(ja)}_{j=0..8}`, and over the `d` preimages `x = (a+k)/d` (`k=0..d-1`) the two point-phases are

> `p-phase = frac(pa/d + pk/d)`,  `q-phase = frac(qa/d + qk/d)`  -- `d` points on a LINE in `T^2`.

So `E[W^i](E) = mean_a mean_k W(AP_9(a); p-phase, q-phase)^i` -- **EXACT** in the point (`k`) direction
(the real `d`-sheet structure), with only the `a`-integral discretized (`Na`-Riemann, converges as
`1/Na`). Being exact in the dangerous point direction, it is immune to the `NG`-grid aliasing that
corrupted the S157 fixed-grid scan at large prim-diam (which reported `0.4464` for a shape whose true
`D3 = 0.4647`).

## The result

- **Decorrelated limit** `D3_inf^{(9)} = 0.524` (margin **+0.19**); the box bound (kps-S89 style: min
  `D3` over the a-priori moment box) covers large scale -- the 2-point phases decorrelate (the rank-2
  discrepancy, S159), so `D3 -> 0.524`.
- **Exact conditional finite check of the binding region** (one near point + one moderate point, per
  S158): binding **min `D3 = 0.467183`** at `block_9(scale 4) u {3,40} = (0,3,4,8,12,16,20,24,28,32,40)`,
  EXACT (the conditional value `0.467172` matches exact Farey `0.467183` to `1e-5`), margin **+0.136 >=
  bar**. Per-`d` minima rise with scale: `0.510 / 0.504 / 0.467 / 0.472 / 0.474 / ... -> 0.524`; the
  min is at small `d` (`d=4`), consistent with the scale-monotonicity (correlation strongest at small
  scale).

Note the exact conditional min (`0.467`) is slightly BELOW S158's reliable-search min (`0.473`) -- the
exact method catches binding configs the random structured search missed, and it is aliasing-immune.
Still `>= bar` with margin `+0.136`.

## The L=9 closure

> **L=9 stratum = [exact conditional finite check of the binding region: min `D3 = 0.467 >= bar`,
> no aliasing] + [box bound: large scale => the 2 points decorrelate => `D3 -> D3_inf^{(9)} = 0.524 >=
> bar`].** With `+0.14` margin throughout (vs the `L=10`-to-bar razor), the closure is comfortable.

This finishes the last stratum with the fleet's exact method: `L=10` (kps-S89, explicit box bound +
conditional check) and `L<=9` (here, same conditional method; the binding `L=9` confirmed exact) both
clear, so with kps-S88's exhaustive `<=30` every primitive 11-set has `D3 >= bar` -- the k=11 tail.

## Honest scope

- The conditional method is EXACT in the point direction (no aliasing) and matches exact Farey to
  `1e-5`; the `a`-integral is `Na=250`-Riemann (converges `1/Na`, `D3` to `~1e-3` -- well inside the
  `+0.14` margin). The binding finite check restricts to the near-point structure (where the min sits,
  S158); the box bound (large scale) rests on the rank-2 decorrelation (S159) + the `+0.19` limit
  margin, at kps-S89's a-priori-constant rigor level.
- NET: the L=9 correlated remainder is closed by the exact conditional method; binding min `0.467 >=
  bar` (exact), large scale `-> 0.524`. Files: `lrc14_L9_conditional_closure_opus_S160.py` (+out);
  exact-verified via `lrc14_d3_exact_verify_klein_S184.D3`. Builds on kps-S89 (conditional method,
  L=10), S159 (rank-2 mechanism), S158 (L<=9 structure).
