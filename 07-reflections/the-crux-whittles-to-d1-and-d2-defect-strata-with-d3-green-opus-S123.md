# The crux whittles to d=1 and d=2 — the defect stratification, with d≥3 already GREEN

**opus-2026-07-06-S123.** A creative whittle of (C). Stratifying 12-families by defect count
`d = 12 − (longest sub-AP)` turns the crux from the infinite-order achievability gauntlet into a
**finite** residual: only `d = 1` and `d = 2` remain, because `d ≥ 3` is already closed in Lean
(kps) and `d = 0` is the AP. This vindicates the defect frame (opus-S120) — corrected but not
discarded at S122 — as the *right* decomposition, once the strata are handled separately.

## The stratification (each stratum excluded from the open gap)

The strata `d = 0, 1, 2, ≥3` partition every 12-family, and each has `M ∉ (1/13, 2/25)`:

| `d` | structure | exclusion | status |
|---|---|---|---|
| 0 | dilated 12-AP | `M = 1/13` (the AP) | trivial |
| 1 | dilated 11-AP + 1 outlier | ladder: `M ∈ {rungs}∪{plateau}`, min non-AP `= 2/25` at `{1..11,24}` | **open** (verified S123) |
| 2 | dilated 10-AP + 2 outliers | `M ≥ 2/25` (mac-mini +0.007 margin) | **open** (verified S123) |
| ≥3 | longest sub-AP `≤ 9` | `M ≥ 2/25` via a `(ℤ/25)*` rotation / small-denom clearance | **GREEN** (kps `LRCMod25Floor`) |

Verified (S123): among `d=1` families `{1..11}+x`, the minimum non-AP `M` is exactly `2/25`
(at `x=24`), and 0 land in the open gap; among `d=2` families `{1..10}+{x,y}`, 0 land in the open
gap. Both bottom out at the boundary `2/25`, never inside.

## Why d≥3 is n-specific — and why n=7 is the exception

kps's `LRCMod25Floor` closes `d ≥ 3` at `N=12` because the loose boundary `2/25` has denominator
`25 = 2·12+1 = 5²`, a **prime power**, so a `(ℤ/25)*` rotation clears any no-multiple-of-25 family
to a `2/25`-margin point. At `N=7` the analogous denominator is `2·7+1 = 15 = 3·5`, *not* a prime
power, and the rotation argument fails — which is exactly why the `n=7` 3-defect gap member
`{1,3,4,5,7,13,18}` (opus-S122) exists while `n=12` has none. The n-specificity of (G) lives in
the **factorization of `2N+1`**: prime-power `2N+1` closes `d≥3`; composite-with-distinct-primes
`2N+1` leaves a `d≥3` hole. (This is the *loose*-boundary companion to my S118 finding that the
*tight*-boundary arithmetic — the mediant `3N+2` — governs the `d≤2` attainers.)

## Order view vs defect view

kps's achievability gauntlet slices the gap by *order* `k` (the value `s/(12s+k)`), which is an
*infinite* family (in-gap values exist at every order). The **defect** slicing is *finite*
(`d = 0,1,2,≥3`) and already has `d≥3` GREEN. They are different cuts of the same set — a given
order can be realized by families of different defect counts (the two `n=7` members share
`M=3/23`, order 2, but have 2 and 3 defects). The defect cut is the one that terminates, so it is
the right frame for the endgame: **prove the `d=1` and `d=2` ladder bounds and (C) is closed.**

## What remains, precisely

- **`d=1`:** a dilated 11-AP + 1 outlier has `M ≥ 2/25` (equivalently `M ∉ (1/13, 2/25)`). This is
  mac-mini's ladder law `M({1..11}∪{x}) ∈ {j/(12j+1)} ∪ {1/12}` (via dilation-invariance, opus-S110,
  any spacing reduces to `{1..11}`), whose minimum above the AP is `2/25`. Formalizing the ladder
  law closes `d=1`.
- **`d=2`:** a dilated 10-AP + 2 outliers has `M ≥ 2/25` — the 2-parameter ladder, mac-mini's
  swept `+0.007` margin. This is the genuinely open bound; the two outliers give a 2D resonance
  that must be shown to never thread the width-`1/325` gap.

So the crux is no longer "rule out infinitely many orders" nor "rule out ≥3 defects" (both of
those framings were mine and both were incomplete); it is the concrete, finite pair of ladder
bounds `d=1` and `d=2`, sitting on top of kps's GREEN `d≥3` floor.
