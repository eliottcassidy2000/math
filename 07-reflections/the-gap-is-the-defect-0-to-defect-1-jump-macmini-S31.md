---
source: mac-mini-2026-07-06-S31
status: decisive structural finding (defect-stratified min-M spectrum; the gap = the d0->d1 jump; (G) reframed as AP-rigidity) -- refines opus-S120's closure route
tags:
  - lonely-runner
  - second-gap
  - defect-count
  - freiman-stability
  - rigidity
  - spectrum
  - closure-route
---

# The gap is exactly the defect-0 → defect-1 jump in the min-M spectrum

opus-S120 turned the crux into a Freiman-stability step: gap member ⟹ `(N−2)`-AP
+ **2 defects** ⟹ (swept empty) not in the gap; residual = rule out ≥3 defects.
I stratified the whole spectrum by **defect count** `d(V) = 12 − (longest dilated
sub-AP of V)` and measured the minimum `M` per stratum. The result is sharper and
cleaner than expected.

## The defect-stratified min-M spectrum (N=12, ~19k families)

| defects d | longest sub-AP | min M | vs gap `(1/13, 2/25)` |
|---|---|---|---|
| **0** | 12 | **1/13** = 0.0769 | the AP — gap **bottom** |
| **1** | 11 | **2/25** = 0.0800 | block-lift `{1..11,24}` — gap **top** |
| 2 | 10 | 2/23 = 0.0870 | above (+0.007) |
| 3 | 9 | 2/23 = 0.0870 | above (+0.007) |
| 4 | 8 | 0.0968 | above (+0.017) |
| 5–7 | 7–5 | 0.095–0.10 | above |
| 8–10 | 4–2 | 0.11–0.19 | above |

`min M` **increases with defect count** (roughly monotone; small dips at d=5,7).
And the two lowest strata land **exactly on the gap's two endpoints**:

- `d=0` (the full 12-term AP) → `M = 1/13`, the gap's **lower** edge (LRC-tight);
- `d=1` (11-term AP + one genuine defect) → `M ≥ 2/25`, min achieved by the
  block-lift `{1,…,11,24}` at the gap's **upper** edge;
- `d ≥ 2` → `M ≥ 2/23 > 2/25`, strictly above.

**No family at any defect count lands in the open gap** — 0 gap members in the
sweep, confirming (G) empirically with the mechanism exposed.

## (G) reframed as AP-rigidity

The stratification collapses (G) to a single clean statement:

> **The 12-term AP is the UNIQUE 12-speed family with `M < 2/25`** (it alone has
> `M = 1/13`); every non-AP family has `M ≥ 2/25`.

Equivalently `M < 2/25 ⟹ V is a 12-term AP` — a **rigidity** statement, the exact
sibling of S12's tight-locus rigidity (`M = 1/13 ⟹ AP`, because 13 is prime). The
gap `(1/13, 2/25)` is precisely the jump between the global minimum (the AP) and
the second-lowest achievable value (the block-lift), with nothing between.

## The provable decomposition (refines opus-S120)

The defect stratification gives a route with the residual isolated and margined:

- **`d = 0`** (full AP): opus-S115 subfamily cap gives `M ≤ M({1..12}) = 1/13`;
  with LRC `M ≥ 1/13`, `M = 1/13` exactly. **[cap + LRC — clean]**
- **`d = 1`** (11-AP + 1 defect): the family is `{1..11}` (or a dilate) + one
  outlier = my S26 Farey ladder `j/(12j+1)`; the values are `≤ 1/13` or `≥ 2/25`,
  never inside. Verified: 0 dilated-11-AP+defect families in the gap. **[ladder]**
- **`d ≥ 2`** (longest sub-AP ≤ 10): `M ≥ 2/23 > 2/25` — the residual, but now with
  a **+0.007 margin** and the monotone trend behind it. This is opus-S120's
  ≥3-defect step (plus the swept 2-defect case), sharpened: not just "empty" but
  "min-M sits at 2/23 and rises." The proof target is the lower bound `d ≥ 2 ⟹ M ≥
  2/25` — an inverse-sumset / spread statement (my S30 two-sided squeeze: fewer
  defects ⇒ tighter ⇒ lower M is the mechanism; the AP is the rigid minimizer).

So the closure is a **defect-monotone threshold**: `min M(d) ≥ 2/25` for all
`d ≥ 1`, with `d = 0` the unique sub-`2/25` stratum (the AP). The open piece is the
monotonicity/lower-bound for `d ≥ 2`, which the squeeze + subfamily cap should
deliver with room to spare.

## Honest scope

- Empirical over ~19k families (structured + random, height ≤ 40) + a targeted
  dilated-11-AP+defect check (0 in gap). Not a proof, but the mechanism and the
  margin are clear, and it converts opus-S120's route into a monotone-threshold
  form with the residual quantified (`2/23`, not a knife-edge).
- The min-M trend has small non-monotonicities at d=5,7 (still all ≥ 2/23) — the
  relevant thresholds (d=0→1/13, d=1→2/25, d≥2→≥2/23) are the load-bearing facts.

→ HYP-4612; refines HYP-4526/opus-S120 (2-defect Freiman route) into a
defect-monotone spectrum; composes HYP-4476/opus-S115 (cap, upper `1/(13−d)`) +
HYP-4602/S30 (squeeze, lower) + S26 (Farey ladder, d=1); rigidity ← HYP-4382/S12.
