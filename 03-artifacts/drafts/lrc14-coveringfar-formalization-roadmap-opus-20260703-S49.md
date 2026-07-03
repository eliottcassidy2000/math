# Formalizing `CoveringFarLonely` — the single open Lean lemma of the LRC(14) endgame

**Author:** opus-2026-07-03-S49 (HYP-4042).
**State of the endgame (klein-S107/S112, kps-S13/S17, mac-mini-S17):** the top theorem
`lrc14_of_covering_far_22 : LRCUpTo13 → CoveringFarLonely 22 → LRC14Statement`
(`LRCWindowData22.lean`) is sorry-free glue. Two of its three inputs are discharged:
- **`LRCUpTo13`** — the LRC(≤13) citation node (owner-sanctioned external hypothesis, by design).
- **the window census** — `hwindow22_closed` closes every family with all `|vᵢ| ≤ 22` by two
  `native_decide` (kernel-gated rows + completeness over `C(22,13)=497420`).

The **only remaining mathematics** is the third input:

```lean
def CoveringFarLonely (W : ℤ) : Prop :=
  ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → CoveringFamily v → 0 < farCountW W v → ∃ t, Lonely 14 v t
```

"every covering family with an entry beyond the window is lonely." Certified so far: an infinite
AP-tail slice (`coveringFar_blockFamily`) and one deep well (`coveringFar_deepWell`). This note gives
the ordered lemma chain to close the general statement, grounded in the landed engines, and records the
one step this session added.

---

## The peel route (single far element), in landed vocabulary

Fix a covering family `v` with a far element. Let `w = vᵢ₀` be one far speed, `B = v ∖ {w}` the other
12 speeds. Everything runs in the rational-`Region` vocabulary (`length`, `inter`, `comb`, `goodRegion`,
`RatIntervals`), whose bridge to loneliness is **landed**:

- `exists_lonely_of_goodRegion_pos` (kps `LRCGoodPipeline`): `0 < length (goodRegion (ofFn v) (1/14)) → ∃ t, Lonely 14 v t`.

The good region of `v` is that of `B` with `w`'s comb subtracted, so:

```
length (goodRegion v)  =  length (goodRegion B)  −  length (goodRegion B  ∩  comb w) .
```

The chain:

| # | Step | Landed? | Where |
|---|------|---------|-------|
| 1 | **Base floor** `0 < length (goodRegion B (1/14))` for any 12 nonzero speeds `B` | **[PROOF]** — the one real remaining lemma | from `LRCUpTo13`: `B` lonely at `1/13` ⟹ an interval of `1/14`-safe times of length `≥ 1/(91·max B)` |
| 2 | **Union = subtract-comb** `length (goodRegion v) = length (goodRegion B) − length (goodRegion B ∩ comb w)` | **[DATA/DECIDE]** diffF/wrap length bookkeeping | `LRCRegionDiff` (`diffF`, `mem_diffF`, `length_diffF`), `wrap` |
| 3 | **Far-peel positivity** `w` large ⟹ `0 < length (goodRegion B) − length (goodRegion B ∩ comb w)` | **[DONE, this session]** `far_peel_length_pos` | `LRCFarPeelCore.lean`, on the landed `length_inter_comb_near_region` |
| 4 | **Positive length ⟹ lonely** | **[DONE]** | `exists_lonely_of_goodRegion_pos` |
| 5 | **Finite window** for `w ≤ threshold(B)` (below where step 3 fires) | **[DATA]** extend the census cut, or dispatch | `LRCWindowData*` band generators |
| 6 | **Multi-far**: nothing extra — step 1 handles *any* 12-speed `B` (far or not) as a citation, so the peel is **single-step**, not an induction | — | — |

**Step 3 is `far_peel_length_pos` (landed this session, kernel-pure):**

```lean
theorem far_peel_length_pos {w : ℕ} (hw : 0 < w) {r φ : ℚ}
    (hr : 0 ≤ r) (hrφ : r ≤ φ) (hφ1 : φ + r ≤ 1)
    (G : Region) (hG : ∀ p ∈ G, 0 ≤ p.1 ∧ p.1 ≤ p.2 ∧ p.2 ≤ 1)
    (hbig : (G.length : ℚ) * (4 * r / w) < (1 - 2 * r) * length G) :
    0 < length G - length (inter G (comb w r φ))
```

with the division-cleared, integer-friendly companion `far_peel_length_pos_of_gt` (hypothesis
`(#comp)·4r < (1−2r)·length G · w`). `#print axioms` = `[propext, Classical.choice, Quot.sound]` — no
`sorryAx`, no `ofReduceBool`. It is the `joint_rate_core_reduced` analogue for the peel: the entire
quantitative content is one inequality on the landed region rate bound.

---

## What remains, precisely

1. **[PROOF — the one genuine lemma] Step 1, the base good-region floor.** Upgrade the LRC(≤13) *point*
   (`∃ t, min_B ‖vᵢt‖ ≥ 1/13`) to a *positive length* `length (goodRegion B (1/14)) > 0`, with an explicit
   floor `μ₀(B) ≥ 1/(91·max B)` (a 1/14-safe interval of half-width `(1/13−1/14)/max B` around the lonely
   point). Elementary (continuity/margin), but not yet written. This floor feeds the `w`-threshold in step 3.
2. **[DATA/DECIDE] Step 2, the length bookkeeping.** `length (goodRegion (B ++ [w])) = length (goodRegion B) −
   length (goodRegion B ∩ wrap (comb w))`, from `diffF`/`wrap`; plus the wrapped-comb form of the rate lemma
   (the landed `length_inter_comb_near_region` is stated for the *non-wrapping* band `r ≤ φ, φ+r ≤ 1`; the
   phase-0 comb wraps at the seam and needs the `wrap`-region version — a one-lemma extension).
3. **[DATA] Step 5, the finite window.** For `w` between the census cut (22) and the step-3 threshold,
   extend the window packs (the generators exist) or route through the dispatch.

With these three, `CoveringFarLonely 22` closes and `theorem lrc14 : LRC14Statement` follows unconditionally
(modulo the `LRCUpTo13` citation, by design). **The sharp F3 / renormalization program (klein-S106 F3-sharp,
opus-S48 (R2), the `joint_rate_core` engine) is what shrinks the step-5 threshold from ~`max(B)²` to `~10³`;
the crude single-step peel above is already valid, just with a wider finite window.**

---

## Build/environment note (for the fleet)

On this Windows box (Lean v4.30.0) `lake build` **segfaults** (`0xC0000005`, exit `3221225477`) inside
`LRCWindowData.lean` — the `native_decide` over the 6084/31471-row packs — matching the S47b "corrupt-olean /
crash forensics." This is machine-specific (klein built the full 8597-job corpus green elsewhere); it blocks
*data-file* builds here but **not analytic files**: `LRCFarPeelCore` builds clean and kernel-pure
(`lake build TournamentH7.LRCFarPeelCore`, 2948 jobs, ✔). Analytic endgame work (steps 1–3 above) is therefore
verifiable on this machine; the data/census steps must be built on a machine that survives the big
`native_decide`.

---

## Status

- **Landed (verified, kernel-pure):** `far_peel_length_pos` + `far_peel_length_pos_of_gt` (`LRCFarPeelCore.lean`,
  registered in the root) — step 3, the measure core of the far-element peel.
- **Route (crisp, named):** `CoveringFarLonely 22` = steps 1–5; the *only* [PROOF] item is step 1 (the base
  good-region floor from LRC(≤13)); steps 2 and 5 are [DATA/DECIDE] on landed machinery.
- Related: HYP-4042 (this); HYP-4017/klein-S112 (`CoveringFarLonely` surface); HYP-3874/mac-mini-S17
  (`joint_rate_core`); HYP-3968/kps-S11 (`rate_core`, `length_inter_comb_near_region`); HYP-3970/kps-S13
  (`goodRegion`, `exists_lonely_of_goodRegion_pos`); HYP-4011/klein-S106 (F3-sharp); HYP-4013/opus-S48 ((R2)).
