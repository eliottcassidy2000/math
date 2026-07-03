/-
  TournamentH7.LRCFarPeelCore  (opus-2026-07-03-S49, HYP-4042)

  THE FAR-ELEMENT PEEL — MEASURE CORE of the covering-far surface (`CoveringFarLonely`,
  the single remaining analytic hypothesis of the LRC(14) endgame, LRC14CoveringFarSurface).

  Removing one far runner `w` from a base region `G ⊆ [0,1)` costs, by the LANDED region
  rate lemma (`length_inter_comb_near_region`, kps LRCFarElementRate),

      length (inter G (comb w r φ))  ≤  2r·length G  +  (#comp)·(4r/w),

  so the PEELED region `G \ comb w` keeps positive length as soon as `w` clears the explicit
  threshold `(#comp)·4r / ((1−2r)·length G)`.  Via `exists_lonely_of_goodRegion_pos` (kps
  LRCGoodPipeline) that positive length is loneliness — so this is the crisp reduction of the
  single-far case of `CoveringFarLonely` to (i) a base good-region FLOOR `length G > 0` and
  (ii) a SIZE THRESHOLD on `w`.  It is the `joint_rate_core_reduced` analogue for the peel:
  the whole quantitative content is one inequality on the LANDED rate bound, sorry-free.

  Purely analytic dependency chain (LRCFarElementRate → RateLemma, LRCSevenTranslate) — it
  does NOT import the native_decide window-data packs, so it builds independently of them.
-/
import TournamentH7.LRCFarElementRate

namespace LonelyRunner
namespace RatIntervals

/-- **Far-peel positivity (measure core).**  For a base region `G ⊆ [0,1)` and a
non-wrapping `w`-comb of radius `r`, if `w` is large enough that the two-partial-wrap loss
`(#comp)·(4r/w)` is strictly below the undamped survival `(1−2r)·length G`, then the peeled
region `G \ comb w` has strictly positive length.  Direct corollary of the landed region
rate lemma. -/
theorem far_peel_length_pos {w : ℕ} (hw : 0 < w) {r φ : ℚ}
    (hr : 0 ≤ r) (hrφ : r ≤ φ) (hφ1 : φ + r ≤ 1)
    (G : Region) (hG : ∀ p ∈ G, 0 ≤ p.1 ∧ p.1 ≤ p.2 ∧ p.2 ≤ 1)
    (hbig : (G.length : ℚ) * (4 * r / w) < (1 - 2 * r) * length G) :
    0 < length G - length (inter G (comb w r φ)) := by
  have hrate := length_inter_comb_near_region hw hr hrφ hφ1 G hG
  have hupper : length (inter G (comb w r φ))
      ≤ 2 * r * length G + (G.length : ℚ) * (4 * r / w) := by
    have h := (abs_le.mp hrate).2
    linarith
  have heq : (1 - 2 * r) * length G = length G - 2 * r * length G := by ring
  rw [heq] at hbig
  linarith

/-- **Far-peel positivity, explicit size threshold** (division cleared — an integer-friendly
`w`-bound for certifying covering-far families).  If the base length is positive, the band is
strict (`2r < 1`) and non-wrapping, and `w` exceeds `(#comp)·4r / ((1−2r)·length G)` — written
without division as `(#comp)·4r < (1−2r)·length G · w` — then the peeled region is positive. -/
theorem far_peel_length_pos_of_gt {w : ℕ} (hw : 0 < w) {r φ : ℚ}
    (hr : 0 ≤ r) (hrφ : r ≤ φ) (hφ1 : φ + r ≤ 1)
    (G : Region) (hG : ∀ p ∈ G, 0 ≤ p.1 ∧ p.1 ≤ p.2 ∧ p.2 ≤ 1)
    (hwbig : (G.length : ℚ) * (4 * r) < (1 - 2 * r) * length G * w) :
    0 < length G - length (inter G (comb w r φ)) := by
  have hwQ : (0 : ℚ) < (w : ℚ) := by exact_mod_cast hw
  refine far_peel_length_pos hw hr hrφ hφ1 G hG ?_
  have h4 : (G.length : ℚ) * (4 * r / w) = (G.length : ℚ) * (4 * r) / w := by ring
  rw [h4, div_lt_iff₀ hwQ]
  linarith

end RatIntervals
end LonelyRunner
