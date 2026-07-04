/-
  TournamentH7.LRCFarPeelDeepWell — THE FAR-PEEL AT ITS ACTUAL (LINEAR) THRESHOLD
  (kind-pasteur-2026-07-03-S39, HYP-4067).

  `far_peel_lonely_of_cite` (LRCFarPeelGood) states an EXPLICIT integer threshold
  `(1+2·ΣB)·(400·ΣB) < 3·w` — quadratic in `ΣB`.  That `V²` is an ARTIFACT of substituting
  two individually-linear bounds (piece count `≤ 1+2·ΣB`, length `≥ 1/(400·ΣB)`) into the
  sharp peel condition.  The BASE lemma `far_peel_lonely` is already sharp: it takes the
  ACTUAL piece count and length via `hbig`, whose true content is `w > p/(3·length)` — LINEAR
  in the actual component count `p`.

  The deep well `{1,…,12,182}` (the covering-min extremizer, `M = 14/183`) is the witness:
  its base `{1,…,12}` has good region with EXACTLY `p = 12` components and length
  `6617/194040`, so `hbig` reads `24/7 < 86021/16170` (TRUE), and `182` clears the true
  threshold `~118` — not the `~1.6·10⁶` the `_of_cite` corollary demands.  So the covering-min
  extremizer is a FAR-PEEL family, closed here by `far_peel_lonely` at its real threshold.

  Uses `native_decide` on the base good-region arithmetic (like the window census); depends on
  `ofReduceBool`.
-/
import TournamentH7.LRCFarPeelGood

namespace LonelyRunner.RatIntervals

/-- The deep well family `{1,2,…,12,182}` (`182 = 13·14`). -/
def deepWell : Fin 13 → ℤ := ![1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 182]

/-- **The deep well is Lonely at the FAR-PEEL's actual threshold.** Peel `182` off the base
`{1,…,12}`: the base good region has `12` components and length `6617/194040`, so the sharp
peel condition `hbig` (`24/7 < 86021/16170`) holds and `far_peel_lonely` closes it.  `182`
clears the true threshold `p/(3·length) ≈ 118`, four orders of magnitude below the `V² ≈
1.6·10⁶` of `far_peel_lonely_of_cite` — demonstrating the `V²` is a bound artifact, not the
mechanism.  (The deep well is the covering-min extremizer `M = 14/183`.) -/
theorem deepWell_far_peel_lonely : ∃ t : ℝ, Lonely 14 deepWell t := by
  refine far_peel_lonely deepWell ?_ ?_
  · decide
  · native_decide

/-- **The far-peel closes the whole deep-well family `{1,…,12, w}` for every `w ≥ 182`** — one
base computation, an INFINITE family.  The base `{1,…,12}` is fixed (`p = 12`, length
`6617/194040`); the sharp threshold is `w > p/(3·length) = 12/(3·6617/194040) ≈ 117`, so every
`w ≥ 182 ≥ 118` clears it.  The covering members are `w = 182k` (`13,14 ∣ 182`), all closed.
This is the linear-in-`p` reach of the far-peel that the `V² ≈ 1.6·10⁶` bound of
`far_peel_lonely_of_cite` hides: a fixed small base peels every large far runner at a bounded
threshold. -/
theorem base12_far_peel (w : ℤ) (hw : 182 ≤ w) :
    ∃ t : ℝ, Lonely 14 ![1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, w] t := by
  refine far_peel_lonely _ ?_ ?_
  · intro i
    fin_cases i <;> simp_all <;> omega
  · have hbase : List.ofFn (Fin.init (![1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, w] : Fin 13 → ℤ))
        = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12] := by rfl
    have hlast : (![1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, w] : Fin 13 → ℤ) (Fin.last 12) = w := by
      rfl
    rw [hbase, hlast]
    have hp : (goodRegion2 [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12] (1 / 14)).length = 12 := by
      native_decide
    have hL : length (goodRegion2 [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12] (1 / 14))
        = 6617 / 194040 := by native_decide
    rw [hp, hL]
    have hwq : (182 : ℚ) ≤ (w.toNat : ℚ) := by
      have : (182 : ℤ) ≤ (w.toNat : ℤ) := by omega
      exact_mod_cast this
    have hcoef : (0 : ℚ) ≤ (1 - 2 * (1 / 14)) * (6617 / 194040) := by norm_num
    calc ((12 : ℕ) : ℚ) * (4 * (1 / 14)) = 24 / 7 := by norm_num
      _ < (1 - 2 * (1 / 14)) * (6617 / 194040) * 182 := by norm_num
      _ ≤ (1 - 2 * (1 / 14)) * (6617 / 194040) * (w.toNat : ℚ) :=
          mul_le_mul_of_nonneg_left hwq hcoef

end LonelyRunner.RatIntervals
