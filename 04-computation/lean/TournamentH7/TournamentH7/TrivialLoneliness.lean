/- TrivialLoneliness.lean — mac-mini-2026-07-16-S129.
   Rung four: the first END-TO-END formal loneliness existence theorem.
   From the built chain (fragmentation → killer_budget → good_floor):
   if the total budget is < 1 on the unit interval, a lonely time exists —
   a point t with lam ≤ |w·t − a| for every speed w ∈ W and every integer a.
   Corollary `trivial_LRC`: every finite speed set W with 6·lam·|W| < 1 and
   lam ≤ 1/2 has a lonely time — the Lonely Runner statement at the trivial
   constant. The remaining formalization program is exactly: improve this
   constant to 1/14 for |W| = 13 along the certified canon chain. -/
import TournamentH7.KillerBudget

open MeasureTheory Set

namespace LRC14

/-- Avoiding the bad set means every integer stays `lam` away from `w*t`. -/
theorem dist_int_ge (w : ℕ) (hw : 0 < w) (lam : ℝ) (hlam : 0 < lam) (t : ℝ)
    (ht : t ∉ badSet w lam) (a : ℤ) :
    lam ≤ |(w : ℝ) * t - a| := by
  have hwR : (0 : ℝ) < (w : ℝ) := by exact_mod_cast hw
  have hwne : (w : ℝ) ≠ 0 := ne_of_gt hwR
  unfold badSet at ht
  rw [mem_iUnion] at ht
  push_neg at ht
  have hta := ht a
  rw [mem_Ioo] at hta
  push_neg at hta
  rcases le_or_gt t ((a : ℝ) / w - lam / w) with hle | hgt
  · -- t on the left of the arc: a - w t ≥ lam
    have : (w : ℝ) * t ≤ (a : ℝ) - lam := by
      have := mul_le_mul_of_nonneg_left hle (le_of_lt hwR)
      have hexp : (w : ℝ) * ((a : ℝ) / w - lam / w) = (a : ℝ) - lam := by field_simp
      linarith [hexp ▸ this]
    rw [abs_sub_comm, abs_of_nonneg (by linarith)]
    linarith
  · have hge := hta hgt
    -- t on the right of the arc: w t - a ≥ lam
    have : (a : ℝ) + lam ≤ (w : ℝ) * t := by
      have := mul_le_mul_of_nonneg_left hge (le_of_lt hwR)
      have hexp : (w : ℝ) * ((a : ℝ) / w + lam / w) = (a : ℝ) + lam := by field_simp
      linarith [hexp ▸ this]
    rw [abs_of_nonneg (by linarith)]
    linarith

/-- **The end-to-end existence theorem.** If the total fragmentation budget on the
    unit interval is `< 1`, a lonely time exists. -/
theorem exists_lonely (W : Finset ℕ) (hW : ∀ w ∈ W, 0 < w) (lam : ℝ) (hlam : 0 < lam)
    (hbudget : ∑ w ∈ W, ((w : ℝ) * 1 + 2 * lam + 1) * (2 * lam / w) < 1) :
    ∃ t ∈ Icc (0 : ℝ) 1, ∀ w ∈ W, ∀ a : ℤ, lam ≤ |(w : ℝ) * t - a| := by
  classical
  have hterm : ∀ w ∈ W, 0 ≤ ((w : ℝ) * 1 + 2 * lam + 1) * (2 * lam / w) := by
    intro w hw
    have hwR : (0 : ℝ) < (w : ℝ) := by exact_mod_cast hW w hw
    positivity
  -- the budget in ENNReal
  have hsum : ∑ w ∈ W, ENNReal.ofReal (((w : ℝ) * 1 + 2 * lam + 1) * (2 * lam / w))
      = ENNReal.ofReal (∑ w ∈ W, ((w : ℝ) * 1 + 2 * lam + 1) * (2 * lam / w)) :=
    (ENNReal.ofReal_sum_of_nonneg hterm).symm
  have hfloor := good_floor W hW 0 1 lam hlam (by norm_num)
  rw [hsum] at hfloor
  have hpos : (0 : ENNReal)
      < volume (Icc (0 : ℝ) (0 + 1) \ ⋃ w ∈ W, badSet w lam) :=
    lt_of_lt_of_le
      (tsub_pos_of_lt ((ENNReal.ofReal_lt_ofReal_iff one_pos).mpr hbudget)) hfloor
  have hne : (Icc (0 : ℝ) (0 + 1) \ ⋃ w ∈ W, badSet w lam).Nonempty := by
    rcases nonempty_of_measure_ne_zero (ne_of_gt hpos) with ⟨t, ht⟩
    exact ⟨t, ht⟩
  rcases hne with ⟨t, htIcc, htBad⟩
  refine ⟨t, by simpa using htIcc, fun w hw a => ?_⟩
  have : t ∉ badSet w lam := fun hmem => htBad (mem_biUnion hw hmem)
  exact dist_int_ge w (hW w hw) lam hlam t this a

/-- Per-speed budget bound: for `w ≥ 1` and `lam ≤ 1/2`, one speed costs at most `6*lam`. -/
theorem per_term_le (w : ℕ) (hw : 0 < w) (lam : ℝ) (hlam : 0 < lam) (hhalf : lam ≤ 1/2) :
    ((w : ℝ) * 1 + 2 * lam + 1) * (2 * lam / w) ≤ 6 * lam := by
  have hwR : (1 : ℝ) ≤ (w : ℝ) := by exact_mod_cast hw
  have hwpos : (0 : ℝ) < (w : ℝ) := by linarith
  have hnum : (w : ℝ) * 1 + 2 * lam + 1 ≤ 3 * (w : ℝ) := by linarith
  have hpos : (0 : ℝ) < 2 * lam / w := by positivity
  calc ((w : ℝ) * 1 + 2 * lam + 1) * (2 * lam / w)
      ≤ 3 * (w : ℝ) * (2 * lam / w) := by
        exact mul_le_mul_of_nonneg_right hnum (le_of_lt hpos)
    _ = 6 * lam := by field_simp; ring

/-- **The Lonely Runner theorem at the trivial constant.** Every finite speed set `W`
    (positive speeds) with `lam ≤ 1/2` and `6 * lam * |W| < 1` admits a time `t` at
    which every speed is `lam`-lonely: `lam ≤ |w·t − a|` for all `w ∈ W, a ∈ ℤ`. -/
theorem trivial_LRC (W : Finset ℕ) (hW : ∀ w ∈ W, 0 < w) (lam : ℝ) (hlam : 0 < lam)
    (hhalf : lam ≤ 1/2) (hsmall : 6 * lam * W.card < 1) :
    ∃ t ∈ Icc (0 : ℝ) 1, ∀ w ∈ W, ∀ a : ℤ, lam ≤ |(w : ℝ) * t - a| := by
  classical
  refine exists_lonely W hW lam hlam (lt_of_le_of_lt ?_ hsmall)
  calc ∑ w ∈ W, ((w : ℝ) * 1 + 2 * lam + 1) * (2 * lam / w)
      ≤ ∑ _w ∈ W, 6 * lam :=
        Finset.sum_le_sum fun w hw => per_term_le w (hW w hw) lam hlam hhalf
    _ = 6 * lam * W.card := by
        rw [Finset.sum_const, nsmul_eq_mul]
        ring

end LRC14
