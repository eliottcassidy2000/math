/-
  TournamentH7.UnitBudgetEndpoint

  Compactness closes the strict endpoint left by the unit-interval measure
  budget.  For thirteen positive integer speeds, the elementary periodic
  argument therefore reaches the exact gap `1/26`, not merely a convenient
  rational gap below it.

  No `sorry`; no `native_decide`.
-/

import TournamentH7.UnitBudget

open MeasureTheory Set Filter Topology

namespace LRC14

noncomputable section

private def endpointLam (n : ℕ) : ℝ :=
  (1 : ℝ) / 26 * (1 - 1 / ((n : ℝ) + 2))

private def endpointGoodSet (W : Finset ℕ) (lam : ℝ) : Set ℝ :=
  Icc (0 : ℝ) 1 ∩ ⋂ w : ℕ, ⋂ (_ : w ∈ W), ⋂ a : ℤ,
    {t : ℝ | lam ≤ |(w : ℝ) * t - a|}

private theorem endpointLam_pos (n : ℕ) : 0 < endpointLam n := by
  unfold endpointLam
  have hn : (0 : ℝ) ≤ n := by positivity
  have hden : (1 : ℝ) < (n : ℝ) + 2 := by linarith
  have hinv : 1 / ((n : ℝ) + 2) < 1 := by
    rw [div_lt_one (by linarith)]
    exact hden
  positivity

private theorem endpointLam_lt_one (n : ℕ) : endpointLam n < 1 := by
  unfold endpointLam
  have hden : 0 < (n : ℝ) + 2 := by positivity
  have hinv : 0 < 1 / ((n : ℝ) + 2) := by positivity
  nlinarith

private theorem endpointLam_budget (n : ℕ) :
    2 * endpointLam n * 13 < 1 := by
  unfold endpointLam
  have hden : 0 < (n : ℝ) + 2 := by positivity
  have hinv : 0 < 1 / ((n : ℝ) + 2) := by positivity
  nlinarith

private theorem endpointLam_mono (n : ℕ) :
    endpointLam n ≤ endpointLam (n + 1) := by
  unfold endpointLam
  have hn : (0 : ℝ) ≤ n := by positivity
  have hden0 : 0 < (n : ℝ) + 2 := by positivity
  have hden1 : 0 < (n : ℝ) + 3 := by positivity
  have hfrac : 1 / ((n : ℝ) + 3) ≤ 1 / ((n : ℝ) + 2) := by
    exact one_div_le_one_div_of_le hden0 (by linarith)
  have hscaled := mul_le_mul_of_nonneg_left
    (sub_le_sub_left hfrac 1) (by norm_num : (0 : ℝ) ≤ 1 / 26)
  convert hscaled using 1 <;> push_cast <;> ring

private theorem endpointGoodSet_closed (W : Finset ℕ) (lam : ℝ) :
    IsClosed (endpointGoodSet W lam) := by
  unfold endpointGoodSet
  apply isClosed_Icc.inter
  apply isClosed_iInter
  intro w
  apply isClosed_iInter
  intro hw
  apply isClosed_iInter
  intro a
  exact isClosed_le continuous_const
    (((continuous_const.mul continuous_id).sub continuous_const).abs)

private theorem endpointGoodSet_antitone (W : Finset ℕ) (n : ℕ) :
    endpointGoodSet W (endpointLam (n + 1)) ⊆
      endpointGoodSet W (endpointLam n) := by
  intro t ht
  simp only [endpointGoodSet, mem_inter_iff, mem_iInter, mem_setOf_eq] at ht ⊢
  exact ⟨ht.1, fun w hw a => (endpointLam_mono n).trans (ht.2 w hw a)⟩

private theorem endpointGoodSet_nonempty
    (W : Finset ℕ) (hW : ∀ w ∈ W, 0 < w) (hcard : W.card = 13)
    (n : ℕ) : (endpointGoodSet W (endpointLam n)).Nonempty := by
  obtain ⟨t, ht, hdist⟩ := exists_lonely_sharp W hW (endpointLam n)
    (endpointLam_pos n) (endpointLam_lt_one n) (by simpa [hcard] using endpointLam_budget n)
  refine ⟨t, ?_⟩
  simp only [endpointGoodSet, mem_inter_iff, mem_iInter, mem_setOf_eq]
  exact ⟨ht, hdist⟩

private theorem endpointLam_tendsto :
    Tendsto endpointLam atTop (nhds ((1 : ℝ) / 26)) := by
  unfold endpointLam
  have hzero : Tendsto (fun n : ℕ => (1 : ℝ) / ((n : ℝ) + 2))
      atTop (nhds 0) := by
    convert ((tendsto_add_atTop_iff_nat 1).2
      (tendsto_one_div_add_atTop_nhds_zero_nat (𝕜 := ℝ))) using 1 <;>
      norm_num [Nat.cast_add, Nat.cast_one] <;> ring
  convert tendsto_const_nhds.mul (tendsto_const_nhds.sub hzero) using 1 <;> norm_num

/-- **Sharp elementary endpoint for thirteen speeds.**  The unit-period bad-set
budget gives every subcritical gap; compactness of `[0,1]` retains a witness at
the limiting gap `1/26`. -/
theorem LRC13speeds_at_gap_26 :
    ∀ W : Finset ℕ, (∀ w ∈ W, 0 < w) → W.card = 13 →
      ∃ t ∈ Icc (0 : ℝ) 1, ∀ w ∈ W, ∀ a : ℤ,
        (1 : ℝ) / 26 ≤ |(w : ℝ) * t - a| := by
  intro W hW hcard
  have hnonempty : (⋂ n : ℕ, endpointGoodSet W (endpointLam n)).Nonempty :=
    IsCompact.nonempty_iInter_of_sequence_nonempty_isCompact_isClosed
      (fun n => endpointGoodSet W (endpointLam n))
      (endpointGoodSet_antitone W)
      (endpointGoodSet_nonempty W hW hcard)
      (isCompact_Icc.of_isClosed_subset
        (endpointGoodSet_closed W (endpointLam 0)) inter_subset_left)
      (fun n => endpointGoodSet_closed W (endpointLam n))
  obtain ⟨t, ht⟩ := hnonempty
  simp only [mem_iInter] at ht
  have ht0 := ht 0
  refine ⟨t, ht0.1, ?_⟩
  intro w hw a
  apply le_of_tendsto' endpointLam_tendsto
  intro n
  have htn := ht n
  simp only [endpointGoodSet, mem_inter_iff, mem_iInter, mem_setOf_eq] at htn
  exact htn.2 w hw a

#print axioms LRC13speeds_at_gap_26

end
end LRC14
