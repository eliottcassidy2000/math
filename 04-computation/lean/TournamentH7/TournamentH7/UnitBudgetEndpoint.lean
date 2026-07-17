/-
  TournamentH7.UnitBudgetEndpoint

  Compactness closes the strict endpoint left by the unit-interval measure
  budget.  For every nonempty finite family of positive integer speeds, the
  elementary periodic argument therefore reaches the exact gap `1/(2|W|)`.
  In particular, thirteen speeds reach `1/26`.

  No `sorry`; no `native_decide`.
-/

import TournamentH7.UnitBudget

open MeasureTheory Set Filter Topology

namespace LRC14

noncomputable section

private def endpointLam (card n : ℕ) : ℝ :=
  (1 : ℝ) / (2 * card) * (1 - 1 / ((n : ℝ) + 2))

private def endpointGoodSet (W : Finset ℕ) (lam : ℝ) : Set ℝ :=
  Icc (0 : ℝ) 1 ∩ ⋂ w : ℕ, ⋂ (_ : w ∈ W), ⋂ a : ℤ,
    {t : ℝ | lam ≤ |(w : ℝ) * t - a|}

private theorem endpointLam_pos (card : ℕ) (hcard : 0 < card) (n : ℕ) :
    0 < endpointLam card n := by
  unfold endpointLam
  have hcardR : (0 : ℝ) < card := by exact_mod_cast hcard
  have hn : (0 : ℝ) ≤ n := by positivity
  have hden : (1 : ℝ) < (n : ℝ) + 2 := by linarith
  have hinv : 1 / ((n : ℝ) + 2) < 1 := by
    rw [div_lt_one (by linarith)]
    exact hden
  positivity

private theorem endpointLam_lt_one (card : ℕ) (hcard : 0 < card) (n : ℕ) :
    endpointLam card n < 1 := by
  unfold endpointLam
  have hcardR : (0 : ℝ) < card := by exact_mod_cast hcard
  have hden : 0 < (n : ℝ) + 2 := by positivity
  have hinv : 0 < 1 / ((n : ℝ) + 2) := by positivity
  have hfactor : 1 - 1 / ((n : ℝ) + 2) < 1 := by linarith
  have htarget : (0 : ℝ) < 1 / (2 * card) := by positivity
  have hcardOne : (1 : ℝ) ≤ card := by
    exact_mod_cast (show 1 ≤ card by omega)
  have htargetLe : (1 : ℝ) / (2 * card) ≤ 1 / 2 := by
    apply (div_le_div_iff₀ (by positivity) (by norm_num)).2
    nlinarith
  calc
    (1 : ℝ) / (2 * card) * (1 - 1 / ((n : ℝ) + 2))
        < (1 : ℝ) / (2 * card) * 1 := mul_lt_mul_of_pos_left hfactor htarget
    _ ≤ 1 / 2 := by simpa using htargetLe
    _ < 1 := by norm_num

private theorem endpointLam_budget (card : ℕ) (hcard : 0 < card) (n : ℕ) :
    2 * endpointLam card n * card < 1 := by
  unfold endpointLam
  have hcardR : (0 : ℝ) < card := by exact_mod_cast hcard
  have hcardNe : (card : ℝ) ≠ 0 := ne_of_gt hcardR
  have hden : 0 < (n : ℝ) + 2 := by positivity
  have hinv : 0 < 1 / ((n : ℝ) + 2) := by positivity
  calc
    2 * ((1 : ℝ) / (2 * card) * (1 - 1 / ((n : ℝ) + 2))) * card
        = 1 - 1 / ((n : ℝ) + 2) := by field_simp
    _ < 1 := by linarith

private theorem endpointLam_mono (card : ℕ) (hcard : 0 < card) (n : ℕ) :
    endpointLam card n ≤ endpointLam card (n + 1) := by
  unfold endpointLam
  have hcardR : (0 : ℝ) < card := by exact_mod_cast hcard
  have hn : (0 : ℝ) ≤ n := by positivity
  have hden0 : 0 < (n : ℝ) + 2 := by positivity
  have hfrac : 1 / (((n + 1 : ℕ) : ℝ) + 2) ≤
      1 / ((n : ℝ) + 2) := by
    apply one_div_le_one_div_of_le hden0
    push_cast
    linarith
  apply mul_le_mul_of_nonneg_left _ (by positivity)
  apply sub_le_sub_left
  exact hfrac

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

private theorem endpointGoodSet_antitone (W : Finset ℕ) (hcard : 0 < W.card)
    (n : ℕ) :
    endpointGoodSet W (endpointLam W.card (n + 1)) ⊆
      endpointGoodSet W (endpointLam W.card n) := by
  intro t ht
  simp only [endpointGoodSet, mem_inter_iff, mem_iInter, mem_setOf_eq] at ht ⊢
  exact ⟨ht.1, fun w hw a =>
    (endpointLam_mono W.card hcard n).trans (ht.2 w hw a)⟩

private theorem endpointGoodSet_nonempty
    (W : Finset ℕ) (hW : ∀ w ∈ W, 0 < w) (hcard : 0 < W.card)
    (n : ℕ) : (endpointGoodSet W (endpointLam W.card n)).Nonempty := by
  obtain ⟨t, ht, hdist⟩ := exists_lonely_sharp W hW (endpointLam W.card n)
    (endpointLam_pos W.card hcard n) (endpointLam_lt_one W.card hcard n)
    (endpointLam_budget W.card hcard n)
  refine ⟨t, ?_⟩
  simp only [endpointGoodSet, mem_inter_iff, mem_iInter, mem_setOf_eq]
  exact ⟨ht, hdist⟩

private theorem endpointLam_tendsto (card : ℕ) :
    Tendsto (endpointLam card) atTop (nhds ((1 : ℝ) / (2 * card))) := by
  unfold endpointLam
  have hzero : Tendsto (fun n : ℕ => (1 : ℝ) / ((n : ℝ) + 2))
      atTop (nhds 0) := by
    convert ((tendsto_add_atTop_iff_nat 1).2
      (tendsto_one_div_add_atTop_nhds_zero_nat (𝕜 := ℝ))) using 1 <;>
      norm_num [Nat.cast_add, Nat.cast_one] <;> ring
  convert tendsto_const_nhds.mul (tendsto_const_nhds.sub hzero) using 1 <;> norm_num

/-- **Sharp elementary endpoint for every finite family.**  The unit-period
bad-set budget gives every subcritical gap; compactness of `[0,1]` retains a
witness at the limiting gap `1/(2|W|)`. -/
theorem exists_lonely_unit_endpoint (W : Finset ℕ) (hW : ∀ w ∈ W, 0 < w)
    (hcard : 0 < W.card) :
    ∃ t ∈ Icc (0 : ℝ) 1, ∀ w ∈ W, ∀ a : ℤ,
      (1 : ℝ) / (2 * W.card) ≤ |(w : ℝ) * t - a| := by
  have hnonempty :
      (⋂ n : ℕ, endpointGoodSet W (endpointLam W.card n)).Nonempty :=
    IsCompact.nonempty_iInter_of_sequence_nonempty_isCompact_isClosed
      (fun n => endpointGoodSet W (endpointLam W.card n))
      (endpointGoodSet_antitone W hcard)
      (endpointGoodSet_nonempty W hW hcard)
      (isCompact_Icc.of_isClosed_subset
        (endpointGoodSet_closed W (endpointLam W.card 0)) inter_subset_left)
      (fun n => endpointGoodSet_closed W (endpointLam W.card n))
  obtain ⟨t, ht⟩ := hnonempty
  simp only [mem_iInter] at ht
  have ht0 := ht 0
  refine ⟨t, ht0.1, ?_⟩
  intro w hw a
  apply le_of_tendsto' (endpointLam_tendsto W.card)
  intro n
  have htn := ht n
  simp only [endpointGoodSet, mem_inter_iff, mem_iInter, mem_setOf_eq] at htn
  exact htn.2 w hw a

/-- Thirteen speeds at the exact elementary endpoint `1/26`. -/
theorem LRC13speeds_at_gap_26 :
    ∀ W : Finset ℕ, (∀ w ∈ W, 0 < w) → W.card = 13 →
      ∃ t ∈ Icc (0 : ℝ) 1, ∀ w ∈ W, ∀ a : ℤ,
        (1 : ℝ) / 26 ≤ |(w : ℝ) * t - a| := by
  intro W hW hcard
  obtain ⟨t, ht, hdist⟩ := exists_lonely_unit_endpoint W hW (by omega)
  refine ⟨t, ht, ?_⟩
  intro w hw a
  have h := hdist w hw a
  norm_num [hcard] at h ⊢
  exact h

/-- Seven speeds at gap `1/14`, the exact endpoint relevant to harmonic
quotient subfamilies after enough repetitions have been removed. -/
theorem LRC7speeds_at_gap_14 :
    ∀ W : Finset ℕ, (∀ w ∈ W, 0 < w) → W.card = 7 →
      ∃ t ∈ Icc (0 : ℝ) 1, ∀ w ∈ W, ∀ a : ℤ,
        (1 : ℝ) / 14 ≤ |(w : ℝ) * t - a| := by
  intro W hW hcard
  obtain ⟨t, ht, hdist⟩ := exists_lonely_unit_endpoint W hW (by omega)
  refine ⟨t, ht, ?_⟩
  intro w hw a
  have h := hdist w hw a
  norm_num [hcard] at h ⊢
  exact h

#print axioms exists_lonely_unit_endpoint
#print axioms LRC13speeds_at_gap_26
#print axioms LRC7speeds_at_gap_14

end
end LRC14
