/-
  TournamentH7.LRCC8Consecutive — THE c = 8 CONSECUTIVE THEOREM, END TO END
  (boxeph-2026-07-17-S77; LEM-044(B) fully in the kernel).

  Every consecutive block of eight runners leaves a good set of positive
  measure at the 1/14 margin: the union bound dies (8 × 1/7 > 1), but the
  seven consecutive-pair credits — each ≥ 1/49 with the k ≡ 3 (mod 7) pair
  STRICTLY above by 9/(49k(k+1)) — cross the wall through the path-Hunter
  skeleton.

  Composition of: `good_pos_of_path_credits` (LRCTreeHunter, S73) +
  `consecutive_credit_closed` (LRCPairOverlapArcs, S75) + `danger_measure_le`
  (S76) + the middle-residue pigeonhole (any seven consecutive k's contain
  one ≡ 3 mod 7).

  Kernel-pure: no `native_decide`, no `sorry`.
-/
import TournamentH7.LRCTreeHunter
import TournamentH7.LRCPairOverlapArcs

open MeasureTheory Set

namespace LonelyRunner.LRC14.Arcs

open LonelyRunner.LRC14.Hunter

noncomputable section

theorem dangerR_isOpen (v : ℕ) : IsOpen (dangerR v) := by
  have h : dangerR v = ⋃ m : ℤ, {t : ℝ | |(v : ℝ) * t - m| < 1/14} := by
    ext t
    simp [dangerR, Set.mem_iUnion, Set.mem_setOf_eq]
  rw [h]
  apply isOpen_iUnion
  intro m
  have hc : Continuous fun t : ℝ => |(v : ℝ) * t - (m : ℝ)| :=
    ((continuous_const.mul continuous_id).sub continuous_const).abs
  exact isOpen_lt hc continuous_const

/-- **The c = 8 consecutive theorem, end to end**: every consecutive block of
eight runners leaves a positive-measure good set on the unit window. -/
theorem c8_consecutive_good_pos (v : ℕ) (hv : 1 ≤ v) :
    0 < (volume.restrict (Ioo (-(1:ℝ)/2) (1/2)))
        ((⋃ i ∈ Finset.range 8, dangerR (v + i))ᶜ) := by
  haveI hprob : IsProbabilityMeasure (volume.restrict (Ioo (-(1:ℝ)/2) (1/2))) := by
    constructor
    rw [Measure.restrict_apply_univ, Real.volume_Ioo]
    norm_num
  have h17 : ENNReal.ofReal (1/7 : ℝ) = 1/7 := by
    rw [ENNReal.ofReal_div_of_pos (by norm_num), ENNReal.ofReal_one,
      ENNReal.ofReal_ofNat]
  set A : ℕ → Set ℝ := fun i => dangerR (v + i) with hA_def
  have hAmeas : ∀ i, MeasurableSet (A i) := fun i =>
    (dangerR_isOpen (v + i)).measurableSet
  have hle : ∀ i ∈ Finset.range 8,
      (volume.restrict (Ioo (-(1:ℝ)/2) (1/2))) (A i) ≤ 1/7 := by
    intro i _
    rw [Measure.restrict_apply (hAmeas i)]
    calc volume (A i ∩ Ioo (-(1:ℝ)/2) (1/2))
        ≤ ENNReal.ofReal (1/7) := danger_measure_le (v + i) (by omega)
      _ = 1/7 := h17
  set c : ℕ → ℝ := fun i =>
    1 / 49 + (((v + i - 1) % 7 : ℕ) : ℝ) * (6 - (((v + i - 1) % 7 : ℕ) : ℝ))
      / (49 * ((v + i - 1 : ℕ) : ℝ) * (((v + i - 1 : ℕ) : ℝ) + 1)) with hc_def
  have hpair : ∀ i ∈ Finset.Ico 1 8, ENNReal.ofReal (c i)
      ≤ (volume.restrict (Ioo (-(1:ℝ)/2) (1/2))) (A i ∩ A (i - 1)) := by
    intro i hi
    rw [Finset.mem_Ico] at hi
    rw [Measure.restrict_apply ((hAmeas i).inter (hAmeas (i - 1)))]
    have hk1 : 1 ≤ v + i - 1 := by omega
    have e1 : A (i - 1) = dangerR (v + i - 1) := by
      simp only [hA_def]
      congr 1
      omega
    have e2 : A i = dangerR ((v + i - 1) + 1) := by
      simp only [hA_def]
      congr 1
      omega
    rw [e1, e2, Set.inter_comm (dangerR ((v + i - 1) + 1)) (dangerR (v + i - 1))]
    exact consecutive_credit_closed (v + i - 1) hk1
  have hc_ge : ∀ i ∈ Finset.Ico 1 8, 1/49 ≤ c i := by
    intro i hi
    rw [Finset.mem_Ico] at hi
    have hk1 : 1 ≤ v + i - 1 := by omega
    have hkR : (0:ℝ) < ((v + i - 1 : ℕ) : ℝ) := by exact_mod_cast hk1
    have hr6 : ((v + i - 1) % 7 : ℕ) ≤ 6 := by omega
    have hr6R : (((v + i - 1) % 7 : ℕ) : ℝ) ≤ 6 := by exact_mod_cast hr6
    have hnn : 0 ≤ (((v + i - 1) % 7 : ℕ) : ℝ) * (6 - (((v + i - 1) % 7 : ℕ) : ℝ))
        / (49 * ((v + i - 1 : ℕ) : ℝ) * (((v + i - 1 : ℕ) : ℝ) + 1)) := by
      apply div_nonneg
      · exact mul_nonneg (Nat.cast_nonneg _) (by linarith)
      · positivity
    simp only [hc_def]
    linarith
  set istar : ℕ := (3 + 7 - v % 7) % 7 + 1 with histar
  have histar_mem : istar ∈ Finset.Ico 1 8 := by
    rw [Finset.mem_Ico]
    omega
  have histar_r : (v + istar - 1) % 7 = 3 := by omega
  have hstrict : 1/49 < c istar := by
    have hk1 : 1 ≤ v + istar - 1 := by omega
    have hkR : (0:ℝ) < ((v + istar - 1 : ℕ) : ℝ) := by exact_mod_cast hk1
    have hD : (0:ℝ) < 49 * ((v + istar - 1 : ℕ) : ℝ) * (((v + istar - 1 : ℕ) : ℝ) + 1) := by
      positivity
    have hex : (0:ℝ) < (((v + istar - 1) % 7 : ℕ) : ℝ) * (6 - (((v + istar - 1) % 7 : ℕ) : ℝ))
        / (49 * ((v + istar - 1 : ℕ) : ℝ) * (((v + istar - 1 : ℕ) : ℝ) + 1)) := by
      rw [histar_r]
      have h9 : ((3:ℕ):ℝ) * (6 - ((3:ℕ):ℝ)) = 9 := by norm_num
      rw [h9]
      exact div_pos (by norm_num) hD
    simp only [hc_def]
    linarith
  have hsum_lt : (1:ℝ)/7 < ∑ i ∈ Finset.Ico 1 8, c i := by
    have h0 : ∑ _i ∈ Finset.Ico 1 8, (1:ℝ)/49 < ∑ i ∈ Finset.Ico 1 8, c i :=
      Finset.sum_lt_sum hc_ge ⟨istar, histar_mem, hstrict⟩
    have h1 : ∑ _i ∈ Finset.Ico 1 8, (1:ℝ)/49 = 7 * (1/49) := by
      rw [Finset.sum_const, Nat.card_Ico, nsmul_eq_mul]
      norm_num
    rw [h1] at h0
    linarith
  have hcred : ((8 : ℕ) : ENNReal) / 7
      < 1 + ∑ i ∈ Finset.Ico 1 8,
          (volume.restrict (Ioo (-(1:ℝ)/2) (1/2))) (A i ∩ A (i - 1)) := by
    have hcnn : ∀ i ∈ Finset.Ico 1 8, 0 ≤ c i := fun i hi =>
      le_trans (by norm_num) (hc_ge i hi)
    have hsum_ge : ENNReal.ofReal (∑ i ∈ Finset.Ico 1 8, c i)
        ≤ ∑ i ∈ Finset.Ico 1 8,
            (volume.restrict (Ioo (-(1:ℝ)/2) (1/2))) (A i ∩ A (i - 1)) := by
      rw [ENNReal.ofReal_sum_of_nonneg hcnn]
      exact Finset.sum_le_sum hpair
    have hstep : (1 : ENNReal)/7 < ENNReal.ofReal (∑ i ∈ Finset.Ico 1 8, c i) := by
      rw [← h17]
      exact (ENNReal.ofReal_lt_ofReal_iff
        (lt_trans (by norm_num) hsum_lt)).mpr hsum_lt
    have h87 : ((8 : ℕ) : ENNReal) / 7 = 1 + 1/7 := by
      rw [show ((8 : ℕ) : ENNReal) = 7 + 1 by norm_num, ENNReal.add_div,
        ENNReal.div_self (by norm_num) (by norm_num)]
    calc ((8 : ℕ) : ENNReal) / 7 = 1 + 1/7 := h87
      _ < 1 + ENNReal.ofReal (∑ i ∈ Finset.Ico 1 8, c i) :=
          ENNReal.add_lt_add_left (by norm_num) hstep
      _ ≤ 1 + ∑ i ∈ Finset.Ico 1 8,
            (volume.restrict (Ioo (-(1:ℝ)/2) (1/2))) (A i ∩ A (i - 1)) :=
          add_le_add (le_refl 1) hsum_ge
  exact good_pos_of_path_credits (volume.restrict (Ioo (-(1:ℝ)/2) (1/2)))
    A hAmeas 8 hle hcred

#print axioms c8_consecutive_good_pos

end

end LonelyRunner.LRC14.Arcs
