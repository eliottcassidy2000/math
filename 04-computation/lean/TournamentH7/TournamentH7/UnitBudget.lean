/- UnitBudget.lean — mac-mini-2026-07-16-S131.
   Rung six: THE PERIODIC SHARPENING. On the unit interval the bad set of modulus w
   costs at most 2*lam exactly: the two boundary arcs (a = 0 and a = w) each
   contribute only half inside [0,1]. Removes the 3x slack of `per_term_le`,
   lifting the ladder's proven loneliness constant to any lam with 2*lam*k < 1 —
   for thirteen speeds, gap 1/27. -/
import TournamentH7.TrivialLoneliness

open MeasureTheory Set

namespace LRC14

/-- On `[0,1]` with `0 < lam < 1`, the live arc indices are exactly `0, …, w`. -/
theorem arcIdx_unit (w : ℕ) (lam : ℝ) (hlam : 0 < lam) (hlam1 : lam < 1) :
    arcIdx w 0 1 lam = Finset.Icc 0 (w : ℤ) := by
  unfold arcIdx
  have hceil : ⌈(w : ℝ) * 0 - lam⌉ = 0 := by
    rw [show (w : ℝ) * 0 - lam = -lam by ring]
    rw [Int.ceil_eq_zero_iff]
    exact Set.mem_Ioc.mpr ⟨by linarith, by linarith⟩
  have hfloor : ⌊(w : ℝ) * 0 + (w : ℝ) * 1 + lam⌋ = (w : ℤ) := by
    rw [show (w : ℝ) * 0 + (w : ℝ) * 1 + lam = (w : ℝ) + lam by ring]
    refine Int.floor_eq_iff.mpr ⟨by push_cast; linarith, by push_cast; linarith⟩
  rw [hceil, hfloor]

/-- **The exact unit budget.** `volume (badSet w lam ∩ [0,1]) ≤ 2*lam`. -/
theorem unit_bad_le (w : ℕ) (hw : 0 < w) (lam : ℝ) (hlam : 0 < lam) (hlam1 : lam < 1) :
    volume (badSet w lam ∩ Icc (0 : ℝ) 1) ≤ ENNReal.ofReal (2 * lam) := by
  classical
  have hwR : (0 : ℝ) < (w : ℝ) := by exact_mod_cast hw
  have hw1 : (1 : ℤ) ≤ (w : ℤ) := by exact_mod_cast hw
  set c : ℝ := lam / w with hc
  have hcpos : 0 < c := by positivity
  -- cover by the arcs indexed 0..w
  have hsub : badSet w lam ∩ Icc (0 : ℝ) 1
      ⊆ ⋃ a ∈ Finset.Icc (0 : ℤ) (w : ℤ),
          (Ioo ((a : ℝ) / w - c) ((a : ℝ) / w + c) ∩ Icc (0 : ℝ) 1) := by
    rintro t ⟨ht, htI⟩
    rcases mem_iUnion.mp ht with ⟨a, hta⟩
    rw [mem_Ioo] at hta
    rw [mem_Icc] at htI
    have hmem : a ∈ arcIdx w 0 1 lam :=
      mem_arcIdx_of_hit w hw 0 1 lam a t hta.1 hta.2 htI.1 (by linarith [htI.2])
    rw [arcIdx_unit w lam hlam hlam1] at hmem
    exact mem_biUnion hmem ⟨by rw [mem_Ioo]; exact hta, by rw [mem_Icc]; exact htI⟩
  -- per-arc bounds
  have hbound : ∀ a ∈ Finset.Icc (0 : ℤ) (w : ℤ),
      volume (Ioo ((a : ℝ) / w - c) ((a : ℝ) / w + c) ∩ Icc (0 : ℝ) 1)
        ≤ ENNReal.ofReal (if a = 0 ∨ a = (w : ℤ) then c else 2 * c) := by
    intro a _
    by_cases h0 : a = 0
    · subst h0
      rw [if_pos (Or.inl rfl)]
      have h00 : ((0 : ℤ) : ℝ) / (w : ℝ) = 0 := by
        norm_num
      have hsub0 : Ioo (((0 : ℤ) : ℝ) / w - c) (((0 : ℤ) : ℝ) / w + c) ∩ Icc (0 : ℝ) 1
          ⊆ Ico (0 : ℝ) c := by
        rintro t ⟨htO, htI⟩
        rw [mem_Ioo, h00] at htO
        rw [mem_Icc] at htI
        exact mem_Ico.mpr ⟨htI.1, by linarith [htO.2]⟩
      calc volume _ ≤ volume (Ico (0 : ℝ) c) := measure_mono hsub0
        _ = ENNReal.ofReal c := by rw [Real.volume_Ico, sub_zero]
    · by_cases hww : a = (w : ℤ)
      · subst hww
        rw [if_pos (Or.inr rfl)]
        have hdiv : (((w : ℤ) : ℝ)) / (w : ℝ) = 1 := by
          push_cast
          field_simp
        have hsubw : Ioo ((((w : ℤ) : ℝ)) / w - c) ((((w : ℤ) : ℝ)) / w + c) ∩ Icc (0 : ℝ) 1
            ⊆ Ioc (1 - c) 1 := by
          rintro t ⟨htO, htI⟩
          rw [mem_Ioo, hdiv] at htO
          rw [mem_Icc] at htI
          exact mem_Ioc.mpr ⟨by linarith [htO.1], htI.2⟩
        calc volume _ ≤ volume (Ioc (1 - c) 1) := measure_mono hsubw
          _ = ENNReal.ofReal c := by rw [Real.volume_Ioc]; congr 1; ring
      · rw [if_neg (by tauto)]
        calc volume _ ≤ volume (Ioo ((a : ℝ) / w - c) ((a : ℝ) / w + c)) :=
              measure_mono inter_subset_left
          _ = ENNReal.ofReal (2 * c) := by rw [Real.volume_Ioo]; congr 1; ring
  -- the peel: Icc 0 w = {0} ∪ {w} ∪ Icc 1 (w-1), valid for every w ≥ 1
  have hpeel : Finset.Icc (0 : ℤ) (w : ℤ)
      = insert 0 (insert (w : ℤ) (Finset.Icc 1 ((w : ℤ) - 1))) := by
    ext a
    simp only [Finset.mem_Icc, Finset.mem_insert]
    omega
  -- the real-valued sum of the bounds is exactly 2*lam
  have hrealsum : (∑ a ∈ Finset.Icc (0 : ℤ) (w : ℤ),
      (if a = 0 ∨ a = (w : ℤ) then c else 2 * c)) = 2 * lam := by
    have hnot1 : (0 : ℤ) ∉ insert (w : ℤ) (Finset.Icc 1 ((w : ℤ) - 1)) := by
      intro h
      rcases Finset.mem_insert.mp h with h | h
      · omega
      · rw [Finset.mem_Icc] at h
        omega
    have hnot2 : (w : ℤ) ∉ Finset.Icc 1 ((w : ℤ) - 1) := by
      intro h
      rw [Finset.mem_Icc] at h
      omega
    rw [hpeel, Finset.sum_insert hnot1, Finset.sum_insert hnot2]
    rw [if_pos (Or.inl rfl), if_pos (Or.inr rfl)]
    have hmid : (∑ a ∈ Finset.Icc (1 : ℤ) ((w : ℤ) - 1),
        (if a = 0 ∨ a = (w : ℤ) then c else 2 * c)) = ((w : ℝ) - 1) * (2 * c) := by
      have hconst : ∀ a ∈ Finset.Icc (1 : ℤ) ((w : ℤ) - 1),
          (if a = 0 ∨ a = (w : ℤ) then c else 2 * c) = 2 * c := by
        intro a ha
        rw [Finset.mem_Icc] at ha
        have hne : ¬(a = 0 ∨ a = (w : ℤ)) := by omega
        rw [if_neg hne]
      rw [Finset.sum_congr rfl hconst, Finset.sum_const, Int.card_Icc]
      have hnat : ((w : ℤ) - 1 + 1 - 1).toNat = w - 1 := by omega
      rw [hnat, nsmul_eq_mul]
      have hcast : ((w - 1 : ℕ) : ℝ) = (w : ℝ) - 1 := by
        have : (1 : ℕ) ≤ w := hw
        push_cast [Nat.cast_sub this]
        ring
      rw [hcast]
    rw [hmid, hc]
    field_simp
    ring
  -- assemble
  calc volume (badSet w lam ∩ Icc (0 : ℝ) 1)
      ≤ volume (⋃ a ∈ Finset.Icc (0 : ℤ) (w : ℤ),
          (Ioo ((a : ℝ) / w - c) ((a : ℝ) / w + c) ∩ Icc (0 : ℝ) 1)) :=
        measure_mono hsub
    _ ≤ ∑ a ∈ Finset.Icc (0 : ℤ) (w : ℤ),
          volume (Ioo ((a : ℝ) / w - c) ((a : ℝ) / w + c) ∩ Icc (0 : ℝ) 1) :=
        measure_biUnion_finset_le _ _
    _ ≤ ∑ a ∈ Finset.Icc (0 : ℤ) (w : ℤ),
          ENNReal.ofReal (if a = 0 ∨ a = (w : ℤ) then c else 2 * c) :=
        Finset.sum_le_sum hbound
    _ = ENNReal.ofReal (∑ a ∈ Finset.Icc (0 : ℤ) (w : ℤ),
          (if a = 0 ∨ a = (w : ℤ) then c else 2 * c)) :=
        (ENNReal.ofReal_sum_of_nonneg (by
          intro a _
          by_cases h : a = 0 ∨ a = (w : ℤ)
          · rw [if_pos h]; positivity
          · rw [if_neg h]; positivity)).symm
    _ = ENNReal.ofReal (2 * lam) := by rw [hrealsum]

end LRC14
