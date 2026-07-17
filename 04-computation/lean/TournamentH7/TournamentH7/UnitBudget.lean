/- UnitBudget.lean — mac-mini-2026-07-16-S130.
   Rung six: THE PERIODIC SHARPENING. On the unit interval the bad set of modulus w
   costs EXACTLY 2*lam, not (w + 2*lam + 1)(2*lam/w): the two boundary arcs (a = 0
   and a = w) each contribute only half inside [0,1]. This removes the 3x slack of
   `per_term_le` and lifts the ladder's proven loneliness constant from 1/(6k+) to
   any lam with 2*lam*k < 1 — for thirteen speeds, gap 1/27. -/
import TournamentH7.TrivialLoneliness

open MeasureTheory Set

namespace LRC14

/-- On `[0,1]` with `0 < lam < 1`, the live arc indices are exactly `0, …, w`. -/
theorem arcIdx_unit (w : ℕ) (hw : 0 < w) (lam : ℝ) (hlam : 0 < lam) (hlam1 : lam < 1) :
    arcIdx w 0 1 lam = Finset.Icc 0 (w : ℤ) := by
  unfold arcIdx
  congr 1
  · -- ⌈w*0 - lam⌉ = 0
    rw [show (w : ℝ) * 0 - lam = -lam by ring]
    rw [Int.ceil_eq_iff (by norm_num : (0:ℤ) ≠ 0 ∨ True) |>.symm.symm]
    · rw [Int.ceil_eq_zero_iff]
      constructor <;> [skip; skip] <;> first
        | (rw [mem_Ioc]; constructor <;> linarith)
        | skip
  · -- ⌊w*0 + w*1 + lam⌋ = w
    rw [show (w : ℝ) * 0 + (w : ℝ) * 1 + lam = (w : ℝ) + lam by ring]
    rw [Int.floor_eq_iff (by positivity)]
    constructor
    · push_cast; linarith
    · push_cast; linarith

/-- **The exact unit budget.** `volume (badSet w lam ∩ [0,1]) ≤ 2*lam`. -/
theorem unit_bad_le (w : ℕ) (hw : 0 < w) (lam : ℝ) (hlam : 0 < lam) (hlam1 : lam < 1) :
    volume (badSet w lam ∩ Icc (0 : ℝ) 1) ≤ ENNReal.ofReal (2 * lam) := by
  classical
  have hwR : (0 : ℝ) < (w : ℝ) := by exact_mod_cast hw
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
    rw [arcIdx_unit w hw lam hlam hlam1] at hmem
    exact mem_biUnion hmem ⟨by rw [mem_Ioo]; exact hta, by rw [mem_Icc]; exact htI⟩
  -- per-arc bounds: boundary arcs cost c, middle arcs cost 2c
  have hbound : ∀ a ∈ Finset.Icc (0 : ℤ) (w : ℤ),
      volume (Ioo ((a : ℝ) / w - c) ((a : ℝ) / w + c) ∩ Icc (0 : ℝ) 1)
        ≤ ENNReal.ofReal (if a = 0 ∨ a = (w : ℤ) then c else 2 * c) := by
    intro a ha
    by_cases h0 : a = 0
    · subst h0
      simp only [if_pos (Or.inl rfl)]
      have : Ioo ((0 : ℝ) / w - c) ((0 : ℝ) / w + c) ∩ Icc (0 : ℝ) 1 ⊆ Ico 0 c := by
        rintro t ⟨htO, htI⟩
        rw [mem_Ioo] at htO
        rw [mem_Icc] at htI
        rw [mem_Ico]
        constructor
        · exact htI.1
        · have : (0 : ℝ) / w + c = c := by rw [zero_div, zero_add]
          linarith [htO.2, this ▸ htO.2]
      calc volume _ ≤ volume (Ico (0 : ℝ) c) := measure_mono this
        _ = ENNReal.ofReal c := by rw [Real.volume_Ico, sub_zero]
    · by_cases hww : a = (w : ℤ)
      · subst hww
        simp only [if_pos (Or.inr rfl)]
        have hdiv : ((w : ℤ) : ℝ) / w = 1 := by
          push_cast
          field_simp
        have : Ioo (((w : ℤ) : ℝ) / w - c) (((w : ℤ) : ℝ) / w + c) ∩ Icc (0 : ℝ) 1
            ⊆ Ioc (1 - c) 1 := by
          rintro t ⟨htO, htI⟩
          rw [mem_Ioo] at htO
          rw [mem_Icc] at htI
          rw [mem_Ioc]
          constructor
          · linarith [hdiv ▸ htO.1]
          · exact htI.2
        calc volume _ ≤ volume (Ioc (1 - c) 1) := measure_mono this
          _ = ENNReal.ofReal c := by rw [Real.volume_Ioc]; congr 1; ring
      · simp only [if_neg (by tauto)]
        calc volume _ ≤ volume (Ioo ((a : ℝ) / w - c) ((a : ℝ) / w + c)) :=
              measure_mono inter_subset_left
          _ = ENNReal.ofReal (2 * c) := by rw [Real.volume_Ioo]; congr 1; ring
  -- assemble: sum ≤ c + c + (w-1)*(2c) = 2*lam
  have hsumbound : ∑ a ∈ Finset.Icc (0 : ℤ) (w : ℤ),
      ENNReal.ofReal (if a = 0 ∨ a = (w : ℤ) then c else 2 * c)
        ≤ ENNReal.ofReal (2 * lam) := by
    rw [← ENNReal.ofReal_sum_of_nonneg (by
      intro a _
      by_cases h : a = 0 ∨ a = (w : ℤ) <;> simp [h] <;> positivity)]
    apply ENNReal.ofReal_le_ofReal
    -- real arithmetic: the sum is 2c + (card - 2) * 2c ≤ 2*lam when w ≥ 1
    have hsplit : ∑ a ∈ Finset.Icc (0 : ℤ) (w : ℤ),
        (if a = 0 ∨ a = (w : ℤ) then c else 2 * c)
        ≤ 2 * c + ((w : ℝ) - 1) * (2 * c) := by
      rcases Nat.eq_or_gt_of_le hw with h1 | h2
      · -- w = 1: the set is {0, 1}, both boundary
        rw [← h1]
        norm_num [Finset.sum_Icc_id_mul_two]
        rw [show Finset.Icc (0 : ℤ) (1 : ℤ) = {0, 1} by decide]
        rw [Finset.sum_insert (by decide), Finset.sum_singleton]
        norm_num
      · -- w ≥ 2: peel 0 and w
        have hw2 : (2 : ℤ) ≤ (w : ℤ) := by exact_mod_cast h2
        have hpeel : Finset.Icc (0 : ℤ) (w : ℤ)
            = insert 0 (insert (w : ℤ) (Finset.Icc 1 ((w : ℤ) - 1))) := by
          ext a
          simp only [Finset.mem_Icc, Finset.mem_insert]
          omega
        rw [hpeel, Finset.sum_insert (by
              simp only [Finset.mem_insert, Finset.mem_Icc]
              omega),
            Finset.sum_insert (by
              simp only [Finset.mem_Icc]
              omega)]
        simp only [if_pos (Or.inl rfl), if_pos (Or.inr rfl)]
        have hmid : ∑ a ∈ Finset.Icc (1 : ℤ) ((w : ℤ) - 1),
            (if a = 0 ∨ a = (w : ℤ) then c else 2 * c) = ((w : ℝ) - 1) * (2 * c) := by
          rw [Finset.sum_congr rfl (fun a ha => by
            rw [Finset.mem_Icc] at ha
            rw [if_neg (by omega)])]
          rw [Finset.sum_const, Int.card_Icc]
          have : ((w : ℤ) - 1 + 1 - 1).toNat = w - 1 := by omega
          rw [this, nsmul_eq_mul]
          have : ((w - 1 : ℕ) : ℝ) = (w : ℝ) - 1 := by
            push_cast [Nat.cast_sub (by omega : 1 ≤ w)]
            ring
          rw [this]
        rw [hmid]
        ring_nf
        linarith [hcpos]
    calc ∑ a ∈ Finset.Icc (0 : ℤ) (w : ℤ),
          (if a = 0 ∨ a = (w : ℤ) then c else 2 * c)
        ≤ 2 * c + ((w : ℝ) - 1) * (2 * c) := hsplit
      _ = (w : ℝ) * (2 * c) := by ring
      _ = 2 * lam := by rw [hc]; field_simp; ring
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
    _ ≤ ENNReal.ofReal (2 * lam) := hsumbound

end LRC14
