/- FragmentationCount.lean — mac-mini-2026-07-16-S127.
   THM-883's kernel: the Fragmentation Lemma, formalized.
   Lemma A (arc count): the arc indices whose arc (a/w ± lam/w) can meet [x, x+L]
   number at most w*L + 2*lam + 1.
   Lemma B (fragmentation): the bad set of modulus w at radius lam meets [x, x+L]
   in measure at most (w*L + 2*lam + 1) * (2*lam/w).
   These are the inequalities behind the j ≤ 5 multi-killer sweep (THM-883) and the
   killer-budget chain; the +2*lam carries the endpoint arcs honestly. -/
import Mathlib.MeasureTheory.Measure.Lebesgue.Basic
import Mathlib.Data.Int.Interval

open MeasureTheory Set

namespace LRC14

/-- Arc indices that can meet `[x, x+L]`: `a ∈ [⌈w*x - lam⌉, ⌊w*x + w*L + lam⌋]`. -/
noncomputable def arcIdx (w : ℕ) (x L lam : ℝ) : Finset ℤ :=
  Finset.Icc ⌈(w : ℝ) * x - lam⌉ ⌊(w : ℝ) * x + (w : ℝ) * L + lam⌋

/-- **Lemma A.** `card (arcIdx) ≤ w*L + 2*lam + 1`. -/
theorem arcIdx_card_le (w : ℕ) (x L lam : ℝ) :
    ((arcIdx w x L lam).card : ℝ) ≤ (w : ℝ) * L + 2 * lam + 1 := by
  classical
  unfold arcIdx
  rcases le_or_lt (⌈(w : ℝ) * x - lam⌉ : ℤ) ⌊(w : ℝ) * x + (w : ℝ) * L + lam⌋ with h | h
  · rw [Int.card_Icc]
    have h1 : ((w : ℝ) * x - lam) ≤ (⌈(w : ℝ) * x - lam⌉ : ℝ) := Int.le_ceil _
    have h2 : ((⌊(w : ℝ) * x + (w : ℝ) * L + lam⌋ : ℤ) : ℝ)
        ≤ (w : ℝ) * x + (w : ℝ) * L + lam := Int.floor_le _
    have hnn : (0 : ℤ) ≤ ⌊(w : ℝ) * x + (w : ℝ) * L + lam⌋ + 1 - ⌈(w : ℝ) * x - lam⌉ := by
      omega
    rw [Int.toNat_of_nonneg hnn |>.symm] at *
    push_cast [Int.toNat_of_nonneg hnn]
    linarith
  · rw [Finset.Icc_eq_empty (not_le.mpr h)]
    have h1 : ((w : ℝ) * x - lam) ≤ (⌈(w : ℝ) * x - lam⌉ : ℝ) := Int.le_ceil _
    have h2 : ((⌊(w : ℝ) * x + (w : ℝ) * L + lam⌋ : ℤ) : ℝ)
        ≤ (w : ℝ) * x + (w : ℝ) * L + lam := Int.floor_le _
    have h3 : ((⌊(w : ℝ) * x + (w : ℝ) * L + lam⌋ : ℝ))
        < ((⌈(w : ℝ) * x - lam⌉ : ℝ)) := by exact_mod_cast h
    simp only [Finset.card_empty, Nat.cast_zero]
    linarith

/-- The bad set of modulus `w` at radius `lam`. -/
def badSet (w : ℕ) (lam : ℝ) : Set ℝ :=
  ⋃ a : ℤ, Ioo ((a : ℝ) / w - lam / w) ((a : ℝ) / w + lam / w)

/-- Membership transfer: a point of `arc_a ∩ [x, x+L]` forces `a ∈ arcIdx`. -/
theorem mem_arcIdx_of_hit (w : ℕ) (hw : 0 < w) (x L lam : ℝ) (a : ℤ) (t : ℝ)
    (h1 : (a : ℝ) / w - lam / w < t) (h2 : t < (a : ℝ) / w + lam / w)
    (h3 : x ≤ t) (h4 : t ≤ x + L) :
    a ∈ arcIdx w x L lam := by
  have hwR : (0 : ℝ) < (w : ℝ) := by exact_mod_cast hw
  unfold arcIdx
  rw [Finset.mem_Icc]
  constructor
  · rw [Int.ceil_le]
    have : (a : ℝ) > (t - lam / w) * 1 := by nlinarith
    have hta : (w : ℝ) * t - lam < (a : ℝ) * 1 := by
      have := mul_lt_mul_of_pos_left h2 hwR
      have hw' : (w : ℝ) ≠ 0 := ne_of_gt hwR
      field_simp at this
      nlinarith
    nlinarith
  · rw [Int.le_floor]
    have hta : ((a : ℝ)) < (w : ℝ) * t + lam := by
      have := mul_lt_mul_of_pos_left h1 hwR
      have hw' : (w : ℝ) ≠ 0 := ne_of_gt hwR
      field_simp at this
      nlinarith
    push_cast
    nlinarith

/-- **Lemma B / THM-883 Fragmentation.** -/
theorem fragmentation (w : ℕ) (hw : 0 < w) (x L lam : ℝ) (hlam : 0 < lam) (hL : 0 ≤ L) :
    volume (badSet w lam ∩ Icc x (x + L))
      ≤ ENNReal.ofReal (((w : ℝ) * L + 2 * lam + 1) * (2 * lam / w)) := by
  classical
  have hwR : (0 : ℝ) < (w : ℝ) := by exact_mod_cast hw
  have hsub : badSet w lam ∩ Icc x (x + L)
      ⊆ ⋃ a ∈ arcIdx w x L lam, Ioo ((a : ℝ) / w - lam / w) ((a : ℝ) / w + lam / w) := by
    rintro t ⟨ht, htI⟩
    rcases mem_iUnion.mp ht with ⟨a, hta⟩
    rw [mem_Ioo] at hta
    rw [mem_Icc] at htI
    exact mem_biUnion (mem_arcIdx_of_hit w hw x L lam a t hta.1 hta.2 htI.1 htI.2)
      (by rw [mem_Ioo]; exact hta)
  calc volume (badSet w lam ∩ Icc x (x + L))
      ≤ volume (⋃ a ∈ arcIdx w x L lam,
          Ioo ((a : ℝ) / w - lam / w) ((a : ℝ) / w + lam / w)) := measure_mono hsub
    _ ≤ ∑ a ∈ arcIdx w x L lam,
          volume (Ioo ((a : ℝ) / w - lam / w) ((a : ℝ) / w + lam / w)) :=
        measure_biUnion_finset_le _ _
    _ = ∑ _a ∈ arcIdx w x L lam, ENNReal.ofReal (2 * (lam / w)) := by
        refine Finset.sum_congr rfl fun a _ => ?_
        rw [Real.volume_Ioo]
        congr 1
        ring
    _ = ((arcIdx w x L lam).card : ℝ≥0∞) * ENNReal.ofReal (2 * (lam / w)) := by
        rw [Finset.sum_const, nsmul_eq_mul]
    _ ≤ ENNReal.ofReal ((w : ℝ) * L + 2 * lam + 1) * ENNReal.ofReal (2 * (lam / w)) := by
        gcongr
        rw [show ((arcIdx w x L lam).card : ℝ≥0∞)
            = ENNReal.ofReal ((arcIdx w x L lam).card : ℝ) by
          simp]
        exact ENNReal.ofReal_le_ofReal (arcIdx_card_le w x L lam)
    _ = ENNReal.ofReal (((w : ℝ) * L + 2 * lam + 1) * (2 * lam / w)) := by
        rw [← ENNReal.ofReal_mul (by positivity)]
        congr 1
        ring

end LRC14
