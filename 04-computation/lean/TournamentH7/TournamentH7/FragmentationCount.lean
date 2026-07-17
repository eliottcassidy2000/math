/- FragmentationCount.lean — mac-mini-2026-07-16-S128 (S127 draft repaired: the S127
   "verified" run trusted `$?` after a pipe — MISTAKE logged; this version builds under
   `lake build`, olean emitted).
   THM-883's kernel: the Fragmentation Lemma.
   Lemma A (arc count): the arc indices whose arc (a/w ± lam/w) can meet [x, x+L]
   number at most w*L + 2*lam + 1.
   Lemma B (fragmentation): the bad set of modulus w at radius lam meets [x, x+L]
   in measure at most (w*L + 2*lam + 1) * (2*lam/w). -/
import Mathlib.MeasureTheory.Measure.Lebesgue.Basic
import Mathlib.Data.Int.Interval

open MeasureTheory Set

namespace LRC14

/-- Arc indices that can meet `[x, x+L]`: `a ∈ [⌈w*x - lam⌉, ⌊w*x + w*L + lam⌋]`. -/
noncomputable def arcIdx (w : ℕ) (x L lam : ℝ) : Finset ℤ :=
  Finset.Icc ⌈(w : ℝ) * x - lam⌉ ⌊(w : ℝ) * x + (w : ℝ) * L + lam⌋

/-- **Lemma A.** `card (arcIdx) ≤ w*L + 2*lam + 1`. -/
theorem arcIdx_card_le (w : ℕ) (x L lam : ℝ) (hlam : 0 < lam) (hL : 0 ≤ L) :
    ((arcIdx w x L lam).card : ℝ) ≤ (w : ℝ) * L + 2 * lam + 1 := by
  classical
  unfold arcIdx
  set c : ℤ := ⌈(w : ℝ) * x - lam⌉ with hc
  set f : ℤ := ⌊(w : ℝ) * x + (w : ℝ) * L + lam⌋ with hf
  have h1 : ((w : ℝ) * x - lam) ≤ (c : ℝ) := hc ▸ Int.le_ceil _
  have h2 : ((f : ℤ) : ℝ) ≤ (w : ℝ) * x + (w : ℝ) * L + lam := hf ▸ Int.floor_le _
  by_cases h : c ≤ f
  · rw [Int.card_Icc]
    have hnn : (0 : ℤ) ≤ f + 1 - c := by omega
    have hcast : (((f + 1 - c).toNat : ℕ) : ℝ) = (f : ℝ) + 1 - (c : ℝ) := by
      have h' : (((f + 1 - c).toNat : ℕ) : ℤ) = f + 1 - c := Int.toNat_of_nonneg hnn
      have := congrArg (fun z : ℤ => (z : ℝ)) h'
      push_cast at this
      linarith [this]
    rw [hcast]
    linarith
  · rw [Finset.Icc_eq_empty h]
    have h3 : (f : ℝ) < (c : ℝ) := by exact_mod_cast lt_of_not_ge h
    simp only [Finset.card_empty, Nat.cast_zero]
    have hw0 : (0 : ℝ) ≤ (w : ℝ) := Nat.cast_nonneg w
    have hwL : (0 : ℝ) ≤ (w : ℝ) * L := mul_nonneg hw0 hL
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
  have hwne : (w : ℝ) ≠ 0 := ne_of_gt hwR
  have key2 : (w : ℝ) * t < (a : ℝ) + lam := by
    have := mul_lt_mul_of_pos_left h2 hwR
    have hexp : (w : ℝ) * ((a : ℝ) / w + lam / w) = (a : ℝ) + lam := by
      field_simp
    linarith [hexp ▸ this]
  have key1 : (a : ℝ) - lam < (w : ℝ) * t := by
    have := mul_lt_mul_of_pos_left h1 hwR
    have hexp : (w : ℝ) * ((a : ℝ) / w - lam / w) = (a : ℝ) - lam := by
      field_simp
    linarith [hexp ▸ this]
  have hxt : (w : ℝ) * x ≤ (w : ℝ) * t :=
    mul_le_mul_of_nonneg_left h3 (le_of_lt hwR)
  have htL : (w : ℝ) * t ≤ (w : ℝ) * x + (w : ℝ) * L := by
    have := mul_le_mul_of_nonneg_left h4 (le_of_lt hwR)
    linarith [this, mul_add (w : ℝ) x L]
  unfold arcIdx
  rw [Finset.mem_Icc]
  constructor
  · rw [Int.ceil_le]
    linarith
  · rw [Int.le_floor]
    linarith

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
  have harc : ∀ a : ℤ,
      volume (Ioo ((a : ℝ) / w - lam / w) ((a : ℝ) / w + lam / w))
        = ENNReal.ofReal (2 * (lam / w)) := by
    intro a
    rw [Real.volume_Ioo]
    congr 1
    ring
  have hcard : ((arcIdx w x L lam).card : ENNReal)
      = ENNReal.ofReal ((arcIdx w x L lam).card : ℝ) := by
    simp
  calc volume (badSet w lam ∩ Icc x (x + L))
      ≤ volume (⋃ a ∈ arcIdx w x L lam,
          Ioo ((a : ℝ) / w - lam / w) ((a : ℝ) / w + lam / w)) := measure_mono hsub
    _ ≤ ∑ a ∈ arcIdx w x L lam,
          volume (Ioo ((a : ℝ) / w - lam / w) ((a : ℝ) / w + lam / w)) :=
        measure_biUnion_finset_le _ _
    _ = ∑ _a ∈ arcIdx w x L lam, ENNReal.ofReal (2 * (lam / w)) :=
        Finset.sum_congr rfl fun a _ => harc a
    _ = ((arcIdx w x L lam).card : ENNReal) * ENNReal.ofReal (2 * (lam / w)) := by
        rw [Finset.sum_const, nsmul_eq_mul]
    _ ≤ ENNReal.ofReal ((w : ℝ) * L + 2 * lam + 1) * ENNReal.ofReal (2 * (lam / w)) := by
        rw [hcard]
        exact mul_le_mul_right'
          (ENNReal.ofReal_le_ofReal (arcIdx_card_le w x L lam hlam hL)) _
    _ = ENNReal.ofReal (((w : ℝ) * L + 2 * lam + 1) * (2 * lam / w)) := by
        rw [← ENNReal.ofReal_mul (by positivity)]
        congr 1
        ring

end LRC14
