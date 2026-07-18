/- LRCArcCounting.lean -- opus-2026-07-17-S353 (HYP-7360 / THM-1012).
   THE COUNTING LEMMA: an interval aligned at HALF-CELLS holding `m` whole
   cells of the `b`-comb contains at least `m` complete arcs, so

     m * (2*lam/b)  ≤  vol (badArcs b lam ∩ Ioo c (c + m/b)),
     c = (j + 1/2)/b.

   This is the engine of THM-1012 (the sharp nesting floor): an `a`-arc of
   length `2*lam/a` swallows ~`2*lam*b/a` whole `b`-cells, each contributing
   a full `b`-arc, which is how the pair overlap reaches the INDEPENDENCE
   constant `4*lam^2` with no sawtooth identity.

   PROOF DEVICE (reused from S351): half-cell alignment.  With the endpoints
   at `(j + 1/2)/b`, no arc straddles a cell boundary, so each cell contains
   exactly one whole arc and the count needs no lattice-point argument —
   induction on `m` with two-set additivity suffices, avoiding any
   `Finset`-indexed disjointness family.
   Kernel-pure target: no sorry, no native_decide. -/
import Mathlib
import TournamentH7.FragmentationLemma

open MeasureTheory

namespace LonelyRunner.LRC14.Hunter

/-- One aligned cell contains one whole arc: `2*lam/b ≤ vol (badArcs ∩ cell)`. -/
theorem one_cell_ge {b : ℕ} {lam : ℝ} (hb : 0 < b) (hlam : 0 < lam)
    (hhalf : 2 * lam ≤ 1) (j : ℤ) :
    ENNReal.ofReal (2 * lam / b)
      ≤ volume (LRC14.badArcs b lam
          ∩ Set.Ioo (((j : ℝ) + 1 / 2) / b) (((j : ℝ) + 1 / 2 + 1) / b)) := by
  have hbR : (0 : ℝ) < b := by exact_mod_cast hb
  -- the arc centred at (j+1)/b lies inside the cell and inside the comb
  have harc : Set.Ioo (((j : ℝ) + 1) / b - lam / b) (((j : ℝ) + 1) / b + lam / b)
      ⊆ LRC14.badArcs b lam
        ∩ Set.Ioo (((j : ℝ) + 1 / 2) / b) (((j : ℝ) + 1 / 2 + 1) / b) := by
    intro t ht
    rw [Set.mem_Ioo] at ht
    constructor
    · rw [LRC14.badArcs, Set.mem_iUnion]
      exact ⟨j + 1, by rw [Set.mem_Ioo]; push_cast; exact ⟨by linarith [ht.1], by linarith [ht.2]⟩⟩
    · rw [Set.mem_Ioo]
      constructor
      · have h : ((j : ℝ) + 1 / 2) / b ≤ ((j : ℝ) + 1) / b - lam / b := by
          rw [div_sub_div_same, div_le_div_iff_of_pos_right hbR]
          linarith
        linarith [ht.1]
      · have h : ((j : ℝ) + 1) / b + lam / b ≤ ((j : ℝ) + 1 / 2 + 1) / b := by
          have hrw : ((j : ℝ) + 1) / b + lam / b = ((j : ℝ) + 1 + lam) / b := by ring
          rw [hrw, div_le_div_iff_of_pos_right hbR]
          linarith
        linarith [ht.2]
  calc ENNReal.ofReal (2 * lam / b)
      = volume (Set.Ioo (((j : ℝ) + 1) / b - lam / b)
          (((j : ℝ) + 1) / b + lam / b)) := by
        rw [Real.volume_Ioo]; congr 1; ring
    _ ≤ _ := measure_mono harc

/-- **THE COUNTING LEMMA**: `m` aligned cells hold `m` whole arcs. -/
theorem aligned_count_ge {b : ℕ} {lam : ℝ} (hb : 0 < b) (hlam : 0 < lam)
    (hhalf : 2 * lam ≤ 1) (j : ℤ) :
    ∀ m : ℕ, ENNReal.ofReal (m * (2 * lam / b))
      ≤ volume (LRC14.badArcs b lam
          ∩ Set.Ioo (((j : ℝ) + 1 / 2) / b) (((j : ℝ) + 1 / 2 + m) / b)) := by
  have hbR : (0 : ℝ) < b := by exact_mod_cast hb
  intro m
  induction m with
  | zero => simp
  | succ m ih =>
      push_cast
      -- split the (m+1)-cell interval at the m-th boundary
      have hsplit : Set.Ioo (((j : ℝ) + 1 / 2) / b) (((j : ℝ) + 1 / 2 + m) / b)
            ∪ Set.Ioo (((j : ℝ) + 1 / 2 + m) / b) (((j : ℝ) + 1 / 2 + (m + 1)) / b)
          ⊆ Set.Ioo (((j : ℝ) + 1 / 2) / b) (((j : ℝ) + 1 / 2 + (m + 1)) / b) := by
        have hmono : (((j : ℝ) + 1 / 2 + m) / b) ≤ (((j : ℝ) + 1 / 2 + (m + 1)) / b) := by
          rw [div_le_div_iff_of_pos_right hbR]; linarith
        rintro t (ht | ht) <;> rw [Set.mem_Ioo] at ht ⊢
        · exact ⟨ht.1, lt_of_lt_of_le ht.2 hmono⟩
        · refine ⟨lt_of_le_of_lt ?_ ht.1, ht.2⟩
          rw [div_le_div_iff_of_pos_right hbR]
          have hm : (0 : ℝ) ≤ m := Nat.cast_nonneg m
          linarith
      -- the two pieces are disjoint
      have hdisj : Disjoint
          (LRC14.badArcs b lam ∩ Set.Ioo (((j : ℝ) + 1 / 2) / b) (((j : ℝ) + 1 / 2 + m) / b))
          (LRC14.badArcs b lam ∩ Set.Ioo (((j : ℝ) + 1 / 2 + m) / b)
            (((j : ℝ) + 1 / 2 + (m + 1)) / b)) := by
        apply Set.disjoint_left.mpr
        rintro t ⟨-, ht1⟩ ⟨-, ht2⟩
        rw [Set.mem_Ioo] at ht1 ht2
        linarith [ht1.2, ht2.1]
      have hmeas1 : MeasurableSet (LRC14.badArcs b lam
          ∩ Set.Ioo (((j : ℝ) + 1 / 2) / b) (((j : ℝ) + 1 / 2 + m) / b)) := by
        apply MeasurableSet.inter _ measurableSet_Ioo
        rw [LRC14.badArcs]
        exact MeasurableSet.iUnion fun _ => measurableSet_Ioo
      -- the last cell, re-indexed as the (j+m)-th aligned cell
      have hcell : ENNReal.ofReal (2 * lam / b)
          ≤ volume (LRC14.badArcs b lam ∩ Set.Ioo (((j : ℝ) + 1 / 2 + m) / b)
              (((j : ℝ) + 1 / 2 + (m + 1)) / b)) := by
        have := one_cell_ge hb hlam hhalf (j + m)
        have hrw1 : ((((j + m : ℤ)) : ℝ) + 1 / 2) / b = (((j : ℝ) + 1 / 2 + m) / b) := by
          push_cast; ring_nf
        have hrw2 : ((((j + m : ℤ)) : ℝ) + 1 / 2 + 1) / b
            = (((j : ℝ) + 1 / 2 + (m + 1)) / b) := by
          push_cast; ring_nf
        rwa [hrw1, hrw2] at this
      calc ENNReal.ofReal (((m : ℝ) + 1) * (2 * lam / b))
          = ENNReal.ofReal ((m : ℝ) * (2 * lam / b))
              + ENNReal.ofReal (2 * lam / b) := by
            rw [← ENNReal.ofReal_add (by positivity) (by positivity)]
            congr 1
            ring
        _ ≤ volume (LRC14.badArcs b lam
              ∩ Set.Ioo (((j : ℝ) + 1 / 2) / b) (((j : ℝ) + 1 / 2 + m) / b))
            + volume (LRC14.badArcs b lam
              ∩ Set.Ioo (((j : ℝ) + 1 / 2 + m) / b)
                (((j : ℝ) + 1 / 2 + (m + 1)) / b)) := add_le_add ih hcell
        _ = volume ((LRC14.badArcs b lam
              ∩ Set.Ioo (((j : ℝ) + 1 / 2) / b) (((j : ℝ) + 1 / 2 + m) / b))
            ∪ (LRC14.badArcs b lam ∩ Set.Ioo (((j : ℝ) + 1 / 2 + m) / b)
                (((j : ℝ) + 1 / 2 + (m + 1)) / b))) :=
            (measure_union hdisj (by
              apply MeasurableSet.inter _ measurableSet_Ioo
              rw [LRC14.badArcs]
              exact MeasurableSet.iUnion fun _ => measurableSet_Ioo)).symm
        _ ≤ _ := by
            apply measure_mono
            intro t ht
            rcases ht with ht | ht
            · exact Set.mem_inter ht.1 (hsplit (Or.inl ht.2))
            · exact Set.mem_inter ht.1 (hsplit (Or.inr ht.2))

/-! ## Axiom audit -/
#print axioms one_cell_ge
#print axioms aligned_count_ge

end LonelyRunner.LRC14.Hunter
