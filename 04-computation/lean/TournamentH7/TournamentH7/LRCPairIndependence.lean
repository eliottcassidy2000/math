/- LRCPairIndependence.lean -- opus-2026-07-17-S354 (HYP-7370 / THM-1012).
   THE ARITHMETIC WRAPPER: from the S353 counting lemma to the sharp nesting
   floor.  Three pieces:

   (1) `consecutive_cells_ge` — THE GENERAL INDUCTION.  If every cell of a
       uniform partition contributes at least `X` to `vol (S ∩ ·)`, then `n`
       consecutive cells contribute at least `n*X`.  This subsumes S353's
       `aligned_count_ge` AND performs the OUTER sum over the `a`-arcs — one
       lemma, both sums, no `Finset`-indexed disjointness family anywhere.

   (2) `aligned_block_in_interval` — THE ALIGNMENT FINDER.  Any interval of
       length `≥ (m+1)/b` contains a half-cell-aligned block of `m` whole
       `b`-cells: take `j = ⌈p*b − 1/2⌉`, so `(j+1/2)/b ≥ p` by `Int.le_ceil`
       and `(j+1/2+m)/b < p + (m+1)/b` by `Int.ceil_lt_add_one`.

   (3) the assembly, giving the INDEPENDENCE constant with a linear defect.

   Kernel-pure target: no sorry, no native_decide. -/
import Mathlib
import TournamentH7.LRCArcCounting

open MeasureTheory

namespace LonelyRunner.LRC14.Hunter

/-- **(1) THE GENERAL CONSECUTIVE-CELLS INDUCTION.** -/
theorem consecutive_cells_ge (S : Set ℝ) (hS : MeasurableSet S)
    (c h X : ℝ) (hh : 0 < h) (hX : 0 ≤ X)
    (hcell : ∀ k : ℕ, ENNReal.ofReal X
      ≤ volume (S ∩ Set.Ioo (c + k * h) (c + (k + 1) * h))) :
    ∀ n : ℕ, ENNReal.ofReal (n * X)
      ≤ volume (S ∩ Set.Ioo c (c + n * h)) := by
  intro n
  induction n with
  | zero => simp
  | succ n ih =>
      push_cast
      have hdisj : Disjoint (S ∩ Set.Ioo c (c + (n : ℝ) * h))
          (S ∩ Set.Ioo (c + (n : ℝ) * h) (c + ((n : ℝ) + 1) * h)) := by
        apply Set.disjoint_left.mpr
        rintro t ⟨-, ht1⟩ ⟨-, ht2⟩
        rw [Set.mem_Ioo] at ht1 ht2
        linarith [ht1.2, ht2.1]
      have hsub : (S ∩ Set.Ioo c (c + (n : ℝ) * h))
            ∪ (S ∩ Set.Ioo (c + (n : ℝ) * h) (c + ((n : ℝ) + 1) * h))
          ⊆ S ∩ Set.Ioo c (c + ((n : ℝ) + 1) * h) := by
        have hnn : (0 : ℝ) ≤ n := Nat.cast_nonneg n
        rintro t (⟨hS', ht⟩ | ⟨hS', ht⟩) <;> rw [Set.mem_Ioo] at ht <;>
          refine ⟨hS', ?_⟩ <;> rw [Set.mem_Ioo]
        · exact ⟨ht.1, by nlinarith [ht.2]⟩
        · exact ⟨by nlinarith [ht.1], ht.2⟩
      calc ENNReal.ofReal (((n : ℝ) + 1) * X)
          = ENNReal.ofReal ((n : ℝ) * X) + ENNReal.ofReal X := by
            rw [← ENNReal.ofReal_add (by positivity) hX]
            congr 1
            ring
        _ ≤ volume (S ∩ Set.Ioo c (c + (n : ℝ) * h))
            + volume (S ∩ Set.Ioo (c + (n : ℝ) * h) (c + ((n : ℝ) + 1) * h)) := by
            refine add_le_add ih ?_
            have := hcell n
            push_cast at this
            exact this
        _ = volume ((S ∩ Set.Ioo c (c + (n : ℝ) * h))
            ∪ (S ∩ Set.Ioo (c + (n : ℝ) * h) (c + ((n : ℝ) + 1) * h))) :=
            (measure_union hdisj (hS.inter measurableSet_Ioo)).symm
        _ ≤ _ := measure_mono hsub

/-- **(2) THE ALIGNMENT FINDER.** -/
theorem aligned_block_in_interval {b : ℕ} (hb : 0 < b) (p q : ℝ) (m : ℕ)
    (hlen : ((m : ℝ) + 1) / b ≤ q - p) :
    ∃ j : ℤ, p ≤ ((j : ℝ) + 1 / 2) / b
      ∧ ((j : ℝ) + 1 / 2 + m) / b ≤ q := by
  have hbR : (0 : ℝ) < b := by exact_mod_cast hb
  refine ⟨⌈p * b - 1 / 2⌉, ?_, ?_⟩
  · rw [le_div_iff₀ hbR]
    have := Int.le_ceil (p * b - 1 / 2)
    linarith
  · rw [div_le_iff₀ hbR]
    have hlt := Int.ceil_lt_add_one (p * b - 1 / 2)
    have hq : ((m : ℝ) + 1) ≤ (q - p) * b := by
      rw [div_le_iff₀ hbR] at hlen
      linarith
    linarith

/-- **(3a) THE PER-`a`-CELL BOUND**: inside one aligned `a`-cell, the `a`-arc
swallows an aligned block of `m` whole `b`-cells, each a full `b`-arc. -/
theorem per_a_cell_ge {a b : ℕ} {lam : ℝ} (ha : 0 < a) (hb : 0 < b)
    (hlam : 0 < lam) (hhalf : 2 * lam ≤ 1) (m : ℕ)
    (hm : ((m : ℝ) + 1) / b ≤ 2 * lam / a) (i : ℤ) (k : ℕ) :
    ENNReal.ofReal ((m : ℝ) * (2 * lam / b))
      ≤ volume ((LRC14.badArcs a lam ∩ LRC14.badArcs b lam)
          ∩ Set.Ioo (((i : ℝ) + 1 / 2) / a + k * (1 / a))
                    (((i : ℝ) + 1 / 2) / a + ((k : ℝ) + 1) * (1 / a))) := by
  have haR : (0 : ℝ) < a := by exact_mod_cast ha
  set p : ℝ := ((i : ℝ) + (k : ℝ) + 1) / a - lam / a with hp
  set q : ℝ := ((i : ℝ) + (k : ℝ) + 1) / a + lam / a with hq
  -- the a-arc has length 2*lam/a, enough for an aligned m-block of b-cells
  have hlen : ((m : ℝ) + 1) / b ≤ q - p := by
    have : q - p = 2 * lam / a := by rw [hp, hq]; ring
    rw [this]; exact hm
  obtain ⟨j, hj1, hj2⟩ := aligned_block_in_interval hb p q m hlen
  -- the aligned block sits inside the a-arc, which sits inside badArcs a and the cell
  have harc_sub : Set.Ioo (((j : ℝ) + 1 / 2) / b) (((j : ℝ) + 1 / 2 + m) / b)
      ⊆ (LRC14.badArcs a lam)
        ∩ Set.Ioo (((i : ℝ) + 1 / 2) / a + k * (1 / a))
                  (((i : ℝ) + 1 / 2) / a + ((k : ℝ) + 1) * (1 / a)) := by
    intro t ht
    rw [Set.mem_Ioo] at ht
    have htp : p < t := lt_of_le_of_lt hj1 ht.1
    have htq : t < q := lt_of_lt_of_le ht.2 hj2
    constructor
    · rw [LRC14.badArcs, Set.mem_iUnion]
      refine ⟨i + k + 1, ?_⟩
      rw [Set.mem_Ioo]
      push_cast
      constructor
      · have : p = ((i : ℝ) + (k : ℝ) + 1) / a - lam / a := hp
        linarith [htp]
      · have : q = ((i : ℝ) + (k : ℝ) + 1) / a + lam / a := hq
        linarith [htq]
    · rw [Set.mem_Ioo]
      have hcellL : ((i : ℝ) + 1 / 2) / a + k * (1 / a) ≤ p := by
        rw [hp]
        have hrw1 : ((i : ℝ) + 1 / 2) / a + (k : ℝ) * (1 / a)
            = ((i : ℝ) + 1 / 2 + k) / a := by field_simp
        have hrw2 : ((i : ℝ) + (k : ℝ) + 1) / a - lam / a
            = ((i : ℝ) + (k : ℝ) + 1 - lam) / a := by ring
        rw [hrw1, hrw2, div_le_div_iff_of_pos_right haR]
        linarith
      have hcellR : q ≤ ((i : ℝ) + 1 / 2) / a + ((k : ℝ) + 1) * (1 / a) := by
        rw [hq]
        have hrw1 : ((i : ℝ) + 1 / 2) / a + ((k : ℝ) + 1) * (1 / a)
            = ((i : ℝ) + 1 / 2 + ((k : ℝ) + 1)) / a := by field_simp
        have hrw2 : ((i : ℝ) + (k : ℝ) + 1) / a + lam / a
            = ((i : ℝ) + (k : ℝ) + 1 + lam) / a := by ring
        rw [hrw1, hrw2, div_le_div_iff_of_pos_right haR]
        linarith
      exact ⟨lt_of_le_of_lt hcellL htp, lt_of_lt_of_le htq hcellR⟩
  calc ENNReal.ofReal ((m : ℝ) * (2 * lam / b))
      ≤ volume (LRC14.badArcs b lam
          ∩ Set.Ioo (((j : ℝ) + 1 / 2) / b) (((j : ℝ) + 1 / 2 + m) / b)) :=
        aligned_count_ge hb hlam hhalf j m
    _ ≤ _ := by
        apply measure_mono
        rintro t ⟨htb, htI⟩
        obtain ⟨hta, htcell⟩ := harc_sub htI
        exact ⟨⟨hta, htb⟩, htcell⟩

/-- **(3b) THM-1012, the sharp nesting floor.**  Over a unit window the pair
overlap is at least `a * m * (2*lam/b)`; with `m ≈ 2*lam*b/a` this is the
INDEPENDENCE constant `4*lam^2` up to a defect linear in `a/b`. -/
theorem pair_overlap_independence {a b : ℕ} {lam : ℝ} (ha : 0 < a) (hb : 0 < b)
    (hlam : 0 < lam) (hhalf : 2 * lam ≤ 1) (m : ℕ)
    (hm : ((m : ℝ) + 1) / b ≤ 2 * lam / a) (i : ℤ) :
    ENNReal.ofReal ((a : ℝ) * ((m : ℝ) * (2 * lam / b)))
      ≤ volume ((LRC14.badArcs a lam ∩ LRC14.badArcs b lam)
          ∩ Set.Ioo (((i : ℝ) + 1 / 2) / a)
                    (((i : ℝ) + 1 / 2) / a + (a : ℝ) * (1 / a))) := by
  have haR : (0 : ℝ) < a := by exact_mod_cast ha
  have hmeas : MeasurableSet (LRC14.badArcs a lam ∩ LRC14.badArcs b lam) := by
    refine MeasurableSet.inter ?_ ?_ <;>
      · rw [LRC14.badArcs]
        exact MeasurableSet.iUnion fun _ => measurableSet_Ioo
  exact consecutive_cells_ge _ hmeas (((i : ℝ) + 1 / 2) / a) (1 / a)
    ((m : ℝ) * (2 * lam / b)) (by positivity) (by positivity)
    (fun k => per_a_cell_ge ha hb hlam hhalf m hm i k) a

/-! ## Axiom audit -/
#print axioms consecutive_cells_ge
#print axioms aligned_block_in_interval
#print axioms per_a_cell_ge
#print axioms pair_overlap_independence

end LonelyRunner.LRC14.Hunter
