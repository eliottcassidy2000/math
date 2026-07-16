/- FragmentationLemma.lean — mac-mini-2026-07-16-S123 (draft), restructured by
   death-star-2026-07-16-S30 with the PERIODICITY architecture.

   THM-883's core counting inequality. IMPORTANT CORRECTION to the draft's proof plan:
   the arc-COUNTING route (step 1 of the draft) is subtly flawed — the number of arcs
   meeting I can reach ⌊Lw⌋ + 2 (two partial edge arcs), overshooting the claimed
   constant (L·w + 1)·(2λ/w). The correct route is PERIODIC WINDOWING:
     (A) badArcs is (1/w)-periodic (index shift);
     (B) any window of length ≤ 1/w meets badArcs in measure ≤ 2λ/w
         (translate into one period; the period's bad measure is exactly 2λ/w for λ ≤ ½);
     (C) chop I into ⌈wL⌉ ≤ wL + 1 windows and add.
   For λ ≥ ½ the target bound is trivial (vol ≤ L ≤ 2λL ≤ RHS).
   The two remaining sorries are (B1) the periodic-window bound and (C1) the chop-and-sum;
   both are localized with exact Mathlib tool pointers. -/
import Mathlib.MeasureTheory.Measure.Lebesgue.Basic
import Mathlib.MeasureTheory.Measure.Lebesgue.EqHaar

open MeasureTheory Set

namespace LRC14

/-- The bad arc set of modulus `w` at radius `lam`, lifted to ℝ. -/
def badArcs (w : ℕ) (lam : ℝ) : Set ℝ :=
  ⋃ a : ℤ, Set.Ioo ((a : ℝ) / w - lam / w) ((a : ℝ) / w + lam / w)

/-- (A) `badArcs` is `(1/w)`-periodic: translation by `1/w` permutes the arcs. -/
theorem badArcs_periodic (w : ℕ) (hw : 1 ≤ w) (lam : ℝ) :
    (fun y => y + 1 / (w : ℝ)) ⁻¹' badArcs w lam = badArcs w lam := by
  ext y
  simp only [badArcs, Set.mem_preimage, Set.mem_iUnion, Set.mem_Ioo]
  constructor
  · rintro ⟨a, h1, h2⟩
    refine ⟨a - 1, ?_, ?_⟩
    · have : ((a - 1 : ℤ) : ℝ) = (a : ℝ) - 1 := by push_cast; ring
      rw [this]
      have hww : (0 : ℝ) < (w : ℝ) := by exact_mod_cast Nat.lt_of_lt_of_le Nat.zero_lt_one hw
      rw [sub_div]
      linarith [h1]
    · have : ((a - 1 : ℤ) : ℝ) = (a : ℝ) - 1 := by push_cast; ring
      rw [this]
      rw [sub_div]
      linarith [h2]
  · rintro ⟨a, h1, h2⟩
    refine ⟨a + 1, ?_, ?_⟩
    · have : ((a + 1 : ℤ) : ℝ) = (a : ℝ) + 1 := by push_cast; ring
      rw [this, add_div]
      linarith [h1]
    · have : ((a + 1 : ℤ) : ℝ) = (a : ℝ) + 1 := by push_cast; ring
      rw [this, add_div]
      linarith [h2]

/-- (B) THE PERIODIC-WINDOW BOUND: a window of length `ℓ ≤ 1/w` meets `badArcs` in
    measure at most `2·lam/w` (for `0 < lam ≤ 1/2`).
    Tools: `badArcs_periodic` + `Real.volume_preserving_add_right` to translate the
    window into `[0, 1/w]`; there only arcs `a = 0, 1` intersect, contributing
    `[0, lam/w) ∪ ((1−lam)/w, 1/w]`, total `2·lam/w` (`Real.volume_Ioo`,
    `measure_union_le`). -/
theorem window_bound (w : ℕ) (hw : 1 ≤ w) (lam : ℝ)
    (hlam : 0 < lam) (hlam2 : lam ≤ 1 / 2) (y ℓ : ℝ) (hℓ0 : 0 ≤ ℓ)
    (hℓ : ℓ ≤ 1 / (w : ℝ)) :
    volume (badArcs w lam ∩ Set.Icc y (y + ℓ))
      ≤ ENNReal.ofReal (2 * lam / w) := by
  sorry

/-- (C) chop-and-sum: `Icc x (x+L)` is covered by `⌈w·L⌉ ≤ w·L + 1` windows of length
    `1/w`; apply `window_bound` to each and sum (`measure_biUnion_finset_le`,
    `Nat.ceil_le`, `ENNReal.ofReal` arithmetic). -/
theorem fragmentation (w : ℕ) (hw : 1 ≤ w) (lam L x : ℝ)
    (hlam : 0 < lam) (hL : 0 ≤ L) :
    volume (badArcs w lam ∩ Set.Icc x (x + L))
      ≤ ENNReal.ofReal ((L * w + 1) * (2 * lam / w)) := by
  have hww : (0 : ℝ) < (w : ℝ) := by exact_mod_cast Nat.lt_of_lt_of_le Nat.zero_lt_one hw
  by_cases hhalf : lam ≤ 1 / 2
  · -- main branch: periodic windowing
    sorry
  · -- trivial branch: lam > 1/2 makes the RHS at least L ≥ vol
    push_neg at hhalf
    have h1 : volume (badArcs w lam ∩ Set.Icc x (x + L)) ≤ volume (Set.Icc x (x + L)) :=
      measure_mono Set.inter_subset_right
    have h2 : volume (Set.Icc x (x + L)) = ENNReal.ofReal L := by
      rw [Real.volume_Icc]; ring_nf
    have h3 : L ≤ (L * w + 1) * (2 * lam / w) := by
      have hlw : 0 ≤ L * w := mul_nonneg hL (le_of_lt hww)
      have e1 : (L * w + 1) * (2 * lam / w) = 2 * lam * L + 2 * lam / w := by
        field_simp
      have e2 : L ≤ 2 * lam * L := by nlinarith
      have e3 : 0 ≤ 2 * lam / w := by positivity
      linarith
    calc volume (badArcs w lam ∩ Set.Icc x (x + L))
        ≤ ENNReal.ofReal L := h1.trans_eq h2
      _ ≤ ENNReal.ofReal ((L * w + 1) * (2 * lam / w)) := ENNReal.ofReal_le_ofReal h3

/-- Killer budget (THM-883 Lemma 2 shape): if the arc-grids of `w₁ … w_j` cover
    `[x, x+L]`, then `L(1 − 2jλ) ≤ 2λ Σ 1/wᵢ`. Follows from `fragmentation` by
    monotonicity + finite subadditivity + algebra. -/
theorem killer_budget (j : ℕ) (ws : Fin j → ℕ) (hws : ∀ i, 1 ≤ ws i)
    (lam L x : ℝ) (hlam : 0 < lam) (hL : 0 ≤ L)
    (hcover : Set.Icc x (x + L) ⊆ ⋃ i, badArcs (ws i) lam) :
    L * (1 - 2 * j * lam) ≤ 2 * lam * ∑ i, (1 : ℝ) / ws i := by
  -- vol(I) ≤ Σ vol(bad_i ∩ I) ≤ Σ (L·wᵢ + 1)(2λ/wᵢ) = 2λ(jL + Σ 1/wᵢ); rearrange.
  sorry

end LRC14
