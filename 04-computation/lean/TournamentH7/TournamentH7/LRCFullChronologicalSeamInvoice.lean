/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Codex (LRC multi-agent project, 2026-07-19)
-/
import Mathlib.Tactic

/-!
# Arithmetic and interval-order core for THM-1253

The paper layer selects a deletion-minimal chain of strict danger teeth.  The
first lemma below is the local redundancy argument that forces
`beta_i <= alpha_(i+2)`; the second records the resulting separation of
consecutive handoff intervals.  The remaining lemmas kernel-check the full
coverage-excess invoice and its lcm and occurrence-count consumers.
-/

namespace LRC14
namespace FullChronologicalSeamInvoice

/-- If the first and third intervals in a monotone three-interval chain
overlap, their union contains the middle interval.  This contradicts deletion
minimality in the paper proof. -/
theorem middle_interval_redundant
    {a₀ a₁ a₂ b₀ b₁ b₂ x : ℝ}
    (ha : a₀ < a₁ ∧ a₁ < a₂)
    (hb : b₀ < b₁ ∧ b₁ < b₂)
    (hoverlap : a₂ < b₀)
    (hx : a₁ < x ∧ x < b₁) :
    (a₀ < x ∧ x < b₀) ∨ (a₂ < x ∧ x < b₂) := by
  by_cases hleft : x < b₀
  · exact Or.inl ⟨lt_trans ha.1 hx.1, hleft⟩
  · right
    exact ⟨lt_of_lt_of_le hoverlap (le_of_not_gt hleft),
      lt_trans hx.2 hb.2⟩

/-- Once `b_i <= a_(i+2)`, adjacent raw handoff intervals have disjoint
interiors. -/
theorem consecutive_handoffs_disjoint
    {a₁ a₂ b₀ b₁ x : ℝ} (hsep : b₀ ≤ a₂) :
    ¬ ((a₁ < x ∧ x < b₀) ∧ (a₂ < x ∧ x < b₁)) := by
  intro hx
  linarith [hx.1.2, hx.2.1]

/-- The whole disjoint chronological handoff family lies in the multiplicity
excess.  Combining that lower bound with the six singleton caps gives the
coefficient `49/6` without a forest or Cayley averaging loss. -/
theorem full_handoff_excess_debt
    {gap invc harmonic totalMass handoffSum : ℝ}
    (hgap : gap = 6 * invc / 7)
    (hsingletons : totalMass ≤ 6 * gap / 7 + 6 * harmonic / 49)
    (hexcess : gap + handoffSum ≤ totalMass) :
    invc + 49 * handoffSum / 6 ≤ harmonic := by
  rw [hgap] at hsingletons hexcess
  linarith

/-- Every raw handoff contributes at least `1/(14*lcm)`, so the full word
receives coefficient `7/12`. -/
theorem full_lcm_occurrence_debt
    {gap invc harmonic totalMass handoffSum lcmSum : ℝ}
    (hgap : gap = 6 * invc / 7)
    (hsingletons : totalMass ≤ 6 * gap / 7 + 6 * harmonic / 49)
    (hexcess : gap + handoffSum ≤ totalMass)
    (hquantum : lcmSum / 14 ≤ handoffSum) :
    invc + 7 * lcmSum / 12 ≤ harmonic := by
  have hfull := full_handoff_excess_debt hgap hsingletons hexcess
  linarith

/-- The direct coefficient is exactly three times the Cayley-averaged
multiplicity coefficient in THM-1250. -/
theorem factor_three_over_tree_average :
    (7 : ℚ) / 12 = 3 * (7 / 36) := by
  norm_num

/-- The same disjoint family upgrades THM-1199's weighted functional drift:
normalization `7c/6`, density floor `3/4`, and tooth quantum `1/14` combine
to the exact coefficient `c/16`. -/
theorem full_weighted_functional_coefficient :
    (7 : ℚ) / 6 * (3 / 4) * (1 / 14) = 1 / 16 := by
  norm_num

/-- Abstract consumer for the full weighted functional seam family. -/
theorem full_weighted_functional_debt
    {functional c seamSum : ℝ}
    (hweighted : c * seamSum / 16 ≤ functional) :
    c * seamSum ≤ 16 * functional := by
  linarith

/-- If every one of `N-1` handoffs has reciprocal lcm at least
`g0/d6^2`, the full invoice has the displayed scalar consumer. -/
theorem occurrence_count_consumer
    {invc harmonic lcmSum g₀ d₆ occurrences : ℝ}
    (hdebt : invc + 7 * lcmSum / 12 ≤ harmonic)
    (hlcm : occurrences * g₀ / d₆ ^ 2 ≤ lcmSum) :
    invc + 7 * g₀ * occurrences / (12 * d₆ ^ 2) ≤ harmonic := by
  have hscaled :
      (7 / 12 : ℝ) * (occurrences * g₀ / d₆ ^ 2) ≤
        (7 / 12 : ℝ) * lcmSum :=
    mul_le_mul_of_nonneg_left hlcm (by norm_num)
  calc
    invc + 7 * g₀ * occurrences / (12 * d₆ ^ 2) =
        invc + (7 / 12 : ℝ) * (occurrences * g₀ / d₆ ^ 2) := by ring
    _ ≤ invc + (7 / 12 : ℝ) * lcmSum := by
      simpa [add_comm] using add_le_add_left hscaled invc
    _ = invc + 7 * lcmSum / 12 := by ring
    _ ≤ harmonic := hdebt

/-- Multiplying the full occurrence invoice by the positive carrier gives
its common-dilate invariant form. -/
theorem scale_covariant_full_invoice
    {c invc harmonic lcmSum : ℝ} (hc : 0 ≤ c)
    (hinv : c * invc = 1)
    (hdebt : invc + 7 * lcmSum / 12 ≤ harmonic) :
    1 + 7 * c * lcmSum / 12 ≤ c * harmonic := by
  nlinarith

#print axioms middle_interval_redundant
#print axioms consecutive_handoffs_disjoint
#print axioms full_handoff_excess_debt
#print axioms full_lcm_occurrence_debt
#print axioms factor_three_over_tree_average
#print axioms full_weighted_functional_coefficient
#print axioms full_weighted_functional_debt
#print axioms occurrence_count_consumer
#print axioms scale_covariant_full_invoice

end FullChronologicalSeamInvoice
end LRC14
