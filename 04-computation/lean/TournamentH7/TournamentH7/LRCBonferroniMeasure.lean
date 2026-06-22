/-
  TournamentH7.LRCBonferroniMeasure -- the Bonferroni / union measure inequality,
  the second elementary input of the witness/p0 UNIFICATION
  (kind-pasteur-2026-06-22-S30, HYP-2832).

  For a probability measure `μ` and sets `A, B`:

      μ(A ∩ B) ≥ μ(A) + μ(B) - 1.

  In the unification, `A = GOOD(E)` (the >1/7 max-gap set), `B = G_P` (the
  small-part safe set), and `μ(A ∩ B) = witnessG2 = rho*_glob`.  Together with the
  inclusion `D ⊆ p0` (`LRCDenseCovers.inner_sector_covered`), this yields
  `rho*_glob ≥ meas(G_P) - p0(E)`, the unification bound.

  Stated and proved generically over any probability measure; sorry-free.
-/

import Mathlib.MeasureTheory.Measure.MeasureSpace
import Mathlib.MeasureTheory.Measure.Typeclasses.Probability
import Mathlib.Tactic

namespace LonelyRunner
namespace BonferroniMeasure

open MeasureTheory

variable {α : Type*} [MeasurableSpace α]

/-- **Bonferroni (ENNReal form).**  For a probability measure `μ` and a measurable
set `B`, `μ A + μ B ≤ μ (A ∩ B) + 1`.  (No subtraction; the clean ENNReal
statement of `μ(A ∩ B) ≥ μA + μB − 1`.) -/
theorem measure_add_le_inter_add_one
    (μ : Measure α) [IsProbabilityMeasure μ] (A B : Set α) (hB : MeasurableSet B) :
    μ A + μ B ≤ μ (A ∩ B) + 1 := by
  have hsum : μ (A ∪ B) + μ (A ∩ B) = μ A + μ B := measure_union_add_inter A hB
  have hle : μ (A ∪ B) ≤ 1 := by
    calc μ (A ∪ B) ≤ μ Set.univ := measure_mono (Set.subset_univ _)
      _ = 1 := measure_univ
  calc μ A + μ B = μ (A ∪ B) + μ (A ∩ B) := hsum.symm
    _ ≤ 1 + μ (A ∩ B) := by gcongr
    _ = μ (A ∩ B) + 1 := by rw [add_comm]

/-- **Bonferroni (real form).**  For a probability measure `μ` and a measurable set
`B`, `(μ A).toReal + (μ B).toReal - 1 ≤ (μ (A ∩ B)).toReal`.  This is the exact
real-valued shape of the unification's `hbonf` hypothesis
`nuShape s + measGP s - 1 ≤ witnessG2 s`, with `nuShape = (μ A).toReal`,
`measGP = (μ B).toReal`, `witnessG2 = (μ (A ∩ B)).toReal`. -/
theorem toReal_bonferroni
    (μ : Measure α) [IsProbabilityMeasure μ] (A B : Set α) (hB : MeasurableSet B) :
    (μ A).toReal + (μ B).toReal - 1 ≤ (μ (A ∩ B)).toReal := by
  have hfin : ∀ s : Set α, μ s ≠ ⊤ := fun s => (measure_ne_top μ s)
  have h := measure_add_le_inter_add_one μ A B hB
  -- push toReal through the ENNReal inequality (all terms finite)
  have h' : (μ A + μ B).toReal ≤ (μ (A ∩ B) + 1).toReal :=
    ENNReal.toReal_mono (by
      simp only [ne_eq, ENNReal.add_eq_top, not_or]
      exact ⟨hfin _, ENNReal.one_ne_top⟩) h
  rw [ENNReal.toReal_add (hfin _) (hfin _),
      ENNReal.toReal_add (hfin _) ENNReal.one_ne_top, ENNReal.toReal_one] at h'
  linarith

/-! ## Axiom audit -/

#print axioms measure_add_le_inter_add_one
#print axioms toReal_bonferroni

end BonferroniMeasure
end LonelyRunner
