/-
  TournamentH7.LRCEventMeasureBridge -- generic event-to-shape handoffs for the
  LRC14 witness/p0 route.

  `LRCWitnessBonferroni` and `LRCWitnessPartA` deliberately carry shape-indexed
  real functions (`nuShape`, `measGP`, `p0Shape`, `witnessG2`) as parameters.
  This file proves the reusable measure-theoretic bridge from actual events in a
  probability space to the two abstract hypotheses used there:

    * Bonferroni: `nuShape s + measGP s - 1 <= witnessG2 s`;
    * monotonicity: a pointwise inclusion of dense-cover events into p0 events
      gives `DShape s <= p0Shape s`, and with `DShape = 1 - nuShape`, the exact
      `hDp0 : 1 - nuShape s <= p0Shape s` hypothesis.

  The file is intentionally generic: it does not define the LRC events.  It
  removes only the boilerplate between future event definitions and the existing
  proof DAG.
-/

import TournamentH7.LRCBonferroniMeasure
import TournamentH7.LRCWitnessBonferroni

namespace LonelyRunner
namespace LRC14
namespace EventMeasureBridge

open MeasureTheory

variable {α : Type*} [MeasurableSpace α]

/-- Real-valued monotonicity for probability measures. -/
theorem measure_toReal_mono_of_subset
    (μ : Measure α) [IsProbabilityMeasure μ] {A B : Set α} (hAB : A ⊆ B) :
    (μ A).toReal ≤ (μ B).toReal := by
  exact ENNReal.toReal_mono (measure_ne_top μ B) (measure_mono hAB)

/-- Shape-indexed Bonferroni handoff from actual events to the abstract
`hbonf` hypothesis used by `LRCWitnessBonferroni` and `LRCWitnessPartA`. -/
theorem shape_bonferroni_handoff
    (μ : Measure α) [IsProbabilityMeasure μ]
    (GOOD GP : Shape → Set α) (hGP : ∀ s, MeasurableSet (GP s))
    (nuShape measGP : Shape → ℝ)
    (hnu : ∀ s, nuShape s = (μ (GOOD s)).toReal)
    (hgp : ∀ s, measGP s = (μ (GP s)).toReal)
    (hwitness : ∀ s, witnessG2 s = (μ ((GOOD s) ∩ (GP s))).toReal) :
    ∀ s, nuShape s + measGP s - 1 ≤ witnessG2 s := by
  intro s
  rw [hnu s, hgp s, hwitness s]
  exact BonferroniMeasure.toReal_bonferroni μ (GOOD s) (GP s) (hGP s)

/-- Shape-indexed measure monotonicity for a pointwise event inclusion.  This is
the measure form of the dense-cover implication once the actual dense and p0
events are defined. -/
theorem shape_measure_mono_handoff
    (μ : Measure α) [IsProbabilityMeasure μ]
    (Dset P0set : Shape → Set α)
    (DShape p0Shape : Shape → ℝ)
    (hsub : ∀ s, Dset s ⊆ P0set s)
    (hD : ∀ s, DShape s = (μ (Dset s)).toReal)
    (hp0 : ∀ s, p0Shape s = (μ (P0set s)).toReal) :
    ∀ s, DShape s ≤ p0Shape s := by
  intro s
  rw [hD s, hp0 s]
  exact measure_toReal_mono_of_subset μ (hsub s)

/-- The exact `hDp0` hypothesis used in the p0-wide-bound route, obtained from
an event inclusion plus the identity `DShape = 1 - nuShape`. -/
theorem shape_D_le_p0_handoff
    (μ : Measure α) [IsProbabilityMeasure μ]
    (Dset P0set : Shape → Set α)
    (nuShape DShape p0Shape : Shape → ℝ)
    (hsub : ∀ s, Dset s ⊆ P0set s)
    (hDmeasure : ∀ s, DShape s = (μ (Dset s)).toReal)
    (hp0measure : ∀ s, p0Shape s = (μ (P0set s)).toReal)
    (hDdef : ∀ s, DShape s = 1 - nuShape s) :
    ∀ s, (1 - nuShape s) ≤ p0Shape s := by
  intro s
  have hmono : DShape s ≤ p0Shape s :=
    shape_measure_mono_handoff μ Dset P0set DShape p0Shape
      hsub hDmeasure hp0measure s
  rw [← hDdef s]
  exact hmono

/-! ## Axiom audit -/

#print axioms measure_toReal_mono_of_subset
#print axioms shape_bonferroni_handoff
#print axioms shape_measure_mono_handoff
#print axioms shape_D_le_p0_handoff

end EventMeasureBridge
end LRC14
end LonelyRunner
