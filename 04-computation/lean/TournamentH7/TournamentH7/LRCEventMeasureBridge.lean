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
import TournamentH7.LRCDenseCovers
import TournamentH7.LRCWitnessFloorConcrete

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

/-! ## Concrete slow-time event specializations -/

/-- Bonferroni handoff specialized to the concrete small-part safe event
`G_P = safeSet P` and the slow-time probability measure on one period.  The
only remaining abstract event is `GOOD`, the max-gap witness event; `safeSet`
measurability is supplied by `LRCDenseCovers`. -/
theorem shape_bonferroni_safeSet_handoff
    (GOOD : Shape → Set ℝ) (Pof : Shape → List ℤ)
    (nuShape measGP : Shape → ℝ)
    (hnu : ∀ s, nuShape s = (DenseCovers.slowμ (GOOD s)).toReal)
    (hgp : ∀ s, measGP s =
      (DenseCovers.slowμ (DenseCovers.safeSet (Pof s))).toReal)
    (hwitness : ∀ s, witnessG2 s =
      (DenseCovers.slowμ ((GOOD s) ∩ DenseCovers.safeSet (Pof s))).toReal) :
    ∀ s, nuShape s + measGP s - 1 ≤ witnessG2 s :=
  shape_bonferroni_handoff DenseCovers.slowμ GOOD
    (fun s => DenseCovers.safeSet (Pof s))
    (fun s => DenseCovers.measurableSet_safeSet (Pof s))
    nuShape measGP hnu hgp hwitness

/-- `D <= p0` handoff specialized to the concrete slow-time events
`D = denseSet E` and `p0 = coverSet E`.  The anchor hypothesis `0 ∈ E` is exactly
what turns the pointwise dense-cover lemma into the set inclusion. -/
theorem shape_D_le_p0_denseCovers_handoff
    (Eof : Shape → List ℤ) (hanchor : ∀ s, (0 : ℤ) ∈ Eof s)
    (nuShape DShape p0Shape : Shape → ℝ)
    (hDmeasure : ∀ s, DShape s =
      (DenseCovers.slowμ (DenseCovers.denseSet (Eof s))).toReal)
    (hp0measure : ∀ s, p0Shape s =
      (DenseCovers.slowμ (DenseCovers.coverSet (Eof s))).toReal)
    (hDdef : ∀ s, DShape s = 1 - nuShape s) :
    ∀ s, (1 - nuShape s) ≤ p0Shape s :=
  shape_D_le_p0_handoff DenseCovers.slowμ
    (fun s => DenseCovers.denseSet (Eof s))
    (fun s => DenseCovers.coverSet (Eof s))
    nuShape DShape p0Shape
    (fun s => DenseCovers.denseSet_subset_coverSet (Eof s) (hanchor s))
    hDmeasure hp0measure hDdef

/-! ## Concrete goodSet witness readout specializations -/

/-- Bonferroni handoff with both slow-time events specialized: the large witness
event is the concrete `GoodSet.goodSet (Eof s)` and the small-part event is the
concrete `safeSet (Pof s)`. -/
theorem shape_bonferroni_goodSet_safeSet_handoff
    (Eof Pof : Shape → List ℤ)
    (nuShape measGP : Shape → ℝ)
    (hnu : ∀ s, nuShape s =
      (DenseCovers.slowμ (TournamentH7.GoodSet.goodSet (Eof s))).toReal)
    (hgp : ∀ s, measGP s =
      (DenseCovers.slowμ (DenseCovers.safeSet (Pof s))).toReal)
    (hwitness : ∀ s, witnessG2 s =
      (DenseCovers.slowμ
        (TournamentH7.GoodSet.goodSet (Eof s) ∩
          DenseCovers.safeSet (Pof s))).toReal) :
    ∀ s, nuShape s + measGP s - 1 ≤ witnessG2 s :=
  shape_bonferroni_safeSet_handoff
    (fun s => TournamentH7.GoodSet.goodSet (Eof s))
    Pof nuShape measGP hnu hgp hwitness

/-- Shape-level strict-cover readout: once `witnessG2` is identified with the
concrete `goodSet ∩ safeSet` slow-time measure, the strict p0 cover bound and
the non-strict cap floor give positive `witnessG2` directly. -/
theorem shape_goodSet_witness_pos_from_strict_cover_bound
    (Eof Pof : Shape → List ℤ) (cap : Shape → ℝ)
    (hanchor : ∀ s, (0 : ℤ) ∈ Eof s)
    (hwitness : ∀ s, witnessG2 s =
      (DenseCovers.slowμ
        (TournamentH7.GoodSet.goodSet (Eof s) ∩
          DenseCovers.safeSet (Pof s))).toReal)
    (hwide : ∀ s,
      (DenseCovers.slowμ (DenseCovers.coverSet (Eof s))).toReal < cap s)
    (hdual : ∀ s, cap s ≤
      (DenseCovers.slowμ (DenseCovers.safeSet (Pof s))).toReal) :
    ∀ s, 0 < witnessG2 s := by
  intro s
  rw [hwitness s]
  exact DenseCovers.goodSet_witness_pos_from_strict_cover_bound
    (Eof s) (Pof s) (cap s) (hanchor s) (hwide s) (hdual s)

/-- Shape-level quantitative readout: the p0 wide margin is transported through
the proved concrete carrier chain into `witnessG2 = μ(goodSet E ∩ safeSet P)`. -/
theorem shape_goodSet_witness_margin_from_wide_bound
    (Eof Pof : Shape → List ℤ) (cap delta : Shape → ℝ)
    (hanchor : ∀ s, (0 : ℤ) ∈ Eof s)
    (hwitness : ∀ s, witnessG2 s =
      (DenseCovers.slowμ
        (TournamentH7.GoodSet.goodSet (Eof s) ∩
          DenseCovers.safeSet (Pof s))).toReal)
    (hwide : ∀ s,
      (DenseCovers.slowμ (DenseCovers.coverSet (Eof s))).toReal ≤
        cap s - delta s)
    (hdual : ∀ s, cap s ≤
      (DenseCovers.slowμ (DenseCovers.safeSet (Pof s))).toReal) :
    ∀ s, delta s ≤ witnessG2 s := by
  intro s
  rw [hwitness s]
  exact DenseCovers.goodSet_witness_margin_from_wide_bound
    (Eof s) (Pof s) (cap s) (delta s) (hanchor s) (hwide s) (hdual s)

/-! ## Axiom audit -/

#print axioms measure_toReal_mono_of_subset
#print axioms shape_bonferroni_handoff
#print axioms shape_measure_mono_handoff
#print axioms shape_D_le_p0_handoff
#print axioms shape_bonferroni_safeSet_handoff
#print axioms shape_D_le_p0_denseCovers_handoff
#print axioms shape_bonferroni_goodSet_safeSet_handoff
#print axioms shape_goodSet_witness_pos_from_strict_cover_bound
#print axioms shape_goodSet_witness_margin_from_wide_bound

end EventMeasureBridge
end LRC14
end LonelyRunner
