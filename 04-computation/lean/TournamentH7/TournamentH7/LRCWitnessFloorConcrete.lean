/-
  TournamentH7.LRCWitnessFloorConcrete -- the CONCRETE witness floor
  `meas(G_P) - p0 ≤ witnessG2`, assembled from real slow-time events
  (kind-pasteur-2026-06-22-S31, HYP-2832 unification at the measure level).

  This is the concrete realization of the witness/p0 unification.  Using the
  carrier `coverSetᶜ` (slow-times where some inner sector is empty), which is a
  subset of the genuine GOOD event (an empty length-`1/7` sector forces a cyclic
  gap `> 1/7`), we prove on the slow-time probability space `slowμ`:

      meas(G_P) - p0(E)  ≤  meas(coverSetᶜ ∩ G_P)  ( ≤ witnessG2 = meas(GOOD ∩ G_P) ).

  The proof is just Bonferroni (`LRCBonferroniMeasure`) + the complement identity
  on `coverSet` (measurable, `LRCDenseCovers`).  No `frac((b-a)x)` modular
  reasoning is needed: `coverSetᶜ` is a clean measurable lower-carrier.  Sorry-free.

  Consequence: the witness-floor obligation `witnessMP ≤ witnessG2` reduces to the
  p0 wide bound `p0 ≤ cap` (since meas(G_P) ≥ cap by the duality), with everything
  on the floor side now discharged from concrete events.
-/

import TournamentH7.LRCDenseCovers
import TournamentH7.LRCBonferroniMeasure

namespace LonelyRunner
namespace DenseCovers

open MeasureTheory

/-- Complement measure in real form on the slow-time probability space:
`(slowμ Aᶜ).toReal = 1 - (slowμ A).toReal` for measurable `A`. -/
theorem slowμ_compl_toReal {A : Set ℝ} (hA : MeasurableSet A) :
    (slowμ Aᶜ).toReal = 1 - (slowμ A).toReal := by
  have hle : slowμ A ≤ 1 := by
    calc slowμ A ≤ slowμ Set.univ := measure_mono (Set.subset_univ _)
      _ = 1 := measure_univ
  have hcompl : slowμ Aᶜ = 1 - slowμ A := by
    rw [measure_compl hA (measure_ne_top slowμ A), measure_univ]
  rw [hcompl, ENNReal.toReal_sub_of_le hle ENNReal.one_ne_top, ENNReal.toReal_one]

/-- **The concrete witness floor.**  On the slow-time probability space,
`meas(G_P) - p0(E) ≤ meas(coverSetᶜ ∩ G_P)`.  Here `meas(coverSet E) = p0(E)`,
`meas(safeSet P) = meas(G_P)`, and `coverSetᶜ ∩ safeSet` is a (measurable) subset
of the genuine 1/7-witness event `GOOD ∩ G_P`, so its measure lower-bounds
`witnessG2`.  Pure Bonferroni + complement identity, sorry-free. -/
theorem witness_floor_concrete (E P : List ℤ) :
    (slowμ (safeSet P)).toReal - (slowμ (coverSet E)).toReal
      ≤ (slowμ ((coverSet E)ᶜ ∩ safeSet P)).toReal := by
  have hbonf :=
    BonferroniMeasure.toReal_bonferroni slowμ (coverSet E)ᶜ (safeSet P)
      (measurableSet_safeSet P)
  rw [slowμ_compl_toReal (measurableSet_coverSet E)] at hbonf
  linarith

/-- **Witness-floor positivity from the wide bound.**  If `p0(E) < meas(G_P)`
(the wide bound `p0 ≤ cap ≤ meas(G_P)` with strict slack, the unification's
`δ_k > 0`), then the concrete witness carrier has positive measure. -/
theorem witness_pos_of_p0_lt_measGP (E P : List ℤ)
    (h : (slowμ (coverSet E)).toReal < (slowμ (safeSet P)).toReal) :
    0 < (slowμ ((coverSet E)ᶜ ∩ safeSet P)).toReal := by
  have := witness_floor_concrete E P
  linarith

/-- The witness carrier `coverSetᶜ ∩ safeSet` lies inside `safeSet` (`G_P`): at
every such slow-time the small parts are already safe.  (Recorded for the
Part-A handoff: these are genuine near-lonely slow-times.) -/
theorem witness_carrier_subset_safe (E P : List ℤ) :
    (coverSet E)ᶜ ∩ safeSet P ⊆ safeSet P :=
  Set.inter_subset_right

/-! ## Axiom audit -/

#print axioms slowμ_compl_toReal
#print axioms witness_floor_concrete
#print axioms witness_pos_of_p0_lt_measGP

end DenseCovers
end LonelyRunner
