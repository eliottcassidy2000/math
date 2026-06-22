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

/-- The concrete lower carrier also lies in the complement of the dense event:
if the anchored phases fail to cover all six inner sectors, they cannot be
1/7-dense.  This is a formal proxy for the max-gap `GOOD` readout that avoids
sorted cyclic-gap infrastructure. -/
theorem witness_carrier_subset_dense_compl (E P : List ℤ) (hE : (0 : ℤ) ∈ E) :
    (coverSet E)ᶜ ∩ safeSet P ⊆ (denseSet E)ᶜ ∩ safeSet P :=
  Set.inter_subset_inter_left _ (coverSet_compl_subset_denseSet_compl E hE)

/-- Measure form of `witness_carrier_subset_dense_compl`: the Bonferroni carrier
lower-bounds the dense-complement witness carrier. -/
theorem witness_carrier_le_dense_compl_measure (E P : List ℤ) (hE : (0 : ℤ) ∈ E) :
    (slowμ ((coverSet E)ᶜ ∩ safeSet P)).toReal ≤
      (slowμ ((denseSet E)ᶜ ∩ safeSet P)).toReal :=
  ENNReal.toReal_mono (measure_ne_top slowμ ((denseSet E)ᶜ ∩ safeSet P))
    (measure_mono (witness_carrier_subset_dense_compl E P hE))

/-- **Witness-floor positivity reduces EXACTLY to the wide bound `p0 ≤ cap`.**
Given the wide bound `p0(E) ≤ cap_k` (`hwide`) and the duality strictness
`cap_k < meas(G_P)` (`hdual`, from `cap_k = min_P meas(G_P)`), the concrete witness
carrier has positive measure.  This isolates the *single* remaining analytic
obligation of the witness route's floor side — the SAME `p0 ≤ cap` the sector
route needs, so the floor adds no new obligation (HYP-2832 unification). -/
theorem witness_pos_from_wide_bound (E P : List ℤ) (capk : ℝ)
    (hwide : (slowμ (coverSet E)).toReal ≤ capk)
    (hdual : capk < (slowμ (safeSet P)).toReal) :
    0 < (slowμ ((coverSet E)ᶜ ∩ safeSet P)).toReal := by
  have hfloor := witness_floor_concrete E P
  linarith

/-- **Witness-floor positivity from a strict p0 cover bound.**  This is the
strictness placement produced by `LRCCoverBound.slowμ_coverSet_lt_cap`: if
`p0(E) < cap_k` and the dual floor only gives `cap_k ≤ meas(G_P)`, then the
concrete witness carrier still has positive measure. -/
theorem witness_pos_from_strict_cover_bound (E P : List ℤ) (capk : ℝ)
    (hwide : (slowμ (coverSet E)).toReal < capk)
    (hdual : capk ≤ (slowμ (safeSet P)).toReal) :
    0 < (slowμ ((coverSet E)ᶜ ∩ safeSet P)).toReal :=
  witness_pos_of_p0_lt_measGP E P (lt_of_lt_of_le hwide hdual)

/-- Strict-cover positivity transferred to the dense-complement witness carrier.
This plugs `LRCCoverBound.slowμ_coverSet_lt_cap` into a concrete max-gap proxy:
`(denseSet E)^c ∩ safeSet P`. -/
theorem dense_compl_witness_pos_from_strict_cover_bound
    (E P : List ℤ) (capk : ℝ) (hE : (0 : ℤ) ∈ E)
    (hwide : (slowμ (coverSet E)).toReal < capk)
    (hdual : capk ≤ (slowμ (safeSet P)).toReal) :
    0 < (slowμ ((denseSet E)ᶜ ∩ safeSet P)).toReal :=
  lt_of_lt_of_le (witness_pos_from_strict_cover_bound E P capk hwide hdual)
    (witness_carrier_le_dense_compl_measure E P hE)

/-- **Concrete margin version of the witness floor.**  The exact top-level p0
route needs a quantitative margin, not only positivity: if the wide bound proves
`p0(E) ≤ cap_k - delta` and the dual cap floor proves `cap_k ≤ meas(G_P)`, then
the concrete witness carrier has measure at least `delta`. -/
theorem witness_margin_from_wide_bound (E P : List ℤ) (capk delta : ℝ)
    (hwide : (slowμ (coverSet E)).toReal ≤ capk - delta)
    (hdual : capk ≤ (slowμ (safeSet P)).toReal) :
    delta ≤ (slowμ ((coverSet E)ᶜ ∩ safeSet P)).toReal := by
  have hfloor := witness_floor_concrete E P
  linarith

/-- Quantitative margin transferred to the dense-complement witness carrier. -/
theorem dense_compl_witness_margin_from_wide_bound
    (E P : List ℤ) (capk delta : ℝ) (hE : (0 : ℤ) ∈ E)
    (hwide : (slowμ (coverSet E)).toReal ≤ capk - delta)
    (hdual : capk ≤ (slowμ (safeSet P)).toReal) :
    delta ≤ (slowμ ((denseSet E)ᶜ ∩ safeSet P)).toReal :=
  le_trans (witness_margin_from_wide_bound E P capk delta hwide hdual)
    (witness_carrier_le_dense_compl_measure E P hE)

/-- Positive-margin corollary of `witness_margin_from_wide_bound`. -/
theorem witness_pos_from_wide_bound_margin (E P : List ℤ) (capk delta : ℝ)
    (hdelta : 0 < delta)
    (hwide : (slowμ (coverSet E)).toReal ≤ capk - delta)
    (hdual : capk ≤ (slowμ (safeSet P)).toReal) :
    0 < (slowμ ((coverSet E)ᶜ ∩ safeSet P)).toReal :=
  lt_of_lt_of_le hdelta
    (witness_margin_from_wide_bound E P capk delta hwide hdual)

/-- Positive-margin corollary for the dense-complement witness carrier. -/
theorem dense_compl_witness_pos_from_wide_bound_margin
    (E P : List ℤ) (capk delta : ℝ) (hE : (0 : ℤ) ∈ E)
    (hdelta : 0 < delta)
    (hwide : (slowμ (coverSet E)).toReal ≤ capk - delta)
    (hdual : capk ≤ (slowμ (safeSet P)).toReal) :
    0 < (slowμ ((denseSet E)ᶜ ∩ safeSet P)).toReal :=
  lt_of_lt_of_le hdelta
    (dense_compl_witness_margin_from_wide_bound E P capk delta hE hwide hdual)

/-! ## Axiom audit -/

#print axioms slowμ_compl_toReal
#print axioms witness_floor_concrete
#print axioms witness_pos_of_p0_lt_measGP
#print axioms witness_carrier_subset_dense_compl
#print axioms witness_carrier_le_dense_compl_measure
#print axioms witness_pos_from_wide_bound
#print axioms witness_pos_from_strict_cover_bound
#print axioms dense_compl_witness_pos_from_strict_cover_bound
#print axioms witness_margin_from_wide_bound
#print axioms dense_compl_witness_margin_from_wide_bound
#print axioms witness_pos_from_wide_bound_margin
#print axioms dense_compl_witness_pos_from_wide_bound_margin

end DenseCovers
end LonelyRunner
