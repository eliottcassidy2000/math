/-
  TournamentH7.LRCC8ConsecutiveWitness

  Proof-facing consumer for `LRCC8Consecutive`.  The upstream measure theorem
  says that eight consecutive danger combs leave positive restricted measure
  in the unit window.  This module extracts an actual point and packages it as
  a literal `Lonely 14` witness for the eight-speed family.

  The faithful carrier remains the ordered path of danger combs with numeric
  overlap credits.  This consumer intentionally adds no graph or tournament
  quotient: those would discard the quantitative Hunter certificate.

  No `sorry`; no `native_decide`.
-/

import TournamentH7.LRCC8Consecutive
import TournamentH7.LonelyRunner

open MeasureTheory Set

namespace LonelyRunner.LRC14.C8ConsecutiveWitness

open Arcs

noncomputable section

/-- Positive restricted good measure gives a point in the unit window that
avoids every danger tooth of all eight consecutive speeds. -/
theorem exists_eight_consecutive_good_time (base : ℕ) (hbase : 1 ≤ base) :
    ∃ time ∈ Ioo (-(1 : ℝ) / 2) (1 / 2),
      ∀ index < 8, ∀ witness : ℤ,
        (1 / 14 : ℝ) ≤
          |((base + index : ℕ) : ℝ) * time - witness| := by
  let dangerUnion : Set ℝ :=
    ⋃ index ∈ Finset.range 8, dangerR (base + index)
  have hmeasurable : MeasurableSet dangerUnion := by
    unfold dangerUnion
    exact Finset.measurableSet_biUnion _ fun index _ =>
      (dangerR_isOpen (base + index)).measurableSet
  have hpositive := c8_consecutive_good_pos base hbase
  have hintersection :
      0 < volume
        (dangerUnionᶜ ∩ Ioo (-(1 : ℝ) / 2) (1 / 2)) := by
    rw [Measure.restrict_apply hmeasurable.compl] at hpositive
    exact hpositive
  have hnonempty :
      (dangerUnionᶜ ∩ Ioo (-(1 : ℝ) / 2) (1 / 2)).Nonempty :=
    nonempty_of_measure_ne_zero (ne_of_gt hintersection)
  obtain ⟨time, htimeGood, htimeWindow⟩ := hnonempty
  refine ⟨time, htimeWindow, ?_⟩
  intro index hindex witness
  by_contra hclose
  push Not at hclose
  apply htimeGood
  simp only [dangerUnion, Set.mem_iUnion, Finset.mem_range]
  exact ⟨index, hindex, witness, hclose⟩

/-- Every positive block of eight consecutive integer speeds has an explicit
`Lonely 14` time. -/
theorem exists_lonely14_eight_consecutive (base : ℕ) (hbase : 1 ≤ base) :
    ∃ time : ℝ,
      LonelyRunner.Lonely 14
        (fun index : Fin 8 => ((base + (index : ℕ) : ℕ) : ℤ)) time := by
  obtain ⟨time, _htimeWindow, hgood⟩ :=
    exists_eight_consecutive_good_time base hbase
  refine ⟨time, ?_⟩
  intro index witness
  simpa only [Nat.cast_add, Nat.cast_ofNat, Int.cast_natCast, Int.cast_add] using
    hgood (index : ℕ) index.isLt witness

#print axioms exists_eight_consecutive_good_time
#print axioms exists_lonely14_eight_consecutive

end

end LonelyRunner.LRC14.C8ConsecutiveWitness
