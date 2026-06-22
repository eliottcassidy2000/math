/-
  TournamentH7.LRCGoodSet  (mac-mini-2026-06-22-S26)

  The CONCRETE "GOOD" event for the LRC(14) witness floor, with measurability —
  the carrier missing from the measure layer (kps built the G_P `safeSet` and
  the sector `coverSet`; this is the max-gap side).

  `GOOD(E) = { x : the phases {frac(e·x) : e ∈ E} have a cyclic gap > 1/7 }`.
  Characterization (clean + measurable): a cyclic gap exceeds 1/7 iff SOME phase
  `frac(a·x)` has the half-open arc `(frac(a·x), frac(a·x)+1/7]` free of phases,
  i.e. for all `b ∈ E`, `frac((b-a)·x) ∉ (0, 1/7]`.  This is a finite union of
  finite intersections of preimages of a Borel set under the (measurable) phase
  maps, hence measurable.  `nuShape s = μ (goodSet (cluster s))` then plugs into
  codex's `shape_bonferroni_handoff`; my pigeonhole (`LRCMaxGapPigeonhole`) gives
  `goodSet = univ` for cluster size ≤ 6 (the `hnu1` node).
-/
import TournamentH7.LRCDenseCovers

open MeasureTheory Set
open scoped Topology

namespace TournamentH7.GoodSet

/-- The concrete GOOD event: some phase `frac(a·x)` leaves the half-open
length-`1/7` arc just after it empty of phases. -/
def goodSet (E : List ℤ) : Set ℝ :=
  ⋃ a ∈ E.toFinset, ⋂ b ∈ E.toFinset,
    {x : ℝ | Int.fract (((b - a : ℤ) : ℝ) * x) ∉ Set.Ioc (0 : ℝ) (1 / 7)}

/-- Each atomic phase-arc condition is a measurable set. -/
lemma measurableSet_arc (c : ℤ) :
    MeasurableSet {x : ℝ | Int.fract ((c : ℝ) * x) ∉ Set.Ioc (0 : ℝ) (1 / 7)} := by
  have hpre : {x : ℝ | Int.fract ((c : ℝ) * x) ∉ Set.Ioc (0 : ℝ) (1 / 7)}
      = (fun x : ℝ => Int.fract ((c : ℝ) * x)) ⁻¹' (Set.Ioc (0 : ℝ) (1 / 7))ᶜ := by
    ext x; simp [Set.mem_preimage]
  rw [hpre]
  exact (LonelyRunner.DenseCovers.measurable_phase c) measurableSet_Ioc.compl

/-- **The GOOD event is measurable.** -/
theorem measurableSet_goodSet (E : List ℤ) : MeasurableSet (goodSet E) := by
  unfold goodSet
  refine Finset.measurableSet_biUnion _ (fun a _ => ?_)
  refine Finset.measurableSet_biInter _ (fun b _ => ?_)
  simpa using measurableSet_arc (b - a)

/-! ## Axiom audit -/

#print axioms measurableSet_arc
#print axioms measurableSet_goodSet

end TournamentH7.GoodSet
