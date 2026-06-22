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

/-- Fractional parts commute with subtraction modulo `1`.  This is the algebraic
identity that converts phase-level empty arcs into speed-difference arcs. -/
theorem fract_sub_eq_fract_fract_sub_fract (u v : ℝ) :
    Int.fract (u - v) = Int.fract (Int.fract u - Int.fract v) := by
  rw [Int.fract_eq_fract]
  refine ⟨⌊u⌋ - ⌊v⌋, ?_⟩
  simp only [Int.fract]
  rw [Int.cast_sub]
  ring

/-- The phase-level empty-arc carrier from `LRCDenseCovers` lies inside the
concrete `goodSet`.  The only content is rewriting phase differences
`fract(b*x)-fract(a*x)` as speed differences `fract((b-a)*x)`. -/
theorem phaseGapSet_subset_goodSet (E : List ℤ) :
    LonelyRunner.DenseCovers.phaseGapSet E ⊆ goodSet E := by
  intro x hx
  rcases hx with ⟨a, haPhase, harc⟩
  rcases LonelyRunner.DenseCovers.mem_phaseFinset.mp haPhase with ⟨e, heE, heq⟩
  unfold goodSet
  simp only [Set.mem_iUnion, Set.mem_iInter, Set.mem_setOf_eq]
  refine ⟨e, by simpa using heE, ?_⟩
  intro b hbE
  have hbList : b ∈ E := by simpa using hbE
  have hbPhase :
      Int.fract ((b : ℝ) * x) ∈ LonelyRunner.DenseCovers.phaseFinset E x :=
    LonelyRunner.DenseCovers.mem_phaseFinset.mpr ⟨b, hbList, rfl⟩
  have h := harc (Int.fract ((b : ℝ) * x)) hbPhase
  rw [← heq] at h
  have harg : (((b - e : ℤ) : ℝ) * x) = (b : ℝ) * x - (e : ℝ) * x := by
    rw [Int.cast_sub]
    ring
  rw [harg, fract_sub_eq_fract_fract_sub_fract]
  exact h

/-- The complement of the dense event is contained in the concrete `GOOD`
carrier. -/
theorem denseSet_compl_subset_goodSet (E : List ℤ) :
    (LonelyRunner.DenseCovers.denseSet E)ᶜ ⊆ goodSet E :=
  subset_trans (LonelyRunner.DenseCovers.denseSet_compl_subset_phaseGapSet E)
    (phaseGapSet_subset_goodSet E)

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
#print axioms fract_sub_eq_fract_fract_sub_fract
#print axioms phaseGapSet_subset_goodSet
#print axioms denseSet_compl_subset_goodSet

end TournamentH7.GoodSet
