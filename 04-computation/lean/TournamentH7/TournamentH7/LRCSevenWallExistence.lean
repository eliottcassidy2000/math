/- LRCSevenWallExistence.lean -- opus-2026-07-17-S349 (HYP-7320).
   THE EXISTENCE CAPSTONE, and an architectural simplification.

   S348 named two bridges as THM-964's remaining Lean surface:
     (1) OverlapMeasureBridge  (muNum = arc-overlap measure, THM-856/865);
     (2) CircleLineReconcile   (the circle probability space vs R-Lebesgue,
         needed to hand a positively-uncovered circle set to the R-based
         window_average / live_window_exists).

   OBSERVATION (S349): bridge (2) is NOT on the critical path for the wall's
   CONCLUSION.  LRC needs a lonely TIME -- a single point at which every
   runner is far -- not a window of positive measure.  Positive uncovered
   measure already gives a NONEMPTY uncovered set, and any point of it is a
   lonely time.  The window machinery is needed only to NEST a further block
   inside (the cascade/gluing architecture), never to conclude.

   So the 7-wall's existence conclusion rests on bridge (1) alone.

   Kernel-pure: no sorry, no native_decide. -/
import Mathlib
import TournamentH7.LRCHunterAssembly

open MeasureTheory
open scoped ENNReal

namespace LonelyRunner.LRC14.Hunter

variable {α : Type*} [MeasurableSpace α]

/-- Positive measure gives a point: the uncovered set is inhabited. -/
theorem uncovered_nonempty_of_pos (μ : Measure α) (S : Set α) (h : 0 < μ S) :
    S.Nonempty :=
  nonempty_of_measure_ne_zero (ne_of_gt h)

/-- **THE EXISTENCE CAPSTONE** (kernel-pure): on a probability space, seven
measurable sets of total measure 1 whose six consecutive path-tree overlaps
carry a strictly positive floor sum leave an INHABITED uncovered set — a point
outside every one of them.  For the wall: `A i = badArcs x_i (1/14)`, and the
produced point is a time at which all seven runners are `1/14`-lonely. -/
theorem seven_block_lonely_point (μ : Measure α) [IsProbabilityMeasure μ]
    (A : ℕ → Set α) (hA : ∀ i, MeasurableSet (A i))
    (hsum : ∑ i ∈ Finset.range 7, μ (A i) = 1)
    (c : ℕ → ℝ≥0∞) (hc : ∀ i ∈ Finset.Ico 1 7, c i ≤ μ (A i ∩ A (i - 1)))
    (hpos : 0 < ∑ i ∈ Finset.Ico 1 7, c i) :
    ∃ t, ∀ i ∈ Finset.range 7, t ∉ A i := by
  obtain ⟨t, ht⟩ := uncovered_nonempty_of_pos μ _
    (seven_block_uncovered_pos μ A hA hsum c hc hpos)
  refine ⟨t, fun i hi => ?_⟩
  intro hmem
  exact ht (Set.mem_biUnion hi hmem)

/-- The same conclusion for a general block size `n`. -/
theorem block_lonely_point (μ : Measure α) [IsProbabilityMeasure μ]
    (A : ℕ → Set α) (hA : ∀ i, MeasurableSet (A i)) (n : ℕ)
    (hsum : ∑ i ∈ Finset.range n, μ (A i) = 1)
    (c : ℕ → ℝ≥0∞) (hc : ∀ i ∈ Finset.Ico 1 n, c i ≤ μ (A i ∩ A (i - 1)))
    (hpos : 0 < ∑ i ∈ Finset.Ico 1 n, c i) :
    ∃ t, ∀ i ∈ Finset.range n, t ∉ A i := by
  obtain ⟨t, ht⟩ := uncovered_nonempty_of_pos μ _
    (uncovered_pos_of_floor_sum μ A hA n hsum c hc hpos)
  refine ⟨t, fun i hi => ?_⟩
  intro hmem
  exact ht (Set.mem_biUnion hi hmem)

/-! ### The `≤ 1` generalization

The assembly's hypothesis `∑ μ (A i) = 1` can be weakened to `≤ 1`: Hunter
gives `μ(⋃) + ∑ overlaps ≤ ∑ μ(A i) ≤ 1`, so `uncovered = 1 - μ(⋃) ≥ ∑
overlaps` all the same.  This matters for the bridge: only an UPPER bound on
the single-comb measure (`μ (badArcs x λ) ≤ 2λ`, the fragmentation direction
already in the corpus) is needed — not the exact value. -/

/-- Hunter path-tree floor under the weakened total-mass hypothesis. -/
theorem uncovered_ge_overlaps_of_sum_le (μ : Measure α) [IsProbabilityMeasure μ]
    (A : ℕ → Set α) (hA : ∀ i, MeasurableSet (A i)) (n : ℕ)
    (hsum : ∑ i ∈ Finset.range n, μ (A i) ≤ 1) :
    ∑ i ∈ Finset.Ico 1 n, μ (A i ∩ A (i - 1))
      ≤ μ (⋃ i ∈ Finset.range n, A i)ᶜ := by
  have htree := tree_hunter_add_le μ A (fun i => i - 1)
    (fun i hi => Nat.sub_lt (by omega) Nat.one_pos) hA n
  have htree' : μ (⋃ i ∈ Finset.range n, A i)
      + ∑ i ∈ Finset.Ico 1 n, μ (A i ∩ A (i - 1)) ≤ 1 := le_trans htree hsum
  have hmeas : MeasurableSet (⋃ i ∈ Finset.range n, A i) :=
    MeasurableSet.biUnion (Set.to_countable _) (fun i _ => hA i)
  rw [measure_compl hmeas (measure_ne_top μ _), measure_univ]
  exact ENNReal.le_sub_of_add_le_left (measure_ne_top μ _) htree'

/-- **THE CAPSTONE UNDER `≤ 1`**: an inhabited uncovered set from a positive
overlap floor sum and total mass at most 1. -/
theorem block_lonely_point_of_sum_le (μ : Measure α) [IsProbabilityMeasure μ]
    (A : ℕ → Set α) (hA : ∀ i, MeasurableSet (A i)) (n : ℕ)
    (hsum : ∑ i ∈ Finset.range n, μ (A i) ≤ 1)
    (c : ℕ → ℝ≥0∞) (hc : ∀ i ∈ Finset.Ico 1 n, c i ≤ μ (A i ∩ A (i - 1)))
    (hpos : 0 < ∑ i ∈ Finset.Ico 1 n, c i) :
    ∃ t, ∀ i ∈ Finset.range n, t ∉ A i := by
  have huncov : 0 < μ (⋃ i ∈ Finset.range n, A i)ᶜ :=
    lt_of_lt_of_le hpos
      (le_trans (Finset.sum_le_sum hc)
        (uncovered_ge_overlaps_of_sum_le μ A hA n hsum))
  obtain ⟨t, ht⟩ := uncovered_nonempty_of_pos μ _ huncov
  refine ⟨t, fun i hi => ?_⟩
  intro hmem
  exact ht (Set.mem_biUnion hi hmem)

/-! ## Axiom audit -/
#print axioms uncovered_nonempty_of_pos
#print axioms seven_block_lonely_point
#print axioms block_lonely_point
#print axioms uncovered_ge_overlaps_of_sum_le
#print axioms block_lonely_point_of_sum_le

end LonelyRunner.LRC14.Hunter
