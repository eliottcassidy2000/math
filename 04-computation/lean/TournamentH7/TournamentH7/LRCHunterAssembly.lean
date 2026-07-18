/- LRCHunterAssembly.lean -- opus-2026-07-17-S348 (HYP-7310 / THM-964 top wire).
   THE HUNTER PATH-TREE ASSEMBLY.  On a probability space, if `n` sets have
   total measure 1, the UNCOVERED measure dominates the sum of consecutive
   (path-tree) overlaps:

     ∑_{i=1}^{n-1} μ(A_i ∩ A_{i-1})  ≤  μ((⋃_i A_i)ᶜ).

   This is boxeph's `tree_hunter_add_le` (parent p i = i−1) + `measure_compl`.
   With per-edge overlap floors `c_i ≤ μ(A_i ∩ A_{i-1})` it gives
   `∑ c_i ≤ μ(uncovered)`, and `LRCWindowAverage.live_window_exists` turns a
   positive uncovered measure into a live window.

   For THM-964 the sets are the 7 danger arcs `badArcs x_i (1/14)` on the unit
   circle: `∑ μ(A_i) = 7·(1/7) = 1` (each arc has measure `2λ = 1/7`) and
   `μ(A_i ∩ A_{i-1}) ≥ floor` from `LRCFloorTable.muNum_lower` THROUGH the
   sawtooth-sum = measure bridge (`muNum a b · g/(14ab) = μ(D_a ∩ D_b)`,
   THM-856/865 — the ONE piece not yet in Lean; named here as the reduction
   target `OverlapMeasureBridge`).  Everything above the bridge is kernel-pure.
   No sorry, no native_decide in this file. -/
import Mathlib
import TournamentH7.LRCTreeHunter

open MeasureTheory
open scoped ENNReal

namespace LonelyRunner.LRC14.Hunter

variable {α : Type*} [MeasurableSpace α]

/-- **THE HUNTER PATH-TREE FLOOR** (kernel-pure): total measure 1 forces the
uncovered set to dominate the consecutive-overlap sum. -/
theorem uncovered_ge_consecutive_overlaps (μ : Measure α) [IsProbabilityMeasure μ]
    (A : ℕ → Set α) (hA : ∀ i, MeasurableSet (A i)) (n : ℕ)
    (hsum : ∑ i ∈ Finset.range n, μ (A i) = 1) :
    ∑ i ∈ Finset.Ico 1 n, μ (A i ∩ A (i - 1))
      ≤ μ (⋃ i ∈ Finset.range n, A i)ᶜ := by
  have htree := tree_hunter_add_le μ A (fun i => i - 1)
    (fun i hi => Nat.sub_lt (by omega) Nat.one_pos) hA n
  rw [hsum] at htree
  have hmeas : MeasurableSet (⋃ i ∈ Finset.range n, A i) :=
    MeasurableSet.biUnion (Set.to_countable _) (fun i _ => hA i)
  rw [measure_compl hmeas (measure_ne_top μ _), measure_univ]
  exact ENNReal.le_sub_of_add_le_left (measure_ne_top μ _) htree

/-- **THE FLOOR-SUM COROLLARY** (kernel-pure): per-edge floors sum to an
uncovered lower bound. -/
theorem uncovered_ge_floor_sum (μ : Measure α) [IsProbabilityMeasure μ]
    (A : ℕ → Set α) (hA : ∀ i, MeasurableSet (A i)) (n : ℕ)
    (hsum : ∑ i ∈ Finset.range n, μ (A i) = 1)
    (c : ℕ → ℝ≥0∞) (hc : ∀ i ∈ Finset.Ico 1 n, c i ≤ μ (A i ∩ A (i - 1))) :
    ∑ i ∈ Finset.Ico 1 n, c i ≤ μ (⋃ i ∈ Finset.range n, A i)ᶜ :=
  le_trans (Finset.sum_le_sum hc) (uncovered_ge_consecutive_overlaps μ A hA n hsum)

/-- **THE POSITIVE-UNCOVERED COROLLARY** (kernel-pure): a strictly positive
floor sum forces positive uncovered measure — the input to `live_window_exists`. -/
theorem uncovered_pos_of_floor_sum (μ : Measure α) [IsProbabilityMeasure μ]
    (A : ℕ → Set α) (hA : ∀ i, MeasurableSet (A i)) (n : ℕ)
    (hsum : ∑ i ∈ Finset.range n, μ (A i) = 1)
    (c : ℕ → ℝ≥0∞) (hc : ∀ i ∈ Finset.Ico 1 n, c i ≤ μ (A i ∩ A (i - 1)))
    (hpos : 0 < ∑ i ∈ Finset.Ico 1 n, c i) :
    0 < μ (⋃ i ∈ Finset.range n, A i)ᶜ :=
  lt_of_lt_of_le hpos (uncovered_ge_floor_sum μ A hA n hsum c hc)

/-- **THE 7-BLOCK CAPSTONE** (kernel-pure): the THM-964 shape.  On a
probability space, seven measurable sets of total measure 1 whose six
consecutive overlaps have a strictly positive floor sum leave a
POSITIVELY-uncovered set.  For the wall, `μ` is the circle's Haar measure,
`A i = badArcs x_i (1/14)`, `∑ μ(A i) = 7·(1/7) = 1`, and the floors are
`LRCFloorTable.muNum_lower` through the overlap-measure bridge. -/
theorem seven_block_uncovered_pos (μ : Measure α) [IsProbabilityMeasure μ]
    (A : ℕ → Set α) (hA : ∀ i, MeasurableSet (A i))
    (hsum : ∑ i ∈ Finset.range 7, μ (A i) = 1)
    (c : ℕ → ℝ≥0∞) (hc : ∀ i ∈ Finset.Ico 1 7, c i ≤ μ (A i ∩ A (i - 1)))
    (hpos : 0 < ∑ i ∈ Finset.Ico 1 7, c i) :
    0 < μ (⋃ i ∈ Finset.range 7, A i)ᶜ :=
  uncovered_pos_of_floor_sum μ A hA 7 hsum c hc hpos

/-!
## The reduction of THM-964 to two named bridges

`seven_block_uncovered_pos` (+ `LRCWindowAverage.live_window_exists`) closes
the 7-wall from these two inputs, the ONLY pieces of THM-964 not yet in Lean:

* **`OverlapMeasureBridge`** — `μ (badArcs a λ ∩ badArcs b λ) = muNum a b · g/(14ab)`
  and `μ (badArcs x λ) = 2λ`, the sawtooth-sum = measure identity (THM-856/865).
  This feeds `LRCFloorTable.muNum_lower` into `hsum`/`hc` above.
* **`CircleLineReconcile`** — the danger sets live on the unit circle
  (a probability space, where the Hunter assembly runs); the window `[x,x+L]`
  and `live_window_exists` run on `ℝ` (Lebesgue).  Identify `[0,1) ⊂ ℝ` with
  the circle so the positively-uncovered circle set maps to a positive-measure
  subset of `ℝ` inside a period.

Every other ingredient — the tree-Hunter inequality (`tree_hunter_add_le`),
the per-edge floor (`muNum_lower`, `overlap_floor_rat`), the allocation law
(`sorted_ratio_pow_le`), this assembly, and the window average
(`window_average`, `live_window_exists`) — is kernel-pure.
-/

/-! ## Axiom audit -/
#print axioms uncovered_ge_consecutive_overlaps
#print axioms uncovered_ge_floor_sum
#print axioms uncovered_pos_of_floor_sum
#print axioms seven_block_uncovered_pos

end LonelyRunner.LRC14.Hunter
