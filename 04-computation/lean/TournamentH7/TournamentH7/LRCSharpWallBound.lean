/- LRCSharpWallBound.lean -- opus-2026-07-17-S355 (HYP-7380).
   THE SHARP WALL BOUND: THM-1012's pair floors wired into the Hunter
   assembly.

   THE BRIDGE.  The assembly (LRCHunterAssembly, LRCSevenWallExistence) runs
   on a PROBABILITY space; the comb estimates (LRCCombUpperBound,
   LRCPairIndependence) live on `ℝ` with Lebesgue measure.  The join is
   `volume.restrict W` for a unit window `W`: `volume W = 1` makes it a
   probability measure, `Measure.restrict_apply` turns its values into
   `volume (· ∩ W)`, and — a bonus — the conclusion then places the lonely
   point INSIDE `W` rather than merely somewhere in `ℝ`.

   `sharp_wall_bound` takes exactly the two families of estimates this
   program produces (per-comb UPPER bounds summing to ≤ 1, per-consecutive-
   pair LOWER bounds) and returns a lonely point of `W`.  The λ = 1/14,
   seven-comb corollary is the wall itself: `7 · 2λ = 1` exactly, so the
   mass hypothesis is tight and any positive pair-floor sum suffices.
   Kernel-pure target: no sorry, no native_decide. -/
import Mathlib
import TournamentH7.LRCSevenWallExistence

open MeasureTheory
open scoped ENNReal

namespace LonelyRunner.LRC14.Hunter

/-- **THE SHARP WALL BOUND** (general form). -/
theorem sharp_wall_bound (W : Set ℝ) (hWmeas : MeasurableSet W)
    (hW : volume W = 1) (A : ℕ → Set ℝ) (hA : ∀ i, MeasurableSet (A i)) (n : ℕ)
    (hsum : ∑ i ∈ Finset.range n, volume (A i ∩ W) ≤ 1)
    (c : ℕ → ℝ≥0∞)
    (hlower : ∀ i ∈ Finset.Ico 1 n, c i ≤ volume (A i ∩ A (i - 1) ∩ W))
    (hpos : 0 < ∑ i ∈ Finset.Ico 1 n, c i) :
    ∃ t ∈ W, ∀ i ∈ Finset.range n, t ∉ A i := by
  haveI : IsProbabilityMeasure (volume.restrict W) :=
    ⟨by rw [Measure.restrict_apply_univ]; exact hW⟩
  set μ := volume.restrict W with hμ
  -- restrict turns the abstract measure into `volume (· ∩ W)`
  have hres : ∀ s : Set ℝ, MeasurableSet s → μ s = volume (s ∩ W) :=
    fun s hs => Measure.restrict_apply hs
  have hsum' : ∑ i ∈ Finset.range n, μ (A i) ≤ 1 := by
    rw [Finset.sum_congr rfl (fun i _ => hres (A i) (hA i))]
    exact hsum
  have hlower' : ∀ i ∈ Finset.Ico 1 n, c i ≤ μ (A i ∩ A (i - 1)) := by
    intro i hi
    rw [hres _ ((hA i).inter (hA (i - 1)))]
    exact hlower i hi
  -- the assembly: positive uncovered measure for μ
  have huncov : 0 < μ (⋃ i ∈ Finset.range n, A i)ᶜ :=
    lt_of_lt_of_le hpos
      (le_trans (Finset.sum_le_sum hlower')
        (uncovered_ge_overlaps_of_sum_le μ A hA n hsum'))
  -- unfold the restriction: the uncovered set meets W in positive measure
  have hmeasU : MeasurableSet (⋃ i ∈ Finset.range n, A i)ᶜ :=
    (MeasurableSet.biUnion (Set.to_countable _) (fun i _ => hA i)).compl
  rw [hres _ hmeasU] at huncov
  obtain ⟨t, htU, htW⟩ := uncovered_nonempty_of_pos volume _ huncov
  refine ⟨t, htW, fun i hi hmem => ?_⟩
  exact htU (Set.mem_biUnion hi hmem)

/-- **THE WALL, at `λ = 1/14` with seven combs.**  `7 · 2λ = 1` exactly, so
the mass hypothesis is tight: any positive pair-floor sum yields a lonely
point of the window. -/
theorem seven_comb_wall (W : Set ℝ) (hWmeas : MeasurableSet W)
    (hW : volume W = 1) (A : ℕ → Set ℝ) (hA : ∀ i, MeasurableSet (A i))
    (hupper : ∀ i ∈ Finset.range 7,
      volume (A i ∩ W) ≤ ENNReal.ofReal (2 * (1 / 14 : ℝ)))
    (c : ℕ → ℝ≥0∞)
    (hlower : ∀ i ∈ Finset.Ico 1 7, c i ≤ volume (A i ∩ A (i - 1) ∩ W))
    (hpos : 0 < ∑ i ∈ Finset.Ico 1 7, c i) :
    ∃ t ∈ W, ∀ i ∈ Finset.range 7, t ∉ A i := by
  refine sharp_wall_bound W hWmeas hW A hA 7 ?_ c hlower hpos
  calc ∑ i ∈ Finset.range 7, volume (A i ∩ W)
      ≤ ∑ _i ∈ Finset.range 7, ENNReal.ofReal (2 * (1 / 14 : ℝ)) :=
        Finset.sum_le_sum hupper
    _ = 7 * ENNReal.ofReal (2 * (1 / 14 : ℝ)) := by
        rw [Finset.sum_const, Finset.card_range, nsmul_eq_mul]
        norm_num
    _ = 1 := by
        rw [show (2 * (1 / 14 : ℝ)) = 1 / 7 by norm_num]
        rw [show ((7 : ENNReal)) = ENNReal.ofReal (7 : ℝ) by simp]
        rw [← ENNReal.ofReal_mul (by norm_num)]
        norm_num

/-! ## Axiom audit -/
#print axioms sharp_wall_bound
#print axioms seven_comb_wall

end LonelyRunner.LRC14.Hunter
