/-
  TournamentH7.LRCWitnessAttainmentBridge

  Connects the general `distZ`/`margin` witness-attainment formulation to the
  existing concrete `nearInt`/`Mreach` LRC(14) path.  This keeps the compactness
  handoff available in either interface without introducing a new analytic
  obligation.
-/

import TournamentH7.LRCMreachConcrete
import TournamentH7.LRCWitnessAttainment

open Set

namespace TournamentH7.LRCWitness

/-- The `infDist` distance to the integers is the same nearest-integer distance
used by the concrete LRC(14) `Mreach` module. -/
theorem distZ_eq_nearInt (x : ℝ) :
    distZ x = LonelyRunner.LRC14Concrete.nearInt x := by
  apply le_antisymm
  · unfold LonelyRunner.LRC14Concrete.nearInt
    apply le_min
    · calc
        distZ x ≤ dist x ((⌊x⌋ : ℤ) : ℝ) := by
          unfold distZ
          exact Metric.infDist_le_dist_of_mem ⟨⌊x⌋, rfl⟩
        _ = |x - ((⌊x⌋ : ℤ) : ℝ)| := Real.dist_eq _ _
        _ = Int.fract x := by
          have h : x - ((⌊x⌋ : ℤ) : ℝ) = Int.fract x := by
            linarith [Int.floor_add_fract x]
          rw [h, abs_of_nonneg (Int.fract_nonneg x)]
    · calc
        distZ x ≤ dist x (((⌊x⌋ + 1 : ℤ) : ℝ)) := by
          unfold distZ
          exact Metric.infDist_le_dist_of_mem ⟨⌊x⌋ + 1, rfl⟩
        _ = |x - (((⌊x⌋ + 1 : ℤ) : ℝ))| := Real.dist_eq _ _
        _ = 1 - Int.fract x := by
          have h : x - (((⌊x⌋ + 1 : ℤ) : ℝ)) = Int.fract x - 1 := by
            push_cast
            linarith [Int.floor_add_fract x]
          rw [h, abs_of_nonpos (by linarith [Int.fract_lt_one x])]
          ring
  · rw [le_distZ_iff]
    intro m
    exact (LonelyRunner.far_iff_fract (LonelyRunner.LRC14Concrete.nearInt x) x).2
      ⟨by
          unfold LonelyRunner.LRC14Concrete.nearInt
          exact min_le_left _ _,
        by
          have hmin : LonelyRunner.LRC14Concrete.nearInt x ≤ 1 - Int.fract x := by
            unfold LonelyRunner.LRC14Concrete.nearInt
            exact min_le_right _ _
          linarith⟩ m

/-- For thirteen runners, the general margin equals the concrete min-reach
function already used in the root LRC(14) skeleton. -/
theorem margin_eq_minReach (v : Fin 13 → ℤ) (t : ℝ) :
    margin v t = LonelyRunner.LRC14Concrete.minReach v t := by
  simp [margin, LonelyRunner.LRC14Concrete.minReach, distZ_eq_nearInt,
    Finset.inf'_univ_eq_ciInf]

/-- The concrete `Mreach` supremum is the supremum of the general margin over
one compact period. -/
theorem Mreach_eq_margin_sSup (v : Fin 13 → ℤ) :
    LonelyRunner.LRC14Concrete.Mreach v =
      sSup (margin v '' Icc (0 : ℝ) 1) := by
  unfold LonelyRunner.LRC14Concrete.Mreach
  congr 1
  ext y
  constructor
  · rintro ⟨t, ht, rfl⟩
    exact ⟨t, ht, margin_eq_minReach v t⟩
  · rintro ⟨t, ht, rfl⟩
    exact ⟨t, ht, (margin_eq_minReach v t).symm⟩

/-- Margin-supremum version of the LRC(14) compactness handoff.  The remaining
proof target is the analytic lower bound on this supremum. -/
theorem exists_lonely_of_margin_sSup_ge (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0)
    (hmargin : (1 : ℝ) / 14 ≤ sSup (margin v '' Icc (0 : ℝ) 1)) :
    ∃ t : ℝ, LonelyRunner.Lonely 14 v t := by
  have hM : (1 : ℝ) / 14 ≤ LonelyRunner.LRC14Concrete.Mreach v := by
    rwa [Mreach_eq_margin_sSup]
  exact LonelyRunner.LRC14Concrete.lonely_of_Mreach_ge v hv hM

/-! ## Axiom audit -/

#print axioms distZ_eq_nearInt
#print axioms margin_eq_minReach
#print axioms Mreach_eq_margin_sSup
#print axioms exists_lonely_of_margin_sSup_ge

end TournamentH7.LRCWitness
