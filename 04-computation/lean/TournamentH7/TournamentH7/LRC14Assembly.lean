/-
  TournamentH7.LRC14Assembly  (mac-mini-2026-07-02-S13)

  THE FINAL ASSEMBLY SURFACE.  One theorem, `lrc14_endgame`, stating LRC(14)
  from the two remaining analytic node-hypotheses of the skeleton, with the
  module map of what already discharges their INGREDIENTS documented inline.
  As modules land, each parameter gets replaced by a proof term and deleted
  from the signature; the file is DONE when the parameter list is empty.

  Current ingredient status (2026-07-02, evening):
    hfloor  (witness floor):  klein's LyWindowEnum (k=8..13 MACHINE-CHECKED
            native_decide) + kps LRCWitnessCert families + opus LRCLadderPack
            multi-cluster certificates — the covering censuses are [LEAN];
            the remaining glue is the case-split bookkeeping.
    hpartA  (G2 ⟹ reach):    CombPatterns (THM-605(i) characterization,
            avoidance + forced, sorry-free) + kps cert_ladder (module 6) +
            opus cert_ladder' (sharp budgets) + RatIntervalsWrap — the
            floors and the recursion glue are [LEAN]; the remaining piece is
            the fuel-checker soundness instantiation (module 6 → skeleton).

  This file is deliberately thin: it pins the TOP of the DAG so every push
  from any agent shrinks a visible, named surface.
-/
import Mathlib.Tactic
import TournamentH7.LRCFourteenSkeleton
import TournamentH7.CombPatterns
import TournamentH7.BonferroniTruncation

namespace TournamentH7.LRC14Assembly

open LonelyRunner.LRC14 LonelyRunner

/-- **The endgame statement.**  LRC(14) follows from the two remaining
skeleton nodes; everything else in the chain is already sorry-free.  The
parameters are EXACTLY the surface left to discharge. -/
theorem lrc14_endgame
    (hfloor : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) →
      (witnessMP : ℝ) ≤ witnessG2 (shapeOf v))
    (hpartA : ∀ v : Fin 13 → ℤ, 0 < witnessG2 (shapeOf v) →
      (1 : ℝ) / 14 ≤ Mreach v) :
    LRC14Statement :=
  lrc14_from_witness_floor hfloor hpartA

/-- Sanity: the endgame consumes the SAME statement shape the sieve layer
already discharges unconditionally on the non-covering branch — recorded here
so the final wiring has its target signature pinned. -/
example (v : Fin 13 → ℤ) (hdiv : ∀ i, ¬ ((14 : ℤ) ∣ v i)) :
    ∃ t : ℝ, LonelyRunner.Lonely 14 v t :=
  ⟨1 / 14, lrc14_no_multiple_of_14 v hdiv⟩

/-! ### The module-6 → Mreach instantiation bridge (the hpartA wiring shape).

A single rational witness time with clearance ≥ 1/14 bounds `Mreach` from
below — the exact shape the fuel checker's soundness theorem must produce.
This is the last interface between the certificate layer and the skeleton. -/

theorem Mreach_ge_of_witness (v : Fin 13 → ℤ) (t : ℝ)
    (ht : t ∈ Set.Icc (0 : ℝ) 1)
    (hcl : (1 : ℝ) / 14 ≤ LonelyRunner.LRC14Concrete.minReach v t) :
    (1 : ℝ) / 14 ≤ LonelyRunner.LRC14Concrete.Mreach v := by
  refine le_trans hcl (le_csSup ?_ ?_)
  · refine ⟨1 / 2, ?_⟩
    rintro y ⟨s, -, rfl⟩
    exact LonelyRunner.LRC14Concrete.minReach_le_half v s
  · exact Set.mem_image_of_mem _ ht

end TournamentH7.LRC14Assembly
