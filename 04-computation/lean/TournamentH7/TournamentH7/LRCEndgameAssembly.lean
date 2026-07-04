/-
  TournamentH7.LRCEndgameAssembly — LRC(14) ⟸ LRC(13) + THE COMPRESSED CASE
  (kind-pasteur-2026-07-04-S6, HYP-4091).

  The covering dispatch (opus `covering_lonely_of_dominant_or_compressed`) splits every covering
  family into DOMINANT (some runner > 13× the rest) and COMPRESSED (no dominant runner).  The
  dominant branch is now discharged by the sharp dominant peel (kps-S5 `hdom_closed_abs`, HYP-4087),
  self-contained from the LRC(13) citation.  Composing with the non-covering sieve
  (`lrc14_of_covering_lonely`) gives the honest remaining surface:

     LRC(14)  ⟸  LRC(13) citation  +  [ every COMPRESSED covering family is lonely ].

  So the entire open obligation of the LRC(14) formalization is now exactly `hcomp` — the compressed
  (census / renormalization / confinement) core.  The far-runner tail is gone.

  Kernel-pure (propext, Classical.choice, Quot.sound); no sorry, no native_decide.
-/
import Mathlib
import TournamentH7.LRC14CertRoute
import TournamentH7.LRCHlargeRoute
import TournamentH7.LRCDominantPeel
import TournamentH7.LRC13Citation

namespace LonelyRunner
namespace LRC14

/-- **hdom discharged.**  The dominant obligation of the covering dispatch is closed by the sharp
dominant peel, from the LRC(13) citation — no extra hypothesis. -/
theorem hdom_discharged (cite : LonelyRunner.LRCUpTo13) :
    ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → CoveringFamily v →
      (∃ i, ∀ j, j ≠ i → 13 * |v j| < |v i|) → ∃ t : ℝ, Lonely 14 v t := by
  intro v hv _ hdomex
  obtain ⟨i, hi⟩ := hdomex
  exact hdom_closed_abs cite v hv i hi

/-- **THE ENDGAME REDUCTION: `LRC(14) ⟸ LRC(13) + compressed`.**  With the dominant branch discharged
by the sharp peel and the non-covering families handled by the denominator sieve, the WHOLE LRC(14)
statement reduces to loneliness of COMPRESSED covering families (every runner within `13×` of some
other).  This is the exact remaining surface — the census / renormalization / confinement core. -/
theorem lrc14_of_compressed (cite : LonelyRunner.LRCUpTo13)
    (hcomp : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → CoveringFamily v →
      (∀ i, ∃ j, j ≠ i ∧ |v i| ≤ 13 * |v j|) → ∃ t : ℝ, Lonely 14 v t) :
    LRC14Statement :=
  lrc14_of_covering_lonely
    (HlargeRoute.covering_lonely_of_dominant_or_compressed (hdom_discharged cite) hcomp)

#print axioms hdom_discharged
#print axioms lrc14_of_compressed

end LRC14
end LonelyRunner
