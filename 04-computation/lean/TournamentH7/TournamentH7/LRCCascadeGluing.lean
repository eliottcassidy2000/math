/-
  TournamentH7.LRCCascadeGluing

  Closed-recurrence addendum to Klein S317's sorry-free formalization of
  Opus S333's cascade and local-density sampling lemmas.  The one-step
  measure bridge lives in `TournamentH7.CascadeGluing`; this module connects
  it to the suffix-debt recurrence proved in
  `TournamentH7.LRCLocalDensityBlockGluing`.
-/

import TournamentH7.CascadeGluing
import TournamentH7.LRCLocalDensityBlockGluing

namespace LRC14.CascadeGluing

open LRC14.LocalDensityBlockGluing

/-- The sharp unrolled composition ledger. At stage `index`, G1 gives
`eta index * (actual index - loss index)`; hence that stage pays
`eta index * loss index` before all later suffix densities. -/
theorem block_gluing_sharp_closed
    (initial : ℝ) (eta loss actual : ℕ → ℝ) (stageCount : ℕ)
    (hinitial : initial ≤ actual 0)
    (hetaNonnegative : ∀ index, 0 ≤ eta index)
    (hstep : ∀ index,
      eta index * (actual index - loss index) ≤ actual (index + 1)) :
    densityProduct eta stageCount * initial
        - ∑ index ∈ Finset.Ico 0 stageCount,
            (eta index * loss index) *
              ∏ later ∈ Finset.Ico (index + 1) stageCount, eta later
      ≤ actual stageCount := by
  rw [← lowerBound_eq_closed initial eta (fun index => eta index * loss index)]
  apply lowerBound_le_actual initial eta (fun index => eta index * loss index)
    actual hinitial hetaNonnegative
  intro index
  simpa [mul_sub] using hstep index

/-- Opus S333's published coarser ledger. If every local density is at most
one and every component loss is nonnegative, discard the leading `eta` on each
loss term in the sharp recurrence. -/
theorem block_gluing_closed
    (initial : ℝ) (eta loss actual : ℕ → ℝ) (stageCount : ℕ)
    (hinitial : initial ≤ actual 0)
    (hetaNonnegative : ∀ index, 0 ≤ eta index)
    (hetaAtMostOne : ∀ index, eta index ≤ 1)
    (hlossNonnegative : ∀ index, 0 ≤ loss index)
    (hstep : ∀ index,
      eta index * (actual index - loss index) ≤ actual (index + 1)) :
    densityProduct eta stageCount * initial
        - ∑ index ∈ Finset.Ico 0 stageCount,
            loss index * ∏ later ∈ Finset.Ico (index + 1) stageCount, eta later
      ≤ actual stageCount := by
  rw [← lowerBound_eq_closed initial eta loss]
  apply lowerBound_le_actual initial eta loss actual hinitial hetaNonnegative
  intro index
  calc
    eta index * actual index - loss index
        ≤ eta index * (actual index - loss index) :=
      fixedScale_weaker_loss (actual index) (loss index) (eta index)
        (hetaAtMostOne index) (hlossNonnegative index)
    _ ≤ actual (index + 1) := hstep index

/-! ## Axiom audit -/

#print axioms block_gluing_sharp_closed
#print axioms block_gluing_closed

end LRC14.CascadeGluing
