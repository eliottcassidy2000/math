/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: opus (LRC multi-agent project, 2026-07-02-S41)
-/
import TournamentH7.RatIntervals
import TournamentH7.RatIntervalsWrap

/-!
# Module 3, ℚ frame: the per-instance commensuration schema (the frame-audit residue)

`LRCCommensuration.lean` (opus-S34) proves the 7-commensuration lemma over ℝ with Haar
measure, for ALL real phases.  The all-ℚ critical path (playbook T1) consumes commensuration
only at RATIONAL phases fixed by certificate data (opus-S38 interface note), so the ℚ-frame
form is a decidable schema, not a measure theorem:

    overlapQ P Q r ψ φ := length (inter (wrap (comb P r ψ)) (wrap (comb Q r φ)))

computes the circle overlap of two wrapped combs as a pure list value; each needed row is one
`native_decide`.  The ℝ theorem `seven_commensuration` remains the mathematical statement
(every phase); this schema is its kernel-checkable trace on the rows module 4 ingests.
Four commensurate rows and a non-commensurate control are provided as the pattern.

The general ∀-phase ℚ theorem (seven-translate exact partition — half-open intervals tile
EXACTLY, no measure theory) remains available as a follow-up via the `Chain` machinery if
module 4 ever quantifies over phases; see the S38/S39 design logs.
-/

namespace LonelyRunner
namespace RatIntervals

/-- The circle overlap of the `P`- and `Q`-combs at radius `r` and phases `ψ, φ`, as an exact
rational list computation (both combs wrapped into `[0,1)`). -/
def overlapQ (P Q : ℕ) (r ψ φ : ℚ) : ℚ :=
  length (inter (wrap (comb P r ψ)) (wrap (comb Q r φ)))

/-! ### Commensurate rows: `7 ∣ Q`, `7 ∤ P` ⟹ overlap `= (2r)² = 1/49` at `r = 1/14`.
Each row is one kernel-checkable evaluation (`native_decide`; ℚ normalization blocks kernel
`decide` via `Nat.gcd`). -/

theorem overlapQ_1_7 : overlapQ 1 7 (1/14) (1/3) (1/5) = 1/49 := by native_decide

theorem overlapQ_2_7 : overlapQ 2 7 (1/14) (2/7) (3/11) = 1/49 := by native_decide

theorem overlapQ_3_14 : overlapQ 3 14 (1/14) 0 (1/2) = 1/49 := by native_decide

theorem overlapQ_5_21 : overlapQ 5 21 (1/14) (4/9) (1/13) = 1/49 := by native_decide

/-- Non-commensurate control: `(P, Q) = (1, 2)` at aligned phases gives `1/14 ≠ 1/49` — the
schema distinguishes, it does not assume. -/
theorem overlapQ_1_2_control : overlapQ 1 2 (1/14) 0 0 = 1/14 := by native_decide

end RatIntervals
end LonelyRunner
