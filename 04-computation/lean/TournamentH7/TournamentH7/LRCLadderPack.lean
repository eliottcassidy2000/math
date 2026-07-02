Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: opus (LRC multi-agent project, 2026-07-02-S39)
-/
import TournamentH7.LRCWitnessCert

/-!
# The first depth-3 ladder pack row (THM-606/607 through module 6)

Instantiates kps's `cert_ladder` (module 6 core, HYP-3963) with the first multi-cluster
13-runner family: `{1, 2} ∪ {V₁ − o : o ∈ {0,1,2}} ∪ {V₂ − o : o ∈ {0,1,3}} ∪
{V₃ − o : o ∈ {0,1,2,4,7}}` at `(V₁, V₂, V₃) = (50, 2200, 100000)`, uniform band
`h + μ = 1/14 + 1/40`, window `[7/20, 3/8]`.  Everything is kernel `decide`.

Parameter note (feedback to module 6): the uniform-`μ` design needs consecutive ratios
`> 1 + 1/μ ≈ 41`; the S36 per-level-`μ` tuple `(50, 2000, 90000)` (ratios 40, 45) fails its
level-2 separation here — per-level budgets (THM-606's `μ_ℓ`, with `μ_d = 0`) are strictly
stronger.  A `cert_ladder'` with per-level `μ` would accept every THM-607 sharp-region tuple.
-/

namespace LonelyRunner
namespace WitnessCert

/-- The three ladder levels of the first depth-3 pack row. -/
def packLevels : List Level :=
  [⟨[0, 1, 2], 27/280, 50⟩, ⟨[0, 1, 3], 31/140, 2200⟩, ⟨[0, 1, 2, 4, 7], 27/280, 100000⟩]

/-- **The first multi-cluster certified 13-runner family instance in Lean**: all thirteen
runners are `1/14`-lonely at a common time. -/
theorem depth3_pack_row :
    ∃ τ : ℝ,
      (∀ s ∈ ([1, 2] : List ℤ), ((1/14 : ℚ) : ℝ) ≤ ‖(((s : ℝ) * τ : ℝ) : UnitAddCircle)‖) ∧
      ∀ L ∈ packLevels, ∀ o ∈ L.offs,
        ((1/14 : ℚ) : ℝ) ≤ ‖((((L.V - o : ℤ) : ℝ) * τ : ℝ) : UnitAddCircle)‖ :=
  cert_ladder (h := 1/14) (μ := 1/40) (by norm_num) (by norm_num)
    (lo := 7/20) (hi := 3/8) [1, 2] packLevels
    (by decide) (by native_decide) (by norm_num) (by native_decide) (by native_decide)

/-! ## The quantified 3-parameter FAMILY (kind-pasteur S7, HYP-3964): variable V-tuples via the symbolic SepChain thresholds. -/


end WitnessCert
end LonelyRunner
