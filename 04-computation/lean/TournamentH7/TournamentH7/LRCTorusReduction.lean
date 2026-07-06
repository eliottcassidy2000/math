/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: opus (LRC multi-agent project, 2026-07-06-S101)
-/
import TournamentH7.LRCTorusProjection

/-!
# (A) ⟸ (C): a loose direction makes the torus loose (HYP-4336)

The (A) leg (no coupled proper 2-subtorus has `M(U) ∈ (1/13, 2/25)`) reduces to the (C)
leg (the 1-D 12-runner Farey gap) plus a pigeonhole rigidity — see
`03-artifacts/drafts/the-A-leg-reduces-to-the-C-leg-opus-S101.md`.  The chain is:

  rank-2  --(pigeonhole rigidity)-->  a full-support direction with a NON-tight projection
          --(1-D gap = hdich/(C))-->  that projection is loose (has a `2/25`-margin point)
          --(projection floor, S99)-->  the torus has a `2/25`-margin point  ⟹  `M(U) ≥ 2/25`.

This file supplies the final arrow as a consumable lemma: **one loose direction ⟹ a loose
torus**, directly from `torus_point_of_projection`.  So (A) is a corollary of (C); the
coupled-torus census (mac-mini HYP-4292, kps `CircleClearFloor`) is EVIDENCE for (A), not a
separate obligation.
-/

namespace LonelyRunner
namespace TorusReduction

open TorusProjection

/-- **A loose direction makes the torus loose.**  If some primitive direction `(a,b)` has a
projection `{a·u_i + b·v_i}` that is loose at radius `2/25` (a point `τ` with every
`‖(a u_i + b v_i)τ‖ ≥ 2/25`), then the coupled 2-torus has a `2/25`-loose point `(t,s)`.
The endpoint of the (A) ⟸ (C) reduction. -/
theorem torus_loose_of_loose_direction (u v : Fin 12 → ℤ)
    (h : ∃ a b : ℤ, ∃ τ : ℝ, ∀ i, ∀ m : ℤ,
        (2 : ℝ) / 25 ≤ |((a * u i + b * v i : ℤ) : ℝ) * τ - m|) :
    ∃ t s : ℝ, ∀ i, ∀ m : ℤ, (2 : ℝ) / 25 ≤ |(u i : ℝ) * t + (v i : ℝ) * s - m| := by
  obtain ⟨a, b, hτ⟩ := h
  exact torus_point_of_projection u v a b (2 / 25) hτ

/-- The same, band-generic: a loose direction at ANY radius `c` makes the torus `c`-loose.
Lets the reduction run at whatever threshold the 1-D gap delivers. -/
theorem torus_loose_of_loose_direction_at (u v : Fin 12 → ℤ) (c : ℝ)
    (h : ∃ a b : ℤ, ∃ τ : ℝ, ∀ i, ∀ m : ℤ,
        c ≤ |((a * u i + b * v i : ℤ) : ℝ) * τ - m|) :
    ∃ t s : ℝ, ∀ i, ∀ m : ℤ, c ≤ |(u i : ℝ) * t + (v i : ℝ) * s - m| := by
  obtain ⟨a, b, hτ⟩ := h
  exact torus_point_of_projection u v a b c hτ

#print axioms torus_loose_of_loose_direction

end TorusReduction
end LonelyRunner
