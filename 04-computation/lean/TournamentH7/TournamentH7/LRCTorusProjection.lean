/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: opus (LRC multi-agent project, 2026-07-06-S99)
-/
import Mathlib

/-!
# The torus projection bridge: settled LRC lands on the 2-torus (HYP-4296)

The (A)-leg's projection floor, formalized. A rank-2 coupled family `(u, v)` restricted to
the sub-circle `(t, s) = (a·τ, b·τ)` collapses to the 1-D family `{a·u i + b·v i}`; a lonely
point of that 1-D family IS a lonely point of the torus family, since
`u i · (a τ) + v i · (b τ) = (a·u i + b·v i)·τ`.  With a FULL-SUPPORT direction the projected
family has 12 nonzero speeds, so the settled `LRC(≤13)` citation gives a point of margin
`≥ 1/13` — hence `M(U) ≥ 1/13` for every coupled 2-subtorus.  This is the formal bridge
that pins the (A)-window's lower edge to the settled LRC floor (reflection
`the-projection-floor-and-the-window-sliver`); the residual "no coupling reaches the sliver
`(1/13, 2/25]`" is the census (mac-mini HYP-4292) + rung (kps HYP-4247).
-/

namespace LonelyRunner
namespace TorusProjection

/-- **The projection bridge.** A margin witness `τ` for the 1-D projected family
`{a·u i + b·v i}` lifts to the torus point `(a·τ, b·τ)` with the SAME margins — restricting
the 2-torus to a rational sub-circle can only lose the max, so the torus dominates every
projection. -/
theorem margin_of_projection {ι : Type*} (u v : ι → ℤ) (a b : ℤ) (τ c : ℝ)
    (h : ∀ i, ∀ m : ℤ, c ≤ |((a * u i + b * v i : ℤ) : ℝ) * τ - m|) :
    ∀ i, ∀ m : ℤ, c ≤ |(u i : ℝ) * (a * τ) + (v i : ℝ) * (b * τ) - m| := by
  intro i m
  have key : (u i : ℝ) * (a * τ) + (v i : ℝ) * (b * τ)
      = ((a * u i + b * v i : ℤ) : ℝ) * τ := by
    push_cast; ring
  rw [key]
  exact h i m

/-- **Existence form**: a lonely point of any 1-D projection is a lonely point of the
coupled torus.  Feeding the full-support projection into `LRC(≤13)` yields the `1/13` floor. -/
theorem torus_point_of_projection {ι : Type*} (u v : ι → ℤ) (a b : ℤ) (c : ℝ)
    (hproj : ∃ τ : ℝ, ∀ i, ∀ m : ℤ, c ≤ |((a * u i + b * v i : ℤ) : ℝ) * τ - m|) :
    ∃ t s : ℝ, ∀ i, ∀ m : ℤ, c ≤ |(u i : ℝ) * t + (v i : ℝ) * s - m| := by
  obtain ⟨τ, hτ⟩ := hproj
  exact ⟨a * τ, b * τ, margin_of_projection u v a b τ c hτ⟩

#print axioms margin_of_projection
#print axioms torus_point_of_projection

end TorusProjection
end LonelyRunner
