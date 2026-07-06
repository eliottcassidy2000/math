/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: opus (LRC multi-agent project, 2026-07-06-S111)
-/
import TournamentH7.LRCWitnessAttainment

/-!
# Leave-one-out nesting: covering IS alignment (HYP, metric side)

mac-mini-S17 (HYP-4452) isolated the *quantitative* content of gap-emptiness — the piece the
structural lenses (Farey ladder, sum-product, roots of unity) cannot see, because the gap is
nonempty at `n = 7, 8` yet empty at `n = 13`.  Its cleanest face is the **leave-one-out
alignment**:

> if a family `S` covers at level `β` (`margin < β` everywhere — no `β`-lonely time), then for
> every runner `j`, the *hole* of the leave-one-out subfamily `S \ {vⱼ}` is contained in the
> *danger arc* of `vⱼ`.

Equivalently, at any time `t` where every *other* runner is `≥ β` from `ℤ`, the dropped runner
`vⱼ` must be `< β` from `ℤ`.  This is the exact, elementary core of the covering-as-nesting
picture: covering forces each thin leave-one-out hole to nest inside `vⱼ`'s harmonic arcs — an
**equidistribution rigidity** at hole-width `1/(n(2n-1))` that (mac-mini's evidence) only the
AP's arc lattice can satisfy at `n = 13`.  This file formalizes the nesting itself; the
*width-vs-budget* infeasibility that makes it decisive at `n = 13` is the remaining analytic
piece.
-/

namespace LonelyRunner
namespace LeaveOneOut

open TournamentH7.LRCWitness

variable {k : ℕ} [Nonempty (Fin k)]

/-- **Leave-one-out nesting (pointwise).**  If the margin is `< β` at `t` (this `t` is *not*
`β`-lonely) and every runner other than `j` is `≥ β` from `ℤ` at `t`, then runner `j` is
`< β` from `ℤ` at `t`.  The dropped runner's danger arc must contain the rest-family's hole. -/
theorem nesting (v : Fin k → ℤ) (β t : ℝ) (j : Fin k)
    (hcov : margin v t < β)
    (hrest : ∀ i, i ≠ j → β ≤ distZ ((v i : ℝ) * t)) :
    distZ ((v j : ℝ) * t) < β := by
  by_contra h
  rw [not_lt] at h
  have hall : ∀ i, β ≤ distZ ((v i : ℝ) * t) := by
    intro i
    rcases eq_or_ne i j with rfl | hij
    · exact h
    · exact hrest i hij
  have hle : β ≤ margin v t := by
    rw [margin, Finset.le_inf'_iff]
    intro i _
    exact hall i
  linarith

/-- **Covering forces nesting (global).**  If `S` covers at level `β` (no time is `β`-lonely,
`margin < β` everywhere), then at every time `t` the leave-one-out hole of `S \ {vⱼ}` sits
inside the danger arc of `vⱼ`: whenever all `i ≠ j` are `β`-far, `vⱼ` is `β`-near. -/
theorem covering_forces_nesting (v : Fin k → ℤ) (β : ℝ)
    (hcov : ∀ t, margin v t < β) (j : Fin k) (t : ℝ)
    (hrest : ∀ i, i ≠ j → β ≤ distZ ((v i : ℝ) * t)) :
    distZ ((v j : ℝ) * t) < β :=
  nesting v β t j (hcov t) hrest

#print axioms nesting
#print axioms covering_forces_nesting

end LeaveOneOut
end LonelyRunner
