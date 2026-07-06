/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: opus (LRC multi-agent project, 2026-07-06-S114)
-/
import TournamentH7.LRCWitnessAttainment

/-!
# The lonely set is open: `safe(S,β) > 0 ⟺ M(S) > β` (HYP-4466)

The fleet's density-floor programme works with the safe-measure
`safe(S,β) = Leb{t : ‖vᵢ t‖ ≥ β ∀i}` and its theta-sum / Fekete / spectral reformulations.
This file records the elementary bridge between that measure and the Lean margin/`M` world:

* the **strictly**-`β`-lonely set `{t : β < margin v t}` is **open** (margin is continuous), so
* if `M(S) > β` (the margin exceeds `β` somewhere) then a whole open interval of times is
  `β`-lonely — hence `safe(S,β) > 0`.

Conversely `M(S) ≤ β` forces the strict lonely set empty, so `safe(S,β) = 0` (measure zero).
Thus **`safe(S,β) > 0 ⟺ M(S) > β`**: the analytic safe-measure routes and the Lean `M`-routes
are interdefinable.  A consequence worth stating plainly (opus-S114): for the *bounded* case,
`safe(S,2/25)=0 ⟹ AP` is a *reformulation* of the gap-emptiness (G), not a reduction of it —
the theta-sum, Freiman, and harmonic-leading-order pictures reconnect (G) to hard theorems but
do not by themselves supply the missing analytic input.  The genuine reduction achieved so far
is the *unbounded* case (safe-equicontinuity compactness); the bounded case still needs either a
height upper bound (making it finite) or a new estimate.
-/

namespace LonelyRunner
namespace LonelyOpen

open TournamentH7.LRCWitness

variable {k : ℕ} [Nonempty (Fin k)]

/-- The strictly-`β`-lonely set `{t : β < margin v t}` is open (the margin is continuous). -/
theorem isOpen_lonely (v : Fin k → ℤ) (β : ℝ) : IsOpen {t : ℝ | β < margin v t} :=
  isOpen_lt continuous_const (continuous_margin v)

/-- **`M > β` yields an open lonely interval (`safe > 0`).**  If the loneliness margin exceeds
`β` at some time `t₀`, a whole open interval `(a,b) ∋ t₀` is strictly `β`-lonely — so the
`β`-lonely set has positive measure.  This is the open-interval form of `safe(S,β) > 0 ⟺
M(S) > β`. -/
theorem exists_Ioo_lonely (v : Fin k → ℤ) (β t0 : ℝ) (h : β < margin v t0) :
    ∃ a b : ℝ, a < t0 ∧ t0 < b ∧ Set.Ioo a b ⊆ {t : ℝ | β < margin v t} := by
  obtain ⟨a, b, ht0, hsub⟩ :=
    mem_nhds_iff_exists_Ioo_subset.mp ((isOpen_lonely v β).mem_nhds h)
  exact ⟨a, b, ht0.1, ht0.2, hsub⟩

/-- Contrapositive packaging: if no open interval around `t₀` is `β`-lonely (e.g. the whole
strict lonely set is empty, `safe = 0`), then the margin never exceeds `β` there — `M ≤ β`. -/
theorem margin_le_of_not_lonely (v : Fin k → ℤ) (β t0 : ℝ)
    (h : ∀ a b : ℝ, a < t0 → t0 < b → ¬ Set.Ioo a b ⊆ {t : ℝ | β < margin v t}) :
    margin v t0 ≤ β := by
  by_contra hlt
  rw [not_le] at hlt
  obtain ⟨a, b, ha, hb, hsub⟩ := exists_Ioo_lonely v β t0 hlt
  exact h a b ha hb hsub

#print axioms isOpen_lonely
#print axioms exists_Ioo_lonely
#print axioms margin_le_of_not_lonely

end LonelyOpen
end LonelyRunner
