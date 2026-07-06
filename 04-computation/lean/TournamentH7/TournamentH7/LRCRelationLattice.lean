/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: opus (LRC multi-agent project, 2026-07-06-S112)
-/
import Mathlib

/-!
# The relation lattice is dilation-invariant (HYP-4446)

The density-floor functional `safe(S, β) = Leb{t : ‖vᵢ t‖ ≥ β ∀i}` has the closed form
`safe(S, β) = Σ_{a ∈ L(S)} ∏ᵢ ĥ_{aᵢ}(β)`, a theta-sum over the **relation lattice**
`L(S) = {a : Σ aᵢ vᵢ = 0}` (opus-S112).  Every property the fleet uses about `safe` in the
height direction — dilated APs all tile (`safe(c·AP) = 0`), `safe` does not degrade with height
(mac-mini S17), `M(c·v) = M(v)` (opus S110) — rests on one algebraic fact: **the relation
lattice is unchanged by dilating every speed.**  Scaling `v ↦ c·v` (`c ≠ 0`) neither creates nor
destroys any relation, so the whole theta-sum, and hence `safe`, is a dilation invariant.
-/

namespace LonelyRunner
namespace RelationLattice

/-- **Dilation invariance of the relation lattice.**  For `c ≠ 0`, an integer combination is a
relation of the dilated speeds `c·v` iff it is a relation of `v`: `Σ aᵢ (c·vᵢ) = 0 ↔ Σ aᵢ vᵢ = 0`.
The algebraic backbone of `safe(c·S) = safe(S)` and `M(c·v) = M(v)` — height is a gauge. -/
theorem sum_smul_eq_zero_iff {k : ℕ} (c : ℤ) (hc : c ≠ 0) (v a : Fin k → ℤ) :
    (∑ i, a i * (c * v i)) = 0 ↔ (∑ i, a i * v i) = 0 := by
  have h : (∑ i, a i * (c * v i)) = c * ∑ i, a i * v i := by
    rw [Finset.mul_sum]
    exact Finset.sum_congr rfl (fun i _ => by ring)
  rw [h, mul_eq_zero]
  constructor
  · rintro (h0 | h0)
    · exact absurd h0 hc
    · exact h0
  · exact fun h0 => Or.inr h0

/-- The relation lattices of `c·v` and `v` coincide as sets (`c ≠ 0`). -/
theorem relationSet_smul {k : ℕ} (c : ℤ) (hc : c ≠ 0) (v : Fin k → ℤ) :
    {a : Fin k → ℤ | ∑ i, a i * (c * v i) = 0} = {a : Fin k → ℤ | ∑ i, a i * v i = 0} := by
  ext a
  exact sum_smul_eq_zero_iff c hc v a

#print axioms sum_smul_eq_zero_iff
#print axioms relationSet_smul

end RelationLattice
end LonelyRunner
