/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Codex (LRC multi-agent project, 2026-07-18)
-/
import Mathlib.Data.Nat.GCD.Basic
import Mathlib.Data.Real.Basic
import Mathlib.Algebra.BigOperators.Group.Finset.Basic
import Mathlib.Algebra.BigOperators.Ring.Finset
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.NormNum
import Mathlib.Tactic.Ring
import Lean.Elab.Tactic.Omega

/-!
# Compact essential-crown arithmetic

This module formalizes two elementary consumers from the proof frame around
THM-1149.

The first is the finite algebraic core of the essential-crown mass identity.
If the multiplicity-one, multiplicity-two, and multiplicity-at-least-three
strata have total mass one and first moment two, then the private
(multiplicity-one) mass is exactly the excess over multiplicity two on the
high strata.

The second is the arithmetic endpoint of the regenerated AP-core branch.  An
analytic/Farey producer may supply `13 * d ∣ v`.  If the AP core
`d * {1, ..., 12}` together with the extra speed `v` is primitive, then
`d = 1`.  If the extra speed also carries modulus 14, it is a positive
multiple of `182`; hence it lies beyond `13 * 12`.  This contradicts the
cleared compactness inequality `v < 13 * (12 * d)`.

Nothing here proves the analytic crown-extraction or Farey-regeneration
hypothesis.  Those inputs remain explicit in the theorem statements.
-/

namespace LonelyRunner
namespace CompactEssentialCrown

/-! ## Finite private-mass balance -/

/-- The finite algebraic form of the essential-crown balance law.

`μ 1` and `μ 2` are the masses of the multiplicity-one and
multiplicity-two strata.  The finite set `high` indexes the remaining strata;
the caller is responsible for establishing that this is the complete
multiplicity-at-least-three decomposition.  The two hypotheses are precisely
total mass `1` and first moment `2`.

The conclusion is

`μ 1 = ∑ k in high, (k - 2) * μ k`.

No positivity hypothesis is needed for this algebraic identity. -/
theorem finite_private_mass_balance
    (high : Finset ℕ) (μ : ℕ → ℝ)
    (total_mass : μ 1 + μ 2 + ∑ k ∈ high, μ k = 1)
    (first_moment :
      μ 1 + 2 * μ 2 + ∑ k ∈ high, (k : ℝ) * μ k = 2) :
    μ 1 = ∑ k ∈ high, ((k : ℝ) - 2) * μ k := by
  have split_high_moment :
      (∑ k ∈ high, (k : ℝ) * μ k) =
        2 * (∑ k ∈ high, μ k) +
          ∑ k ∈ high, ((k : ℝ) - 2) * μ k := by
    rw [Finset.mul_sum, ← Finset.sum_add_distrib]
    apply Finset.sum_congr rfl
    intro k hk
    ring
  rw [split_high_moment] at first_moment
  linarith

/-! ## Regenerated AP-core endpoint -/

/-- If the regenerated-tooth conclusion `13 * d ∣ v` holds and the AP scale
`d` is coprime to the extra speed, then the scale is primitive: `d = 1`.

For the family `d * {1, ..., 12} ∪ {v}`, the coprimality hypothesis is exactly
the arithmetic content of primitivity. -/
theorem ap_scale_eq_one_of_regenerated_and_primitive
    {d v : ℕ} (regenerated : 13 * d ∣ v) (primitive : d.Coprime v) :
    d = 1 := by
  have d_dvd_thirteen_mul_d : d ∣ 13 * d := by
    exact ⟨13, by omega⟩
  have d_dvd_v : d ∣ v := d_dvd_thirteen_mul_d.trans regenerated
  exact Nat.eq_one_of_dvd_coprimes primitive dvd_rfl d_dvd_v

/-- Simultaneous divisibility by the coprime moduli `13` and `14` gives the
deep-well modulus `182`. -/
theorem one_eighty_two_dvd_of_thirteen_dvd_of_fourteen_dvd
    {v : ℕ} (h13 : 13 ∣ v) (h14 : 14 ∣ v) : 182 ∣ v := by
  have h182 : 13 * 14 ∣ v :=
    (by decide : Nat.Coprime 13 14).mul_dvd_of_dvd_of_dvd h13 h14
  norm_num at h182 ⊢
  exact h182

/-- The complete post-extraction arithmetic conclusion.  A positive extra
speed which kills all regenerated denominator-13 teeth, is primitive with the
AP scale, and carries modulus 14 lies strictly beyond the compactness wall
`13 * (12 * d)`.

The numerical endpoint is `13 * 12 = 156 < 182 ≤ v` after primitivity forces
`d = 1`. -/
theorem regenerated_primitive_cover14_forces_beyond_compact_wall
    {d v : ℕ} (hv : 0 < v) (regenerated : 13 * d ∣ v)
    (primitive : d.Coprime v) (cover14 : 14 ∣ v) :
    13 * (12 * d) < v := by
  have hd : d = 1 :=
    ap_scale_eq_one_of_regenerated_and_primitive regenerated primitive
  subst d
  have h13 : 13 ∣ v := by
    simpa using regenerated
  have h182 : 182 ∣ v :=
    one_eighty_two_dvd_of_thirteen_dvd_of_fourteen_dvd h13 cover14
  have hv182 : 182 ≤ v := Nat.le_of_dvd hv h182
  norm_num
  omega

/-- Corollary in the exact shape used by the compact AP-core branch: the
cleared ratio inequality `v < 13 * (12 * d)` is incompatible with the
regenerated-tooth, primitive, and Cover14 hypotheses.

This theorem deliberately starts *after* a tight AP-core has been extracted;
it does not assert that every essential crown contains such a core. -/
theorem no_compact_regenerated_primitive_cover14_branch
    {d v : ℕ} (hv : 0 < v) (regenerated : 13 * d ∣ v)
    (primitive : d.Coprime v) (cover14 : 14 ∣ v)
    (compact : v < 13 * (12 * d)) : False := by
  have beyond : 13 * (12 * d) < v :=
    regenerated_primitive_cover14_forces_beyond_compact_wall
      hv regenerated primitive cover14
  omega

end CompactEssentialCrown
end LonelyRunner

#print axioms LonelyRunner.CompactEssentialCrown.finite_private_mass_balance
#print axioms LonelyRunner.CompactEssentialCrown.no_compact_regenerated_primitive_cover14_branch
