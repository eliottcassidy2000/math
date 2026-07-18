/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: codex-S67 (LRC multi-agent project, 2026-07-18)
-/
import Mathlib

/-!
# Rational-scale and cyclic-coset guardrails

Two elementary facts prevent misleading reductions in the LRC14 rational-witness
language.

* For a positive denominator `q`, the upper inequality `val / q < 1 / 13` says
  exactly `13 * val < q`.  The distinct inequality `q < 14 * val` says exactly
  `1 / 14 < val / q`; it is not a consequence of the upper inequality.
* An integer affine lift `u + q * k` can be made divisible by `d` exactly when
  `gcd q d` divides `u`.  This is the precise compatibility condition between
  two modular charts; coprimality is only the special case in which every
  integer `u` lifts.
* If `val` is coprime to `q`, multiplication by `val` permutes `ZMod q`.
  Consequently the set of all multiples, and every additive translate of it, is
  the whole cyclic group.  Such a translate is not a proper coset quotient.

The results are deliberately representation-level guardrails.  They make no claim
about which rational point maximizes an LRC margin.
-/

namespace LonelyRunner
namespace RationalScaleGuardrails

/-- The upper rational-scale inequality, with its denominator cleared exactly. -/
theorem quotient_lt_one_thirteenth_iff {α : Type*} [Field α] [LinearOrder α]
    [IsStrictOrderedRing α]
    {val q : α} (hq : 0 < q) :
    val / q < 1 / 13 ↔ 13 * val < q := by
  rw [div_lt_iff₀ hq]
  constructor <;> intro h <;> nlinarith

/-- The lower LRC14 inequality, with its denominator cleared exactly. -/
theorem one_fourteenth_lt_quotient_iff {α : Type*} [Field α] [LinearOrder α]
    [IsStrictOrderedRing α]
    {val q : α} (hq : 0 < q) :
    1 / 14 < val / q ↔ q < 14 * val := by
  rw [lt_div_iff₀ hq]
  constructor <;> intro h <;> nlinarith

/-- The open critical window requires both independent cleared inequalities. -/
theorem quotient_in_critical_window_iff {α : Type*} [Field α] [LinearOrder α]
    [IsStrictOrderedRing α]
    {val q : α} (hq : 0 < q) :
    1 / 14 < val / q ∧ val / q < 1 / 13 ↔
      q < 14 * val ∧ 13 * val < q := by
  rw [one_fourteenth_lt_quotient_iff hq, quotient_lt_one_thirteenth_iff hq]

/-- A concrete counterexample to inferring the lower window inequality from the
upper one: `1/15 < 1/13`, but `1/15` is not above `1/14`. -/
theorem upper_scale_bound_does_not_force_lower_scale_bound :
    (1 : ℚ) / 15 < 1 / 13 ∧ ¬ (1 : ℚ) / 14 < 1 / 15 := by
  norm_num

/-- Exact integer affine-lift compatibility.  There is an integer translate
`u + q * k` divisible by `d` if and only if the common divisor `gcd q d`
already divides `u`.

No positivity hypotheses are needed: the equivalence also covers `q = 0`
and/or `d = 0`.  The reverse implication is the explicit Bézout lift. -/
theorem exists_integer_affine_lift_iff_gcd_dvd (q d : ℕ) (u : ℤ) :
    (∃ k : ℤ, (d : ℤ) ∣ u + (q : ℤ) * k) ↔
      (Nat.gcd q d : ℤ) ∣ u := by
  constructor
  · rintro ⟨k, hk⟩
    have hgq : (Nat.gcd q d : ℤ) ∣ (q : ℤ) := by
      exact_mod_cast Nat.gcd_dvd_left q d
    have hgd : (Nat.gcd q d : ℤ) ∣ (d : ℤ) := by
      exact_mod_cast Nat.gcd_dvd_right q d
    have hsum : (Nat.gcd q d : ℤ) ∣ u + (q : ℤ) * k := hgd.trans hk
    have hmul : (Nat.gcd q d : ℤ) ∣ (q : ℤ) * k :=
      dvd_mul_of_dvd_left hgq k
    simpa using Int.dvd_sub hsum hmul
  · rintro ⟨m, hm⟩
    refine ⟨-(Nat.gcdA q d * m), ?_⟩
    refine ⟨Nat.gcdB q d * m, ?_⟩
    rw [hm, Nat.gcd_eq_gcd_ab]
    ring

/-- If `val` is coprime to `q`, every additive translate of the multiples of
`val` modulo `q` is all of `ZMod q`. -/
theorem coprime_affine_range_eq_univ (val q : ℕ) (shift : ZMod q)
    (hcop : val.Coprime q) :
    Set.range (fun x : ZMod q => shift + (val : ZMod q) * x) = Set.univ := by
  apply Set.eq_univ_of_forall
  intro y
  let u : (ZMod q)ˣ := ZMod.unitOfCoprime val hcop
  refine ⟨(↑(u⁻¹) : ZMod q) * (y - shift), ?_⟩
  change shift + (u : ZMod q) * ((↑(u⁻¹) : ZMod q) * (y - shift)) = y
  rw [← mul_assoc]
  change shift + (u.val * u.inv) * (y - shift) = y
  rw [u.val_inv, one_mul]
  simpa only [add_comm] using sub_add_cancel y shift

/-- In particular, the `val`-multiple set itself is the whole cyclic group. -/
theorem coprime_multiple_range_eq_univ (val q : ℕ) (hcop : val.Coprime q) :
    Set.range (fun x : ZMod q => (val : ZMod q) * x) = Set.univ := by
  simpa using coprime_affine_range_eq_univ val q 0 hcop

end RationalScaleGuardrails
end LonelyRunner

#print axioms LonelyRunner.RationalScaleGuardrails.quotient_in_critical_window_iff
#print axioms LonelyRunner.RationalScaleGuardrails.exists_integer_affine_lift_iff_gcd_dvd
#print axioms LonelyRunner.RationalScaleGuardrails.coprime_affine_range_eq_univ
