/-
  TournamentH7.LRCUnitGrid14 -- exact denominator-14 unit-grid sieve for LRC(14).

  This module packages the `q = 14` specialization of `LonelyRunner.sieve_frac`
  in the form used by the LRC14 tightness/apex-denominator reductions:

      a ∈ (Z/14Z)^×  ==>  Lonely 14 v (a/14)  iff  no speed is divisible by 14.

  It is the Lean side of the THM-523 denominator witness at the apex grid, and
  the `Lonely`-predicate companion to `LRCApex7Floor`: without a multiple of
  `14`, every unit apex point is a witness; with a multiple of `14`, no apex
  point can be lonely because that runner sits exactly on the observer.
-/

import TournamentH7.LonelyRunner

namespace LonelyRunner
namespace UnitGrid14

variable {ι : Type*}

/-- Any reduced denominator-14 point is a lonely time if no speed is divisible by
`14`.  This is `sieve_frac` with `n = q = 14`. -/
theorem denom14_unit_lonely (a : ℤ) (hcop : IsCoprime (14 : ℤ) a)
    (v : ι → ℤ) (hdiv : ∀ i, ¬ ((14 : ℤ) ∣ v i)) :
    Lonely 14 v ((a : ℝ) / 14) := by
  simpa using sieve_frac 14 14 a v (le_refl 14) (by norm_num) hcop hdiv

/-- A multiple-of-14 speed blocks every denominator-14 point in the concrete
`Lonely` predicate: the runner is exactly at an integer position. -/
theorem not_denom14_lonely_of_dvd (a : ℤ) (v : ι → ℤ) {i : ι}
    (hdiv : (14 : ℤ) ∣ v i) :
    ¬ Lonely 14 v ((a : ℝ) / 14) := by
  intro h
  obtain ⟨m, hm⟩ := hdiv
  have hfar := h i (m * a)
  have hzero :
      (v i : ℝ) * ((a : ℝ) / 14) - ((m * a : ℤ) : ℝ) = 0 := by
    rw [hm]
    push_cast
    ring
  rw [hzero, abs_zero] at hfar
  norm_num at hfar

/-- Exact apex-grid criterion.  At any unit denominator-14 point, loneliness is
equivalent to the absence of a multiple-of-14 speed. -/
theorem denom14_unit_lonely_iff_no_multiple (a : ℤ) (hcop : IsCoprime (14 : ℤ) a)
    (v : ι → ℤ) :
    Lonely 14 v ((a : ℝ) / 14) ↔ ∀ i, ¬ ((14 : ℤ) ∣ v i) := by
  constructor
  · intro h i hdiv
    exact not_denom14_lonely_of_dvd a v hdiv h
  · intro hdiv
    exact denom14_unit_lonely a hcop v hdiv

theorem denom14_one_lonely_iff_no_multiple (v : ι → ℤ) :
    Lonely 14 v ((1 : ℝ) / 14) ↔ ∀ i, ¬ ((14 : ℤ) ∣ v i) :=
  by simpa using (denom14_unit_lonely_iff_no_multiple (1 : ℤ) (by norm_num) v)

theorem denom14_three_lonely_iff_no_multiple (v : ι → ℤ) :
    Lonely 14 v ((3 : ℝ) / 14) ↔ ∀ i, ¬ ((14 : ℤ) ∣ v i) :=
  by simpa using (denom14_unit_lonely_iff_no_multiple (3 : ℤ) (by norm_num) v)

theorem denom14_five_lonely_iff_no_multiple (v : ι → ℤ) :
    Lonely 14 v ((5 : ℝ) / 14) ↔ ∀ i, ¬ ((14 : ℤ) ∣ v i) :=
  by simpa using (denom14_unit_lonely_iff_no_multiple (5 : ℤ) (by norm_num) v)

theorem denom14_nine_lonely_iff_no_multiple (v : ι → ℤ) :
    Lonely 14 v ((9 : ℝ) / 14) ↔ ∀ i, ¬ ((14 : ℤ) ∣ v i) :=
  by simpa using (denom14_unit_lonely_iff_no_multiple (9 : ℤ) (by norm_num) v)

theorem denom14_eleven_lonely_iff_no_multiple (v : ι → ℤ) :
    Lonely 14 v ((11 : ℝ) / 14) ↔ ∀ i, ¬ ((14 : ℤ) ∣ v i) :=
  by simpa using (denom14_unit_lonely_iff_no_multiple (11 : ℤ) (by norm_num) v)

theorem denom14_thirteen_lonely_iff_no_multiple (v : ι → ℤ) :
    Lonely 14 v ((13 : ℝ) / 14) ↔ ∀ i, ¬ ((14 : ℤ) ∣ v i) :=
  by simpa using (denom14_unit_lonely_iff_no_multiple (13 : ℤ) (by norm_num) v)

/-- Specializing the general q-sieve counterexample condition: any family with
no `14`-lonely time must contain a multiple of `14`. -/
theorem counterexample_needs_multiple14 (v : ι → ℤ)
    (hcex : ∀ t : ℝ, ¬ Lonely 14 v t) :
    ∃ i, (14 : ℤ) ∣ v i :=
  counterexample_needs_all_divisors 14 v hcex 14 (by norm_num) (le_refl 14)

/-! ## Axiom audit -/

#print axioms denom14_unit_lonely_iff_no_multiple
#print axioms counterexample_needs_multiple14

end UnitGrid14
end LonelyRunner
