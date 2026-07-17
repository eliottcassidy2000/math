/-
  Cofinal clean moduli for the normalized B5 ledger. The explicit ruler is
  1 mod 14, coprime to every nonzero speed, and converts every modular relation
  in a chosen coefficient-height box into an exact integer relation. This is
  the finite-height bridge to THM-939's dense-core traps. No \`sorry\`; no
  \`native_decide\`.
-/

import Mathlib

namespace LonelyRunner
namespace LRCB5CleanModulus

open Finset

/-- A nonzero multiple of a positive modulus cannot lie strictly between its
negative and positive copies. -/
theorem eq_zero_of_dvd_of_abs_lt {q z : ℤ} (hq : 0 < q)
    (hdiv : q ∣ z) (hsmall : |z| < q) : z = 0 := by
  obtain ⟨multiple, rfl⟩ := hdiv
  by_cases hm : multiple = 0
  · simp [hm]
  · have hone : (1 : ℤ) ≤ |multiple| := Int.one_le_abs hm
    rw [abs_mul, abs_of_pos hq] at hsmall
    nlinarith [abs_nonneg multiple]

/-- Height times total speed mass bounds an integer linear form. -/
theorem abs_relation_le_height_mass {n : ℕ} (coeff speed : Fin n → ℤ)
    (height : ℤ) (_hheight : 0 ≤ height)
    (hcoeff : ∀ i, |coeff i| ≤ height) :
    |∑ i, coeff i * speed i| ≤ height * ∑ i, |speed i| := by
  calc
    |∑ i, coeff i * speed i| ≤ ∑ i, |coeff i * speed i| :=
      Finset.abs_sum_le_sum_abs _ _
    _ = ∑ i, |coeff i| * |speed i| := by
      apply Finset.sum_congr rfl
      intro i _
      rw [abs_mul]
    _ ≤ ∑ i, height * |speed i| := by
      apply Finset.sum_le_sum
      intro i _
      exact mul_le_mul_of_nonneg_right (hcoeff i) (abs_nonneg _)
    _ = height * ∑ i, |speed i| := by rw [Finset.mul_sum]

/-- Low-height modular relations are genuine integer relations once the ruler
is larger than `height * ∑ |vᵢ|`.  This is the exact finite bridge needed
before applying THM-939's relation traps. -/
theorem modular_relation_is_exact {n : ℕ} (coeff speed : Fin n → ℤ)
    (height q : ℤ) (hheight : 0 ≤ height) (hq : 0 < q)
    (hcoeff : ∀ i, |coeff i| ≤ height)
    (hqLarge : height * ∑ i, |speed i| < q)
    (hmod : q ∣ ∑ i, coeff i * speed i) :
    ∑ i, coeff i * speed i = 0 := by
  apply eq_zero_of_dvd_of_abs_lt hq hmod
  exact lt_of_le_of_lt
    (abs_relation_le_height_mass coeff speed height hheight hcoeff) hqLarge

/-- A concrete cofinal ruler family.  It is `1 mod 14`, exceeds the chosen
height/mass scale, and is coprime to every nonzero speed. -/
def cleanModulus {n : ℕ} (speed : Fin n → ℤ) (height : ℕ) : ℕ :=
  14 * (height * ∑ i, (speed i).natAbs + 1) * (∏ i, (speed i).natAbs) + 1

theorem cleanModulus_mod_fourteen {n : ℕ} (speed : Fin n → ℤ) (height : ℕ) :
    cleanModulus speed height % 14 = 1 := by
  simp [cleanModulus, Nat.add_mod, Nat.mul_assoc]

theorem one_lt_cleanModulus {n : ℕ} (speed : Fin n → ℤ) (height : ℕ)
    (hspeed : ∀ i, speed i ≠ 0) :
    1 < cleanModulus speed height := by
  have hprod : 0 < ∏ i, (speed i).natAbs := by
    exact Finset.prod_pos fun i _ => Int.natAbs_pos.mpr (hspeed i)
  unfold cleanModulus
  nlinarith [Nat.zero_le (height * ∑ i, (speed i).natAbs)]

theorem fourteen_le_cleanModulus {n : ℕ} (speed : Fin n → ℤ) (height : ℕ)
    (hspeed : ∀ i, speed i ≠ 0) :
    14 ≤ cleanModulus speed height := by
  have hprod : 0 < ∏ i, (speed i).natAbs := by
    exact Finset.prod_pos fun i _ => Int.natAbs_pos.mpr (hspeed i)
  unfold cleanModulus
  nlinarith [Nat.zero_le (height * ∑ i, (speed i).natAbs)]

theorem height_mass_lt_cleanModulus {n : ℕ} (speed : Fin n → ℤ)
    (height : ℕ) (hspeed : ∀ i, speed i ≠ 0) :
    height * ∑ i, (speed i).natAbs < cleanModulus speed height := by
  have hprod : 0 < ∏ i, (speed i).natAbs := by
    exact Finset.prod_pos fun i _ => Int.natAbs_pos.mpr (hspeed i)
  unfold cleanModulus
  nlinarith [Nat.zero_le (height * ∑ i, (speed i).natAbs)]

/-- A clean modulus is not merely larger than the requested relation-height
box: for a nonzero family it is larger than fourteen times every individual
speed.  This useful scale fact also exposes why the first multiplier is a
near-zero obstruction to any depth-six coverage cap at this ruler. -/
theorem fourteen_natAbs_lt_cleanModulus {n : ℕ} (speed : Fin n → ℤ)
    (height : ℕ) (hspeed : ∀ i, speed i ≠ 0) (i : Fin n) :
    14 * (speed i).natAbs < cleanModulus speed height := by
  have hprod : 0 < ∏ j, (speed j).natAbs := by
    exact Finset.prod_pos fun j _ => Int.natAbs_pos.mpr (hspeed j)
  have hdiv : (speed i).natAbs ∣ ∏ j, (speed j).natAbs := by
    exact Finset.dvd_prod_of_mem _ (Finset.mem_univ i)
  have hle : (speed i).natAbs ≤ ∏ j, (speed j).natAbs :=
    Nat.le_of_dvd hprod hdiv
  unfold cleanModulus
  nlinarith [Nat.zero_le (height * ∑ j, (speed j).natAbs)]

theorem speed_coprime_cleanModulus {n : ℕ} (speed : Fin n → ℤ)
    (height : ℕ) (i : Fin n) :
    Nat.Coprime (speed i).natAbs (cleanModulus speed height) := by
  have hdiv : (speed i).natAbs ∣ ∏ j, (speed j).natAbs := by
    exact Finset.dvd_prod_of_mem _ (Finset.mem_univ i)
  obtain ⟨factor, hfactor⟩ := hdiv
  unfold cleanModulus
  rw [hfactor]
  rw [show 14 * (height * ∑ j, (speed j).natAbs + 1) *
      ((speed i).natAbs * factor) + 1 =
      (14 * (height * ∑ j, (speed j).natAbs + 1) * factor) *
        (speed i).natAbs + 1 by ring]
  exact (Nat.coprime_mul_right_add_right ((speed i).natAbs) 1
    (14 * (height * ∑ j, (speed j).natAbs + 1) * factor)).mpr
      (Nat.coprime_one_right _)


/-- Integer-gcd form consumed by the singleton-deviation theorem. -/
theorem int_gcd_speed_cleanModulus_eq_one {n : ℕ} (speed : Fin n → ℤ)
    (height : ℕ) (i : Fin n) :
    Int.gcd (speed i) (cleanModulus speed height : ℤ) = 1 := by
  rw [Int.gcd_eq_natAbs]
  simpa using (Nat.coprime_iff_gcd_eq_one.mp
    (speed_coprime_cleanModulus speed height i))

/-! ## Axiom audit -/

#print axioms modular_relation_is_exact
#print axioms cleanModulus_mod_fourteen
#print axioms fourteen_le_cleanModulus
#print axioms height_mass_lt_cleanModulus
#print axioms fourteen_natAbs_lt_cleanModulus
#print axioms int_gcd_speed_cleanModulus_eq_one

end LRCB5CleanModulus
end LonelyRunner
