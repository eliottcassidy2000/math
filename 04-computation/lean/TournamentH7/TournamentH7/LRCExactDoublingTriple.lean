/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: codex-S59 (LRC multi-agent project, 2026-07-17)
-/
import Mathlib

set_option maxHeartbeats 800000
set_option maxRecDepth 100000

/-!
# Exact Bernoulli residue law for doubling triples

This file formalizes the finite algebraic core of THM-948's exact support law for
`{a, 2a, b}`.  If `d = gcd(a,b)`, `A = a/d`, and `B = b/d`, the analytic argument
in THM-948 reduces the centered support mass to

`C(A mod 14, B mod 28) / (686 A B)`.

Here `C` is defined directly from the periodic Bernoulli `B₂` formula.  Kernel
reduction proves all 392 residue cells at once:

* the Bernoulli expression is integral;
* `C = 0` exactly on `A = 0,7 (mod 14)` or `B = 0,14 (mod 28)`;
* its sign is the product of the two centered-residue signs;
* reflection in either residue circle reverses `C`.

The finite result is then lifted to arbitrary natural residues, giving
`C(A,B)=0 ↔ 7 ∣ A ∨ 14 ∣ B` without any coprimality assumption.  The separate
Fourier/integration step identifying a concrete support mass with the Bernoulli
formula is deliberately not asserted here.

No `sorry`; no `native_decide`.  The 392-cell proofs use kernel `decide`.
-/

namespace LonelyRunner
namespace ExactDoublingTriple

/--
If `r = n mod m`, this is `6m² B₂(r/m)`.  Keeping the periodic Bernoulli
evaluation in integers is much lighter for kernel reduction than normalizing
392 rational floors.
-/
def bernoulliTwoNumerator (n : ℤ) (m : ℕ) : ℤ :=
  let M : ℤ := m
  let r := n % M
  6 * r * r - 6 * r * M + M * M

/-- Exact rational implementation of `B₂({n/m})` for the positive moduli used below. -/
def periodicBernoulliTwoAt (n : ℤ) (m : ℕ) : ℚ :=
  (bernoulliTwoNumerator n m : ℚ) / (6 * (m : ℚ) ^ 2)

/-- The numerator of `J(A/7,B/14)`, whose fixed denominator is `6*28²`. -/
def sinePair28Numerator (A B : ℕ) : ℤ :=
  bernoulliTwoNumerator (2 * (A : ℤ) - B) 28 -
    bernoulliTwoNumerator (2 * (A : ℤ) + B) 28

/-- The numerator of `J(A/7,B/7)`, whose fixed denominator is `6*14²`. -/
def sinePair14Numerator (A B : ℕ) : ℤ :=
  bernoulliTwoNumerator ((A : ℤ) - B) 14 -
    bernoulliTwoNumerator ((A : ℤ) + B) 14

/-- The numerator of `J(2A/7,B/7)`, whose fixed denominator is `6*14²`. -/
def sinePair14DoubleNumerator (A B : ℕ) : ℤ :=
  bernoulliTwoNumerator (2 * (A : ℤ) - B) 14 -
    bernoulliTwoNumerator (2 * (A : ℤ) + B) 14

/-- Exact `J(A/7,B/14)` Bernoulli difference. -/
def sinePair28 (A B : ℕ) : ℚ :=
  (sinePair28Numerator A B : ℚ) / (6 * (28 : ℚ) ^ 2)

/-- Exact `J(A/7,B/7)` Bernoulli difference. -/
def sinePair14 (A B : ℕ) : ℚ :=
  (sinePair14Numerator A B : ℚ) / (6 * (14 : ℚ) ^ 2)

/-- Exact `J(2A/7,B/7)` Bernoulli difference. -/
def sinePair14Double (A B : ℕ) : ℚ :=
  (sinePair14DoubleNumerator A B : ℚ) / (6 * (14 : ℚ) ^ 2)

/--
The normalized Fourier/Bernoulli expression `Q`, so that the closed support
formula is `M({A,2A,B}) = Q(A mod 14,B mod 28)/(A B)`.
-/
def normalizedFourierMass (aResidue bResidue : ℕ) : ℚ :=
  if bResidue % 2 = 0 then
    (5 / 7) * sinePair28 aResidue bResidue -
      (1 / 7) * sinePair14 aResidue bResidue
  else
    sinePair28 aResidue bResidue -
      (1 / 14) * sinePair14Double aResidue bResidue -
      (1 / 7) * sinePair14 aResidue bResidue

/--
After clearing the fixed Bernoulli denominators, `C` is this integer divided by 48.
The even and odd formulas are respectively
`5 N₂₈ - 4 N₁₄` and `7 N₂₈ - 2 N₁₄' - 4 N₁₄`.
-/
def coefficientNumerator (aResidue bResidue : ℕ) : ℤ :=
  if bResidue % 2 = 0 then
    5 * sinePair28Numerator aResidue bResidue -
      4 * sinePair14Numerator aResidue bResidue
  else
    7 * sinePair28Numerator aResidue bResidue -
      2 * sinePair14DoubleNumerator aResidue bResidue -
      4 * sinePair14Numerator aResidue bResidue

/-- The Bernoulli expression scaled by 686; it will be proved integral. -/
def cellCoefficientRat (aResidue bResidue : ℕ) : ℚ :=
  686 * normalizedFourierMass aResidue bResidue

/-- Integer numerator of the scaled Bernoulli expression. -/
def cellCoefficient (aResidue bResidue : ℕ) : ℤ :=
  coefficientNumerator aResidue bResidue / 48

/-- Centered sign on the 14-cycle. -/
def centeredSign14 (a : ℕ) : ℤ :=
  if a % 14 = 0 ∨ a % 14 = 7 then 0 else if a % 14 < 7 then 1 else -1

/-- Centered sign on the 28-cycle. -/
def centeredSign28 (b : ℕ) : ℤ :=
  if b % 28 = 0 ∨ b % 28 = 14 then 0 else if b % 28 < 14 then 1 else -1

/-- The integer coefficient on arbitrary natural inputs, reduced to its residue cell. -/
def residueCoefficient (A B : ℕ) : ℤ :=
  cellCoefficient (A % 14) (B % 28)

/-- The rational Bernoulli expression on arbitrary natural inputs. -/
def residueCoefficientRat (A B : ℕ) : ℚ :=
  cellCoefficientRat (A % 14) (B % 28)

/--
The closed rational expression supplied to the analytic support-mass interface.
For positive coprime `A,B`, THM-948 identifies the concrete mass with this value.
-/
def closedResidueMass (A B : ℕ) : ℚ :=
  (residueCoefficient A B : ℚ) / ((686 : ℚ) * A * B)

/-! ## The denominator-clearing identity -/

/-- Clearing the three fixed Bernoulli denominators leaves exactly `/48`. -/
theorem cellCoefficientRat_eq_numerator_div (a b : ℕ) :
    cellCoefficientRat a b = (coefficientNumerator a b : ℚ) / 48 := by
  by_cases h : b % 2 = 0
  · simp [cellCoefficientRat, normalizedFourierMass, coefficientNumerator,
      sinePair28, sinePair14, h]
    ring
  · simp [cellCoefficientRat, normalizedFourierMass, coefficientNumerator,
      sinePair28, sinePair14, sinePair14Double, h]
    ring

/-! ## The complete 392-cell kernel audit -/

/--
One shared kernel computation audits integrality, the zero/sign laws, and both
residue reflections.  The named theorems below are projections of this result,
so the Bernoulli table is reduced only once during a build.
-/
theorem cell_audit :
    ∀ a : Fin 14, ∀ b : Fin 28,
      (let c := cellCoefficient a b
       coefficientNumerator a b = 48 * c ∧
       (c = 0 ↔
         ((a : ℕ) = 0 ∨ (a : ℕ) = 7) ∨ ((b : ℕ) = 0 ∨ (b : ℕ) = 14)) ∧
       c.sign = centeredSign14 a * centeredSign28 b ∧
       cellCoefficient ((14 - (a : ℕ)) % 14) b = -c ∧
       cellCoefficient a ((28 - (b : ℕ)) % 28) = -c) := by
  intro a b
  fin_cases a <;> fin_cases b <;> decide

/-- Every one of the 392 scaled Bernoulli values is the advertised integer. -/
theorem cell_coefficient_integral :
    ∀ a : Fin 14, ∀ b : Fin 28,
      cellCoefficientRat a b = (cellCoefficient a b : ℚ) := by
  intro a b
  rw [cellCoefficientRat_eq_numerator_div, (cell_audit a b).1]
  push_cast
  ring

/-- The exact finite zero locus. -/
theorem cell_coefficient_zero_iff :
    ∀ a : Fin 14, ∀ b : Fin 28,
      cellCoefficient a b = 0 ↔
        ((a : ℕ) = 0 ∨ (a : ℕ) = 7) ∨ ((b : ℕ) = 0 ∨ (b : ℕ) = 14) := by
  intro a b
  exact (cell_audit a b).2.1

/-- The complete sign law, including its zero cells. -/
theorem cell_coefficient_sign :
    ∀ a : Fin 14, ∀ b : Fin 28,
      (cellCoefficient a b).sign = centeredSign14 a * centeredSign28 b := by
  intro a b
  exact (cell_audit a b).2.2.1

/-- Reflecting the 14-residue reverses the coefficient. -/
theorem cell_reflect_fourteen :
    ∀ a : Fin 14, ∀ b : Fin 28,
      cellCoefficient ((14 - (a : ℕ)) % 14) b = -cellCoefficient a b := by
  intro a b
  exact (cell_audit a b).2.2.2.1

/-- Reflecting the 28-residue reverses the coefficient. -/
theorem cell_reflect_twenty_eight :
    ∀ a : Fin 14, ∀ b : Fin 28,
      cellCoefficient a ((28 - (b : ℕ)) % 28) = -cellCoefficient a b := by
  intro a b
  exact (cell_audit a b).2.2.2.2

/-! ## Lifting the finite audit to arbitrary residues -/

/-- The arbitrary-residue Bernoulli expression is always the integer coefficient. -/
theorem residue_coefficient_integral (A B : ℕ) :
    residueCoefficientRat A B = (residueCoefficient A B : ℚ) := by
  let a : Fin 14 := ⟨A % 14, Nat.mod_lt _ (by omega)⟩
  let b : Fin 28 := ⟨B % 28, Nat.mod_lt _ (by omega)⟩
  simpa [residueCoefficientRat, residueCoefficient, a, b] using
    cell_coefficient_integral a b

/-- Exact arbitrary-residue form of the finite zero locus. -/
theorem residue_coefficient_zero_iff_residues (A B : ℕ) :
    residueCoefficient A B = 0 ↔
      (A % 14 = 0 ∨ A % 14 = 7) ∨ (B % 28 = 0 ∨ B % 28 = 14) := by
  let a : Fin 14 := ⟨A % 14, Nat.mod_lt _ (by omega)⟩
  let b : Fin 28 := ⟨B % 28, Nat.mod_lt _ (by omega)⟩
  simpa [residueCoefficient, a, b] using cell_coefficient_zero_iff a b

/-- Exact arbitrary-residue centered-sign law. -/
theorem residue_coefficient_sign (A B : ℕ) :
    (residueCoefficient A B).sign = centeredSign14 A * centeredSign28 B := by
  let a : Fin 14 := ⟨A % 14, Nat.mod_lt _ (by omega)⟩
  let b : Fin 28 := ⟨B % 28, Nat.mod_lt _ (by omega)⟩
  simpa [residueCoefficient, centeredSign14, centeredSign28, a, b,
    Nat.mod_mod] using cell_coefficient_sign a b

/-- Multiples of seven are exactly the two zero classes on the 14-cycle. -/
theorem mod_fourteen_zero_or_seven_iff (A : ℕ) :
    A % 14 = 0 ∨ A % 14 = 7 ↔ 7 ∣ A := by
  rw [Nat.dvd_iff_mod_eq_zero,
    ← Nat.mod_mod_of_dvd A (by norm_num : 7 ∣ 14)]
  omega

/-- Multiples of fourteen are exactly the two zero classes on the 28-cycle. -/
theorem mod_twenty_eight_zero_or_fourteen_iff (B : ℕ) :
    B % 28 = 0 ∨ B % 28 = 14 ↔ 14 ∣ B := by
  rw [Nat.dvd_iff_mod_eq_zero,
    ← Nat.mod_mod_of_dvd B (by norm_num : 14 ∣ 28)]
  omega

/--
The full zero theorem: the Bernoulli coefficient vanishes exactly on the two
Fourier divisor branches `7 ∣ A` or `14 ∣ B`.
-/
theorem residue_coefficient_zero_iff_dvd (A B : ℕ) :
    residueCoefficient A B = 0 ↔ 7 ∣ A ∨ 14 ∣ B := by
  rw [residue_coefficient_zero_iff_residues,
    mod_fourteen_zero_or_seven_iff, mod_twenty_eight_zero_or_fourteen_iff]

/-- The same full zero locus at the closed rational mass interface. -/
theorem closed_residue_mass_zero_iff_dvd {A B : ℕ} (hA : 0 < A) (hB : 0 < B) :
    closedResidueMass A B = 0 ↔ 7 ∣ A ∨ 14 ∣ B := by
  simp [closedResidueMass, residue_coefficient_zero_iff_dvd A B,
    Nat.ne_of_gt hA, Nat.ne_of_gt hB]

/-- The original negative support-three witness is an exact coefficient cell. -/
theorem coefficient_one_fifteen : residueCoefficient 1 15 = -24 := by
  decide

/-- The entire `15 mod 28` family stays on that same negative coefficient cell. -/
theorem coefficient_one_fifteen_add_twenty_eight_mul (m : ℕ) :
    residueCoefficient 1 (15 + 28 * m) = -24 := by
  calc
    residueCoefficient 1 (15 + 28 * m) = residueCoefficient 1 15 := by
      simp [residueCoefficient, Nat.add_mod]
    _ = -24 := coefficient_one_fifteen

/--
The infinite exact negative family from THM-948:
`M({1,2,15+28m}) = -12/(343(15+28m))` at the closed-form interface.
-/
theorem closed_residue_mass_one_fifteen_add_twenty_eight_mul (m : ℕ) :
    closedResidueMass 1 (15 + 28 * m) =
      -12 / (343 * (15 + 28 * m) : ℚ) := by
  rw [closedResidueMass, coefficient_one_fifteen_add_twenty_eight_mul]
  push_cast
  have hm : (15 : ℚ) + 28 * (m : ℚ) ≠ 0 := by positivity
  field_simp [hm]
  norm_num

/-- The coefficient is 14-periodic in its first input. -/
@[simp] theorem residue_coefficient_add_fourteen (A B : ℕ) :
    residueCoefficient (A + 14) B = residueCoefficient A B := by
  simp [residueCoefficient]

/-- The coefficient is 28-periodic in its second input. -/
@[simp] theorem residue_coefficient_add_twenty_eight (A B : ℕ) :
    residueCoefficient A (B + 28) = residueCoefficient A B := by
  simp [residueCoefficient]

/-- Reflection in the first residue circle reverses the arbitrary-residue coefficient. -/
theorem residue_coefficient_reflect_fourteen (A B : ℕ) :
    residueCoefficient ((14 - A % 14) % 14) B = -residueCoefficient A B := by
  let a : Fin 14 := ⟨A % 14, Nat.mod_lt _ (by omega)⟩
  let b : Fin 28 := ⟨B % 28, Nat.mod_lt _ (by omega)⟩
  simpa [residueCoefficient, a, b, Nat.mod_mod] using cell_reflect_fourteen a b

/-- Reflection in the second residue circle reverses the arbitrary-residue coefficient. -/
theorem residue_coefficient_reflect_twenty_eight (A B : ℕ) :
    residueCoefficient A ((28 - B % 28) % 28) = -residueCoefficient A B := by
  let a : Fin 14 := ⟨A % 14, Nat.mod_lt _ (by omega)⟩
  let b : Fin 28 := ⟨B % 28, Nat.mod_lt _ (by omega)⟩
  simpa [residueCoefficient, a, b, Nat.mod_mod] using cell_reflect_twenty_eight a b

#print axioms cell_coefficient_integral
#print axioms cell_coefficient_zero_iff
#print axioms cell_coefficient_sign
#print axioms residue_coefficient_zero_iff_dvd
#print axioms closed_residue_mass_one_fifteen_add_twenty_eight_mul

end ExactDoublingTriple
end LonelyRunner
