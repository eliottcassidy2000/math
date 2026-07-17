import Mathlib

/-!
# Continuous relation lock

At an arbitrary real phase, integer witnesses inside the open `1/14` danger
bands inherit every integer speed relation of coefficient weight at most
fourteen.  This is the continuous counterpart of `LRCRelationLock`, whose
sampled statement is phrased at `p/q`.

The proof carrier is a relation circuit with real errors and integer witness
labels.  A runner tournament forgets the coefficients and cannot preserve the
strict total error budget.  Relation circuits, rather than runners, are the
useful vertices for the selected-witness exact-frequency branches.
-/

namespace LonelyRunner
namespace LRCRealRelationLock

open Finset
open scoped BigOperators

/-- Integer witnesses in open danger bands inherit every speed relation whose
total absolute coefficient weight is at most fourteen. -/
theorem real_relation_lock
    {ι : Type*} (support : Finset ι)
    (coefficient speed witness : ι → ℤ) (phase : ℝ)
    (hrelation : ∑ index ∈ support,
      coefficient index * speed index = 0)
    (hweight : ∑ index ∈ support, |coefficient index| ≤ 14)
    (hbad : ∀ index ∈ support,
      |(speed index : ℝ) * phase - witness index| < 1 / 14) :
    ∑ index ∈ support, coefficient index * witness index = 0 := by
  by_cases hnonzero : ∃ index ∈ support, coefficient index ≠ 0
  · have htermLe : ∀ index ∈ support,
        |(coefficient index : ℝ) *
            ((speed index : ℝ) * phase - witness index)| ≤
          (|coefficient index| : ℝ) * (1 / 14) := by
      intro index hindex
      rw [abs_mul]
      exact mul_le_mul_of_nonneg_left (le_of_lt (hbad index hindex))
        (abs_nonneg _)
    have htermStrict : ∃ index ∈ support,
        |(coefficient index : ℝ) *
            ((speed index : ℝ) * phase - witness index)| <
          (|coefficient index| : ℝ) * (1 / 14) := by
      obtain ⟨index, hindex, hcoefficient⟩ := hnonzero
      refine ⟨index, hindex, ?_⟩
      rw [abs_mul]
      exact mul_lt_mul_of_pos_left (hbad index hindex)
        (by exact_mod_cast abs_pos.mpr hcoefficient)
    have hsumStrict :
        (∑ index ∈ support,
          |(coefficient index : ℝ) *
            ((speed index : ℝ) * phase - witness index)|) <
        ∑ index ∈ support,
          (|coefficient index| : ℝ) * (1 / 14) :=
      Finset.sum_lt_sum htermLe htermStrict
    have hright :
        (∑ index ∈ support,
          (|coefficient index| : ℝ) * (1 / 14)) ≤ 1 := by
      rw [← Finset.sum_mul]
      have habsCast :
          (∑ index ∈ support, |(coefficient index : ℝ)|) =
            ((∑ index ∈ support, |coefficient index| : ℤ) : ℝ) := by
        rw [Int.cast_sum]
        apply Finset.sum_congr rfl
        intro index _hindex
        exact Int.cast_abs.symm
      rw [habsCast]
      have hweightReal :
          ((∑ index ∈ support, |coefficient index| : ℤ) : ℝ) ≤ 14 := by
        exact_mod_cast hweight
      norm_num
      nlinarith
    have hidentity :
        ∑ index ∈ support,
          (coefficient index : ℝ) *
            ((speed index : ℝ) * phase - witness index) =
          -((∑ index ∈ support,
            coefficient index * witness index : ℤ) : ℝ) := by
      calc
        (∑ index ∈ support,
          (coefficient index : ℝ) *
            ((speed index : ℝ) * phase - witness index)) =
            ((∑ index ∈ support,
              coefficient index * speed index : ℤ) : ℝ) * phase -
            ((∑ index ∈ support,
              coefficient index * witness index : ℤ) : ℝ) := by
          push_cast
          rw [Finset.sum_mul, ← Finset.sum_sub_distrib]
          apply Finset.sum_congr rfl
          intro index _hindex
          ring
        _ = -((∑ index ∈ support,
            coefficient index * witness index : ℤ) : ℝ) := by
          rw [hrelation]
          norm_num
    have habsReal :
        |((∑ index ∈ support,
          coefficient index * witness index : ℤ) : ℝ)| < 1 := by
      rw [← abs_neg,
        ← hidentity]
      exact lt_of_le_of_lt
        (Finset.abs_sum_le_sum_abs
          (fun index =>
            (coefficient index : ℝ) *
              ((speed index : ℝ) * phase - witness index)) support)
        (lt_of_lt_of_le hsumStrict hright)
    have habsInteger :
        |∑ index ∈ support,
          coefficient index * witness index| < 1 := by
      exact_mod_cast habsReal
    exact Int.abs_lt_one_iff.mp habsInteger
  · push Not at hnonzero
    have hcoeffZero : ∀ index ∈ support, coefficient index = 0 :=
      hnonzero
    apply Finset.sum_eq_zero
    intro index hindex
    rw [hcoeffZero index hindex, zero_mul]

/-- Three-term workhorse for selected-witness scalar walls. -/
theorem real_relation_lock3
    (speed₁ speed₂ speed₃ witness₁ witness₂ witness₃
      coefficient₁ coefficient₂ coefficient₃ : ℤ)
    (phase : ℝ)
    (hrelation : coefficient₁ * speed₁ + coefficient₂ * speed₂ +
      coefficient₃ * speed₃ = 0)
    (hweight : |coefficient₁| + |coefficient₂| + |coefficient₃| ≤ 14)
    (hbad₁ : |(speed₁ : ℝ) * phase - witness₁| < 1 / 14)
    (hbad₂ : |(speed₂ : ℝ) * phase - witness₂| < 1 / 14)
    (hbad₃ : |(speed₃ : ℝ) * phase - witness₃| < 1 / 14) :
    coefficient₁ * witness₁ + coefficient₂ * witness₂ +
      coefficient₃ * witness₃ = 0 := by
  have hsum := real_relation_lock
    ({0, 1, 2} : Finset (Fin 3))
    (fun index => ![coefficient₁, coefficient₂, coefficient₃] index)
    (fun index => ![speed₁, speed₂, speed₃] index)
    (fun index => ![witness₁, witness₂, witness₃] index)
    phase
  simpa [Fin.sum_univ_three, add_assoc] using hsum
    (by simpa [Fin.sum_univ_three, add_assoc] using hrelation)
    (by simpa [Fin.sum_univ_three, add_assoc] using hweight)
    (by
      intro index hindex
      fin_cases index <;> simp_all)

#print axioms real_relation_lock
#print axioms real_relation_lock3

end LRCRealRelationLock
end LonelyRunner
