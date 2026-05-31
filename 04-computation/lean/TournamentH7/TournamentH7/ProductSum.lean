/-
  TournamentH7.ProductSum -- additive/multiplicative defect normal form.

  This module formalizes the list-level core of THM-361.  Removing all
  entries equal to 1 preserves the product and records exactly how much
  additive slack was removed from the sum.
-/

import Mathlib.Tactic

namespace Tournament
namespace ProductSum

/-- The number of entries equal to 1. -/
def ones (xs : List Nat) : Nat :=
  (xs.filter (fun x => x = 1)).length

/-- The product-sum core obtained by deleting all entries equal to 1. -/
def core (xs : List Nat) : List Nat :=
  xs.filter (fun x => x ≠ 1)

/-- A finite list is product-sum when its product equals its sum. -/
def IsProductSum (xs : List Nat) : Prop :=
  xs.prod = xs.sum

/-- Removing ones does not change the product. -/
theorem prod_eq_core_prod (xs : List Nat) :
    xs.prod = (core xs).prod := by
  induction xs with
  | nil =>
      simp [core]
  | cons x xs ih =>
      by_cases hx : x = 1
      · simp [core, hx, ih]
      · simp [core, hx, ih]

/-- Removing ones subtracts exactly `ones xs` from the sum. -/
theorem sum_eq_core_sum_add_ones (xs : List Nat) :
    xs.sum = (core xs).sum + ones xs := by
  induction xs with
  | nil =>
      simp [core, ones]
  | cons x xs ih =>
      by_cases hx : x = 1
      · subst x
        simp [core, ones, ih, Nat.add_assoc, Nat.add_comm]
      · simp [core, ones, hx, ih, Nat.add_comm, Nat.add_left_comm]

/--
Product-sum normal form: a list is product-sum iff its core has product equal
to core sum plus the number of removed ones.
-/
theorem product_sum_iff_core (xs : List Nat) :
    IsProductSum xs ↔ (core xs).prod = (core xs).sum + ones xs := by
  unfold IsProductSum
  rw [prod_eq_core_prod xs, sum_eq_core_sum_add_ones xs]

/-- Padding a core by ones preserves product. -/
theorem prod_replicate_one_append (d : Nat) (c : List Nat) :
    ((List.replicate d 1) ++ c).prod = c.prod := by
  simp

/-- Padding a core by `d` ones adds exactly `d` to the sum. -/
theorem sum_replicate_one_append (d : Nat) (c : List Nat) :
    ((List.replicate d 1) ++ c).sum = d + c.sum := by
  simp

/-- Any nonnegative defect can be repaired by adjoining that many ones. -/
theorem pad_core_product_sum {d : Nat} {c : List Nat}
    (h : c.prod = c.sum + d) :
    IsProductSum ((List.replicate d 1) ++ c) := by
  unfold IsProductSum
  rw [prod_replicate_one_append d c, sum_replicate_one_append d c]
  simpa [Nat.add_comm] using h

/-- Product-sum lists have core defect equal to their number of ones. -/
theorem core_defect_eq_ones_of_product_sum {xs : List Nat}
    (h : IsProductSum xs) :
    (core xs).prod - (core xs).sum = ones xs := by
  have hcore := (product_sum_iff_core xs).mp h
  omega

/--
The subtraction form of the normal form, with the usual natural-number side
condition that the core sum is at most the core product.
-/
theorem product_sum_of_core_defect_eq_ones {xs : List Nat}
    (hle : (core xs).sum ≤ (core xs).prod)
    (h : (core xs).prod - (core xs).sum = ones xs) :
    IsProductSum xs := by
  apply (product_sum_iff_core xs).mpr
  omega

/-- The ordered positive two-entry product-sum resonance is only `(2,2)`. -/
theorem two_entry_product_sum {a b : Nat} (ha : 0 < a) (hb : 0 < b)
    (h : IsProductSum [a, b]) :
    a = 2 ∧ b = 2 := by
  unfold IsProductSum at h
  simp at h
  rcases Nat.exists_eq_succ_of_ne_zero (Nat.ne_of_gt ha) with ⟨a0, rfl⟩
  rcases Nat.exists_eq_succ_of_ne_zero (Nat.ne_of_gt hb) with ⟨b0, rfl⟩
  simp at h
  have hprod : a0 * b0 = 1 := by
    nlinarith
  have ha0 : a0 = 1 := Nat.eq_one_of_dvd_one ⟨b0, hprod.symm⟩
  have hb0 : b0 = 1 := Nat.eq_one_of_dvd_one ⟨a0, by
    rw [Nat.mul_comm]
    exact hprod.symm⟩
  constructor <;> omega

/-- Positive ordered pairs are product-sum iff they are `(2,2)`. -/
theorem two_entry_product_sum_iff {a b : Nat} (ha : 0 < a) (hb : 0 < b) :
    IsProductSum [a, b] ↔ a = 2 ∧ b = 2 := by
  constructor
  · exact two_entry_product_sum ha hb
  · rintro ⟨rfl, rfl⟩
    simp [IsProductSum]

end ProductSum
end Tournament
