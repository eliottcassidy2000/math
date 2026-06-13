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

/-- Three product-sum entries cannot all be at least `2`. -/
theorem no_three_ge_two_product_sum {a b c : Nat}
    (ha : 2 ≤ a) (hb : 2 ≤ b) (hc : 2 ≤ c)
    (h : IsProductSum [a, b, c]) : False := by
  unfold IsProductSum at h
  simp at h
  have h' : a * b * c = a + b + c := by
    simpa [Nat.mul_assoc, Nat.add_assoc] using h
  have h4a : 4 * a ≤ a * b * c := by
    have hbc4 : 4 ≤ b * c := by
      have := Nat.mul_le_mul hb hc
      simpa using this
    calc
      4 * a = a * 4 := by rw [Nat.mul_comm]
      _ ≤ a * (b * c) := Nat.mul_le_mul_left a hbc4
      _ = a * b * c := by rw [Nat.mul_assoc]
  have h4b : 4 * b ≤ a * b * c := by
    have hac4 : 4 ≤ a * c := by
      have := Nat.mul_le_mul ha hc
      simpa using this
    calc
      4 * b = b * 4 := by rw [Nat.mul_comm]
      _ ≤ b * (a * c) := Nat.mul_le_mul_left b hac4
      _ = a * b * c := by ring
  have h4c : 4 * c ≤ a * b * c := by
    have hab4 : 4 ≤ a * b := by
      have := Nat.mul_le_mul ha hb
      simpa using this
    calc
      4 * c = c * 4 := by rw [Nat.mul_comm]
      _ ≤ c * (a * b) := Nat.mul_le_mul_left c hab4
      _ = a * b * c := by ring
  rw [h'] at h4a h4b h4c
  omega

/-- With one explicit unit, the remaining positive nonunit pair is `(2,3)`
or `(3,2)`. -/
theorem one_cons_two_entry_product_sum {a b : Nat}
    (ha : 0 < a) (hb : 0 < b) (ha1 : a ≠ 1) (hb1 : b ≠ 1)
    (h : IsProductSum [1, a, b]) :
    (a = 2 ∧ b = 3) ∨ (a = 3 ∧ b = 2) := by
  have ha2 : 2 ≤ a := by
    by_contra hlt
    have hle : a ≤ 1 := by omega
    interval_cases a; simp_all
  have hb2 : 2 ≤ b := by
    by_contra hlt
    have hle : b ≤ 1 := by omega
    interval_cases b; simp_all
  let a0 := a - 2
  let b0 := b - 2
  have ha_eq : a = a0 + 2 := by omega
  have hb_eq : b = b0 + 2 := by omega
  unfold IsProductSum at h
  rw [ha_eq, hb_eq] at h
  simp at h
  have hab : (a0 + 1) * (b0 + 1) = 2 := by
    nlinarith
  have hb0_pos : 0 < b0 + 1 := by omega
  have ha0_pos : 0 < a0 + 1 := by omega
  have ha0_le1 : a0 ≤ 1 := by
    have hle : a0 + 1 ≤ (a0 + 1) * (b0 + 1) :=
      Nat.le_mul_of_pos_right (a0 + 1) hb0_pos
    omega
  have hb0_le1 : b0 ≤ 1 := by
    have hle : b0 + 1 ≤ (a0 + 1) * (b0 + 1) :=
      Nat.le_mul_of_pos_left (b0 + 1) ha0_pos
    omega
  interval_cases a0 <;> interval_cases b0 <;> simp_all

/-- The only ordered positive, pairwise-distinct, three-entry product-sum
lists are permutations of `[1, 2, 3]`. -/
theorem three_entry_distinct_product_sum {a b c : Nat}
    (ha : 0 < a) (hb : 0 < b) (hc : 0 < c)
    (hab : a ≠ b) (hac : a ≠ c) (hbc : b ≠ c)
    (h : IsProductSum [a, b, c]) :
    List.Perm [a, b, c] [1, 2, 3] := by
  by_cases ha1 : a = 1
  · subst a
    have hb1 : b ≠ 1 := by
      intro hb_eq
      exact hab hb_eq.symm
    have hc1 : c ≠ 1 := by
      intro hc_eq
      exact hac hc_eq.symm
    have h' : IsProductSum [1, b, c] := by simpa using h
    rcases one_cons_two_entry_product_sum hb hc hb1 hc1 h' with h23 | h32
    · rcases h23 with ⟨rfl, rfl⟩
      decide
    · rcases h32 with ⟨rfl, rfl⟩
      decide
  by_cases hb1 : b = 1
  · subst b
    have ha1' : a ≠ 1 := ha1
    have hc1 : c ≠ 1 := by
      intro hc_eq
      exact hbc hc_eq.symm
    have h' : IsProductSum [1, a, c] := by
      simpa [IsProductSum, Nat.add_assoc, Nat.add_comm, Nat.add_left_comm,
        Nat.mul_assoc, Nat.mul_comm, Nat.mul_left_comm] using h
    rcases one_cons_two_entry_product_sum ha hc ha1' hc1 h' with h23 | h32
    · rcases h23 with ⟨rfl, rfl⟩
      decide
    · rcases h32 with ⟨rfl, rfl⟩
      decide
  by_cases hc1 : c = 1
  · subst c
    have ha1' : a ≠ 1 := ha1
    have hb1' : b ≠ 1 := hb1
    have h' : IsProductSum [1, a, b] := by
      simpa [IsProductSum, Nat.add_assoc, Nat.add_comm, Nat.add_left_comm,
        Nat.mul_assoc, Nat.mul_comm, Nat.mul_left_comm] using h
    rcases one_cons_two_entry_product_sum ha hb ha1' hb1' h' with h23 | h32
    · rcases h23 with ⟨rfl, rfl⟩
      decide
    · rcases h32 with ⟨rfl, rfl⟩
      decide
  exfalso
  have ha2 : 2 ≤ a := by
    by_contra hlt
    have hle : a ≤ 1 := by omega
    interval_cases a; simp_all
  have hb2 : 2 ≤ b := by
    by_contra hlt
    have hle : b ≤ 1 := by omega
    interval_cases b; simp_all
  have hc2 : 2 ≤ c := by
    by_contra hlt
    have hle : c ≤ 1 := by omega
    interval_cases c; simp_all
  exact no_three_ge_two_product_sum ha2 hb2 hc2 h

end ProductSum
end Tournament
