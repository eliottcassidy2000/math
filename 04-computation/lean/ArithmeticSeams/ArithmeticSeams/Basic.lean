import Std

set_option autoImplicit false

namespace ArithmeticSeams

/-- The shift identifies natural-number labels with positive cardinalities. -/
def shift (a : Nat) : Nat := a + 1

def augAdd (a b : Nat) : Nat := a + b + 1

def augMul (a b : Nat) : Nat := a * b + a + b

theorem shift_injective {a b : Nat} (h : shift a = shift b) : a = b := by
  unfold shift at h
  omega

theorem shift_positive (a : Nat) : 0 < shift a := by
  unfold shift
  omega

theorem shift_hits_exactly_positive (m : Nat) :
    (∃ a, shift a = m) ↔ 0 < m := by
  constructor
  · intro ⟨a, h⟩
    rw [← h]
    exact shift_positive a
  · intro h
    exact ⟨m - 1, by unfold shift; omega⟩

theorem shift_augAdd (a b : Nat) :
    shift (augAdd a b) = shift a + shift b := by
  unfold shift augAdd
  omega

theorem shift_augMul (a b : Nat) :
    shift (augMul a b) = shift a * shift b := by
  unfold shift augMul
  simp only [Nat.add_mul, Nat.mul_add, Nat.one_mul, Nat.mul_one]
  omega

theorem augAdd_comm (a b : Nat) : augAdd a b = augAdd b a := by
  apply shift_injective
  rw [shift_augAdd, shift_augAdd, Nat.add_comm]

theorem augAdd_assoc (a b c : Nat) :
    augAdd (augAdd a b) c = augAdd a (augAdd b c) := by
  apply shift_injective
  simp only [shift_augAdd]
  exact Nat.add_assoc _ _ _

theorem augMul_comm (a b : Nat) : augMul a b = augMul b a := by
  apply shift_injective
  rw [shift_augMul, shift_augMul, Nat.mul_comm]

theorem augMul_assoc (a b c : Nat) :
    augMul (augMul a b) c = augMul a (augMul b c) := by
  apply shift_injective
  simp only [shift_augMul]
  exact Nat.mul_assoc _ _ _

theorem augMul_left_distrib (a b c : Nat) :
    augMul a (augAdd b c) = augAdd (augMul a b) (augMul a c) := by
  apply shift_injective
  simp only [shift_augMul, shift_augAdd]
  exact Nat.mul_add _ _ _

theorem augMul_right_distrib (a b c : Nat) :
    augMul (augAdd a b) c = augAdd (augMul a c) (augMul b c) := by
  apply shift_injective
  simp only [shift_augMul, shift_augAdd]
  exact Nat.add_mul _ _ _

theorem augMul_zero_left (a : Nat) : augMul 0 a = a := by
  simp only [augMul, Nat.zero_mul, Nat.zero_add]

theorem augMul_zero_right (a : Nat) : augMul a 0 = a := by
  simp only [augMul, Nat.mul_zero, Nat.zero_add, Nat.add_zero]

/-- On Nat, augmented addition has no identity; one must adjoin the label -1. -/
theorem augAdd_has_no_identity : ¬ ∃ e : Nat, ∀ a, augAdd e a = a := by
  intro ⟨e, h⟩
  have hzero := h 0
  unfold augAdd at hzero
  omega

theorem augAdd_diagonal (a : Nat) : augAdd a a = 2 * a + 1 := by
  unfold augAdd
  omega

theorem augMul_diagonal (a : Nat) : augMul a a = a ^ 2 + 2 * a := by
  simp only [augMul, Nat.pow_succ, Nat.pow_zero, Nat.one_mul]
  omega

/-- Zeroth power uses the multiplicative identity, whose augmented label is 0. -/
def augPow (a : Nat) : Nat → Nat
  | 0 => 0
  | n + 1 => augMul a (augPow a n)

theorem shift_augPow (a n : Nat) : shift (augPow a n) = (a + 1) ^ n := by
  induction n with
  | zero => rfl
  | succ n ih =>
    change shift (augMul a (augPow a n)) = (a + 1) ^ (n + 1)
    rw [shift_augMul, ih]
    change (a + 1) * (a + 1) ^ n = (a + 1) ^ (n + 1)
    rw [Nat.pow_succ, Nat.mul_comm]

theorem augPow_add (a m n : Nat) :
    augPow a (m + n) = augMul (augPow a m) (augPow a n) := by
  apply shift_injective
  rw [shift_augMul, shift_augPow, shift_augPow, shift_augPow]
  exact Nat.pow_add _ _ _

theorem augPow_zero_base (n : Nat) : augPow 0 n = 0 := by
  induction n with
  | zero => rfl
  | succ n ih =>
    simp only [augPow, augMul_zero_left, ih]

theorem pointed_product_control :
    augMul 2 3 = 11 ∧ shift (augMul 2 3) = 12 ∧
    augPow 1 6 = 63 ∧ augPow 1 16 = 65535 := by
  decide

end ArithmeticSeams
