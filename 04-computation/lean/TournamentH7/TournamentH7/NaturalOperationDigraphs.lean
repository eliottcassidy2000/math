/-
  TournamentH7.NaturalOperationDigraphs -- natural-number operation shadows

  This module formalizes the S365 observation that the simple one-input
  shadow of the additive gate `{x,y} -> x+y` is just the strict order on
  positive naturals, while the corresponding multiplication shadow is the
  divisibility DAG.

  Product-sum equations enter as critical pairs between the two operation
  systems.  The two-nonunit stratum has the clean divisor form

      r + (a+1) + (b+1) = (a+1)(b+1)  iff  a*b = r+1.
-/

import Mathlib.Tactic

namespace NatOperation

/-- Additive one-input shadow: `x` points to `z` if some positive `y`
    satisfies `x+y=z`. -/
def AddShadow (x z : ℕ) : Prop :=
  ∃ y : ℕ, 1 ≤ y ∧ x + y = z

/-- Multiplicative one-input shadow with the unit allowed. -/
def MulUnitShadow (x z : ℕ) : Prop :=
  ∃ y : ℕ, 1 ≤ y ∧ x * y = z

/-- Multiplicative one-input shadow using a nonunit second input. -/
def MulShadow (x z : ℕ) : Prop :=
  ∃ y : ℕ, 2 ≤ y ∧ x * y = z

/-- The additive shadow is exactly the strict natural-number order. -/
theorem addShadow_iff_lt (x z : ℕ) :
    AddShadow x z ↔ x < z := by
  constructor
  · intro h
    rcases h with ⟨y, hy, hxy⟩
    omega
  · intro hx
    refine ⟨z - x, ?_, ?_⟩
    · omega
    · omega

/-- With the multiplicative unit allowed, the multiplication shadow is
    divisibility, on positive targets. -/
theorem mulUnitShadow_iff_dvd {x z : ℕ} (hz : 1 ≤ z) :
    MulUnitShadow x z ↔ x ∣ z := by
  constructor
  · intro h
    rcases h with ⟨y, _hy, hxy⟩
    exact ⟨y, hxy.symm⟩
  · intro h
    rcases h with ⟨y, hzy⟩
    cases y with
    | zero =>
        simp at hzy
        omega
    | succ y =>
        refine ⟨Nat.succ y, by omega, ?_⟩
        exact hzy.symm

/-- With a nonunit multiplier, the multiplication shadow is proper
    divisibility on positive sources. -/
theorem mulShadow_iff_dvd_and_lt {x z : ℕ} (hx : 1 ≤ x) :
    MulShadow x z ↔ x ∣ z ∧ x < z := by
  constructor
  · intro h
    rcases h with ⟨y, hy, hxy⟩
    constructor
    · exact ⟨y, hxy.symm⟩
    · have hy1 : 1 < y := by omega
      have hxpos : 0 < x := by omega
      have hlt : x * 1 < x * y := Nat.mul_lt_mul_of_pos_left hy1 hxpos
      rw [Nat.mul_one] at hlt
      rw [hxy] at hlt
      exact hlt
  · intro h
    rcases h with ⟨hdvd, hlt⟩
    rcases hdvd with ⟨y, hzy⟩
    cases y with
    | zero =>
        simp at hzy
        omega
    | succ y =>
        cases y with
        | zero =>
            simp at hzy
            omega
        | succ y =>
            refine ⟨Nat.succ (Nat.succ y), by omega, ?_⟩
            exact hzy.symm

/-- Binary product-sum collision in shifted positive variables.

    Setting `x=a+1`, `y=b+1`, the equation `x+y=xy` is equivalent to
    `a*b=1`. -/
theorem shifted_binary_collision_iff (a b : ℕ) :
    (a + 1) + (b + 1) = (a + 1) * (b + 1) ↔ a * b = 1 := by
  constructor <;> intro h <;> nlinarith

/-- The positive binary product-sum collision is unique. -/
theorem shifted_binary_collision_unique {a b : ℕ}
    (h : (a + 1) + (b + 1) = (a + 1) * (b + 1)) :
    a = 1 ∧ b = 1 := by
  have hab : a * b = 1 := (shifted_binary_collision_iff a b).mp h
  have ha_pos : 1 ≤ a := by
    cases a with
    | zero =>
        simp at hab
    | succ a =>
        omega
  have hb_pos : 1 ≤ b := by
    cases b with
    | zero =>
        simp at hab
    | succ b =>
        omega
  constructor
  · by_contra hne
    have ha2 : 2 ≤ a := by omega
    have hge : 2 ≤ a * b := by
      simpa using Nat.mul_le_mul ha2 hb_pos
    omega
  · by_contra hne
    have hb2 : 2 ≤ b := by omega
    have hge : 2 ≤ a * b := by
      simpa using Nat.mul_le_mul ha_pos hb2
    omega

/-- The two-nonunit product-sum layer is a divisor equation.

    Here `r` is the number of unit paddings and the nonunit factors are
    `a+1` and `b+1`. -/
theorem twoFactor_productSum_iff (r a b : ℕ) :
    r + (a + 1) + (b + 1) = (a + 1) * (b + 1) ↔ a * b = r + 1 := by
  constructor <;> intro h <;> nlinarith

/-- The universal binary witness: `r` ones plus factors `2` and `r+2`
    always solve `sum=product`.  In arity `k=r+2`, this gives the
    computational bound `m(k) <= 2k`. -/
theorem trivial_twoFactor_productSum (r : ℕ) :
    r + (1 + 1) + ((r + 1) + 1) = (1 + 1) * ((r + 1) + 1) := by
  nlinarith

/-- Product-minus-sum defect: if a nonunit seed has product at least its sum,
    the defect is exactly the number of ones needed to repair the additive
    fold to the multiplicative endpoint. -/
def SeedDefect (L : List ℕ) : ℕ :=
  L.prod - L.sum

theorem seedDefect_add_sum {L : List ℕ} (h : L.sum ≤ L.prod) :
    SeedDefect L + L.sum = L.prod := by
  unfold SeedDefect
  exact Nat.sub_add_cancel h

theorem sum_replicate_seedDefect_ones_append {L : List ℕ} (h : L.sum ≤ L.prod) :
    ((List.replicate (SeedDefect L) 1) ++ L).sum = L.prod := by
  simp [SeedDefect, Nat.sub_add_cancel h]

theorem prod_replicate_seedDefect_ones_append (L : List ℕ) :
    ((List.replicate (SeedDefect L) 1) ++ L).prod = L.prod := by
  simp [SeedDefect]

end NatOperation
