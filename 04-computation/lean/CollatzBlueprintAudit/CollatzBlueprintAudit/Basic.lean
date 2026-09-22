import Std

set_option autoImplicit false

namespace CollatzBlueprintAudit

/-! ## The blueprint's impossible all-integer residue premise -/

def AllIntegersHaveOddResidue : Prop :=
  ∀ n : Int, n % 6 = 1 ∨ n % 6 = 3 ∨ n % 6 = 5

theorem not_allIntegersHaveOddResidue : ¬ AllIntegersHaveOddResidue := by
  intro h
  have hzero := h 0
  omega

/-- This retains only the problematic field, without the template's syntax errors. -/
structure BlueprintResiduePremise where
  monodromy_mod6 : AllIntegersHaveOddResidue

theorem blueprintResiduePremise_uninhabited : ¬ Nonempty BlueprintResiduePremise := by
  intro ⟨m⟩
  exact not_allIntegersHaveOddResidue m.monodromy_mod6

/-- The corrected statement restricts the quantifier to odd integers. -/
theorem odd_integer_residues (n : Int) (hn : n % 2 = 1) :
    n % 6 = 1 ∨ n % 6 = 3 ∨ n % 6 = 5 := by
  omega

/-! ## An elementary iteration API and a general descent criterion -/

def iterate {α : Type} (f : α → α) : Nat → α → α
  | 0, x => x
  | k + 1, x => f (iterate f k x)

theorem iterate_add {α : Type} (f : α → α) (a b : Nat) (x : α) :
    iterate f (a + b) x = iterate f b (iterate f a x) := by
  induction b with
  | zero => rfl
  | succ b ih => exact congrArg f ih

def PositivePreserving (f : Nat → Nat) : Prop :=
  ∀ n, 0 < n → 0 < f n

theorem iterate_positive (f : Nat → Nat) (hf : PositivePreserving f)
    (k n : Nat) (hn : 0 < n) : 0 < iterate f k n := by
  induction k with
  | zero => exact hn
  | succ k ih => exact hf (iterate f k n) ih

def EventuallyReachesOne (f : Nat → Nat) : Prop :=
  ∀ n, 0 < n → ∃ k, iterate f k n = 1

def EventuallyStrictlyDescends (f : Nat → Nat) : Prop :=
  ∀ n, 1 < n → ∃ k, 0 < k ∧ iterate f k n < n

/-- Strict descent is required at some positive time depending on the start.
The proof uses strong induction, not a density or average-drift assertion. -/
theorem reachesOne_iff_strictDescent (f : Nat → Nat) (hf : PositivePreserving f) :
    EventuallyReachesOne f ↔ EventuallyStrictlyDescends f := by
  constructor
  · intro h n hn
    obtain ⟨k, hk⟩ := h n (by omega)
    refine ⟨k, ?_, ?_⟩
    · cases k with
      | zero => simp only [iterate] at hk; omega
      | succ k => omega
    · rw [hk]
      exact hn
  · intro h n
    induction n using Nat.strongRecOn with
    | ind n ih =>
      intro hn
      by_cases h1 : n = 1
      · exact ⟨0, h1⟩
      · obtain ⟨k, _hk, hsmall⟩ := h n (by omega)
        have hpos := iterate_positive f hf k n hn
        obtain ⟨j, hj⟩ := ih (iterate f k n) hsmall hpos
        exact ⟨k + j, by rw [iterate_add]; exact hj⟩

/-! ## Exact affine-word certificates: retain the carry term -/

/-- The orbit identity is an explicit hypothesis; no word realization is assumed. -/
theorem descent_iff_affine_inequality (n y K L B : Nat)
    (hidentity : 2 ^ K * y = 3 ^ L * n + B) :
    y < n ↔ 3 ^ L * n + B < 2 ^ K * n := by
  constructor
  · intro h
    rw [← hidentity]
    exact Nat.mul_lt_mul_of_pos_left h (Nat.two_pow_pos K)
  · intro h
    by_cases hy : y < n
    · exact hy
    · have hny : n ≤ y := by omega
      have hle := Nat.mul_le_mul_left (2 ^ K) hny
      rw [hidentity] at hle
      omega

/-- A future certificate must provide the true endpoint AND the full inequality. -/
theorem strictDescent_of_affineCertificates (f : Nat → Nat)
    (hcert : ∀ n, 1 < n → ∃ k K L B,
      0 < k ∧
      2 ^ K * iterate f k n = 3 ^ L * n + B ∧
      3 ^ L * n + B < 2 ^ K * n) :
    EventuallyStrictlyDescends f := by
  intro n hn
  obtain ⟨k, K, L, B, hk, hid, hlt⟩ := hcert n hn
  exact ⟨k, hk, (descent_iff_affine_inequality n (iterate f k n) K L B hid).mpr hlt⟩

/-! ## Standard unaccelerated Collatz on natural numbers -/

def collatz (n : Nat) : Nat :=
  if n % 2 = 0 then n / 2 else 3 * n + 1

theorem collatz_positive : PositivePreserving collatz := by
  intro n hn
  unfold collatz
  split <;> omega

/-- A proved equivalence; neither side is proved here. -/
theorem collatz_reachesOne_iff_strictDescent :
    EventuallyReachesOne collatz ↔ EventuallyStrictlyDescends collatz :=
  reachesOne_iff_strictDescent collatz collatz_positive

/-- Explicit conditional interface for future work on the still-open descent side. -/
theorem collatz_reachesOne_of_strictDescent
    (h : EventuallyStrictlyDescends collatz) : EventuallyReachesOne collatz :=
  collatz_reachesOne_iff_strictDescent.mpr h

/-! ## Positive and hostile controls -/

/-- The generic positivity and reachability hypotheses have a concrete model. -/
theorem constantOne_control :
    PositivePreserving (fun _ : Nat => 1) ∧
    EventuallyReachesOne (fun _ : Nat => 1) ∧
    EventuallyStrictlyDescends (fun _ : Nat => 1) := by
  have hp : PositivePreserving (fun _ : Nat => 1) := by
    intro n hn
    exact Nat.zero_lt_succ 0
  have hr : EventuallyReachesOne (fun _ : Nat => 1) := by
    intro n hn
    exact ⟨1, rfl⟩
  exact ⟨hp, hr, (reachesOne_iff_strictDescent _ hp).mp hr⟩

/-- Without positive preservation, descent can mean falling into zero. -/
theorem positivity_hypothesis_is_necessary :
    EventuallyStrictlyDescends (fun _ : Nat => 0) ∧
    ¬ EventuallyReachesOne (fun _ : Nat => 0) := by
  constructor
  · intro n hn
    refine ⟨1, Nat.zero_lt_succ 0, ?_⟩
    change 0 < n
    omega
  · intro h
    obtain ⟨k, hk⟩ := h 2 (by decide)
    cases k with
    | zero => change 2 = 1 at hk; contradiction
    | succ k => change 0 = 1 at hk; contradiction

theorem collatz_seven_reachesOne : iterate collatz 16 7 = 1 := by decide

/-- On the odd accelerated orbit, this is one step; on the standard map, two. -/
theorem collatz_nineteen_growth_control :
    iterate collatz 2 19 = 29 ∧ (3 * 19 + 1) % 4 = 2 ∧ 19 < (29 : Nat) := by
  decide

theorem collatz_zero_fixed (k : Nat) : iterate collatz k 0 = 0 := by
  induction k with
  | zero => rfl
  | succ k ih => simp only [iterate, ih, collatz]; rfl

/-- The restriction to positive starts cannot simply be deleted. -/
theorem collatz_zero_never_reachesOne : ¬ ∃ k, iterate collatz k 0 = 1 := by
  intro ⟨k, hk⟩
  rw [collatz_zero_fixed] at hk
  contradiction

/-! ## Exact signed return examples, without an exhaustion claim -/

def signedCollatz (n : Int) : Int :=
  if n % 2 = 0 then n / 2 else 3 * n + 1

theorem signed_negOne_returns : iterate signedCollatz 2 (-1) = -1 := by decide

theorem signed_negFive_returns : iterate signedCollatz 5 (-5) = -5 := by decide

theorem signed_negSeventeen_returns : iterate signedCollatz 18 (-17) = -17 := by decide

end CollatzBlueprintAudit
