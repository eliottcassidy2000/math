import CollatzBlueprintAudit.Basic

set_option autoImplicit false

namespace CollatzBlueprintAudit

/-!
The standard map has two kinds of step. Total time is not the number of
tripling operations. These counters follow the actual orbit, rather than
being freely chosen existential certificate data.
-/

def evenSteps (n : Nat) : Nat → Nat
  | 0 => 0
  | t + 1 =>
    if iterate collatz t n % 2 = 0 then evenSteps n t + 1 else evenSteps n t

def oddSteps (n : Nat) : Nat → Nat
  | 0 => 0
  | t + 1 =>
    if iterate collatz t n % 2 = 0 then oddSteps n t else oddSteps n t + 1

def standardCarry (n : Nat) : Nat → Nat
  | 0 => 0
  | t + 1 =>
    if iterate collatz t n % 2 = 0 then standardCarry n t
    else 3 * standardCarry n t + 2 ^ evenSteps n t

theorem standard_clock_partition (n t : Nat) : evenSteps n t + oddSteps n t = t := by
  induction t with
  | zero => rfl
  | succ t ih =>
    by_cases he : iterate collatz t n % 2 = 0
    · simp only [evenSteps, oddSteps, he, if_true]
      omega
    · simp only [evenSteps, oddSteps, he, if_false]
      omega

/-- Exact identity for every standard-map prefix, including starts at zero. -/
theorem standard_affine_tally (n t : Nat) :
    2 ^ evenSteps n t * iterate collatz t n =
    3 ^ oddSteps n t * n + standardCarry n t := by
  induction t with
  | zero => simp only [evenSteps, oddSteps, standardCarry, iterate, Nat.pow_zero,
                       Nat.one_mul, Nat.add_zero]
  | succ t ih =>
    by_cases he : iterate collatz t n % 2 = 0
    · have hdiv : 2 * (iterate collatz t n / 2) = iterate collatz t n := by
        have hd := Nat.mod_add_div (iterate collatz t n) 2
        rw [he, Nat.zero_add] at hd
        exact hd
      simp only [evenSteps, oddSteps, standardCarry, iterate, collatz, he, if_true]
      calc
        2 ^ (evenSteps n t + 1) * (iterate collatz t n / 2) =
            2 ^ evenSteps n t * (2 * (iterate collatz t n / 2)) := by
          rw [Nat.pow_succ, Nat.mul_assoc]
        _ = 2 ^ evenSteps n t * iterate collatz t n := by rw [hdiv]
        _ = 3 ^ oddSteps n t * n + standardCarry n t := ih
    · simp only [evenSteps, oddSteps, standardCarry, iterate, collatz, he, if_false]
      calc
        2 ^ evenSteps n t * (3 * iterate collatz t n + 1) =
            3 * (2 ^ evenSteps n t * iterate collatz t n) + 2 ^ evenSteps n t := by
          rw [Nat.mul_add, Nat.mul_one]
          ac_rfl
        _ = 3 * (3 ^ oddSteps n t * n + standardCarry n t) + 2 ^ evenSteps n t := by
          rw [ih]
        _ = 3 ^ (oddSteps n t + 1) * n +
            (3 * standardCarry n t + 2 ^ evenSteps n t) := by
          rw [Nat.pow_succ, Nat.mul_add]
          ac_rfl

/-- The certificate now uses counters and carry computed from the actual orbit. -/
theorem standard_descent_iff_tally_margin (n t : Nat) :
    iterate collatz t n < n ↔
    3 ^ oddSteps n t * n + standardCarry n t < 2 ^ evenSteps n t * n :=
  descent_iff_affine_inequality n (iterate collatz t n)
    (evenSteps n t) (oddSteps n t) (standardCarry n t) (standard_affine_tally n t)

def GlobalTallyMargin : Prop :=
  ∀ n, 1 < n → ∃ t, 0 < t ∧
    3 ^ oddSteps n t * n + standardCarry n t < 2 ^ evenSteps n t * n

/-- A precise global reduction; the right-hand side is still unproved. -/
theorem collatz_reachesOne_iff_globalTallyMargin :
    EventuallyReachesOne collatz ↔ GlobalTallyMargin := by
  rw [collatz_reachesOne_iff_strictDescent]
  constructor
  · intro h n hn
    obtain ⟨t, ht, hd⟩ := h n hn
    exact ⟨t, ht, (standard_descent_iff_tally_margin n t).mp hd⟩
  · intro h n hn
    obtain ⟨t, ht, hm⟩ := h n hn
    exact ⟨t, ht, (standard_descent_iff_tally_margin n t).mpr hm⟩

/-! The attachment instead leaves K and B free and uses total time in 3^t. -/

def RawLocalCertificate (f : Nat → Nat) (n : Nat) : Prop :=
  ∃ t K B, 2 ^ K * iterate f t n = 3 ^ t * n + B ∧
    3 ^ t * n + B < 2 ^ K * n

theorem rawLocalCertificate_implies_localDescent (f : Nat → Nat) (n : Nat)
    (h : RawLocalCertificate f n) :
    ∃ t, 0 < t ∧ iterate f t n < n := by
  obtain ⟨t, K, B, hid, hm⟩ := h
  have hd := (descent_iff_affine_inequality n (iterate f t n) K t B hid).mpr hm
  refine ⟨t, ?_, hd⟩
  cases t with
  | zero => change n < n at hd; omega
  | succ t => exact Nat.zero_lt_succ t

/-- Freely inflated metadata make this certificate true at n=2. -/
theorem rawLocalCertificate_two : RawLocalCertificate collatz 2 := by
  exact ⟨1, 3, 2, by decide⟩

theorem true_tally_two :
    iterate collatz 1 2 = 1 ∧ evenSteps 2 1 = 1 ∧
    oddSteps 2 1 = 0 ∧ standardCarry 2 1 = 0 := by decide

/-- Retaining the true halving count exposes the wrong total-time exponent. -/
theorem wrong_tripling_clock_two :
    ¬ ∃ B : Nat, 2 ^ evenSteps 2 1 * iterate collatz 1 2 = 3 ^ 1 * 2 + B := by
  intro ⟨B, h⟩
  change 2 = 6 + B at h
  omega

/-! A positive-preserving map separates a local witness from a global claim. -/

def localTrap (n : Nat) : Nat := if n = 2 then 1 else n

theorem localTrap_positive : PositivePreserving localTrap := by
  intro n hn
  unfold localTrap
  split
  · exact Nat.zero_lt_succ 0
  · exact hn

theorem localTrap_three_fixed (t : Nat) : iterate localTrap t 3 = 3 := by
  induction t with
  | zero => rfl
  | succ t ih => simp only [iterate, ih, localTrap]; rfl

/-- This refutes the generic inference, not the still-open Collatz-specific iff. -/
theorem local_certificate_does_not_imply_global :
    RawLocalCertificate localTrap 2 ∧ ¬ EventuallyStrictlyDescends localTrap := by
  constructor
  · exact ⟨1, 3, 2, by decide⟩
  · intro h
    obtain ⟨t, _ht, hd⟩ := h 3 (by decide)
    rw [localTrap_three_fixed] at hd
    omega

/-! Exact paired controls: same residue and carry, different halving totals. -/

theorem twentyThree_tally_control :
    iterate collatz 10 23 = 5 ∧ evenSteps 23 10 = 7 ∧
    oddSteps 23 10 = 3 ∧ standardCarry 23 10 = 19 := by decide

theorem ninetyFive_tally_control :
    iterate collatz 6 95 = 323 ∧ evenSteps 95 6 = 3 ∧
    oddSteps 95 6 = 3 ∧ standardCarry 95 6 = 19 := by decide

theorem residue_and_carry_do_not_force_descent :
    23 % 9 = 5 ∧ 95 % 9 = 5 ∧
    standardCarry 23 10 = standardCarry 95 6 ∧
    iterate collatz 10 23 < 23 ∧ 95 < iterate collatz 6 95 := by decide

/-- Unlike a positive-start formulation, the attachment's all-Nat conclusion is false. -/
theorem not_all_naturals_reachOne : ¬ ∀ n, ∃ t, iterate collatz t n = 1 := by
  intro h
  exact collatz_zero_never_reachesOne (h 0)

end CollatzBlueprintAudit
