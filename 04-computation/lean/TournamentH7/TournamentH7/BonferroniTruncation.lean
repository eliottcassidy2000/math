/-
  TournamentH7.BonferroniTruncation  (mac-mini-2026-07-01-S101)

  THM-599's combinatorial core, Lean-checked: the partial alternating binomial
  sum identity and the odd-depth Bonferroni lower bound.

    Σ_{d ≤ D} (−1)^d C(c, d)  =  (−1)^D C(c−1, D)      (c ≥ 1),

  hence for ODD D the partial sum is ≤ 0 = 1_{c = 0} on c ≥ 1, and equals 1 at
  c = 0: the truncated inclusion–exclusion at odd depth is a pointwise lower
  bound for the uncovered indicator — the pointwise engine of the quintic
  closure (THM-599, kappa_13 = 2052/7^5).  Sorry-free.
-/
import Mathlib.Tactic

namespace TournamentH7.BonferroniTruncation

open Finset

/-- Partial alternating sum of binomial coefficients:
`Σ_{d=0}^{D} (−1)^d C(c, d) = (−1)^D C(c−1, D)` for `1 ≤ c`. -/
theorem partial_alternating_choose (c D : ℕ) (hc : 1 ≤ c) :
    ∑ d ∈ range (D + 1), (-1 : ℤ) ^ d * (c.choose d) = (-1) ^ D * ((c - 1).choose D) := by
  induction D with
  | zero => simp
  | succ D ih =>
    rw [sum_range_succ, ih]
    have hpascal : (c.choose (D + 1) : ℤ)
        = ((c - 1).choose D : ℤ) + ((c - 1).choose (D + 1) : ℤ) := by
      have hc1 : c - 1 + 1 = c := Nat.sub_add_cancel hc
      have h : (c - 1 + 1).choose (D + 1) = (c - 1).choose D + (c - 1).choose (D + 1) :=
        Nat.choose_succ_succ (c - 1) D
      rw [hc1] at h
      exact_mod_cast h
    rw [hpascal, pow_succ]
    ring

/-- **The odd-depth Bonferroni lower bound (pointwise).**  For odd `D` and any
`c : ℕ`: `Σ_{d ≤ D} (−1)^d C(c,d) ≤ (if c = 0 then 1 else 0)`. -/
theorem odd_truncation_le_uncovered (c D : ℕ) (hD : D % 2 = 1) :
    ∑ d ∈ range (D + 1), (-1 : ℤ) ^ d * (c.choose d) ≤ (if c = 0 then 1 else 0) := by
  rcases Nat.eq_zero_or_pos c with hc | hc
  · subst hc
    simp only [if_pos rfl]
    have hsum : ∑ d ∈ range (D + 1), (-1 : ℤ) ^ d * ((0 : ℕ).choose d)
        = ∑ d ∈ range (D + 1), (if d = 0 then (1 : ℤ) else 0) := by
      apply sum_congr rfl
      intro d _
      match d with
      | 0 => simp
      | (e + 1) =>
        rw [Nat.choose_eq_zero_of_lt (Nat.succ_pos e)]
        simp
    rw [hsum, sum_ite_eq' (range (D + 1)) 0 (fun _ => (1 : ℤ))]
    simp
  · rw [if_neg (by omega : c ≠ 0)]
    rw [partial_alternating_choose c D hc]
    have hodd : Odd D := Nat.odd_iff.mpr hD
    rw [hodd.neg_one_pow]
    have hpos : (0 : ℤ) ≤ ((c - 1).choose D : ℤ) := by positivity
    linarith

end TournamentH7.BonferroniTruncation
