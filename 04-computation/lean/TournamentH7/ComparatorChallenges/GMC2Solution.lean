/-
Copyright (c) 2026 The multi-agent math project contributors.
Released under Apache 2.0.
Authors: klein-S428

# Comparator solution: GMC(2)

Discharges `GMC2Challenge.gmc2` from the repo's `GMC2.gmc2`.

The only mathematical content here is the **bridge**: every functional satisfying
`GMC2Challenge.WickFunctional` is equal to `GMC2.E`. That is what upgrades the
existing certificate from "kernel-pure relative to our definition of `E`" to
"kernel-pure relative to the Gaussian moment law", which is the thing a reader
can actually check.

The bridge is three steps:
  1. `E 0 = 0`, from additivity alone;
  2. `E (monomial s c) = c * GMC2.wt s`, by writing `monomial s c` as
     `C c * (X 0 ^ s 0 * X 1 ^ s 1)` and applying the moment law -- this is where
     the Wick weight `wt` is *derived* rather than assumed;
  3. `E = GMC2.E`, since a polynomial is the finite sum of its monomials.
-/
import TournamentH7.GMC2Main
import ComparatorChallenges.GMC2Challenge

noncomputable section

open MvPolynomial Finset

namespace GMC2Challenge

variable {E : MvPolynomial (Fin 2) ℂ → ℂ}

/-- Additivity alone forces `E 0 = 0`. -/
lemma wick_zero (hE : WickFunctional E) : E 0 = 0 := by
  have h := hE.add 0 0
  rw [add_zero] at h
  exact self_eq_add_left.mp h

/-- `E` commutes with finite sums. -/
lemma wick_finsetSum (hE : WickFunctional E) {α : Type*} (s : Finset α)
    (f : α → MvPolynomial (Fin 2) ℂ) :
    E (∑ i ∈ s, f i) = ∑ i ∈ s, E (f i) := by
  classical
  induction s using Finset.induction with
  | empty => simpa using wick_zero hE
  | insert a s ha ih =>
      rw [Finset.sum_insert ha, Finset.sum_insert ha, hE.add, ih]

/-- **The Wick weight is derived, not assumed.** The moment law on `X 0 ^ a * X 1 ^ b`
already determines `E` on every monomial, and the value is exactly `GMC2.wt`. -/
lemma wick_monomial (hE : WickFunctional E) (s : Fin 2 →₀ ℕ) (c : ℂ) :
    E (monomial s c) = c * GMC2.wt s := by
  have hs : (monomial s c : MvPolynomial (Fin 2) ℂ)
      = C c * (X 0 ^ s 0 * X 1 ^ s 1) := by
    rw [MvPolynomial.monomial_eq]
    congr 1
    rw [Finsupp.prod_fintype _ _ (fun i => pow_zero _)]
    simp [Fin.prod_univ_two]
  rw [hs, hE.smul, hE.moment, GMC2.wt]

/-- **The bridge.** Any Wick functional is the repo's Gaussian expectation. -/
theorem wick_eq_E (hE : WickFunctional E) (P : MvPolynomial (Fin 2) ℂ) :
    E P = GMC2.E P := by
  conv_lhs => rw [P.as_sum]
  rw [wick_finsetSum hE]
  simp only [wick_monomial hE]
  rfl

/-- `GMC2.E` itself satisfies the challenge's hypothesis, so the hypothesis is
not vacuous -- there is at least one Wick functional. -/
theorem wickFunctional_E : WickFunctional GMC2.E where
  add := GMC2.E_add
  smul := by
    intro c P
    conv_lhs => rw [P.as_sum, Finset.mul_sum]
    conv_rhs => rw [P.as_sum]
    rw [GMC2.E_finset_sum, GMC2.E_finset_sum, Finset.mul_sum]
    refine Finset.sum_congr rfl (fun s _ => ?_)
    rw [C_mul_monomial, GMC2.E_monomial, GMC2.E_monomial]
    ring
  moment := by
    intro a b
    have hs : (X 0 ^ a * X 1 ^ b : MvPolynomial (Fin 2) ℂ)
        = monomial (Finsupp.single (0 : Fin 2) a + Finsupp.single (1 : Fin 2) b) 1 := by
      rw [X_pow_eq_monomial, X_pow_eq_monomial, monomial_mul, one_mul]
    rw [hs, GMC2.E_monomial, GMC2.wt]
    simp

/-- **GMC(2), against the independently stated Wick hypothesis.** -/
theorem gmc2' (E : MvPolynomial (Fin 2) ℂ → ℂ) (hE : WickFunctional E)
    (P Q : MvPolynomial (Fin 2) ℂ)
    (hnull : ∀ m : ℕ, 1 ≤ m → E (P ^ m) = 0) :
    ∃ N : ℕ, ∀ m ≥ N, E (Q * P ^ m) = 0 := by
  have hnull' : ∀ m : ℕ, 1 ≤ m → GMC2.E (P ^ m) = 0 := by
    intro m hm; rw [← wick_eq_E hE]; exact hnull m hm
  obtain ⟨N, hN⟩ := GMC2.gmc2 P Q hnull'
  exact ⟨N, fun m hm => by rw [wick_eq_E hE]; exact hN m hm⟩

end GMC2Challenge

end
