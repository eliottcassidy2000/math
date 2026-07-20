/-
Copyright (c) 2026 The multi-agent math project contributors.
Released under Apache 2.0.
Authors: death-star-2026-07-20-S61h

The combinatorial reduction NC2 => GMC(2), formalized (kernel-pure, no sorry/native_decide).

Two real Gaussians = one complex Gaussian Z; the Gaussian expectation E on C[Z,W]
(Z = var 0, W = var 1) has E[Z^a W^b] = a! * [a=b] (Wick).  The CHARGE of a monomial
Z^a W^b is a - b; E kills every nonzero charge.  The 2-D Nullcone Conjecture (NC2) says a
nullcone member is charge one-sided.  This file proves the step everyone uses:

  ONE-SIDED CHARGE SUPPORT  =>  MATHIEU-ZHAO  (i.e. GMC(2) for that P),

by pure charge arithmetic: if every monomial of P has charge >= 1, then for any Q the
product Q * P^m has all charges >= m - deg_W Q, which is > 0 for m large, so E(Q P^m) = 0.

This is the analysis-free spine (klein-S351 "three lines"; mac-mini/kp THM-1540 corollary).
The analytic input -- that NC2 itself holds (klein Gamma Bridge + Duistermaat-van der Kallen
n=1 for the toral layer + boxeph's Radial Lemma) -- is NOT formalized here; it enters as the
hypothesis `hP` (all charges >= 1), exactly the one-sided conclusion of NC2.
-/
import Mathlib

open MvPolynomial Finset

namespace GMC2

/-- Charge of a monomial (exponent vector) in two variables Z = var 0, W = var 1. -/
def charge (s : Fin 2 →₀ ℕ) : ℤ := (s 0 : ℤ) - (s 1 : ℤ)

/-- Wick weight of a monomial: `a!` on the diagonal `Z^a W^a` (charge 0), else 0. -/
noncomputable def wt (s : Fin 2 →₀ ℕ) : ℂ := if s 0 = s 1 then (Nat.factorial (s 0) : ℂ) else 0

/-- The Gaussian expectation `E[Z^a W^b] = a! [a=b]`, extended linearly. -/
noncomputable def E (P : MvPolynomial (Fin 2) ℂ) : ℂ := ∑ s ∈ P.support, P.coeff s * wt s

@[simp] lemma charge_add (s t : Fin 2 →₀ ℕ) : charge (s + t) = charge s + charge t := by
  simp only [charge, Finsupp.add_apply]; push_cast; ring

/-- Charge 0 is exactly the diagonal; off-diagonal the Wick weight is 0. -/
lemma wt_of_charge_ne {s : Fin 2 →₀ ℕ} (h : charge s ≠ 0) : wt s = 0 := by
  unfold wt
  split_ifs with hh
  · exact absurd (by simp [charge, hh]) h
  · rfl

/-- If every monomial of `P` has nonzero charge, then `E P = 0`. -/
lemma E_eq_zero_of_charges_ne {P : MvPolynomial (Fin 2) ℂ}
    (h : ∀ s ∈ P.support, charge s ≠ 0) : E P = 0 := by
  unfold E
  refine Finset.sum_eq_zero (fun s hs => ?_)
  rw [wt_of_charge_ne (h s hs), mul_zero]

/-- Charges add and are bounded below on powers: if every monomial of `P` has charge `≥ 1`,
then every monomial of `P^m` has charge `≥ m`. -/
lemma le_charge_of_mem_support_pow {P : MvPolynomial (Fin 2) ℂ}
    (hP : ∀ s ∈ P.support, (1 : ℤ) ≤ charge s) :
    ∀ (m : ℕ), ∀ s ∈ (P ^ m).support, (m : ℤ) ≤ charge s := by
  intro m
  induction m with
  | zero =>
    intro s hs
    simp only [pow_zero] at hs
    have : s = 0 := by
      have := MvPolynomial.support_one (R := ℂ) (σ := Fin 2)
      rw [this] at hs; exact Finset.mem_singleton.mp hs
    simp [this, charge]
  | succ n ih =>
    intro s hs
    rw [pow_succ] at hs
    have hmem := MvPolynomial.support_mul _ _ hs
    rw [Finset.mem_add] at hmem
    obtain ⟨a, ha, b, hb, rfl⟩ := hmem
    rw [charge_add]
    have h1 := ih a ha
    have h2 := hP b hb
    push_cast
    linarith

/-- **The reduction NC2 ⇒ GMC(2).** If every monomial of `P` has charge `≥ 1` (the one-sided
conclusion of the 2-D nullcone conjecture), then `ker E` is Mathieu–Zhao at `P`: for every `Q`
there is a threshold `N` beyond which `E (Q * P^m) = 0`.  Explicit threshold: `N = deg_W Q + 1`. -/
theorem mathieuZhao_of_charge_pos (P Q : MvPolynomial (Fin 2) ℂ)
    (hP : ∀ s ∈ P.support, (1 : ℤ) ≤ charge s) :
    ∃ N : ℕ, ∀ m ≥ N, E (Q * P ^ m) = 0 := by
  refine ⟨(Q.support.sup fun s => s 1) + 1, fun m hm => ?_⟩
  apply E_eq_zero_of_charges_ne
  intro s hs
  have hmem := MvPolynomial.support_mul _ _ hs
  rw [Finset.mem_add] at hmem
  obtain ⟨c, hc, b, hb, rfl⟩ := hmem
  -- charge b ≥ m (from the power), charge c ≥ -(c 1) ≥ -(sup)
  have hb2 := le_charge_of_mem_support_pow hP m b hb
  have hc1 : (c 1 : ℤ) ≤ ((Q.support.sup fun s => s 1 : ℕ) : ℤ) := by
    exact_mod_cast Finset.le_sup (f := fun s => s 1) hc
  have hcc : -(c 1 : ℤ) ≤ charge c := by
    have : (0 : ℤ) ≤ (c 0 : ℤ) := by positivity
    simp only [charge]; linarith
  have hmN : ((Q.support.sup fun s => s 1 : ℕ) : ℤ) + 1 ≤ (m : ℤ) := by exact_mod_cast hm
  -- charge (c + b) ≥ 1 > 0
  rw [charge_add]
  have : (1 : ℤ) ≤ charge c + charge b := by linarith
  linarith

end GMC2

-- Axiom audit: should be only propext / Classical.choice / Quot.sound (Mathlib standard).
#print axioms GMC2.mathieuZhao_of_charge_pos
