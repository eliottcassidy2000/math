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

/-- The negative-charge mirror of `le_charge_of_mem_support_pow`: if every monomial of
`P` has charge at most `-1`, every monomial of `P^m` has charge at most `-m`. -/
lemma charge_le_neg_of_mem_support_pow {P : MvPolynomial (Fin 2) ℂ}
    (hP : ∀ s ∈ P.support, charge s ≤ (-1 : ℤ)) :
    ∀ (m : ℕ), ∀ s ∈ (P ^ m).support, charge s ≤ -(m : ℤ) := by
  intro m
  induction m with
  | zero =>
    intro s hs
    simp only [pow_zero] at hs
    have : s = 0 := by
      have := MvPolynomial.support_one (R := ℂ) (σ := Fin 2)
      rw [this] at hs
      exact Finset.mem_singleton.mp hs
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

/-- Strictly positive charge support is in the moment nullcone: every positive power has
zero Gaussian expectation. -/
theorem moments_zero_of_charge_pos (P : MvPolynomial (Fin 2) ℂ)
    (hP : ∀ s ∈ P.support, (1 : ℤ) ≤ charge s) :
    ∀ m : ℕ, 1 ≤ m → E (P ^ m) = 0 := by
  intro m hm
  apply E_eq_zero_of_charges_ne
  intro s hs hzero
  have hcharge := le_charge_of_mem_support_pow hP m s hs
  rw [hzero] at hcharge
  exact (not_le_of_gt (by exact_mod_cast hm)) hcharge

/-- Strictly negative charge support is in the moment nullcone: every positive power has
zero Gaussian expectation. -/
theorem moments_zero_of_charge_neg (P : MvPolynomial (Fin 2) ℂ)
    (hP : ∀ s ∈ P.support, charge s ≤ (-1 : ℤ)) :
    ∀ m : ℕ, 1 ≤ m → E (P ^ m) = 0 := by
  intro m hm
  apply E_eq_zero_of_charges_ne
  intro s hs hzero
  have hcharge := charge_le_neg_of_mem_support_pow hP m s hs
  rw [hzero] at hcharge
  have hmz : (0 : ℤ) < (m : ℤ) := by exact_mod_cast hm
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

/-- Negative-charge branch of the NC2-to-GMC(2) reduction. The explicit threshold is one
more than the largest `Z` exponent in `Q`. -/
theorem mathieuZhao_of_charge_neg (P Q : MvPolynomial (Fin 2) ℂ)
    (hP : ∀ s ∈ P.support, charge s ≤ (-1 : ℤ)) :
    ∃ N : ℕ, ∀ m ≥ N, E (Q * P ^ m) = 0 := by
  refine ⟨(Q.support.sup fun s => s 0) + 1, fun m hm => ?_⟩
  apply E_eq_zero_of_charges_ne
  intro s hs
  have hmem := MvPolynomial.support_mul _ _ hs
  rw [Finset.mem_add] at hmem
  obtain ⟨c, hc, b, hb, rfl⟩ := hmem
  have hb2 := charge_le_neg_of_mem_support_pow hP m b hb
  have hc0 : (c 0 : ℤ) ≤ ((Q.support.sup fun s => s 0 : ℕ) : ℤ) := by
    exact_mod_cast Finset.le_sup (f := fun s => s 0) hc
  have hcc : charge c ≤ (c 0 : ℤ) := by
    have : (0 : ℤ) ≤ (c 1 : ℤ) := by positivity
    simp only [charge]
    linarith
  have hmN : ((Q.support.sup fun s => s 0 : ℕ) : ℤ) + 1 ≤ (m : ℤ) := by
    exact_mod_cast hm
  rw [charge_add]
  have : charge c + charge b ≤ (-1 : ℤ) := by linarith
  linarith

/-- The charge-one-sided conclusion of NC2, with both orientations retained explicitly. -/
def ChargeOneSided (P : MvPolynomial (Fin 2) ℂ) : Prop :=
  (∀ s ∈ P.support, (1 : ℤ) ≤ charge s) ∨
  (∀ s ∈ P.support, charge s ≤ (-1 : ℤ))

/-- **Entry point of the THM-2022 contrapositive.** `P` is *not* charge-one-sided exactly when its
support carries a monomial of charge `≤ 0` and one of charge `≥ 0` — i.e. `0 ∈ conv(charges)`, the
hypothesis from which the lowest balanced face is built. -/
lemma not_chargeOneSided_iff (P : MvPolynomial (Fin 2) ℂ) :
    ¬ ChargeOneSided P ↔
      (∃ s ∈ P.support, charge s ≤ 0) ∧ (∃ t ∈ P.support, 0 ≤ charge t) := by
  unfold ChargeOneSided
  push_neg
  constructor
  · rintro ⟨⟨s, hs, hs1⟩, ⟨t, ht, ht1⟩⟩
    exact ⟨⟨s, hs, by omega⟩, ⟨t, ht, by omega⟩⟩
  · rintro ⟨⟨s, hs, hs0⟩, ⟨t, ht, ht0⟩⟩
    exact ⟨⟨s, hs, by omega⟩, ⟨t, ht, by omega⟩⟩

/-- A pointwise formulation of the two-dimensional nullcone statement for one polynomial. -/
def NC2At (P : MvPolynomial (Fin 2) ℂ) : Prop :=
  (∀ m : ℕ, 1 ≤ m → E (P ^ m) = 0) → ChargeOneSided P

/-- Both strict one-sided loci are moment-null. This is the easy converse in NC2. -/
theorem moments_zero_of_charge_oneSided (P : MvPolynomial (Fin 2) ℂ)
    (hP : ChargeOneSided P) : ∀ m : ℕ, 1 ≤ m → E (P ^ m) = 0 := by
  rcases hP with hpos | hneg
  · exact moments_zero_of_charge_pos P hpos
  · exact moments_zero_of_charge_neg P hneg

/-- **NC2 implies GMC(2), formal interface.** Once null moments classify `P` as strictly
charge-one-sided, the elementary charge bounds prove eventual vanishing of `E (Q * P^m)`
for every multiplier `Q`. No tournament ordering, tie-breaker, or Paley identification enters. -/
theorem mathieuZhao_of_nc2At (P Q : MvPolynomial (Fin 2) ℂ)
    (hNC2 : NC2At P) (hnull : ∀ m : ℕ, 1 ≤ m → E (P ^ m) = 0) :
    ∃ N : ℕ, ∀ m ≥ N, E (Q * P ^ m) = 0 := by
  rcases hNC2 hnull with hpos | hneg
  · exact mathieuZhao_of_charge_pos P Q hpos
  · exact mathieuZhao_of_charge_neg P Q hneg

/-- The full 2-D nullcone conjecture NC2 (both directions), as a single statement over all `P`.
The easy direction (`←`, one-sided ⇒ null) is `moments_zero_of_charge_oneSided`; the hard
direction (`→`, null ⇒ one-sided) is the content of THM-2022. -/
def NC2 : Prop := ∀ P : MvPolynomial (Fin 2) ℂ, NC2At P

/-- **GMC(2) is a corollary of NC2, with no analytic gap.** Given NC2 (the hard classification),
for every `P` whose positive moments all vanish and every multiplier `Q`, the products
`E (Q * P^m)` vanish for all large `m`.  This is the clean logical reduction: the entire GMC(2)
problem rests on the single theorem `NC2`. -/
theorem gmc2_of_nc2 (hNC2 : NC2) (P Q : MvPolynomial (Fin 2) ℂ)
    (hnull : ∀ m : ℕ, 1 ≤ m → E (P ^ m) = 0) :
    ∃ N : ℕ, ∀ m ≥ N, E (Q * P ^ m) = 0 :=
  mathieuZhao_of_nc2At P Q (hNC2 P) hnull

/-! ### THM-2022 §4 — Kummer/Lucas isolation of the dilated face layer

At a good prime `p > m0`, dividing the moment of order `p·m0` by `(p·A0)!` and reducing mod `p`
kills every channel that is not `p`-dilated (Kummer's carries) and, on the surviving `p·s`
channels, replaces the multinomial coefficient by its undilated value (Lucas).  The multinomial
Lucas below is the exact form needed; Mathlib had the binomial case but not this. -/

/-- **Multinomial Lucas.** Dilating every part of a multiplicity vector by a prime `p` fixes the
multinomial coefficient modulo `p`: `multinomial S (p • k) ≡ multinomial S k [MOD p]`.  Assembled
from the binomial Lucas `Choose.choose_mul_mul_modEq_choose_nat` along the `multinomial_insert`
recursion (THM-2022 (14)). -/
theorem multinomial_dilate_modEq {α : Type*} [DecidableEq α] (p : ℕ) [Fact p.Prime]
    (S : Finset α) (f : α → ℕ) :
    Nat.multinomial S (fun i => p * f i) ≡ Nat.multinomial S f [MOD p] := by
  classical
  induction S using Finset.induction with
  | empty => rw [Nat.multinomial_empty, Nat.multinomial_empty]
  | insert a s ha ih =>
      rw [Nat.multinomial_insert ha (fun i => p * f i), Nat.multinomial_insert ha f]
      have hsum : ∑ i ∈ s, p * f i = p * ∑ i ∈ s, f i := by rw [Finset.mul_sum]
      simp only [hsum, ← Nat.mul_add]
      exact Nat.ModEq.mul Choose.choose_mul_mul_modEq_choose_nat ih

/-- If `p` is prime, `p ∣ n`, and `p ∤ k ≤ n`, then `p ∣ C(n,k)`.  (Absorption identity
`k · C(n,k) = n · C(n-1,k-1)` plus primality; the elementary half of Kummer.) -/
lemma dvd_choose_of_dvd {p : ℕ} (hp : p.Prime) {n k : ℕ} (hn : p ∣ n) (hk : ¬ p ∣ k)
    (hkn : k ≤ n) : p ∣ n.choose k := by
  have hk1 : 1 ≤ k := Nat.one_le_iff_ne_zero.mpr (fun h => hk (h ▸ dvd_zero p))
  have hn1 : 1 ≤ n := le_trans hk1 hkn
  have habs : n * (n - 1).choose (k - 1) = n.choose k * k := by
    have h := Nat.succ_mul_choose_eq (n - 1) (k - 1)
    simp only [Nat.succ_eq_add_one, Nat.sub_add_cancel hn1, Nat.sub_add_cancel hk1] at h
    linarith [h]
  have hdvd : p ∣ n.choose k * k := habs ▸ (dvd_mul_of_dvd_left hn _)
  exact (hp.dvd_mul.mp hdvd).resolve_right hk

/-- **THM-2022 §4 no-carry channel survival (11)–(12).** If `p` is prime, `p ∣ ∑ r`, and some part
`r i` is *not* divisible by `p`, then `p ∣ multinomial S r`.  Thus at mass `p·m0` the only channels
surviving mod `p` are the `p`-dilated ones — a non-dilated allocation crosses a carry wall. -/
lemma multinomial_dvd_of_exists_not_dvd {α : Type*} [DecidableEq α] {p : ℕ} (hp : p.Prime)
    (S : Finset α) (r : α → ℕ) (hsum : p ∣ ∑ i ∈ S, r i)
    (hex : ∃ i ∈ S, ¬ p ∣ r i) : p ∣ Nat.multinomial S r := by
  classical
  revert hsum hex
  induction S using Finset.induction with
  | empty => intro _ hex; obtain ⟨i, hi, _⟩ := hex; exact absurd hi (Finset.notMem_empty i)
  | insert a s ha ih =>
      intro hsum hex
      rw [Nat.multinomial_insert ha]
      rw [Finset.sum_insert ha] at hsum
      by_cases hra : p ∣ r a
      · have hsum_s : p ∣ ∑ i ∈ s, r i := (Nat.dvd_add_right hra).mp hsum
        have hex_s : ∃ i ∈ s, ¬ p ∣ r i := by
          obtain ⟨i, hi, hpi⟩ := hex
          rcases Finset.mem_insert.mp hi with rfl | his
          · exact absurd hra hpi
          · exact ⟨i, his, hpi⟩
        exact dvd_mul_of_dvd_right (ih hsum_s hex_s) _
      · exact dvd_mul_of_dvd_left (dvd_choose_of_dvd hp hsum hra (Nat.le_add_right _ _)) _

/-- **THM-2022 §4 off-face vanishing (case 2).** An *off-face* dilated channel (radial height
`A' > A0`) has its Wick factorial ratio killed by `p`: `p · (p·A0)! ∣ (p·A')!`.  After normalizing
the moment by `(p·A0)!`, its quotient `(p·A')!/(p·A0)!` is divisible by `p` — it contains the factor
`p·(A0+1)`.  So the surviving residue layer is exactly the *on-face* dilated channels. -/
lemma factorial_dilate_dvd {p : ℕ} (hp : 1 ≤ p) {A0 A' : ℕ} (h : A0 < A') :
    p * (p * A0).factorial ∣ (p * A').factorial := by
  have h1 : p * (A0 + 1) ≤ p * A' := mul_le_mul_left' (by omega) p
  have hdvd1 : (p * (A0 + 1)).factorial ∣ (p * A').factorial := Nat.factorial_dvd_factorial h1
  have heq : p * (A0 + 1) = p * A0 + p := by ring
  rw [heq] at hdvd1
  have hpdvd : p ∣ (p * A0 + 1).ascFactorial p :=
    dvd_trans (Nat.dvd_factorial hp le_rfl) (Nat.factorial_dvd_ascFactorial _ _)
  calc p * (p * A0).factorial
      ∣ (p * A0).factorial * (p * A0 + 1).ascFactorial p := by
        rw [mul_comm p]; exact mul_dvd_mul_left _ hpdvd
    _ = (p * A0 + p).factorial := Nat.factorial_mul_ascFactorial _ _
    _ ∣ (p * A').factorial := hdvd1

/-! ### THM-2022 §5 — Frobenius non-cancellation of the lowest balanced face

In the residue field at a good prime `p`, the natural-number Wick/multinomial weights are
Frobenius-fixed, so the normalized moment layer collapses to an exact `p`-th power `Q̄^p`.
This is the algebraic engine of THM-2022: an entire *tied* balanced face survives as one
Frobenius power and therefore cannot cancel.  Proved here in full generality over any
commutative ring of characteristic `p`. -/

/-- **Frobenius collapse of a natCast-weighted channel sum.** Over a commutative ring of
characteristic `p`, `(∑ w_s g_s)^p = ∑ w_s g_s^p` — the weights `w_s : ℕ` are fixed by the
Frobenius endomorphism.  With `g_s` the channel monomial value and `w_s` the multinomial
coefficient, the right side is the `p`-dilated face layer and the left is `Q̄^p` (THM-2022 (15)). -/
theorem sum_natCast_mul_pow_char {R : Type*} [CommRing R] (p : ℕ) [ExpChar R p]
    {ι : Type*} (S : Finset ι) (w : ι → ℕ) (g : ι → R) :
    (∑ s ∈ S, (w s : R) * g s) ^ p = ∑ s ∈ S, (w s : R) * (g s) ^ p := by
  rw [sum_pow_char]
  refine Finset.sum_congr rfl (fun s _ => ?_)
  rw [mul_pow, ← frobenius_def, map_natCast]

/-- **THM-2022 (15): the tied balanced face survives Frobenius as `Q̄^p`.** In the residue field of
characteristic `p`, the `p`-dilated balanced-face channel sum `∑_t multinomial(S, p·s_t) · (h_t)^p`
equals the `p`-th power of the face constant term `Q̄ = ∑_t multinomial(S, s_t) · h_t`.  Combines the
multinomial Lucas `multinomial_dilate_modEq` (§4, weights unchanged mod `p`) with the Frobenius
collapse `sum_natCast_mul_pow_char` (§5).  This is the exact non-cancellation identity: `Q̄ ≠ 0`
forces the dilated moment layer `≠ 0`. -/
theorem face_sum_frobenius {F : Type*} [CommRing F] (p : ℕ) [Fact p.Prime] [CharP F p]
    {α : Type*} [DecidableEq α] (S : Finset α) {ι : Type*} (T : Finset ι)
    (s : ι → α → ℕ) (h : ι → F) :
    ∑ t ∈ T, (Nat.multinomial S (fun i => p * s t i) : F) * (h t) ^ p
      = (∑ t ∈ T, (Nat.multinomial S (s t) : F) * h t) ^ p := by
  haveI : ExpChar F p := ExpChar.prime Fact.out
  rw [sum_natCast_mul_pow_char]
  refine Finset.sum_congr rfl (fun t _ => ?_)
  congr 1
  exact (CharP.cast_eq_iff_mod_eq F p).mpr (multinomial_dilate_modEq p S (s t))

/-- **Non-cancellation of the balanced face.** Over a *field* of characteristic `p`, a nonzero
face constant term `Q̄` forces the dilated balanced-face channel sum to be nonzero — the exact
statement THM-2022 needs at the good prime: a whole tied face survives as `Q̄^p ≠ 0` and cannot
cancel.  (`Q̄^p ≠ 0` is `pow_ne_zero` in a field/domain.) -/
theorem face_sum_ne_zero {F : Type*} [Field F] (p : ℕ) [Fact p.Prime] [CharP F p]
    {α : Type*} [DecidableEq α] (S : Finset α) {ι : Type*} (T : Finset ι)
    (s : ι → α → ℕ) (h : ι → F)
    (hQ : (∑ t ∈ T, (Nat.multinomial S (s t) : F) * h t) ≠ 0) :
    ∑ t ∈ T, (Nat.multinomial S (fun i => p * s t i) : F) * (h t) ^ p ≠ 0 := by
  rw [face_sum_frobenius]
  exact pow_ne_zero p hQ

/-! ### THM-2022 §1 — the exact Wick channel expansion of the Gaussian moment

`E` is a linear functional with `E (monomial v c) = c * wt v`.  Expanding `P^m` by the
multinomial theorem and evaluating term by term gives the moment as a sum over multiplicity
vectors (channels): the Wick factorial `wt (radial)` vanishes off the charge-0 diagonal, so only
balanced channels survive.  This is THM-2022 (1), the bridge from the abstract expectation to the
channel arithmetic on which the face geometry and the Frobenius layer act. -/

/-- `E P` may be summed over any finite superset of `P.support` (extra coefficients are `0`). -/
lemma E_eq_sum_of_subset {P : MvPolynomial (Fin 2) ℂ} {D : Finset (Fin 2 →₀ ℕ)}
    (h : P.support ⊆ D) : E P = ∑ s ∈ D, P.coeff s * wt s :=
  Finset.sum_subset h (fun s _ hs => by rw [MvPolynomial.notMem_support_iff.mp hs, zero_mul])

/-- The Gaussian expectation is additive. -/
lemma E_add (P Q : MvPolynomial (Fin 2) ℂ) : E (P + Q) = E P + E Q := by
  rw [E_eq_sum_of_subset (P := P) (D := P.support ∪ Q.support) Finset.subset_union_left,
      E_eq_sum_of_subset (P := Q) (D := P.support ∪ Q.support) Finset.subset_union_right,
      E_eq_sum_of_subset (P := P + Q) (D := P.support ∪ Q.support) MvPolynomial.support_add,
      ← Finset.sum_add_distrib]
  exact Finset.sum_congr rfl (fun s _ => by rw [MvPolynomial.coeff_add]; ring)

/-- `E` on a single monomial is its coefficient times the Wick weight. -/
lemma E_monomial (v : Fin 2 →₀ ℕ) (c : ℂ) : E (monomial v c) = c * wt v := by
  rw [E_eq_sum_of_subset (P := monomial v c) (D := {v}) MvPolynomial.support_monomial_subset,
      Finset.sum_singleton]
  simp [MvPolynomial.coeff_monomial]

/-- `E` commutes with finite sums. -/
lemma E_sum {ι : Type*} (s : Finset ι) (f : ι → MvPolynomial (Fin 2) ℂ) :
    E (∑ i ∈ s, f i) = ∑ i ∈ s, E (f i) := by
  classical
  induction s using Finset.induction with
  | empty => simp [show E 0 = 0 by simp [E]]
  | insert a s ha ih => rw [Finset.sum_insert ha, E_add, ih, Finset.sum_insert ha]

/-- A finite product of monomials is the monomial of the summed exponents and multiplied
coefficients. -/
lemma prod_monomial {ι : Type*} (s : Finset ι) (e : ι → (Fin 2 →₀ ℕ)) (c : ι → ℂ) :
    ∏ i ∈ s, monomial (e i) (c i) = monomial (∑ i ∈ s, e i) (∏ i ∈ s, c i) := by
  classical
  induction s using Finset.induction with
  | empty => simp
  | insert a s ha ih =>
      rw [Finset.prod_insert ha, ih, monomial_mul, Finset.sum_insert ha, Finset.prod_insert ha]

/-- **THM-2022 (1): the exact Wick channel expansion of the Gaussian moment.** The `m`-th moment
`E (P^m)` equals the sum, over multiplicity vectors `k` of mass `m` supported on `P.support`, of
the multinomial coefficient `multinomial P.support k` times the coefficient monomial
`∏ v, (P.coeff v)^{k v}` times the Wick factorial `wt` of the radial exponent `∑ v, k v • v`.
Since `wt` is `(radial₀)!` on charge-0 (balanced) channels and `0` off them, this is exactly the
balanced-channel sum `M_m(c)` of THM-2022. -/
theorem wick_expansion (P : MvPolynomial (Fin 2) ℂ) (m : ℕ) :
    E (P ^ m) = ∑ k ∈ P.support.piAntidiag m,
      (Nat.multinomial P.support k : ℂ) * (∏ v ∈ P.support, P.coeff v ^ (k v))
        * wt (∑ v ∈ P.support, k v • v) := by
  have key : P ^ m = ∑ k ∈ P.support.piAntidiag m,
      (Nat.multinomial P.support k : MvPolynomial (Fin 2) ℂ)
        * ∏ v ∈ P.support, (monomial v (P.coeff v)) ^ (k v) := by
    nth_rewrite 1 [MvPolynomial.as_sum P]
    exact Finset.sum_pow_eq_sum_piAntidiag _ _ _
  have key2 : P ^ m = ∑ k ∈ P.support.piAntidiag m,
      monomial (∑ v ∈ P.support, k v • v)
        ((Nat.multinomial P.support k : ℂ) * ∏ v ∈ P.support, P.coeff v ^ (k v)) := by
    rw [key]
    refine Finset.sum_congr rfl (fun k _ => ?_)
    have hp : (∏ v ∈ P.support, (monomial v (P.coeff v)) ^ (k v))
        = monomial (∑ v ∈ P.support, k v • v) (∏ v ∈ P.support, P.coeff v ^ (k v)) := by
      simp_rw [MvPolynomial.monomial_pow]
      exact prod_monomial _ _ _
    rw [hp, show (↑(Nat.multinomial P.support k) : MvPolynomial (Fin 2) ℂ)
          = C ((Nat.multinomial P.support k : ℂ)) from (map_natCast C _).symm,
        MvPolynomial.C_mul_monomial]
  rw [key2, E_sum]
  exact Finset.sum_congr rfl (fun k _ => by rw [E_monomial])

/-- The charge of a channel's radial exponent `∑ᵢ kᵢ • vᵢ` is `∑ᵢ kᵢ · charge vᵢ`.  A channel is
**balanced** (`R_m` in THM-2022) exactly when this vanishes; `wt` is nonzero only there, so
`wick_expansion` is effectively the balanced-channel sum `M_m`. -/
lemma charge_radial {ι : Type*} (S : Finset ι) (k : ι → ℕ) (v : ι → (Fin 2 →₀ ℕ)) :
    charge (∑ i ∈ S, k i • v i) = ∑ i ∈ S, (k i : ℤ) * charge (v i) := by
  classical
  induction S using Finset.induction with
  | empty => simp [charge]
  | insert a s ha ih =>
      rw [Finset.sum_insert ha, Finset.sum_insert ha, ← ih]
      simp only [charge, Finsupp.add_apply, Finsupp.smul_apply, smul_eq_mul]
      push_cast; ring

/-- **The moment is the balanced-channel sum `M_m`.** Only *balanced* channels — those whose radial
exponent `∑ v, k v • v` has charge `0` (the set `R_m` of THM-2022) — contribute to `E (P^m)`, since
the Wick weight `wt` vanishes off the charge-0 diagonal.  This is `wick_expansion` restricted to
`R_m`. -/
theorem wick_expansion_balanced (P : MvPolynomial (Fin 2) ℂ) (m : ℕ) :
    E (P ^ m) = ∑ k ∈ (P.support.piAntidiag m).filter
        (fun k => charge (∑ v ∈ P.support, k v • v) = 0),
      (Nat.multinomial P.support k : ℂ) * (∏ v ∈ P.support, P.coeff v ^ (k v))
        * wt (∑ v ∈ P.support, k v • v) := by
  classical
  rw [wick_expansion]
  refine (Finset.sum_subset (Finset.filter_subset _ _) (fun k hkw hkf => ?_)).symm
  rw [Finset.mem_filter, not_and] at hkf
  rw [wt_of_charge_ne (hkf hkw), mul_zero]

end GMC2

-- Axiom audit: should be only propext / Classical.choice / Quot.sound (Mathlib standard).
#print axioms GMC2.mathieuZhao_of_charge_pos
#print axioms GMC2.mathieuZhao_of_charge_neg
#print axioms GMC2.moments_zero_of_charge_oneSided
#print axioms GMC2.mathieuZhao_of_nc2At
#print axioms GMC2.not_chargeOneSided_iff
#print axioms GMC2.gmc2_of_nc2
#print axioms GMC2.multinomial_dilate_modEq
#print axioms GMC2.dvd_choose_of_dvd
#print axioms GMC2.multinomial_dvd_of_exists_not_dvd
#print axioms GMC2.factorial_dilate_dvd
#print axioms GMC2.sum_natCast_mul_pow_char
#print axioms GMC2.face_sum_frobenius
#print axioms GMC2.face_sum_ne_zero
#print axioms GMC2.wick_expansion
#print axioms GMC2.charge_radial
#print axioms GMC2.wick_expansion_balanced
