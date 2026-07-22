import TournamentH7.GMC2ConstantTermRelations

/-!
# Positive-coefficient DvdK: no cancellation, any support

For strictly positive real coefficients the constant term of `f^m` is a sum of *positive* terms, so it is
strictly positive as soon as ONE balanced composition of size `m` exists — cancellation cannot occur. This
is the general (any-charge) positive case of the DvdK existence input, DvdK-premise-free. The hard part of
the complex-coefficient theorem is exactly the cancellation this rules out; codex THM-2067 (Galois orbit
product) handles that.
-/

open MvPolynomial Finset

namespace GMC2DvdKPositive

/-- **Positive-coefficient constant-term positivity.** With `c i > 0` and a balanced composition `r0` of
size `m`, `CT(f^m) > 0`. -/
theorem ct_pos_of_balanced {ι : Type*} [Fintype ι] [DecidableEq ι]
    (q : ι → ℤ) (c : ι → ℝ) (hc : ∀ i, 0 < c i) (m : ℕ)
    (r0 : ι → ℕ)
    (hr0 : r0 ∈ Finset.piAntidiag (Finset.univ : Finset ι) m)
    (hbal : GMC2ConstantTermRelations.totalCharge q r0 = 0) :
    0 < MvPolynomial.aeval c (GMC2ConstantTermRelations.constantTermRelation q m) := by
  rw [GMC2ConstantTermRelations.aeval_constantTermRelation]
  apply Finset.sum_pos'
  · intro r _
    split_ifs with h
    · exact mul_nonneg (Nat.cast_nonneg _)
        (Finset.prod_nonneg (fun i _ => (pow_pos (hc i) _).le))
    · exact le_refl 0
  · refine ⟨r0, hr0, ?_⟩
    rw [if_pos hbal]
    apply mul_pos
    · exact_mod_cast Nat.multinomial_pos _ _
    · exact Finset.prod_pos (fun i _ => pow_pos (hc i) _)

/-- Consequently, given a balanced composition of positive size, `CT(f^m) ≠ 0`. -/
theorem ct_ne_zero_of_balanced {ι : Type*} [Fintype ι] [DecidableEq ι]
    (q : ι → ℤ) (c : ι → ℝ) (hc : ∀ i, 0 < c i) (m : ℕ)
    (r0 : ι → ℕ)
    (hr0 : r0 ∈ Finset.piAntidiag (Finset.univ : Finset ι) m)
    (hbal : GMC2ConstantTermRelations.totalCharge q r0 = 0) :
    MvPolynomial.aeval c (GMC2ConstantTermRelations.constantTermRelation q m) ≠ 0 :=
  (ct_pos_of_balanced q c hc m r0 hr0 hbal).ne'

/-- **Feasibility.** A two-sided support (a `+` charge at `i`, a `-` charge at `j`) always admits a balanced
composition of positive size: `|q j|` copies of the `+` charge and `|q i|` copies of the `-` charge. -/
theorem exists_balanced_of_twosided {ι : Type*} [Fintype ι] [DecidableEq ι]
    (q : ι → ℤ) (i j : ι) (hi : 0 < q i) (hj : q j < 0) :
    ∃ (m : ℕ) (r : ι → ℕ), 1 ≤ m ∧
      r ∈ Finset.piAntidiag (Finset.univ : Finset ι) m ∧
      GMC2ConstantTermRelations.totalCharge q r = 0 := by
  set r : ι → ℕ := fun k =>
    (if k = i then (q j).natAbs else 0) + (if k = j then (q i).natAbs else 0) with hr
  have hsum : ∑ k, r k = (q j).natAbs + (q i).natAbs := by
    simp only [hr, Finset.sum_add_distrib, Finset.sum_ite_eq' Finset.univ i,
      Finset.sum_ite_eq' Finset.univ j, Finset.mem_univ, if_true]
  refine ⟨(q j).natAbs + (q i).natAbs, r, ?_, ?_, ?_⟩
  · have : 0 < (q j).natAbs := Int.natAbs_pos.mpr hj.ne
    omega
  · rw [Finset.mem_piAntidiag]
    exact ⟨hsum, by intro k _; exact Finset.mem_univ k⟩
  · -- totalCharge = |q j| * q i + |q i| * q j = -q j * q i + q i * q j = 0
    simp only [GMC2ConstantTermRelations.totalCharge, hr]
    push_cast
    simp only [add_mul, Finset.sum_add_distrib, ite_mul, zero_mul,
      Finset.sum_ite_eq' Finset.univ, Finset.mem_univ, if_true]
    rw [abs_of_nonneg hi.le, abs_of_nonpos hj.le]; ring

/-- **Positive-coefficient DvdK1.** For any two-sided support and strictly positive real coefficients, some
positive power has a nonzero (indeed positive) constant term. No cancellation, DvdK-premise-free. -/
theorem dvdk1_positive {ι : Type*} [Fintype ι] [DecidableEq ι]
    (q : ι → ℤ) (c : ι → ℝ) (hc : ∀ i, 0 < c i)
    (i j : ι) (hi : 0 < q i) (hj : q j < 0) :
    ∃ m, 1 ≤ m ∧
      0 < MvPolynomial.aeval c (GMC2ConstantTermRelations.constantTermRelation q m) := by
  obtain ⟨m, r0, hm, hr0, hbal⟩ := exists_balanced_of_twosided q i j hi hj
  exact ⟨m, hm, ct_pos_of_balanced q c hc m r0 hr0 hbal⟩

end GMC2DvdKPositive

#print axioms GMC2DvdKPositive.ct_pos_of_balanced
#print axioms GMC2DvdKPositive.ct_ne_zero_of_balanced
#print axioms GMC2DvdKPositive.dvdk1_positive
