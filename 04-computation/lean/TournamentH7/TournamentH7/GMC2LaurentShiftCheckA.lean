import TournamentH7.GMC2ConstantTermRelations
import TournamentH7.GMC2OrbitProduct

/-!
# Check A for THM-2067: a Laurent constant term is an ordinary coefficient

If `R` is an ordinary polynomial and

`Λ(T) = T⁻ᴹ R(T)`,

then the constant coefficient of `Λ^m` is the coefficient of `X^(M*m)` in
`R^m`.  This is the coefficient-shift identity called **Check A** in the
THM-2067 / DvdK formalization map.  Keeping it over an arbitrary commutative
semiring makes the bridge independent of the later complex/Galois argument.
-/

open Polynomial
open scoped BigOperators

namespace GMC2LaurentShiftCheckA

open LaurentPolynomial

variable {A : Type*} [CommSemiring A]

/-- Multiplying a Laurent polynomial by `T^s` translates its coefficient at
`k` to the old coefficient at `k-s`. -/
theorem coeff_mul_T (f : LaurentPolynomial A) (s k : ℤ) :
    (f * T s) k = f (k - s) := by
  simp [T, sub_eq_add_neg]

/-- **Check A.**  For `Λ = R(T) T^{-M}`, `CT(Λ^m) = [X^{Mm}]R^m`.

The left side uses evaluation at exponent `0` because a Laurent polynomial is
implemented as a finitely supported function `ℤ → A`; the right side is the
ordinary polynomial coefficient API consumed by the THM-2067 generating
polynomial `Φ(X) = X^M - t R(X)`. -/
theorem constantCoeff_shifted_pow_eq_coeff_pow (R : A[X]) (M m : ℕ) :
    ((Polynomial.toLaurent R * T (-(M : ℤ))) ^ m) 0 =
      (R ^ m).coeff (M * m) := by
  rw [mul_pow, map_pow, T_pow, coeff_mul_T]
  simp [Polynomial.toLaurent_apply, Nat.mul_comm]

/-- A product of Laurent monomials retains exactly two coordinates: the
product of coefficients and the sum of exponents. -/
theorem prod_monomial_eq (s : Finset ι) (q : ι → ℤ) (c : ι → A) (r : ι → ℕ) :
    ∏ i ∈ s, (C (c i) * T (q i)) ^ r i =
      C (∏ i ∈ s, c i ^ r i) * T (∑ i ∈ s, (r i : ℤ) * q i) := by
  classical
  induction s using Finset.induction_on with
  | empty => simp
  | @insert a s ha ih =>
      simp only [Finset.prod_insert ha, Finset.sum_insert ha, ih, mul_pow, map_pow, T_pow,
        map_mul]
      rw [← T_add]
      ring

/-- The universal relation in `GMC2ConstantTermRelations` really is the
constant coefficient of the corresponding finite Laurent polynomial power.
This is the semantic half of Check A and fixes the exact interface to
`GMC2DvdKInterface.DvdK1`. -/
theorem constantCoeff_pow_eq_aeval_constantTermRelation
    {ι : Type*} [Fintype ι] [DecidableEq ι] [Algebra ℚ A]
    (q : ι → ℤ) (c : ι → A) (m : ℕ) :
    ((∑ i : ι, C (c i) * T (q i)) ^ m) 0 =
      MvPolynomial.aeval c
        (GMC2ConstantTermRelations.constantTermRelation q m) := by
  classical
  rw [Finset.sum_pow_eq_sum_piAntidiag]
  rw [GMC2ConstantTermRelations.aeval_constantTermRelation]
  simp only [Finsupp.sum_apply]
  apply Finset.sum_congr rfl
  intro r hr
  rw [prod_monomial_eq]
  by_cases hq : GMC2ConstantTermRelations.totalCharge q r = 0
  · simp [GMC2ConstantTermRelations.totalCharge, hq]
  · simp [GMC2ConstantTermRelations.totalCharge, hq]

/-- Shift every Laurent exponent upward by `M` and regard the result as an
ordinary polynomial.  A lower-bound hypothesis below guarantees that
`Int.toNat` discards no exponent information. -/
noncomputable def shiftedPolynomial
    {ι : Type*} [Fintype ι] (q : ι → ℤ) (c : ι → A) (M : ℕ) : A[X] :=
  ∑ i : ι, Polynomial.monomial (q i + (M : ℤ)).toNat (c i)

/-- Shifting `shiftedPolynomial` back down by `M` recovers the original
Laurent polynomial, provided every shifted exponent is nonnegative. -/
theorem toLaurent_shiftedPolynomial_mul_T
    {ι : Type*} [Fintype ι] [DecidableEq ι]
    (q : ι → ℤ) (c : ι → A) (M : ℕ)
    (hM : ∀ i, -(M : ℤ) ≤ q i) :
    Polynomial.toLaurent (shiftedPolynomial q c M) * T (-(M : ℤ)) =
      ∑ i : ι, C (c i) * T (q i) := by
  classical
  simp only [shiftedPolynomial, map_sum, Finset.sum_mul]
  apply Finset.sum_congr rfl
  intro i hi
  rw [Polynomial.toLaurent_C_mul_T, ← mul_assoc, ← T_add]
  congr 2
  have hnonneg : 0 ≤ q i + (M : ℤ) := by omega
  rw [Int.toNat_of_nonneg hnonneg]
  omega

/-- **Repository-facing Check A.**  If `R` is obtained by shifting the finite
charge support upward by any common `M`, then its central coefficient at every
power is exactly the evaluated universal constant-term relation.

No injectivity, nonzero-coefficient, or straddling hypothesis is needed for
this identity; those hypotheses enter only in the later nonvanishing theorem.
-/
theorem coeff_shiftedPolynomial_pow_eq_aeval_constantTermRelation
    {ι : Type*} [Fintype ι] [DecidableEq ι] [Algebra ℚ A]
    (q : ι → ℤ) (c : ι → A) (M m : ℕ)
    (hM : ∀ i, -(M : ℤ) ≤ q i) :
    ((shiftedPolynomial q c M) ^ m).coeff (M * m) =
      MvPolynomial.aeval c
        (GMC2ConstantTermRelations.constantTermRelation q m) := by
  rw [← constantCoeff_shifted_pow_eq_coeff_pow]
  rw [toLaurent_shiftedPolynomial_mul_T q c M hM]
  exact constantCoeff_pow_eq_aeval_constantTermRelation q c m

end GMC2LaurentShiftCheckA

#print axioms GMC2LaurentShiftCheckA.coeff_mul_T
#print axioms GMC2LaurentShiftCheckA.constantCoeff_shifted_pow_eq_coeff_pow
#print axioms GMC2LaurentShiftCheckA.constantCoeff_pow_eq_aeval_constantTermRelation
#print axioms GMC2LaurentShiftCheckA.coeff_shiftedPolynomial_pow_eq_aeval_constantTermRelation

namespace GMC2AdditiveOrbitSum

open Finset MulAction

variable {G Ω B : Type*} [Group G] [Fintype G] [MulAction G Ω]
  [Fintype Ω] [DecidableEq Ω]

/-- Additive companion to `GMC2OrbitProduct.prod_smul_eq_prod_pow_card_stabilizer`.
For a transitive finite action, summing a function along one group orbit
counts every point with the cardinality of a stabilizer. -/
theorem sum_smul_eq_card_stabilizer_nsmul [IsPretransitive G Ω] [AddCommMonoid B]
    (f : Ω → B) (x : Ω) :
    ∑ g : G, f (g • x) =
      Fintype.card (stabilizer G x) • ∑ α : Ω, f α := by
  rw [← Fintype.sum_fiberwise (fun g : G ↦ g • x) (fun g ↦ f (g • x))]
  have key : ∀ y : Ω,
      (∑ i : {g : G // g • x = y}, f ((i : G) • x)) =
        Fintype.card (stabilizer G x) • f y := by
    intro y
    rw [Finset.sum_congr rfl (fun i _ ↦ congrArg f i.2)]
    simp [GMC2OrbitProduct.card_fiber_smul_eq_card_stabilizer x y]
  rw [Finset.sum_congr rfl (fun y _ ↦ key y), Finset.smul_sum]

/-- Uniform incidence identity for a subset all of whose group translates
have the same weighted sum `a`.  This is the additive Galois-orbit mechanism:
the left side counts translates, while the right side counts root incidence.
-/
theorem card_nsmul_translateSum_eq [IsPretransitive G Ω] [AddCommMonoid B]
    (f : Ω → B) (S : Finset Ω) (x : Ω) (a : B)
    (htranslate : ∀ g : G, ∑ β ∈ S, f (g • β) = a) :
    Fintype.card G • a =
      (S.card * Fintype.card (stabilizer G x)) • ∑ α : Ω, f α := by
  calc
    Fintype.card G • a = ∑ _g : G, a := by simp
    _ = ∑ g : G, ∑ β ∈ S, f (g • β) :=
      Finset.sum_congr rfl (fun g _ ↦ (htranslate g).symm)
    _ = ∑ β ∈ S, ∑ g : G, f (g • β) := by rw [Finset.sum_comm]
    _ = ∑ β ∈ S,
        Fintype.card (stabilizer G x) • ∑ α : Ω, f α := by
      apply Finset.sum_congr rfl
      intro β hβ
      rw [sum_smul_eq_card_stabilizer_nsmul f β,
        GMC2OrbitProduct.card_stabilizer_eq_card_stabilizer β x]
    _ = (S.card * Fintype.card (stabilizer G x)) • ∑ α : Ω, f α := by
      simp [mul_nsmul]

/-- **Additive orbit contradiction.**  Over a characteristic-zero field, no
subset can have every translate of its weighted sum equal to `1` while the
full-orbit weighted sum is `0`.

In the proposed THM-2067 route, `f(α)=α^(M-1)/Φ'(α)`, `S` is the small-root
subset, the germ identity gives every translated sum as `1`, and Lagrange
interpolation gives the full-root sum as `0`.
-/
theorem translateSum_one_ne_fullSum_zero
    [IsPretransitive G Ω] {K : Type*} [Field K] [CharZero K]
    (f : Ω → K) (S : Finset Ω) (x : Ω)
    (htranslate : ∀ g : G, ∑ β ∈ S, f (g • β) = 1)
    (hfull : ∑ α : Ω, f α = 0) : False := by
  have h := card_nsmul_translateSum_eq f S x 1 htranslate
  rw [hfull, nsmul_zero] at h
  have hcast : (Fintype.card G : K) = 0 := by simpa using h
  exact (Nat.cast_ne_zero.mpr Fintype.card_ne_zero) hcast

end GMC2AdditiveOrbitSum

#print axioms GMC2AdditiveOrbitSum.card_nsmul_translateSum_eq
#print axioms GMC2AdditiveOrbitSum.translateSum_one_ne_fullSum_zero

namespace GMC2FullRootLagrange

open Polynomial
open scoped BigOperators

/-- Lagrange's top-coefficient identity specialized below the top degree:
for `k+1 < |s|`, the sum of `α^k` divided by the nodal derivative is zero.
This is the residue-at-infinity identity needed by the additive THM-2067
route, stated directly for a finite set of distinct roots. -/
theorem sum_pow_div_derivative_nodal_eq_zero
    {K : Type*} [Field K] (s : Finset K) (k : ℕ)
    (hk : k + 1 < s.card) :
    ∑ α ∈ s, α ^ k /
        (Polynomial.derivative (Lagrange.nodal s id)).eval α = 0 := by
  classical
  have hdeg : ((Polynomial.X : K[X]) ^ k).degree < s.card := by
    rw [Polynomial.degree_X_pow]
    exact_mod_cast (lt_trans (Nat.lt_succ_self k) hk)
  have hlag := Lagrange.coeff_eq_sum (s := s) (v := id)
    (P := (Polynomial.X : K[X]) ^ k) Function.injective_id.injOn hdeg
  have hcoeff : ((Polynomial.X : K[X]) ^ k).coeff (s.card - 1) = 0 := by
    rw [Polynomial.coeff_X_pow]
    simp only [ite_eq_right_iff]
    intro heq
    omega
  rw [hcoeff] at hlag
  rw [← hlag]
  apply Finset.sum_congr rfl
  intro α hα
  rw [Lagrange.eval_nodal_derivative_eval_node_eq hα, Lagrange.eval_nodal]
  simp

/-- The exact exponent used by THM-2067: if `0 < M < |s|`, then the full
root sum of `α^(M-1)/Φ'(α)` vanishes for the monic nodal polynomial `Φ`.
-/
theorem sum_pow_pred_div_derivative_nodal_eq_zero
    {K : Type*} [Field K] (s : Finset K) (M : ℕ)
    (hM0 : 0 < M) (hMcard : M < s.card) :
    ∑ α ∈ s, α ^ (M - 1) /
        (Polynomial.derivative (Lagrange.nodal s id)).eval α = 0 := by
  apply sum_pow_div_derivative_nodal_eq_zero s (M - 1)
  omega

end GMC2FullRootLagrange

#print axioms GMC2FullRootLagrange.sum_pow_div_derivative_nodal_eq_zero
#print axioms GMC2FullRootLagrange.sum_pow_pred_div_derivative_nodal_eq_zero
