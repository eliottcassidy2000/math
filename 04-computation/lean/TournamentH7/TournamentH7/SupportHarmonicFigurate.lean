/-
  TournamentH7.SupportHarmonicFigurate — elementary algebraic kernel for
  THM-2000 (support-harmonic Abel--Dini law and the figurate mass surface)
  and THM-2005 (support-Dirichlet automatic/tournament atlas).

  This module deliberately formalizes only the finite, algebraic claims used
  by the analytic proof.  It introduces no axioms and contains no placeholders.
  In particular it certifies:

    * the master-figurate binomial factorization and reciprocal factorization;
    * the concrete ordinary-polygonal factorization and reciprocal split;
    * triangular and triangular-square reciprocal decompositions;
    * the cross-multiplied square-pyramidal numerator cancellation;
    * the odd/even maximum-cyclic-triple denominator algebra and reciprocal
      decompositions, including the two condensation-ratio defects (but not
      the tournament extremality theorem itself);
    * the `k^2 - 1` ladder decomposition; and
    * the finite multiplicative-block reciprocal sandwich, including its
      dyadic specialization.

  Infinite sums, beta-integral interchange, Abel--Stieltjes integration, and
  the score-theoretic proof of the maximum-cyclic-triple formulas remain in
  the paper proof; none is postulated here.
-/

import Mathlib

namespace SupportHarmonicFigurate

/-! ## Master figurate denominator -/

/-- The algebraic denominator underlying the master figurate array.  In the
application, `top = C(m+d-2,d)` and `side = C(m+d-2,d-1)`. -/
def masterDenom (a _d _m top side : ℚ) : ℚ := a * top + side

/-- The binomial ratio `d * top = (m-1) * side` turns the master denominator
into the factor appearing in THM-2000, equation (16). -/
theorem master_denom_factor {a d m top side : ℚ}
    (hbinom : d * top = (m - 1) * side) :
    d * masterDenom a d m top side = (d + a * (m - 1)) * side := by
  unfold masterDenom
  calc
    d * (a * top + side) = a * (d * top) + d * side := by ring
    _ = a * ((m - 1) * side) + d * side := by rw [hbinom]
    _ = (d + a * (m - 1)) * side := by ring

/-- Reciprocal form of `master_denom_factor`.  The explicit nonvanishing
hypotheses are exactly the denominators appearing in the formula. -/
theorem master_reciprocal_factor {a d m top side : ℚ}
    (hbinom : d * top = (m - 1) * side)
    (hside : side ≠ 0)
    (hlinear : d + a * (m - 1) ≠ 0)
    (hmaster : masterDenom a d m top side ≠ 0) :
    1 / masterDenom a d m top side =
      (1 / side) * (d / (d + a * (m - 1))) := by
  have hfactor := master_denom_factor (a := a) hbinom
  field_simp [hside, hlinear, hmaster]
  linarith

/-! ## Ordinary polygonal denominators -/

/-- Polynomial continuation of the ordinary `h`-gonal denominator
`P_h(n) = ((h-2)n^2-(h-4)n)/2`. -/
def ordinaryPolygonalDenom (h n : ℚ) : ℚ :=
  ((h - 2) * n ^ 2 - (h - 4) * n) / 2

/-- Concrete linear factorization of the ordinary polygonal denominator. -/
theorem ordinary_polygonal_denom_factor (h n : ℚ) :
    ordinaryPolygonalDenom h n = n * ((h - 2) * (n - 1) + 2) / 2 := by
  unfold ordinaryPolygonalDenom
  ring

/-- The finite rational identity underlying THM-2005, equation (13).  The
hypotheses guarantee that every denominator in the displayed identity is
nonzero. -/
theorem ordinary_polygonal_reciprocal_split {h n : ℚ}
    (hh2 : h - 2 ≠ 0) (hh4 : h - 4 ≠ 0) (hn : n ≠ 0)
    (hlinear : (h - 2) * (n - 1) + 2 ≠ 0) :
    1 / ordinaryPolygonalDenom h n =
      2 / (h - 4) * (1 / (n - 1 + 2 / (h - 2)) - 1 / n) := by
  rw [ordinary_polygonal_denom_factor]
  field_simp [hh2, hh4, hn, hlinear]
  ring

/-! ## Faulhaber term identities -/

/-- Triangular reciprocal telescoping, the termwise core of mass `2`. -/
theorem triangular_reciprocal_split {n : ℚ}
    (hn : n ≠ 0) (hn1 : n + 1 ≠ 0) :
    2 / (n * (n + 1)) = 2 / n - 2 / (n + 1) := by
  field_simp [hn, hn1]
  ring

/-- Cross-multiplied numerator cancellation behind the square-pyramidal
partial fractions in THM-2000, equation (35).  Together with nonzero factors,
this is precisely the finite algebraic content of that rational identity. -/
theorem square_pyramidal_numerator_cancel (n : ℚ) :
    6 * (n + 1) * (2 * n + 1) + 6 * n * (2 * n + 1) -
      24 * n * (n + 1) = 6 := by
  ring

/-! ## Maximum cyclic-triple denominator algebra -/

/-- The rational odd-order branch obtained by substituting `2m+1` into
`n(n^2-1)/24`.  This definition records the extremal formula algebraically;
it does not assert the tournament extremality theorem. -/
def oddMaxC3Denom (m : ℚ) : ℚ :=
  (2 * m + 1) * ((2 * m + 1) ^ 2 - 1) / 24

/-- The rational even-order branch obtained by substituting `2m` into
`n(n^2-4)/24`. -/
def evenMaxC3Denom (m : ℚ) : ℚ :=
  (2 * m) * ((2 * m) ^ 2 - 4) / 24

/-- The odd maximum-`c3` denominator is the square-pyramidal/Faulhaber
denominator `m(m+1)(2m+1)/6`. -/
theorem odd_max_c3_denom_factor (m : ℚ) :
    oddMaxC3Denom m = m * (m + 1) * (2 * m + 1) / 6 := by
  unfold oddMaxC3Denom
  ring

/-- The even maximum-`c3` denominator factors as `m(m-1)(m+1)/3`. -/
theorem even_max_c3_denom_factor (m : ℚ) :
    evenMaxC3Denom m = m * (m - 1) * (m + 1) / 3 := by
  unfold evenMaxC3Denom
  ring

/-- Odd-branch partial fractions.  Finite summation of this identity is the
algebraic input for the odd parity mass in THM-2005. -/
theorem odd_max_c3_reciprocal_split {m : ℚ}
    (hm : m ≠ 0) (hmp1 : m + 1 ≠ 0) (h2mp1 : 2 * m + 1 ≠ 0) :
    1 / oddMaxC3Denom m =
      6 / m + 6 / (m + 1) - 24 / (2 * m + 1) := by
  rw [odd_max_c3_denom_factor]
  have hcomm : 2 * m + 1 = 1 + m * 2 := by ring
  rw [hcomm] at h2mp1 ⊢
  field_simp [hm, hmp1, h2mp1]
  ring

/-- Even-branch second-difference decomposition.  Its finite partial sums
telescope to boundary terms, giving the algebraic input for the even parity
mass in THM-2005. -/
theorem even_max_c3_reciprocal_split {m : ℚ}
    (hmm1 : m - 1 ≠ 0) (hm : m ≠ 0) (hmp1 : m + 1 ≠ 0) :
    1 / evenMaxC3Denom m =
      (3 / 2) * (1 / (m - 1) - 2 / m + 1 / (m + 1)) := by
  rw [even_max_c3_denom_factor]
  field_simp [hmm1, hm, hmp1]
  ring

/-- At even tournament order `2m`, the ratio of consecutive algebraic
maximum-`c3` denominators is `(2m-1)/(2m+2)`.  THM-2016 supplies the separate
tournament statement that this ratio is the reducibility temperature. -/
theorem even_condensation_ratio {m : ℚ}
    (hmm1 : m - 1 ≠ 0) (hm : m ≠ 0) (hmp1 : m + 1 ≠ 0) :
    oddMaxC3Denom (m - 1) / evenMaxC3Denom m =
      (2 * m - 1) / (2 * m + 2) := by
  rw [odd_max_c3_denom_factor, even_max_c3_denom_factor]
  field_simp [hmm1, hm, hmp1]
  ring

/-- At odd tournament order `2m+1`, the consecutive denominator ratio is
`(2m-2)/(2m+1)`. -/
theorem odd_condensation_ratio {m : ℚ}
    (hm : m ≠ 0) (hmp1 : m + 1 ≠ 0) :
    evenMaxC3Denom m / oddMaxC3Denom m =
      (2 * m - 2) / (2 * m + 1) := by
  rw [even_max_c3_denom_factor, odd_max_c3_denom_factor]
  field_simp [hm, hmp1]
  ring

/-- The even condensation defect is the harmonic term `3/(2m+2)`. -/
theorem even_condensation_defect {m : ℚ}
    (hmm1 : m - 1 ≠ 0) (hm : m ≠ 0) (hmp1 : m + 1 ≠ 0) :
    1 - oddMaxC3Denom (m - 1) / evenMaxC3Denom m =
      3 / (2 * m + 2) := by
  rw [even_condensation_ratio hmm1 hm hmp1]
  field_simp [hmp1]
  ring

/-- The odd condensation defect is the harmonic term `3/(2m+1)`. -/
theorem odd_condensation_defect {m : ℚ}
    (hm : m ≠ 0) (hmp1 : m + 1 ≠ 0) (h2mp1 : 2 * m + 1 ≠ 0) :
    1 - evenMaxC3Denom m / oddMaxC3Denom m =
      3 / (2 * m + 1) := by
  rw [odd_condensation_ratio hm hmp1]
  field_simp [h2mp1]
  ring

/-- Polynomial continuation of the square of the `n`-th triangular number. -/
def triangularSquare (n : ℚ) : ℚ := (n * (n + 1) / 2) ^ 2

/-- The reciprocal decomposition used for the cubic Faulhaber mass
`4*pi^2/3 - 12` (THM-2000, equation (36), including the outer factor `4`). -/
theorem triangular_square_reciprocal_split {n : ℚ}
    (hn : n ≠ 0) (hn1 : n + 1 ≠ 0) :
    1 / triangularSquare n =
      4 * (1 / n ^ 2 + 1 / (n + 1) ^ 2 - 2 / (n * (n + 1))) := by
  unfold triangularSquare
  field_simp [hn, hn1]
  ring

/-- Partial-sum ladder decomposition for `k^2-1`; summing this identity gives
the reciprocal mass `3/4` after the two boundary terms are evaluated. -/
theorem ladder_sum_reciprocal_split {k : ℚ}
    (hkm1 : k - 1 ≠ 0) (hkp1 : k + 1 ≠ 0) :
    1 / (k ^ 2 - 1) = (1 / 2) * (1 / (k - 1) - 1 / (k + 1)) := by
  have hsq : k ^ 2 - 1 ≠ 0 := by
    intro h
    have hprod : (k - 1) * (k + 1) = 0 := by nlinarith
    exact (mul_ne_zero hkm1 hkp1) hprod
  field_simp [hkm1, hkp1, hsq]
  ring

/-! ## Finite multiplicative-block sandwich -/

/-- If every positive value in a finite block lies between `lo` and `hi`, its
reciprocal mass lies between `card / hi` and `card / lo`.  This is the finite
pointwise core of the multiplicative-block Dini criterion. -/
theorem reciprocal_block_bounds {ι : Type*} (S : Finset ι) (x : ι → ℝ)
    (lo hi : ℝ) (hlo : 0 < lo)
    (hlower : ∀ i ∈ S, lo ≤ x i)
    (hupper : ∀ i ∈ S, x i ≤ hi) :
    (S.card : ℝ) / hi ≤ ∑ i ∈ S, 1 / x i ∧
      ∑ i ∈ S, 1 / x i ≤ (S.card : ℝ) / lo := by
  by_cases hS : S.Nonempty
  · have hxi : ∀ i ∈ S, 0 < x i := fun i hiS => lt_of_lt_of_le hlo (hlower i hiS)
    have hhi : 0 < hi := by
      obtain ⟨i, hiS⟩ := hS
      exact lt_of_lt_of_le (hxi i hiS) (hupper i hiS)
    constructor
    · calc
        (S.card : ℝ) / hi = ∑ i ∈ S, 1 / hi := by
          simp [div_eq_mul_inv]
        _ ≤ ∑ i ∈ S, 1 / x i := by
          exact Finset.sum_le_sum fun i hiS =>
            one_div_le_one_div_of_le (hxi i hiS) (hupper i hiS)
    · calc
        ∑ i ∈ S, 1 / x i ≤ ∑ i ∈ S, 1 / lo := by
          exact Finset.sum_le_sum fun i hiS =>
            one_div_le_one_div_of_le hlo (hlower i hiS)
        _ = (S.card : ℝ) / lo := by
          simp [div_eq_mul_inv]
  · rw [Finset.not_nonempty_iff_eq_empty] at hS
    subst S
    simp

/-- Dyadic specialization: a block in `[2^k,2^(k+1)]` contributes between
`card/2^(k+1)` and `card/2^k`. -/
theorem dyadic_reciprocal_block_bounds {ι : Type*} (S : Finset ι)
    (x : ι → ℝ) (k : ℕ)
    (hlower : ∀ i ∈ S, (2 : ℝ) ^ k ≤ x i)
    (hupper : ∀ i ∈ S, x i ≤ (2 : ℝ) ^ (k + 1)) :
    (S.card : ℝ) / (2 : ℝ) ^ (k + 1) ≤ ∑ i ∈ S, 1 / x i ∧
      ∑ i ∈ S, 1 / x i ≤ (S.card : ℝ) / (2 : ℝ) ^ k := by
  exact reciprocal_block_bounds S x ((2 : ℝ) ^ k) ((2 : ℝ) ^ (k + 1))
    (pow_pos (by norm_num) k) hlower hupper

#print axioms master_denom_factor
#print axioms master_reciprocal_factor
#print axioms ordinary_polygonal_denom_factor
#print axioms ordinary_polygonal_reciprocal_split
#print axioms triangular_reciprocal_split
#print axioms square_pyramidal_numerator_cancel
#print axioms odd_max_c3_denom_factor
#print axioms even_max_c3_denom_factor
#print axioms odd_max_c3_reciprocal_split
#print axioms even_max_c3_reciprocal_split
#print axioms even_condensation_ratio
#print axioms odd_condensation_ratio
#print axioms even_condensation_defect
#print axioms odd_condensation_defect
#print axioms triangular_square_reciprocal_split
#print axioms ladder_sum_reciprocal_split
#print axioms reciprocal_block_bounds
#print axioms dyadic_reciprocal_block_bounds

end SupportHarmonicFigurate
