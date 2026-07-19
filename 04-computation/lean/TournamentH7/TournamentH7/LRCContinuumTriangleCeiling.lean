/-
  TournamentH7.LRCContinuumTriangleCeiling -- THM-1210 finite triangle core

  The continuum bad-set argument reduces to the additive three-band event

      J(p,q) = meas {u : {p*u}, {q*u}, {(p+q)*u} all lie in A},
      A = (1/7,2/7) union (3/7,4/7) union (5/7,6/7).

  This file contains the elementary arithmetic tail cutoff, the selection of
  a non-arithmetic triangle from four ordered integer frequencies, and an
  exact rational carry evaluator for the 99-pair reduced core
  `p < q`, `p+q <= 25`.

  The evaluator is deliberately a finite arithmetic object.  Its formula is
  the one-dimensional carry quotient of the six triangular cells of the
  three-band event.  Identifying that formula with Haar measure, proving the
  shifted-grid discrepancy bound, and proving the BAD-to-three-band
  implication belong to the analytic bridge, not to the kernel `decide`
  certificate below.
-/

import Mathlib

namespace LonelyRunner
namespace LRC14
namespace ContinuumTriangleCeiling

set_option maxHeartbeats 1000000

/-! ## The uniform tail arithmetic -/

/-- Sharper analytic tail after conditioning on the total frequency
`N = p+q`: the area `3/49` plus the total-frequency discrepancy `6/(7N)`
fits in the `2/21` budget as soon as `N >= 26`. -/
theorem sum_tail_cutoff (N : ℚ) (hN : 26 ≤ N) :
    (3 : ℚ) / 49 + 6 / (7 * N) < 2 / 21 := by
  have hden : (0 : ℚ) < 7 * N := by positivity
  have hterm : (6 : ℚ) / (7 * N) ≤ 3 / 91 := by
    apply (div_le_iff₀ hden).2
    nlinarith
  calc
    (3 : ℚ) / 49 + 6 / (7 * N) ≤ 3 / 49 + 3 / 91 :=
      add_le_add_right hterm _
    _ < 2 / 21 := by norm_num

/-- At the first tail value `N=26`, the upper bound is `60/637`, leaving
the exact positive margin `2/1911` below `2/21`. -/
theorem sum_tail_cutoff_endpoint :
    (3 : ℚ) / 49 + 6 / (7 * 26) = 60 / 637 ∧
      (2 : ℚ) / 21 - 60 / 637 = 2 / 1911 := by
  norm_num

/-- The exact arithmetic split used by the sharper `N = p+q` conditioning:
every pair is either in the finite range `N <= 25` or in the analytic tail
`N >= 26`. -/
theorem sum_core_or_tail (p q : ℕ) : p + q ≤ 25 ∨ 26 ≤ p + q := by
  omega

/-! ## Selecting a non-arithmetic triangle -/

/-- Three integer frequencies form an arithmetic progression when their two
successive gaps agree. -/
def IsAP3 (a b c : ℤ) : Prop := b - a = c - b

/-- Four ordered frequencies form an arithmetic progression when both
consecutive triples do. -/
def IsAP4 (a b c d : ℤ) : Prop := IsAP3 a b c ∧ IsAP3 b c d

/-- A non-AP quadruple has a non-AP consecutive triple.  This is the minimal
logical form useful when the surrounding proof already distinguishes the
four-term arithmetic-progression equality case. -/
theorem nonAP_consecutive_triangle {a b c d : ℤ}
    (h : ¬ IsAP4 a b c d) :
    ¬ IsAP3 a b c ∨ ¬ IsAP3 b c d := by
  simpa [IsAP4] using not_and_or.mp h

/-- In fact every four strictly ordered frequencies contain a non-AP triple,
even when the four frequencies themselves form a four-term AP.  One of the
two endpoint triples `(a,b,d)` and `(a,c,d)` must have unequal successive
gaps: equality for both would force `b = c`. -/
theorem sorted_four_has_nonAP_triangle {a b c d : ℤ}
    (hab : a < b) (hbc : b < c) (hcd : c < d) :
    ∃ x y z : ℤ,
      x < y ∧ y < z ∧
      ((x = a ∧ y = b ∧ z = d) ∨ (x = a ∧ y = c ∧ z = d)) ∧
      ¬ IsAP3 x y z := by
  by_cases h : b - a = d - b
  · refine ⟨a, c, d, lt_trans hab hbc, hcd, Or.inr ⟨rfl, rfl, rfl⟩, ?_⟩
    intro hacd
    simp only [IsAP3] at hacd
    omega
  · exact ⟨a, b, d, hab, lt_trans hbc hcd,
      Or.inl ⟨rfl, rfl, rfl⟩, by simpa [IsAP3] using h⟩

/-! ## Exact finite-core carry evaluator -/

/-- The three residue classes contributed by the upper triangular cells
`(1,1)`, `(1,3)`, and `(3,1)`.  The reflected lower cells account for the
factor two in `triangleMeasure`. -/
def triangleResidues (p q : ℕ) : List ℤ :=
  [2 * (p : ℤ) - q, 4 * (p : ℤ) - q, 2 * (p : ℤ) - 3 * q].map (· % 7)

/-- Multiplicity of the carry index `k mod 7` among the three upper-cell
residue classes. -/
def residueMultiplicity (p q k : ℕ) : ℕ :=
  (triangleResidues p q).count ((k : ℤ) % 7)

/-- The denominator-cleared length of the `k`-th diagonal slice in the
`p` by `q` carry rectangle.  Dividing this weight by `p*q` gives the actual
rational slice length. -/
def carryWeight (p q k : ℕ) : ℕ :=
  if k ≤ q then k * p else (p + q - k) * q

/-- The integer numerator in the carry formula. -/
def triangleNumerator (p q : ℕ) : ℕ :=
  ∑ k ∈ Finset.Ioo 0 (p + q), residueMultiplicity p q k * carryWeight p q k

/-- Exact carry formula for the additive three-band event `J(p,q)`:

`J(p,q) = 2*S/(7(p+q)pq)`, where
`S = sum_{0<k<p+q} c_k W_k`.

Here `c_k` is `residueMultiplicity` and `W_k` is `carryWeight`.  This
denominator-cleared presentation lets the finite core use kernel `decide` on
natural-number inequalities rather than evaluating rational comparisons. -/
def triangleMeasure (p q : ℕ) : ℚ :=
  ((2 * triangleNumerator p q : ℕ) : ℚ) /
    ((7 * (p + q) * p * q : ℕ) : ℚ)

/-- The cleared form of `triangleMeasure p q <= 2/21`. -/
def triangleBoundOK (p q : ℕ) : Prop :=
  3 * triangleNumerator p q ≤ (p + q) * p * q

instance (p q : ℕ) : Decidable (triangleBoundOK p q) :=
  inferInstanceAs (Decidable (3 * triangleNumerator p q ≤ (p + q) * p * q))

/-- The cleared form of equality `triangleMeasure p q = 2/21`. -/
def triangleEquality (p q : ℕ) : Prop :=
  3 * triangleNumerator p q = (p + q) * p * q

instance (p q : ℕ) : Decidable (triangleEquality p q) :=
  inferInstanceAs (Decidable (3 * triangleNumerator p q = (p + q) * p * q))

private theorem cleared_le_iff (numerator denominator : ℕ)
    (hdenominator : 0 < denominator) :
    3 * numerator ≤ denominator ↔
      (((2 * numerator : ℕ) : ℚ) / ((7 * denominator : ℕ) : ℚ) ≤ 2 / 21) := by
  have hdenNat : 0 < 7 * denominator := Nat.mul_pos (by norm_num) hdenominator
  have hden : (0 : ℚ) < ((7 * denominator : ℕ) : ℚ) := by
    exact_mod_cast hdenNat
  rw [div_le_iff₀ hden]
  constructor
  · intro h
    have hcast : (((3 * numerator : ℕ) : ℚ)) ≤ (denominator : ℚ) := by
      exact_mod_cast h
    push_cast at hcast ⊢
    norm_num at hcast ⊢
    linarith
  · intro h
    have hcast : (((3 * numerator : ℕ) : ℚ)) ≤ (denominator : ℚ) := by
      push_cast at h ⊢
      norm_num at h ⊢
      linarith
    exact_mod_cast hcast

private theorem cleared_eq_iff (numerator denominator : ℕ)
    (hdenominator : 0 < denominator) :
    3 * numerator = denominator ↔
      (((2 * numerator : ℕ) : ℚ) / ((7 * denominator : ℕ) : ℚ) = 2 / 21) := by
  have hdenNat : 0 < 7 * denominator := Nat.mul_pos (by norm_num) hdenominator
  have hden : (0 : ℚ) < ((7 * denominator : ℕ) : ℚ) := by
    exact_mod_cast hdenNat
  rw [div_eq_iff hden.ne']
  constructor
  · intro h
    have hcast : (((3 * numerator : ℕ) : ℚ)) = (denominator : ℚ) := by
      exact_mod_cast h
    push_cast at hcast ⊢
    norm_num at hcast ⊢
    linarith
  · intro h
    have hcast : (((3 * numerator : ℕ) : ℚ)) = (denominator : ℚ) := by
      push_cast at h ⊢
      norm_num at h ⊢
      linarith
    exact_mod_cast hcast

/-- For positive gaps, the natural-number predicate is exactly the rational
`2/21` ceiling for the carry evaluator. -/
theorem triangleBoundOK_iff_measure_le (p q : ℕ) (hp : 0 < p) (hq : 0 < q) :
    triangleBoundOK p q ↔ triangleMeasure p q ≤ 2 / 21 := by
  have hden : 0 < (p + q) * p * q := by positivity
  simpa only [triangleBoundOK, triangleMeasure, Nat.mul_assoc] using
    cleared_le_iff (triangleNumerator p q) ((p + q) * p * q) hden

/-- For positive gaps, the natural-number equality predicate is exactly
equality in the rational `2/21` ceiling. -/
theorem triangleEquality_iff_measure_eq (p q : ℕ) (hp : 0 < p) (hq : 0 < q) :
    triangleEquality p q ↔ triangleMeasure p q = 2 / 21 := by
  have hden : 0 < (p + q) * p * q := by positivity
  simpa only [triangleEquality, triangleMeasure, Nat.mul_assoc] using
    cleared_eq_iff (triangleNumerator p q) ((p + q) * p * q) hden

/-- The 99 reduced pairs left after the sharper `p+q >= 26` analytic tail. -/
def corePairs : Finset (ℕ × ℕ) :=
  ((Finset.Icc 1 24).product (Finset.Icc 2 24)).filter fun pair =>
    pair.1 < pair.2 ∧ pair.1 + pair.2 ≤ 25 ∧ Nat.Coprime pair.1 pair.2

def coreViolators : Finset (ℕ × ℕ) :=
  corePairs.filter fun pair => ¬ triangleBoundOK pair.1 pair.2

def coreEqualityPairs : Finset (ℕ × ℕ) :=
  corePairs.filter fun pair => triangleEquality pair.1 pair.2

/-- Exact checksum for the reduced finite universe. -/
theorem corePairs_card : corePairs.card = 99 := by
  set_option maxRecDepth 100000 in
    decide

/-- No reduced pair below the cutoff exceeds `2/21`. -/
theorem coreViolators_empty : coreViolators = ∅ := by
  set_option maxRecDepth 100000 in
    decide

/-- The unique reduced equality pair is `(1,2)`. -/
theorem coreEqualityPairs_exact : coreEqualityPairs = {(1, 2)} := by
  set_option maxRecDepth 100000 in
    decide

/-- Proof-facing form of the finite certificate. -/
theorem triangleBoundOK_of_core_mem {p q : ℕ}
    (hmem : (p, q) ∈ corePairs) :
    triangleBoundOK p q := by
  by_contra hbound
  have hbad : (p, q) ∈ coreViolators := by
    simp only [coreViolators, Finset.mem_filter]
    exact ⟨hmem, hbound⟩
  rw [coreViolators_empty] at hbad
  simp at hbad

/-- Convenient quantified wrapper: every positive coprime `p < q` with
`p+q <= 25`
obeys the exact finite-core ceiling. -/
theorem triangleBoundOK_core (p q : ℕ)
    (hp : 0 < p) (hpq : p < q) (hsum : p + q ≤ 25)
    (hcop : Nat.Coprime p q) :
    triangleBoundOK p q := by
  apply triangleBoundOK_of_core_mem
  rw [corePairs]
  apply Finset.mem_filter.mpr
  constructor
  · apply Finset.mem_product.mpr
    constructor
    · exact Finset.mem_Icc.mpr ⟨by omega, by omega⟩
    · exact Finset.mem_Icc.mpr ⟨by omega, by omega⟩
  · exact ⟨hpq, hsum, hcop⟩

/-- Rational proof-facing wrapper for every positive coprime
`p < q` with `p+q <= 25`. -/
theorem triangleMeasure_le_core (p q : ℕ)
    (hp : 0 < p) (hpq : p < q) (hsum : p + q ≤ 25)
    (hcop : Nat.Coprime p q) :
    triangleMeasure p q ≤ 2 / 21 := by
  have hqpos : 0 < q := by omega
  exact (triangleBoundOK_iff_measure_le p q hp hqpos).mp
    (triangleBoundOK_core p q hp hpq hsum hcop)

/-- Equality classification inside the reduced core, in cleared form. -/
theorem triangleEquality_iff_of_core_mem {p q : ℕ}
    (hmem : (p, q) ∈ corePairs) :
    triangleEquality p q ↔ (p, q) = (1, 2) := by
  constructor
  · intro heq
    have heqmem : (p, q) ∈ coreEqualityPairs := by
      exact Finset.mem_filter.mpr ⟨hmem, heq⟩
    rw [coreEqualityPairs_exact] at heqmem
    simpa using heqmem
  · intro hpq
    have heqmem : (p, q) ∈ coreEqualityPairs := by
      rw [coreEqualityPairs_exact]
      simpa [hpq]
    exact (Finset.mem_filter.mp heqmem).2

/-- Rational equality classification inside the reduced core. -/
theorem triangleMeasure_eq_top_iff_of_core_mem {p q : ℕ}
    (hmem : (p, q) ∈ corePairs) :
    triangleMeasure p q = 2 / 21 ↔ (p, q) = (1, 2) := by
  have hshape := hmem
  rw [corePairs] at hshape
  have hproduct := (Finset.mem_filter.mp hshape).1
  have hpBounds := Finset.mem_Icc.mp (Finset.mem_product.mp hproduct).1
  have hqBounds := Finset.mem_Icc.mp (Finset.mem_product.mp hproduct).2
  have hppos : 0 < p := by omega
  have hqpos : 0 < q := by omega
  exact (triangleEquality_iff_measure_eq p q hppos hqpos).symm.trans
    (triangleEquality_iff_of_core_mem hmem)

/-! ## Equality-locus rigidity -/

/-- Two positive gaps have the reduced equality ratio `(1,2)`, allowing
either orientation, or are equal. -/
def Ratio12 (x y : ℕ) : Prop := x = y ∨ y = 2 * x ∨ x = 2 * y

/-- If all four deletion triangles of the gap triple `(a,b,c)` lie in the
additive-triangle equality locus, then all three gaps agree.  This is the
finite combinatorial step that upgrades triangle equality to a four-term
arithmetic progression. -/
theorem ratio12_four_triangles_rigid (a b c : ℕ)
    (ha : 0 < a) (hb : 0 < b) (hc : 0 < c)
    (hab : Ratio12 a b) (hbc : Ratio12 b c)
    (ha_bc : Ratio12 a (b + c)) (hab_c : Ratio12 (a + b) c) :
    a = b ∧ b = c := by
  simp only [Ratio12] at hab hbc ha_bc hab_c
  rcases hab with hab | hab | hab <;>
    rcases hbc with hbc | hbc | hbc <;>
      rcases ha_bc with ha_bc | ha_bc | ha_bc <;>
        rcases hab_c with hab_c | hab_c | hab_c <;>
          omega

#print axioms sum_tail_cutoff
#print axioms sum_tail_cutoff_endpoint
#print axioms sorted_four_has_nonAP_triangle
#print axioms coreViolators_empty
#print axioms coreEqualityPairs_exact
#print axioms triangleBoundOK_core
#print axioms triangleMeasure_le_core
#print axioms triangleMeasure_eq_top_iff_of_core_mem
#print axioms ratio12_four_triangles_rigid

end ContinuumTriangleCeiling
end LRC14
end LonelyRunner
