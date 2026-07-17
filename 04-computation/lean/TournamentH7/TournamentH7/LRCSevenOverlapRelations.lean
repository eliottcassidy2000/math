/-
  The exact algebraic bridge from THM-947 seven-overlap witnesses to sparse
  integer relations.  Pair cross-products are small, and every bad triple
  satisfies a Plucker relation with those cross-products as coefficients.

  Tournament-analysis audit: vertices are bad runners, while an edge is
  colored by its integer cross-product.  Forgetting the color destroys the
  telescoping triangle law, so an uncolored tournament or support graph is not
  a faithful quotient.
-/

import TournamentH7.LRCArcWire

namespace LonelyRunner
namespace LRC14Concrete

open scoped Classical

/-- Cross-product color carried by a pair of failing runners. -/
def overlapDet (v : Fin 13 → ℤ) (q : ℕ) (p : ℕ)
    (i j : Fin 13) : ℤ :=
  v i * failWitness v q p j - v j * failWitness v q p i

/-- THM-947 with the common positive modulus cancelled. -/
theorem overlapDet_bound (v : Fin 13 → ℤ) (q p : ℕ) (i j : Fin 13)
    (hq : 0 < q) (hvi : v i ≠ 0) (hvj : v j ≠ 0)
    (hbadi : ¬ inBand v q p i) (hbadj : ¬ inBand v q p j) :
    14 * |overlapDet v q p i j| < |v i| + |v j| := by
  have h := seven_overlap_pair_constraint
    v q p i j hq hvi hvj hbadi hbadj
  have hqZ : (0 : ℤ) < (q : ℤ) := by exact_mod_cast hq
  have hscaled :
      (q : ℤ) * (14 * |overlapDet v q p i j|) <
        (q : ℤ) * (|v i| + |v j|) := by
    simpa [overlapDet, mul_assoc, mul_left_comm, mul_comm] using h
  by_contra hnot
  push Not at hnot
  have hproduct : 0 ≤ (q : ℤ) *
      (14 * |overlapDet v q p i j| - (|v i| + |v j|)) :=
    mul_nonneg hqZ.le (sub_nonneg.mpr hnot)
  nlinarith

theorem overlapDet_skew (v : Fin 13 → ℤ) (q p : ℕ) (i j : Fin 13) :
    overlapDet v q p j i = -overlapDet v q p i j := by
  simp [overlapDet]

/-- The edge colors on every triangle obey the exact Plucker relation. -/
theorem overlapDet_plucker (v : Fin 13 → ℤ) (q p : ℕ)
    (i j k : Fin 13) :
    overlapDet v q p j k * v i - overlapDet v q p i k * v j +
        overlapDet v q p i j * v k = 0 := by
  simp [overlapDet]
  ring

/-- Sparse coefficient vector exported to the dense-core relation traps. -/
def overlapTripleCoeff (v : Fin 13 → ℤ) (q p : ℕ)
    (i j k index : Fin 13) : ℤ :=
  (if index = i then overlapDet v q p j k else 0) +
  (if index = j then -overlapDet v q p i k else 0) +
  (if index = k then overlapDet v q p i j else 0)

/-- The sparse coefficient vector is an exact integer relation. -/
theorem overlapTripleCoeff_sum_eq_zero (v : Fin 13 → ℤ) (q p : ℕ)
    (i j k : Fin 13) :
    ∑ index, overlapTripleCoeff v q p i j k index * v index = 0 := by
  classical
  simp_rw [overlapTripleCoeff, add_mul]
  rw [Finset.sum_add_distrib, Finset.sum_add_distrib]
  simp only [ite_mul, zero_mul, Finset.sum_ite_eq', Finset.mem_univ, if_true]
  simpa [sub_eq_add_neg, neg_mul] using overlapDet_plucker v q p i j k

/-- Each of the three nonzero slots inherits the corresponding THM-947
cross-product bound whenever the whole triple is bad. -/
theorem overlapTripleCoeff_slot_bounds
    (v : Fin 13 → ℤ) (q p : ℕ) (i j k : Fin 13)
    (hq : 0 < q) (hvi : v i ≠ 0) (hvj : v j ≠ 0) (hvk : v k ≠ 0)
    (hbadi : ¬ inBand v q p i) (hbadj : ¬ inBand v q p j)
    (hbadk : ¬ inBand v q p k) :
    14 * |overlapDet v q p j k| < |v j| + |v k| ∧
    14 * |overlapDet v q p i k| < |v i| + |v k| ∧
    14 * |overlapDet v q p i j| < |v i| + |v j| := by
  exact ⟨overlapDet_bound v q p j k hq hvj hvk hbadj hbadk,
    overlapDet_bound v q p i k hq hvi hvk hbadi hbadk,
    overlapDet_bound v q p i j hq hvi hvj hbadi hbadj⟩

/-! ## Axiom audit -/

#print axioms overlapDet_bound
#print axioms overlapDet_plucker
#print axioms overlapTripleCoeff_sum_eq_zero
#print axioms overlapTripleCoeff_slot_bounds

end LRC14Concrete
end LonelyRunner
