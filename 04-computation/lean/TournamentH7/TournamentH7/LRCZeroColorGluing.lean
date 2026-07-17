import TournamentH7.LRCSevenOverlapRelations

/-!
# Zero-color witness gluing

At fixed `(q,p)`, a zero determinant says that two witness-speed vectors have
the same rational slope.  For nonzero speeds this relation is reflexive,
symmetric, and transitive.  Hence every connected zero-color stalk is already
a zero-color clique and carries one common rational witness parameter.

This closes the local gluing primitive left open by `LRCAlignedResonance`.
It does not bound how many multiplier events share the same parameter.  The
faithful vertices are witness-speed vectors; the zero/nonzero determinant is
the binary switch and connected components are the tie classes.  A tournament
orientation of runners preserves only an order between classes and destroys
the common rational parameter inside each class.
-/

namespace LonelyRunner
namespace LRC14Concrete

open scoped Classical

/-- Rational nearest-integer witness slope at one multiplier event. -/
def witnessSlope (v : Fin 13 → ℤ) (q p : ℕ) (index : Fin 13) : ℚ :=
  (failWitness v q p index : ℚ) / (v index : ℚ)

/-- Zero determinant is exactly equality of rational witness slopes. -/
theorem overlapDet_eq_zero_iff_witnessSlope_eq
    (v : Fin 13 → ℤ) (q p : ℕ) (first second : Fin 13)
    (hfirst : v first ≠ 0) (hsecond : v second ≠ 0) :
    overlapDet v q p first second = 0 ↔
      witnessSlope v q p first = witnessSlope v q p second := by
  have hfirstQ : (v first : ℚ) ≠ 0 := by exact_mod_cast hfirst
  have hsecondQ : (v second : ℚ) ≠ 0 := by exact_mod_cast hsecond
  constructor
  · intro hzero
    rw [witnessSlope, witnessSlope,
      div_eq_div_iff hfirstQ hsecondQ]
    have hcross :
        v first * failWitness v q p second =
          v second * failWitness v q p first := by
      exact sub_eq_zero.mp hzero
    exact_mod_cast (show
      failWitness v q p first * v second =
        failWitness v q p second * v first by
      simpa [mul_comm] using hcross.symm)
  · intro hslope
    rw [witnessSlope, witnessSlope,
      div_eq_div_iff hfirstQ hsecondQ] at hslope
    have hcross :
        failWitness v q p first * v second =
          failWitness v q p second * v first := by
      exact_mod_cast hslope
    exact sub_eq_zero.mpr (by simpa [mul_comm] using hcross.symm)

theorem overlapDet_zero_refl
    (v : Fin 13 → ℤ) (q p : ℕ) (index : Fin 13) :
    overlapDet v q p index index = 0 := by
  simp [overlapDet]

theorem overlapDet_zero_symm
    (v : Fin 13 → ℤ) (q p : ℕ) (first second : Fin 13)
    (hzero : overlapDet v q p first second = 0) :
    overlapDet v q p second first = 0 := by
  rw [overlapDet_skew, hzero, neg_zero]

/-- Transitivity is the Plücker identity with the middle speed cancelled. -/
theorem overlapDet_zero_trans
    (v : Fin 13 → ℤ) (q p : ℕ) (first middle last : Fin 13)
    (hmiddle : v middle ≠ 0)
    (hfirstMiddle : overlapDet v q p first middle = 0)
    (hmiddleLast : overlapDet v q p middle last = 0) :
    overlapDet v q p first last = 0 := by
  have hidentity :
      v middle * overlapDet v q p first last =
        v first * overlapDet v q p middle last +
          v last * overlapDet v q p first middle := by
    simp [overlapDet]
    ring
  rw [hfirstMiddle, hmiddleLast, mul_zero, mul_zero, add_zero] at hidentity
  exact (mul_eq_zero.mp hidentity).resolve_left hmiddle

/-- The zero-color relation on a nonzero speed family is a setoid. -/
def overlapZeroSetoid
    (v : Fin 13 → ℤ) (hv : ∀ index, v index ≠ 0) (q p : ℕ) :
    Setoid (Fin 13) where
  r first second := overlapDet v q p first second = 0
  iseqv := ⟨
    overlapDet_zero_refl v q p,
    fun hzero => overlapDet_zero_symm v q p _ _ hzero,
    fun hfirstMiddle hmiddleLast =>
      overlapDet_zero_trans v q p _ _ _ (hv _) hfirstMiddle hmiddleLast⟩

/-- A path of zero-color edges has zero-color endpoints.  Thus every connected
component of the zero-color graph is a clique. -/
theorem overlapDet_zero_of_reflTransGen
    (v : Fin 13 → ℤ) (hv : ∀ index, v index ≠ 0) (q p : ℕ)
    (first last : Fin 13)
    (hpath : Relation.ReflTransGen
      (fun left right => overlapDet v q p left right = 0) first last) :
    overlapDet v q p first last = 0 := by
  exact Relation.ReflTransGen.head_induction_on hpath
    (overlapDet_zero_refl v q p last)
    (fun hfirstMiddle _hpath hinduction =>
      overlapDet_zero_trans v q p _ _ last
        (hv _) hfirstMiddle hinduction)

/-- Pairwise zero color glues a nonempty stalk to one rational parameter. -/
theorem exists_common_witnessSlope_of_pairwise_zero
    (v : Fin 13 → ℤ) (hv : ∀ index, v index ≠ 0) (q p : ℕ)
    (stalk : Finset (Fin 13)) (hstalk : stalk.Nonempty)
    (hzero : Set.Pairwise (stalk : Set (Fin 13))
      fun first second => overlapDet v q p first second = 0) :
    ∃ slope : ℚ, ∀ index ∈ stalk,
      witnessSlope v q p index = slope := by
  obtain ⟨root, hroot⟩ := hstalk
  refine ⟨witnessSlope v q p root, ?_⟩
  intro index hindex
  by_cases hindexRoot : index = root
  · subst index
    rfl
  · apply (overlapDet_eq_zero_iff_witnessSlope_eq
      v q p index root (hv index) (hv root)).mp
    exact hzero hindex hroot hindexRoot

/-- Connected zero color suffices: transitivity first upgrades every root path
to a root edge, after which all vertices share the root slope. -/
theorem exists_common_witnessSlope_of_zero_connected
    (v : Fin 13 → ℤ) (hv : ∀ index, v index ≠ 0) (q p : ℕ)
    (stalk : Finset (Fin 13)) (root : Fin 13) (_hroot : root ∈ stalk)
    (hconnected : ∀ index ∈ stalk,
      Relation.ReflTransGen
        (fun left right => overlapDet v q p left right = 0) root index) :
    ∃ slope : ℚ, ∀ index ∈ stalk,
      witnessSlope v q p index = slope := by
  refine ⟨witnessSlope v q p root, ?_⟩
  intro index hindex
  have hzero := overlapDet_zero_of_reflTransGen
    v hv q p root index (hconnected index hindex)
  exact (overlapDet_eq_zero_iff_witnessSlope_eq
    v q p root index (hv root) (hv index)).mp hzero |>.symm

#print axioms overlapDet_eq_zero_iff_witnessSlope_eq
#print axioms overlapDet_zero_trans
#print axioms overlapZeroSetoid
#print axioms overlapDet_zero_of_reflTransGen
#print axioms exists_common_witnessSlope_of_pairwise_zero
#print axioms exists_common_witnessSlope_of_zero_connected

end LRC14Concrete
end LonelyRunner
