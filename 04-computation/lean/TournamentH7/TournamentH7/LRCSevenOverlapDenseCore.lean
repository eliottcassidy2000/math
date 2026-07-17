/-
  Dense-core consequences of the colored seven-overlap Plucker relation.
  High nonzero base colors force explicit determinant mass on incident edges;
  two unit spokes force a zero-colored base. No `sorry`; no `native_decide`.
-/

import TournamentH7.LRCSevenOverlapRelations
import TournamentH7.LRCDenseCoreRelationTrap

namespace LonelyRunner
namespace LRC14Grand

open LonelyRunner
open LRC14Concrete
open scoped Classical

/-- Plucker coefficients transported to the absolute-speed ordering. -/
def signedPermutedOverlapCoeff
    (v : Fin 13 → ℤ) (q p : ℕ) (σ : Equiv.Perm (Fin 13))
    (i j k index : Fin 13) : ℤ :=
  Int.sign (v (σ index)) *
    overlapTripleCoeff v q p (σ i) (σ j) (σ k) (σ index)

theorem signedPermutedOverlapCoeff_sum_abs_eq_zero
    (v : Fin 13 → ℤ) (q p : ℕ) (σ : Equiv.Perm (Fin 13))
    (i j k : Fin 13) :
    ∑ index, signedPermutedOverlapCoeff v q p σ i j k index * |v (σ index)| = 0 := by
  have h := overlapTripleCoeff_sum_eq_zero v q p (σ i) (σ j) (σ k)
  rw [← Equiv.sum_comp σ] at h
  calc
    ∑ index, signedPermutedOverlapCoeff v q p σ i j k index * |v (σ index)| =
        ∑ index, overlapTripleCoeff v q p (σ i) (σ j) (σ k) (σ index) *
          v (σ index) := by
      apply Finset.sum_congr rfl
      intro index _
      rw [signedPermutedOverlapCoeff]
      calc
        Int.sign (v (σ index)) *
              overlapTripleCoeff v q p (σ i) (σ j) (σ k) (σ index) *
              |v (σ index)| =
            overlapTripleCoeff v q p (σ i) (σ j) (σ k) (σ index) *
              (Int.sign (v (σ index)) * |v (σ index)|) := by ring
        _ = overlapTripleCoeff v q p (σ i) (σ j) (σ k) (σ index) *
              v (σ index) := by rw [Int.sign_mul_abs]
    _ = 0 := h

/-- Any Plucker vector whose top lies above the dense pair spends at least
three units of coefficient mass below that top. -/
theorem signedPermutedOverlapCoeff_below_mass_ge_three
    (v : Fin 13 → ℤ) (hv : ∀ index, v index ≠ 0)
    (q p : ℕ) (σ : Equiv.Perm (Fin 13))
    (hmono : Monotone (fun index => |v (σ index)|))
    (dense : Fin 12)
    (hladder : ∀ r : Fin 12, dense < r →
      3 * |v (σ r.castSucc)| ≤ |v (σ r.succ)|)
    (i j k top : Fin 13)
    (htopIndex : (dense : ℕ) + 2 ≤ (top : ℕ))
    (htop : signedPermutedOverlapCoeff v q p σ i j k top ≠ 0)
    (hhigh : ∀ index : Fin 13, (top : ℕ) < (index : ℕ) →
      signedPermutedOverlapCoeff v q p σ i j k index = 0) :
    3 ≤ ∑ index ∈ Finset.univ.filter
      (fun index : Fin 13 => (index : ℕ) < (top : ℕ)),
      |signedPermutedOverlapCoeff v q p σ i j k index| := by
  let w : Fin 13 → ℤ := fun index => |v (σ index)|
  let coeff : Fin 13 → ℤ := signedPermutedOverlapCoeff v q p σ i j k
  have hpos : ∀ index, 0 < w index := by
    intro index
    exact abs_pos.mpr (hv (σ index))
  have hrelation : ∑ index, coeff index * w index = 0 := by
    exact signedPermutedOverlapCoeff_sum_abs_eq_zero v q p σ i j k
  by_contra hnot
  push Not at hnot
  have hmassRaw : (∑ index ∈ Finset.univ.filter
      (fun index : Fin 13 => (index : ℕ) < (top : ℕ)),
      |signedPermutedOverlapCoeff v q p σ i j k index|) ≤ 2 := by omega
  have hmass : (∑ index ∈ Finset.univ.filter
      (fun index : Fin 13 => (index : ℕ) < (top : ℕ)), |coeff index|) ≤ 2 := by
    simpa [coeff] using hmassRaw
  exact (no_low_mass_relation_above_pair w hpos hmono dense hladder
    coeff top htopIndex htop hhigh hmass) hrelation

/-- Triangle specialization: a nonzero color below a high vertex forces at
least three units of determinant mass on the two incident colors. -/
theorem overlapDet_incident_mass_ge_three_of_high_base_ne_zero
    (v : Fin 13 → ℤ) (hv : ∀ index, v index ≠ 0)
    (q p : ℕ) (σ : Equiv.Perm (Fin 13))
    (hmono : Monotone (fun index => |v (σ index)|))
    (dense : Fin 12)
    (hladder : ∀ r : Fin 12, dense < r →
      3 * |v (σ r.castSucc)| ≤ |v (σ r.succ)|)
    (left right top : Fin 13)
    (hleft : (left : ℕ) < (top : ℕ))
    (hright : (right : ℕ) < (top : ℕ))
    (htopIndex : (dense : ℕ) + 2 ≤ (top : ℕ))
    (hbase : overlapDet v q p (σ left) (σ right) ≠ 0) :
    3 ≤ |overlapDet v q p (σ right) (σ top)| +
      |overlapDet v q p (σ left) (σ top)| := by
  have htop13 : (top : ℕ) < 13 := top.isLt
  let previous : Fin 13 := ⟨(top : ℕ) - 1, by omega⟩
  have hpreviousPos : 0 < |v (σ previous)| := abs_pos.mpr (hv (σ previous))
  have hstep : 3 * |v (σ previous)| ≤ |v (σ top)| := by
    have h := ladder_top_step (fun index => |v (σ index)|) dense hladder
      (top : ℕ) htopIndex htop13
    have htopEq : (⟨(top : ℕ), htop13⟩ : Fin 13) = top := Fin.ext rfl
    simpa [previous, htopEq] using h
  have hleftLe : |v (σ left)| ≤ |v (σ previous)| := by
    apply hmono
    show (left : ℕ) ≤ (top : ℕ) - 1
    omega
  have hrightLe : |v (σ right)| ≤ |v (σ previous)| := by
    apply hmono
    show (right : ℕ) ≤ (top : ℕ) - 1
    omega
  have hrelation := overlapDet_plucker v q p (σ left) (σ right) (σ top)
  have heq :
      overlapDet v q p (σ left) (σ right) * v (σ top) =
        -(overlapDet v q p (σ right) (σ top) * v (σ left) -
          overlapDet v q p (σ left) (σ top) * v (σ right)) := by
    linarith
  have htriangle :
      |overlapDet v q p (σ left) (σ right)| * |v (σ top)| ≤
        |overlapDet v q p (σ right) (σ top)| * |v (σ left)| +
          |overlapDet v q p (σ left) (σ top)| * |v (σ right)| := by
    calc
      |overlapDet v q p (σ left) (σ right)| * |v (σ top)| =
          |overlapDet v q p (σ left) (σ right) * v (σ top)| := by
            rw [abs_mul]
      _ = |overlapDet v q p (σ right) (σ top) * v (σ left) -
            overlapDet v q p (σ left) (σ top) * v (σ right)| := by
            rw [heq, abs_neg]
      _ ≤ |overlapDet v q p (σ right) (σ top) * v (σ left)| +
            |overlapDet v q p (σ left) (σ top) * v (σ right)| :=
          abs_sub _ _
      _ = |overlapDet v q p (σ right) (σ top)| * |v (σ left)| +
            |overlapDet v q p (σ left) (σ top)| * |v (σ right)| := by
          rw [abs_mul, abs_mul]
  have hbaseOne : (1 : ℤ) ≤ |overlapDet v q p (σ left) (σ right)| :=
    Int.one_le_abs hbase
  have htopPos : 0 < |v (σ top)| := abs_pos.mpr (hv (σ top))
  have hincidentNonnegA : 0 ≤ |overlapDet v q p (σ right) (σ top)| :=
    abs_nonneg _
  have hincidentNonnegB : 0 ≤ |overlapDet v q p (σ left) (σ top)| :=
    abs_nonneg _
  nlinarith

/-- A high vertex cannot have two unit-color incident edges above a nonzero
base edge.  Equivalently, its unit-color neighborhood is a zero-color clique. -/
theorem overlapDet_base_eq_zero_of_two_incident_unit
    (v : Fin 13 → ℤ) (hv : ∀ index, v index ≠ 0)
    (q p : ℕ) (σ : Equiv.Perm (Fin 13))
    (hmono : Monotone (fun index => |v (σ index)|))
    (dense : Fin 12)
    (hladder : ∀ r : Fin 12, dense < r →
      3 * |v (σ r.castSucc)| ≤ |v (σ r.succ)|)
    (left right top : Fin 13)
    (hleft : (left : ℕ) < (top : ℕ))
    (hright : (right : ℕ) < (top : ℕ))
    (htopIndex : (dense : ℕ) + 2 ≤ (top : ℕ))
    (hleftUnit : |overlapDet v q p (σ left) (σ top)| ≤ 1)
    (hrightUnit : |overlapDet v q p (σ right) (σ top)| ≤ 1) :
    overlapDet v q p (σ left) (σ right) = 0 := by
  by_contra hbase
  have hmass := overlapDet_incident_mass_ge_three_of_high_base_ne_zero
    v hv q p σ hmono dense hladder left right top hleft hright htopIndex hbase
  omega

/-- A nonzero-color triangle below one high vertex costs at least five units
of determinant mass on the three spokes to that vertex. -/
theorem overlapDet_three_spoke_mass_ge_five_of_high_nonzero_triangle
    (v : Fin 13 → ℤ) (hv : ∀ index, v index ≠ 0)
    (q p : ℕ) (σ : Equiv.Perm (Fin 13))
    (hmono : Monotone (fun index => |v (σ index)|))
    (dense : Fin 12)
    (hladder : ∀ r : Fin 12, dense < r →
      3 * |v (σ r.castSucc)| ≤ |v (σ r.succ)|)
    (a b c top : Fin 13)
    (ha : (a : ℕ) < (top : ℕ))
    (hb : (b : ℕ) < (top : ℕ))
    (hc : (c : ℕ) < (top : ℕ))
    (htopIndex : (dense : ℕ) + 2 ≤ (top : ℕ))
    (hab : overlapDet v q p (σ a) (σ b) ≠ 0)
    (hac : overlapDet v q p (σ a) (σ c) ≠ 0)
    (hbc : overlapDet v q p (σ b) (σ c) ≠ 0) :
    5 ≤ |overlapDet v q p (σ a) (σ top)| +
      |overlapDet v q p (σ b) (σ top)| +
      |overlapDet v q p (σ c) (σ top)| := by
  have hAB := overlapDet_incident_mass_ge_three_of_high_base_ne_zero
    v hv q p σ hmono dense hladder a b top ha hb htopIndex hab
  have hAC := overlapDet_incident_mass_ge_three_of_high_base_ne_zero
    v hv q p σ hmono dense hladder a c top ha hc htopIndex hac
  have hBC := overlapDet_incident_mass_ge_three_of_high_base_ne_zero
    v hv q p σ hmono dense hladder b c top hb hc htopIndex hbc
  omega

/-! ## Axiom audit -/

#print axioms signedPermutedOverlapCoeff_sum_abs_eq_zero
#print axioms signedPermutedOverlapCoeff_below_mass_ge_three
#print axioms overlapDet_incident_mass_ge_three_of_high_base_ne_zero
#print axioms overlapDet_base_eq_zero_of_two_incident_unit
#print axioms overlapDet_three_spoke_mass_ge_five_of_high_nonzero_triangle

end LRC14Grand
end LonelyRunner
