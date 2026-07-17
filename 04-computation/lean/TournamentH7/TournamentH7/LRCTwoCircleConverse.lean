/-
  TournamentH7.LRCTwoCircleConverse — THE TWO-CIRCLE DEEP CERTIFICATE,
  CONVERSE PART II, THE SPEED-ONE CASE.

  `LRCTwoCircle` proves that both resonance circles lie in the canonical
  depth-six set.  This module starts the converse case split and closes its
  first genuine branch: if the canonical phase is deep and speed one fails,
  then the phase lies on the integer circle.

  The proof keeps the witnesses.  Among six distinct failing speeds one has
  speed at least six.  The coefficient-weight-14 relation lock forces that
  runner's witness to be the corresponding multiple of the speed-one
  witness.  The latter is either zero or one, and the large runner's strict
  band inequality sharpens the speed-one `1/14` neighborhood to the required
  `1/84` resonance circle.

  This is exactly the `k₀ = 1` branch of the structured Part II plan.  It does
  not claim the remaining `k₀ = 2,...,13` cases or the full converse.

  Kernel-pure: no `sorry`, no `native_decide`, no new axioms.
-/
import Mathlib
import TournamentH7.LRCTwoCircle

namespace LonelyRunner
namespace LRC14Concrete

open Finset

/-- If speed one and one speed `k >= 6` both fail, relation locking places
their witnesses on the same integer ray and the phase lies on Circle I. -/
theorem circleI_of_one_and_large_witness
    (q p k : ℕ) (wOne wLarge : ℤ)
    (hq : 0 < q) (hp : 0 < p) (hpq : p < q)
    (hk6 : 6 ≤ k) (hk13 : k ≤ 13)
    (hOne : 14 * |(p : ℤ) - wOne * q| < q)
    (hLarge : 14 * |(k : ℤ) * (p : ℤ) - wLarge * q| < q) :
    84 * (p : ℤ) < q ∨ 84 * ((q : ℤ) - p) < q := by
  have hk1 : 1 ≤ k := by omega
  have hlock : (k : ℤ) * wOne = wLarge := by
    have h := rational_lock_weight14
      1 (k : ℤ) wOne wLarge 1 k q p
      (by omega) hk1 (by omega) (by norm_num) (by simpa using hOne) hLarge
    simpa using h
  have hqZ : (0 : ℤ) < (q : ℤ) := by exact_mod_cast hq
  have hpZ : (0 : ℤ) < (p : ℤ) := by exact_mod_cast hp
  have hpqZ : (p : ℤ) < q := by exact_mod_cast hpq
  have hwLower : 0 ≤ wOne := by
    by_contra hneg
    have hw : wOne ≤ -1 := by omega
    have hdist : (q : ℤ) ≤ |(p : ℤ) - wOne * q| := by
      have hpos : (0 : ℤ) < (p : ℤ) - wOne * q := by nlinarith
      rw [abs_of_pos hpos]
      nlinarith
    nlinarith [hOne]
  have hwUpper : wOne ≤ 1 := by
    by_contra hhigh
    have hw : 2 ≤ wOne := by omega
    have hdist : (q : ℤ) ≤ |(p : ℤ) - wOne * q| := by
      have hneg : (p : ℤ) - wOne * q < 0 := by nlinarith
      rw [abs_of_neg hneg]
      nlinarith
    nlinarith [hOne]
  rcases (show wOne = 0 ∨ wOne = 1 by omega) with rfl | rfl
  · left
    have hwLarge : wLarge = 0 := by simpa using hlock.symm
    rw [hwLarge, zero_mul, sub_zero, abs_of_pos (mul_pos (by positivity) hpZ)] at hLarge
    have hkp : 84 * (p : ℤ) ≤ 14 * (k : ℤ) * (p : ℤ) := by nlinarith
    exact lt_of_le_of_lt hkp hLarge
  · right
    have hwLarge : wLarge = (k : ℤ) := by simpa using hlock.symm
    rw [hwLarge] at hLarge
    have hrewrite : (k : ℤ) * (p : ℤ) - (k : ℤ) * (q : ℤ)
        = -((k : ℤ) * ((q : ℤ) - p)) := by ring
    rw [hrewrite, abs_neg, abs_of_pos (mul_pos (by positivity) (by linarith))] at hLarge
    have hkqp : 84 * ((q : ℤ) - p)
        ≤ 14 * (k : ℤ) * ((q : ℤ) - p) := by nlinarith
    exact lt_of_le_of_lt hkqp hLarge

/-- **Part II, `k₀ = 1`.**  A canonical depth-six phase at which speed one
fails lies on the integer resonance circle. -/
theorem circleI_of_deep_and_speedOne_bad
    (q p : ℕ) (hq : 0 < q) (hp : 0 < p) (hpq : p < q)
    (hdeep : 6 ≤ bandCount canonical q p)
    (hOneBad : ¬ inBand canonical q p (0 : Fin 13)) :
    84 * (p : ℤ) < q ∨ 84 * ((q : ℤ) - p) < q := by
  let bad : Finset (Fin 13) :=
    Finset.univ.filter fun i : Fin 13 => ¬ inBand canonical q p i
  have hbadCard : 6 ≤ bad.card := by
    simpa [bad, bandCount] using hdeep
  have hlarge : ∃ i : Fin 13, i ∈ bad ∧ 5 ≤ (i : ℕ) := by
    by_contra hnone
    push_neg at hnone
    have hsub : bad ⊆ Finset.univ.filter fun i : Fin 13 => (i : ℕ) < 5 := by
      intro i hi
      rw [Finset.mem_filter]
      exact ⟨Finset.mem_univ i, by omega⟩
    have hsmallCard :
        (Finset.univ.filter fun i : Fin 13 => (i : ℕ) < 5).card = 5 := by
      decide
    have hcardLe := Finset.card_le_card hsub
    rw [hsmallCard] at hcardLe
    omega
  obtain ⟨i, hiBad, hi5⟩ := hlarge
  have hiBad' : ¬ inBand canonical q p i := by
    exact (Finset.mem_filter.mp hiBad).2
  let k : ℕ := (i : ℕ) + 1
  have hk6 : 6 ≤ k := by omega
  have hk13 : k ≤ 13 := by
    exact Nat.succ_le_iff.mpr i.isLt
  have hOne := bad_at_witness canonical q p (0 : Fin 13) hq hOneBad
  have hLarge := bad_at_witness canonical q p i hq hiBad'
  have hOne' :
      14 * |(p : ℤ) - failWitness canonical q p (0 : Fin 13) * q| < q := by
    simpa [canonical] using hOne
  have hLarge' :
      14 * |(k : ℤ) * (p : ℤ) - failWitness canonical q p i * q| < q := by
    simpa [canonical, k] using hLarge
  exact circleI_of_one_and_large_witness q p k
    (failWitness canonical q p (0 : Fin 13))
    (failWitness canonical q p i) hq hp hpq hk6 hk13 hOne' hLarge'

/-- A depth-six canonical phase has a failing speed among `1,...,8`.  This is
the finite endpoint of the Part II smallest-failing-speed split: `k₀ >= 9`
cannot occur because only five canonical speeds remain. -/
theorem exists_bad_speed_le_eight_of_deep
    (q p : ℕ) (hdeep : 6 ≤ bandCount canonical q p) :
    ∃ i : Fin 13, ¬ inBand canonical q p i ∧ (i : ℕ) < 8 := by
  let bad : Finset (Fin 13) :=
    Finset.univ.filter fun i : Fin 13 => ¬ inBand canonical q p i
  have hbadCard : 6 ≤ bad.card := by
    simpa [bad, bandCount] using hdeep
  by_contra hnone
  push_neg at hnone
  have hsub : bad ⊆ Finset.univ.filter fun i : Fin 13 => 8 ≤ (i : ℕ) := by
    intro i hi
    rw [Finset.mem_filter]
    exact ⟨Finset.mem_univ i, hnone i (Finset.mem_filter.mp hi).2⟩
  have hhighCard :
      (Finset.univ.filter fun i : Fin 13 => 8 ≤ (i : ℕ)).card = 5 := by
    decide
  have hcardLe := Finset.card_le_card hsub
  rw [hhighCard] at hcardLe
  omega

/-! ## Axiom audit -/
#print axioms circleI_of_one_and_large_witness
#print axioms circleI_of_deep_and_speedOne_bad
#print axioms exists_bad_speed_le_eight_of_deep

end LRC14Concrete
end LonelyRunner
