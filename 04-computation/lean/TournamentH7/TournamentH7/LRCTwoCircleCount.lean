/-
  TournamentH7.LRCTwoCircleCount — EXACT CANONICAL DEEP-MULTIPLIER COUNT.

  `LRCTwoCircleII.deep_iff_circles` identifies the canonical depth-six set
  with two endpoint circles and one half circle.  This module turns those
  three strict integer inequalities into literal natural-number intervals.

  With

      B = (q - 1) / 84,
      L = (q - B + 1) / 2,
      U = (q + B) / 2,

  the low, high, and half circles are respectively

      Icc 1 B,   Icc (q-B) (q-1),   Icc L U.

  The first two intervals each have `B` points.  All three are pairwise
  disjoint.  The half interval has

      q even:  1 + 2 * (B / 2),
      q odd:   2 * ((B + 1) / 2)

  points.  Hence their sum is the exact number of `p ∈ Ioo 0 q` with
  canonical `bandCount ≥ 6`.

  Kernel-pure: no `sorry`, no `native_decide`, no new axioms.
-/
import Mathlib
import TournamentH7.LRCTwoCircleII

namespace LonelyRunner
namespace LRC14Concrete

open Finset

/-- Radius, in integer multiplier steps, of either endpoint circle. -/
def circleBlock (q : ℕ) : ℕ := (q - 1) / 84

/-- First multiplier in the half circle. -/
def halfLower (q : ℕ) : ℕ := (q - circleBlock q + 1) / 2

/-- Last multiplier in the half circle. -/
def halfUpper (q : ℕ) : ℕ := (q + circleBlock q) / 2

def lowCircleMultipliers (q : ℕ) : Finset ℕ :=
  Finset.Icc 1 (circleBlock q)

def highCircleMultipliers (q : ℕ) : Finset ℕ :=
  Finset.Icc (q - circleBlock q) (q - 1)

def halfCircleMultipliers (q : ℕ) : Finset ℕ :=
  Finset.Icc (halfLower q) (halfUpper q)

def canonicalCircleMultipliers (q : ℕ) : Finset ℕ :=
  (lowCircleMultipliers q ∪ highCircleMultipliers q) ∪
    halfCircleMultipliers q

def canonicalDeepMultipliers (q : ℕ) : Finset ℕ :=
  (Finset.Ioo 0 q).filter fun p => 6 ≤ bandCount canonical q p

/-- The defining quotient gives the useful non-division inequality. -/
theorem circleBlock_mul_le (q : ℕ) : circleBlock q * 84 ≤ q - 1 := by
  exact (Nat.le_div_iff_mul_le (by norm_num : 0 < 84)).mp le_rfl

theorem circleBlock_lt_q (q : ℕ) (hq : 0 < q) : circleBlock q < q := by
  have h := circleBlock_mul_le q
  omega

theorem two_circleBlock_lt_q (q : ℕ) (hq : 0 < q) :
    2 * circleBlock q < q := by
  have h := circleBlock_mul_le q
  omega

theorem three_circleBlock_lt_q (q : ℕ) (hq : 0 < q) :
    3 * circleBlock q < q := by
  have h := circleBlock_mul_le q
  omega

/-- Strict `q/84` proximity is exactly membership in the integer block. -/
theorem eightyFour_mul_lt_iff_le_circleBlock (q x : ℕ) (hq : 0 < q) :
    (84 : ℤ) * (x : ℤ) < q ↔ x ≤ circleBlock q := by
  have hdiv : x ≤ circleBlock q ↔ x * 84 ≤ q - 1 := by
    exact Nat.le_div_iff_mul_le (by norm_num : 0 < 84)
  rw [hdiv]
  constructor
  · intro h
    exact_mod_cast (show (x : ℤ) * 84 ≤ (q - 1 : ℕ) by omega)
  · intro h
    exact_mod_cast (show 84 * (x : ℕ) < q by omega)

/-- The low endpoint circle is the literal interval `Icc 1 B`. -/
theorem mem_lowCircleMultipliers_iff (q p : ℕ) (hq : 0 < q) :
    p ∈ lowCircleMultipliers q ↔
      0 < p ∧ p < q ∧ (84 : ℤ) * (p : ℤ) < q := by
  rw [lowCircleMultipliers, Finset.mem_Icc]
  constructor
  · rintro ⟨hp1, hpB⟩
    refine ⟨by omega, hpB.trans_lt (circleBlock_lt_q q hq), ?_⟩
    exact (eightyFour_mul_lt_iff_le_circleBlock q p hq).mpr hpB
  · rintro ⟨hp0, _hpq, hnear⟩
    exact ⟨by omega, (eightyFour_mul_lt_iff_le_circleBlock q p hq).mp hnear⟩

/-- The high endpoint circle is the reflected interval `Icc (q-B) (q-1)`. -/
theorem mem_highCircleMultipliers_iff (q p : ℕ) (hq : 0 < q) :
    p ∈ highCircleMultipliers q ↔
      0 < p ∧ p < q ∧ (84 : ℤ) * ((q : ℤ) - p) < q := by
  rw [highCircleMultipliers, Finset.mem_Icc]
  constructor
  · rintro ⟨hlow, hupp⟩
    have hpq : p < q := by omega
    have hdiff : q - p ≤ circleBlock q := by omega
    have hnear :=
      (eightyFour_mul_lt_iff_le_circleBlock q (q - p) hq).mpr hdiff
    have hcast : ((q - p : ℕ) : ℤ) = (q : ℤ) - p := by omega
    rw [hcast] at hnear
    exact ⟨by omega, hpq, hnear⟩
  · rintro ⟨hp0, hpq, hnear⟩
    have hcast : ((q - p : ℕ) : ℤ) = (q : ℤ) - p := by omega
    have hnear' : (84 : ℤ) * ((q - p : ℕ) : ℤ) < q := by
      rwa [hcast]
    have hdiff :=
      (eightyFour_mul_lt_iff_le_circleBlock q (q - p) hq).mp hnear'
    exact ⟨by omega, by omega⟩

/-- The half-resonance inequality is the literal interval `Icc L U`. -/
theorem mem_halfCircleMultipliers_iff (q p : ℕ) (hq : 0 < q) :
    p ∈ halfCircleMultipliers q ↔
      0 < p ∧ p < q ∧
        84 * |2 * (p : ℤ) - q| < (q : ℤ) := by
  rw [halfCircleMultipliers, Finset.mem_Icc]
  have hBq := circleBlock_lt_q q hq
  constructor
  · rintro ⟨hlo, hup⟩
    have hleft : (q : ℤ) - circleBlock q ≤ 2 * (p : ℤ) := by
      dsimp [halfLower] at hlo
      omega
    have hright : 2 * (p : ℤ) ≤ (q : ℤ) + circleBlock q := by
      dsimp [halfUpper] at hup
      omega
    have habs : |2 * (p : ℤ) - q| ≤ (circleBlock q : ℤ) :=
      (abs_le.mpr (by constructor <;> linarith))
    have hp0 : 0 < p := by
      dsimp [halfLower] at hlo
      omega
    have hpq : p < q := by
      dsimp [halfUpper] at hup
      omega
    have hnatAbs : Int.natAbs (2 * (p : ℤ) - q) ≤ circleBlock q := by
      have hcast : ((Int.natAbs (2 * (p : ℤ) - q) : ℕ) : ℤ) ≤
          (circleBlock q : ℤ) := by
        simpa only [Int.natCast_natAbs] using habs
      exact_mod_cast hcast
    have hnear := (eightyFour_mul_lt_iff_le_circleBlock q
      (Int.natAbs (2 * (p : ℤ) - q)) hq).mpr hnatAbs
    refine ⟨hp0, hpq, ?_⟩
    simpa only [Int.natCast_natAbs] using hnear
  · rintro ⟨hp0, hpq, hnear⟩
    have hbound : |2 * (p : ℤ) - q| ≤ (circleBlock q : ℤ) := by
      have hnat : Int.natAbs (2 * (p : ℤ) - q) ≤ circleBlock q := by
        apply (eightyFour_mul_lt_iff_le_circleBlock q
          (Int.natAbs (2 * (p : ℤ) - q)) hq).mp
        rw [Int.natCast_natAbs]
        exact hnear
      have hcast : ((Int.natAbs (2 * (p : ℤ) - q) : ℕ) : ℤ) ≤
          (circleBlock q : ℤ) := by
        exact_mod_cast hnat
      simpa only [Int.natCast_natAbs] using hcast
    have hsides := abs_le.mp hbound
    constructor
    · dsimp [halfLower]
      omega
    · dsimp [halfUpper]
      omega

theorem lowCircleMultipliers_card (q : ℕ) :
    (lowCircleMultipliers q).card = circleBlock q := by
  rw [lowCircleMultipliers, Nat.card_Icc]
  omega

theorem highCircleMultipliers_card (q : ℕ) (hq : 0 < q) :
    (highCircleMultipliers q).card = circleBlock q := by
  rw [highCircleMultipliers, Nat.card_Icc]
  have hBq := circleBlock_lt_q q hq
  omega

/-- Parity form of the half-circle count. -/
theorem halfCircleMultipliers_card (q : ℕ) (hq : 0 < q) :
    (halfCircleMultipliers q).card =
      if q % 2 = 0 then 1 + 2 * (circleBlock q / 2)
      else 2 * ((circleBlock q + 1) / 2) := by
  rw [halfCircleMultipliers, Nat.card_Icc]
  have hBq := circleBlock_lt_q q hq
  rcases Nat.mod_two_eq_zero_or_one q with hqEven | hqOdd
  · rw [if_pos hqEven]
    dsimp [halfLower, halfUpper]
    rcases Nat.mod_two_eq_zero_or_one (circleBlock q) with hBEven | hBOdd <;>
      omega
  · rw [if_neg (by omega)]
    dsimp [halfLower, halfUpper]
    rcases Nat.mod_two_eq_zero_or_one (circleBlock q) with hBEven | hBOdd <;>
      omega

/-- Compact parity form used by the exact-count ledger. -/
theorem halfCircleMultipliers_card_compact (q : ℕ) (hq : 0 < q) :
    (halfCircleMultipliers q).card =
      circleBlock q + 1 - ((circleBlock q + q) % 2) := by
  rw [halfCircleMultipliers_card q hq]
  rcases Nat.mod_two_eq_zero_or_one q with hqEven | hqOdd
  · rw [if_pos hqEven]
    rcases Nat.mod_two_eq_zero_or_one (circleBlock q) with hBEven | hBOdd <;>
      omega
  · rw [if_neg (by omega)]
    rcases Nat.mod_two_eq_zero_or_one (circleBlock q) with hBEven | hBOdd <;>
      omega

theorem lowCircle_disjoint_highCircle (q : ℕ) (hq : 0 < q) :
    Disjoint (lowCircleMultipliers q) (highCircleMultipliers q) := by
  rw [Finset.disjoint_left]
  intro p hpLow hpHigh
  rw [lowCircleMultipliers, Finset.mem_Icc] at hpLow
  rw [highCircleMultipliers, Finset.mem_Icc] at hpHigh
  have htwo := two_circleBlock_lt_q q hq
  omega

theorem lowCircle_disjoint_halfCircle (q : ℕ) (hq : 0 < q) :
    Disjoint (lowCircleMultipliers q) (halfCircleMultipliers q) := by
  rw [Finset.disjoint_left]
  intro p hpLow hpHalf
  rw [lowCircleMultipliers, Finset.mem_Icc] at hpLow
  rw [halfCircleMultipliers, Finset.mem_Icc] at hpHalf
  have hthree := three_circleBlock_lt_q q hq
  dsimp [halfLower] at hpHalf
  omega

theorem highCircle_disjoint_halfCircle (q : ℕ) (hq : 0 < q) :
    Disjoint (highCircleMultipliers q) (halfCircleMultipliers q) := by
  rw [Finset.disjoint_left]
  intro p hpHigh hpHalf
  rw [highCircleMultipliers, Finset.mem_Icc] at hpHigh
  rw [halfCircleMultipliers, Finset.mem_Icc] at hpHalf
  have hthree := three_circleBlock_lt_q q hq
  dsimp [halfUpper] at hpHalf
  omega

theorem endpointCircles_disjoint_halfCircle (q : ℕ) (hq : 0 < q) :
    Disjoint (lowCircleMultipliers q ∪ highCircleMultipliers q)
      (halfCircleMultipliers q) := by
  rw [Finset.disjoint_union_left]
  exact ⟨lowCircle_disjoint_halfCircle q hq,
    highCircle_disjoint_halfCircle q hq⟩

/-- The in-kernel two-circle theorem as an equality of finite multiplier sets. -/
theorem canonicalDeepMultipliers_eq_circleMultipliers (q : ℕ) (hq : 0 < q) :
    canonicalDeepMultipliers q = canonicalCircleMultipliers q := by
  ext p
  rw [canonicalDeepMultipliers, Finset.mem_filter, Finset.mem_Ioo]
  rw [canonicalCircleMultipliers, Finset.mem_union, Finset.mem_union]
  rw [mem_lowCircleMultipliers_iff q p hq,
    mem_highCircleMultipliers_iff q p hq,
    mem_halfCircleMultipliers_iff q p hq]
  constructor
  · rintro ⟨⟨hp, hpq⟩, hdeep⟩
    rcases (deep_iff_circles q p hq hp hpq).mp hdeep with hendpoint | hhalf
    · rcases hendpoint with hlow | hhigh
      · exact Or.inl (Or.inl ⟨hp, hpq, hlow⟩)
      · exact Or.inl (Or.inr ⟨hp, hpq, hhigh⟩)
    · exact Or.inr ⟨hp, hpq, hhalf⟩
  · intro hcircles
    rcases hcircles with hendpoint | hhalf
    · rcases hendpoint with hlow | hhigh
      · exact ⟨⟨hlow.1, hlow.2.1⟩,
          (deep_iff_circles q p hq hlow.1 hlow.2.1).mpr
            (Or.inl (Or.inl hlow.2.2))⟩
      · exact ⟨⟨hhigh.1, hhigh.2.1⟩,
          (deep_iff_circles q p hq hhigh.1 hhigh.2.1).mpr
            (Or.inl (Or.inr hhigh.2.2))⟩
    · exact ⟨⟨hhalf.1, hhalf.2.1⟩,
        (deep_iff_circles q p hq hhalf.1 hhalf.2.1).mpr
          (Or.inr hhalf.2.2)⟩

/-- Exact cardinality of the canonical deep multiplier set. -/
theorem canonicalDeepMultipliers_card (q : ℕ) (hq : 0 < q) :
    (canonicalDeepMultipliers q).card =
      2 * circleBlock q +
        if q % 2 = 0 then 1 + 2 * (circleBlock q / 2)
        else 2 * ((circleBlock q + 1) / 2) := by
  rw [canonicalDeepMultipliers_eq_circleMultipliers q hq]
  rw [canonicalCircleMultipliers]
  rw [Finset.card_union_of_disjoint
    (endpointCircles_disjoint_halfCircle q hq)]
  rw [Finset.card_union_of_disjoint (lowCircle_disjoint_highCircle q hq)]
  rw [lowCircleMultipliers_card, highCircleMultipliers_card q hq,
    halfCircleMultipliers_card q hq]
  omega

/-- Exact count in the compact `2B + (B+1-(B+q)%2)` form. -/
theorem canonicalDeepMultipliers_card_compact (q : ℕ) (hq : 0 < q) :
    (canonicalDeepMultipliers q).card =
      2 * circleBlock q +
        (circleBlock q + 1 - ((circleBlock q + q) % 2)) := by
  rw [canonicalDeepMultipliers_card q hq]
  rw [← halfCircleMultipliers_card q hq]
  rw [halfCircleMultipliers_card_compact q hq]

/-- At the LRC-14 target modulus, the unique deep multiplier is the half point. -/
theorem canonicalDeepMultipliers_fourteen :
    canonicalDeepMultipliers 14 = {7} := by
  rw [canonicalDeepMultipliers_eq_circleMultipliers 14 (by norm_num)]
  decide

theorem canonicalDeepMultipliers_card_fourteen :
    (canonicalDeepMultipliers 14).card = 1 := by
  rw [canonicalDeepMultipliers_fourteen]
  simp

/-! ## Axiom audit -/
#print axioms mem_lowCircleMultipliers_iff
#print axioms mem_highCircleMultipliers_iff
#print axioms mem_halfCircleMultipliers_iff
#print axioms halfCircleMultipliers_card
#print axioms halfCircleMultipliers_card_compact
#print axioms endpointCircles_disjoint_halfCircle
#print axioms canonicalDeepMultipliers_eq_circleMultipliers
#print axioms canonicalDeepMultipliers_card
#print axioms canonicalDeepMultipliers_card_compact
#print axioms canonicalDeepMultipliers_fourteen
#print axioms canonicalDeepMultipliers_card_fourteen

end LRC14Concrete
end LonelyRunner
