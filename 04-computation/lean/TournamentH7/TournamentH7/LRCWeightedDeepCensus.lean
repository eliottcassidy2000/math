/-
  TournamentH7.LRCWeightedDeepCensus

  The exact THM-950 census ledger.  A multiplier of depth `c` does not in
  general cost the worst-case constant `792`; its precise depth-five
  Bonferroni debt is `choose (c-1) 5`.  Consequently

    B5 = liveCount - sum_p choose (bandCount(p)-1) 5.

  This makes the intended seven-overlap bridge visible.  Rooting a seven-set
  at one bad runner gives `choose (c-1) 6` events, and for every `8 ≤ c ≤ 13`
  the exact debt is at most three times this count.  Thus the three-unit
  colored-spoke charge can in principle pay all depth-at-least-eight debt;
  depths six and seven, and rank-one aligned events, are the honest residue.

  No `sorry`; no `native_decide`.
-/

import TournamentH7.LRCDeepCertificate

namespace LonelyRunner
namespace LRC14Concrete

open Finset

/-- Exact depth-five leverage cost over all sampled multipliers.  Depths zero
through five contribute zero automatically. -/
def weightedDeepCost (v : Fin 13 → ℤ) (q : ℕ) : ℕ :=
  ∑ p ∈ Finset.Ioo 0 q, (bandCount v q p - 1).choose 5

/-- Pointwise depth-five Bonferroni identity on the only possible coverage
values `0,...,13`. -/
theorem bled_eq_live_indicator_sub_weightedCost :
    ∀ c ∈ Finset.range 14,
      (∑ d ∈ range 6, (-1 : ℤ) ^ d * Nat.choose c d) =
        (if c = 0 then 1 else 0) - (c - 1).choose 5 := by
  decide

/-- **Exact weighted census identity.**  This is the sharp form behind
THM-950's uniform `792` lower bound. -/
theorem B5_eq_live_sub_weightedDeepCost (v : Fin 13 → ℤ) (q : ℕ) :
    B5 v q = (liveCount v q : ℤ) - weightedDeepCost v q := by
  have hswap : B5 v q =
      ∑ p ∈ Finset.Ioo 0 q,
        ∑ d ∈ range 6, (-1 : ℤ) ^ d * Nat.choose (bandCount v q p) d := by
    unfold B5 momentS
    push_cast
    have hstep : ∀ d ∈ range 6,
        ((-1 : ℤ)) ^ d * (∑ p ∈ Finset.Ioo 0 q,
          (Nat.choose (bandCount v q p) d : ℤ)) =
        ∑ p ∈ Finset.Ioo 0 q,
          (-1 : ℤ) ^ d * (Nat.choose (bandCount v q p) d : ℤ) := by
      intro d _
      rw [Finset.mul_sum]
    rw [Finset.sum_congr rfl hstep, Finset.sum_comm]
  rw [hswap]
  calc
    (∑ p ∈ Finset.Ioo 0 q,
        ∑ d ∈ range 6, (-1 : ℤ) ^ d * Nat.choose (bandCount v q p) d) =
        ∑ p ∈ Finset.Ioo 0 q,
          ((if bandCount v q p = 0 then (1 : ℤ) else 0) -
            ((bandCount v q p - 1).choose 5 : ℤ)) := by
      apply Finset.sum_congr rfl
      intro p _hp
      exact bled_eq_live_indicator_sub_weightedCost
        (bandCount v q p)
        (Finset.mem_range.mpr (by
          have := bandCount_le_thirteen v q p
          omega))
    _ = (liveCount v q : ℤ) - weightedDeepCost v q := by
      rw [Finset.sum_sub_distrib]
      congr 1
      · rw [Finset.sum_boole]
        unfold liveCount
        norm_cast
      · unfold weightedDeepCost
        push_cast
        rfl

/-- The exact finite weighted race is equivalent to positive `B5`; unlike an
abstract certificate, neither direction loses information. -/
theorem weightedDeepCost_lt_liveCount_iff_B5_pos
    (v : Fin 13 → ℤ) (q : ℕ) :
    weightedDeepCost v q < liveCount v q ↔ 0 < B5 v q := by
  rw [B5_eq_live_sub_weightedDeepCost]
  omega

/-- End-to-end loneliness from the exact weighted census. -/
theorem lonely_of_exact_weighted_census
    (v : Fin 13 → ℤ) (q : ℕ) (hq : 0 < q)
    (hv : ∀ i, v i ≠ 0)
    (hrace : weightedDeepCost v q < liveCount v q) :
    ∃ t : ℝ, Lonely 14 v t := by
  exact lonely_of_Mreach_ge v hv
    (mreach_ge_of_B5_pos v q hq
      ((weightedDeepCost_lt_liveCount_iff_B5_pos v q).mp hrace))

/-- Number of sampled multipliers at one exact depth. -/
def exactDepthCount (v : Fin 13 → ℤ) (q depth : ℕ) : ℕ :=
  ((Finset.Ioo 0 q).filter fun p => bandCount v q p = depth).card

/-- Total number of rooted seven-subsets, counted with multiplicity over
multiplier events.  At depth `c` this is exactly `choose (c-1) 6`. -/
def rootedSevenActivity (v : Fin 13 → ℤ) (q : ℕ) : ℕ :=
  ∑ p ∈ Finset.Ioo 0 q, (bandCount v q p - 1).choose 6

/-- Pointwise form of the exact `6/7/colored` split. -/
theorem weightedDepthCost_le_six_threeSeven_threeRooted :
    ∀ c ∈ Finset.range 14,
      (c - 1).choose 5 ≤
        (if c = 6 then 1 else 0) +
          3 * (if c = 7 then 1 else 0) +
          3 * (c - 1).choose 6 := by
  decide

/-- **Global exact-residue bound.**  All weighted B5 debt is bounded by one
unit per depth-six event, three extra units per depth-seven event, and three
units per rooted seven-overlap event.  This is the precise consumer shape for
the activity-weighted colored graph. -/
theorem weightedDeepCost_le_six_threeSeven_threeRooted
    (v : Fin 13 → ℤ) (q : ℕ) :
    weightedDeepCost v q ≤
      exactDepthCount v q 6 + 3 * exactDepthCount v q 7 +
        3 * rootedSevenActivity v q := by
  unfold weightedDeepCost exactDepthCount rootedSevenActivity
  calc
    (∑ p ∈ Finset.Ioo 0 q, (bandCount v q p - 1).choose 5) ≤
        ∑ p ∈ Finset.Ioo 0 q,
          ((if bandCount v q p = 6 then 1 else 0) +
            3 * (if bandCount v q p = 7 then 1 else 0) +
            3 * (bandCount v q p - 1).choose 6) := by
      apply Finset.sum_le_sum
      intro p _hp
      exact weightedDepthCost_le_six_threeSeven_threeRooted
        (bandCount v q p)
        (Finset.mem_range.mpr (by
          have := bandCount_le_thirteen v q p
          omega))
    _ = ((Finset.Ioo 0 q).filter fun p => bandCount v q p = 6).card +
          3 * ((Finset.Ioo 0 q).filter fun p => bandCount v q p = 7).card +
          3 * ∑ p ∈ Finset.Ioo 0 q,
            (bandCount v q p - 1).choose 6 := by
      rw [Finset.sum_add_distrib, Finset.sum_add_distrib,
        ← Finset.mul_sum, ← Finset.mul_sum,
        Finset.sum_boole, Finset.sum_boole]
      norm_cast

/-- Exact pointwise costs for depths six through thirteen. -/
theorem weighted_depth_cost_table :
    ((Nat.choose 5 5, Nat.choose 6 5, Nat.choose 7 5, Nat.choose 8 5,
      Nat.choose 9 5, Nat.choose 10 5, Nat.choose 11 5, Nat.choose 12 5)) =
      (1, 6, 21, 56, 126, 252, 462, 792) := by
  norm_num [Nat.choose]

/-- **Rooted seven-set payment inequality.**  For all possible depths at least
eight, three units per rooted seven-overlap event cover the entire exact B5
debt.  Equality occurs at depth eight (`21 = 3*7`). -/
theorem weightedDepthCost_le_three_mul_rootedSeven
    (c : ℕ) (hc8 : 8 ≤ c) (hc13 : c ≤ 13) :
    (c - 1).choose 5 ≤ 3 * (c - 1).choose 6 := by
  interval_cases c <;> decide

/-- The sole high-depth equality case in the rooted-seven payment inequality. -/
theorem weightedDepthCost_eq_three_mul_rootedSeven_iff
    (c : ℕ) (hc8 : 8 ≤ c) (hc13 : c ≤ 13) :
    (c - 1).choose 5 = 3 * (c - 1).choose 6 ↔ c = 8 := by
  interval_cases c <;> decide

/-! ## Axiom audit -/

#print axioms bled_eq_live_indicator_sub_weightedCost
#print axioms B5_eq_live_sub_weightedDeepCost
#print axioms weightedDeepCost_lt_liveCount_iff_B5_pos
#print axioms lonely_of_exact_weighted_census
#print axioms weightedDeepCost_le_six_threeSeven_threeRooted
#print axioms weightedDepthCost_le_three_mul_rootedSeven
#print axioms weightedDepthCost_eq_three_mul_rootedSeven_iff

end LRC14Concrete
end LonelyRunner
