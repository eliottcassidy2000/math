/- LRCFloorAllocation.lean -- opus-2026-07-17-S342 (HYP-7260).
   THM-964's (P) link: THE ALLOCATION LAW.  A comparable 7-block's six
   sorted consecutive ratios multiply to ≤ 13, so the j-th largest ratio
   satisfies r_j^(j+1) ≤ 13 — each Hunter path-tree edge lands in an
   explicit sawtooth-floor class, and the six class floors sum to the
   global μ(U₇) lower bound (THM-964; floors via boxeph's muNum, LEM-042).

   MINING NOTE: the folded factor r·(6 − r) in the all-k closed form
   (consecutive_closed_form, LRCTreeHunter) is the same reflection as the
   staircase triangle's leg product k·(n−1−k) from the project's original
   tournament geometry — the fold at modulus 7 = the λ = 1/14 triangle's
   hypotenuse fold.  One geometry, both projects.

   Kernel-pure: no sorry, no native_decide. -/
import Mathlib

namespace LonelyRunner.LRC14.Hunter

/-- **The sorted-ratio power law**: if `r 0 ≥ r 1 ≥ … ≥ r (n−1) ≥ 1` and the
product is ≤ `C`, then the `j`-th ratio satisfies `(r j)^(j+1) ≤ C`. -/
theorem sorted_ratio_pow_le {n : ℕ} (r : Fin n → ℚ) (C : ℚ)
    (hmono : ∀ i j : Fin n, i ≤ j → r j ≤ r i)
    (h1 : ∀ i, 1 ≤ r i)
    (hprod : ∏ i, r i ≤ C) (j : Fin n) :
    (r j) ^ ((j : ℕ) + 1) ≤ C := by
  have hnonneg : ∀ i, (0 : ℚ) ≤ r i := fun i => le_trans zero_le_one (h1 i)
  -- (r j)^(j+1) ≤ ∏_{i ≤ j} r i ≤ ∏ all r i ≤ C
  have hstep : (r j) ^ ((j : ℕ) + 1) ≤ ∏ i ∈ Finset.univ.filter (· ≤ j), r i := by
    have hcard : (Finset.univ.filter (· ≤ j)).card = (j : ℕ) + 1 := by
      have : Finset.univ.filter (· ≤ j) = Finset.Iic j := by
        ext x; simp
      rw [this, Fin.card_Iic]
    calc (r j) ^ ((j : ℕ) + 1)
        = ∏ _i ∈ Finset.univ.filter (· ≤ j), r j := by
          rw [Finset.prod_const, hcard]
      _ ≤ ∏ i ∈ Finset.univ.filter (· ≤ j), r i := by
          apply Finset.prod_le_prod
          · intro i _; exact hnonneg j
          · intro i hi
            exact hmono i j (by simpa using (Finset.mem_filter.mp hi).2)
  have hrest : ∏ i ∈ Finset.univ.filter (· ≤ j), r i ≤ ∏ i, r i := by
    have hsplit := Finset.prod_filter_mul_prod_filter_not
      Finset.univ (· ≤ j) r
    have hP1 : (0 : ℚ) ≤ ∏ i ∈ Finset.univ.filter (· ≤ j), r i :=
      Finset.prod_nonneg fun i _ => hnonneg i
    have hP2 : (1 : ℚ) ≤ ∏ i ∈ Finset.univ.filter (fun i => ¬ i ≤ j), r i := by
      refine Finset.prod_induction r (fun x => (1 : ℚ) ≤ x)
        (fun a b ha hb => ?_) le_rfl (fun i _ => h1 i)
      nlinarith
    calc ∏ i ∈ Finset.univ.filter (· ≤ j), r i
        ≤ (∏ i ∈ Finset.univ.filter (· ≤ j), r i)
          * ∏ i ∈ Finset.univ.filter (fun i => ¬ i ≤ j), r i :=
          by
            have := mul_le_mul_of_nonneg_left hP2 hP1
            simpa using this
      _ = ∏ i, r i := hsplit
  exact le_trans hstep (le_trans hrest hprod)

/-- **The floor combination** (THM-964 (P), abstract form): with per-edge
floors chosen by ratio class, the Hunter tree sum dominates the sum of the
class floors.  The class floors themselves are muNum minima (LEM-042);
their numeric values are the THM-964 table. -/
theorem floor_combination {n : ℕ} (edgeOverlap classFloor : Fin n → ℚ)
    (hfloor : ∀ j, classFloor j ≤ edgeOverlap j) :
    ∑ j, classFloor j ≤ ∑ j, edgeOverlap j :=
  Finset.sum_le_sum fun j _ => hfloor j

/-! ## Axiom audit -/
#print axioms sorted_ratio_pow_le
#print axioms floor_combination

end LonelyRunner.LRC14.Hunter
