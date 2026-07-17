import TournamentH7.LRCB5NormalizedBridge
import TournamentH7.LRCPairRatioLayerArithmetic

namespace LonelyRunner
namespace LRCPairContinuumBridge

open Finset
open LRCB5CleanModulus
open LRCB5NormalizedBridge
open LRCPairRatioLayerArithmetic

noncomputable section

/-- The elementary error term expected from comparing all pair intersections
on the full circle with the nonzero `q`-grid.  The geometric discrepancy
theorem will instantiate this budget. -/
def pairGridErrorBudget (v : Fin 13 → ℤ) (q : ℕ) : ℚ :=
  (24 * (∑ i, (v i).natAbs : ℚ) + 78) / ((q : ℚ) - 1)

theorem speedMass_ge_thirteen (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0) :
    13 ≤ ∑ i, (v i).natAbs := by
  calc
    13 = ∑ _i : Fin 13, 1 := by simp
    _ ≤ ∑ i, (v i).natAbs := by
      apply Finset.sum_le_sum
      intro i _
      exact Int.natAbs_pos.mpr (hv i)

private theorem pair_error_arithmetic (speedMass speedProduct : ℚ)
    (hmass : 13 ≤ speedMass) (hproduct : 1 ≤ speedProduct) :
    (24 * speedMass + 78) /
        (14 * (534 * speedMass + 1) * speedProduct) ≤ 5 / 1246 := by
  have hfactor :
      534 * speedMass ≤ (534 * speedMass + 1) * speedProduct := by
    calc
      534 * speedMass ≤ 534 * speedMass + 1 := by linarith
      _ ≤ (534 * speedMass + 1) * speedProduct :=
        le_mul_of_one_le_right (by positivity) hproduct
  have hdenLower :
      14 * 534 * speedMass ≤ 14 * (534 * speedMass + 1) * speedProduct := by
    calc
      14 * 534 * speedMass = 14 * (534 * speedMass) := by ring
      _ ≤ 14 * ((534 * speedMass + 1) * speedProduct) :=
        mul_le_mul_of_nonneg_left hfactor (by norm_num)
      _ = 14 * (534 * speedMass + 1) * speedProduct := by ring
  have hnum : 24 * speedMass + 78 ≤ 30 * speedMass := by linarith
  have hdenPos : 0 < 14 * (534 * speedMass + 1) * speedProduct := by
    have : 0 < speedMass := lt_of_lt_of_le (by norm_num) hmass
    positivity
  apply (div_le_iff₀ hdenPos).2
  calc
    24 * speedMass + 78 ≤ 30 * speedMass := hnum
    _ = (5 / 1246 : ℚ) * (14 * 534 * speedMass) := by ring
    _ ≤ (5 / 1246 : ℚ) *
        (14 * (534 * speedMass + 1) * speedProduct) :=
      mul_le_mul_of_nonneg_left hdenLower (by norm_num)

/-- Height `534` makes the proposed clean-grid error smaller than `5/1246`.
This is the exact arithmetic half of the continuum-to-grid bridge. -/
theorem pairGridErrorBudget_cleanModulus_534_le
    (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0) :
    pairGridErrorBudget v (cleanModulus v 534) ≤ 5 / 1246 := by
  let speedMass : ℕ := ∑ i, (v i).natAbs
  let speedProduct : ℕ := ∏ i, (v i).natAbs
  have hmassNat : 13 ≤ speedMass := by
    simpa [speedMass] using speedMass_ge_thirteen v hv
  have hproductNat : 1 ≤ speedProduct := by
    have hpos : 0 < speedProduct := by
      dsimp [speedProduct]
      exact Finset.prod_pos fun i _ => Int.natAbs_pos.mpr (hv i)
    omega
  have hmass : (13 : ℚ) ≤ speedMass := by exact_mod_cast hmassNat
  have hproduct : (1 : ℚ) ≤ speedProduct := by exact_mod_cast hproductNat
  have hdenEq :
      ((cleanModulus v 534 : ℕ) : ℚ) - 1 =
        14 * (534 * (speedMass : ℚ) + 1) * speedProduct := by
    simp only [cleanModulus, speedMass, speedProduct]
    push_cast
    ring
  rw [pairGridErrorBudget, hdenEq]
  have hmassCast : (speedMass : ℚ) = ∑ i, ((v i).natAbs : ℚ) := by
    simp [speedMass]
  rw [← hmassCast]
  exact pair_error_arithmetic speedMass speedProduct hmass hproduct

/-- The explicit clean-grid budget fits strictly inside THM-954's conservative
two-path margin. -/
theorem cleanModulus_534_error_lt_path_margin
    (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0) :
    pairGridErrorBudget v (cleanModulus v 534) <
      13 / 30 - negativePairTierBoundPathOnly := by
  calc
    pairGridErrorBudget v (cleanModulus v 534) ≤ 5 / 1246 :=
      pairGridErrorBudget_cleanModulus_534_le v hv
    _ < 8270807 / 2058376320 := by norm_num
    _ = 13 / 30 - negativePairTierBoundPathOnly :=
      negativePairTierBoundPathOnly_margin.symm

#print axioms speedMass_ge_thirteen
#print axioms pairGridErrorBudget_cleanModulus_534_le
#print axioms cleanModulus_534_error_lt_path_margin

end
end LRCPairContinuumBridge
end LonelyRunner
