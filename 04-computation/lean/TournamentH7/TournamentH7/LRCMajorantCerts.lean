/-
  TournamentH7.LRCMajorantCerts — THM-705/711/712's LP-vertex certificates
  (mac-mini cont.38): the optimal majorants, their contact patterns, the
  product identity, and the exact thresholds — all as decidable rational facts.
-/
import TournamentH7.LRCFourteenSkeleton

namespace LonelyRunner
namespace LRC14
namespace MajorantCerts

/-- The target g(N) = 1_{N=0} + (1/3)·1_{N=1} on the cells. -/
def gQ : Fin 8 → ℚ
  | 0 => 1
  | 1 => 1/3
  | _ => 0

/-- THM-705's optimal quadratic majorant q₂(N) = 1 − N/2 + N(N−1)/12. -/
def q2 (N : ℚ) : ℚ := 1 - N/2 + N*(N-1)/12

/-- THM-712's optimal cubic majorant. -/
def q3 (N : ℚ) : ℚ := 1 - 2*N/3 + 47*N*(N-1)/252 - 5*N*(N-1)*(N-2)/252

/-- **THM-705 certificate**: q₂ majorizes g on the cells. -/
theorem q2_majorizes : ∀ N : Fin 8, gQ N ≤ q2 (N : ℚ) := by
  intro N
  fin_cases N <;> norm_num [gQ, q2]

/-- **THM-712 certificate**: q₃ majorizes g on the cells. -/
theorem q3_majorizes : ∀ N : Fin 8, gQ N ≤ q3 (N : ℚ) := by
  intro N
  fin_cases N <;> norm_num [gQ, q3]

/-- **THM-711's product identity**: 1 − q₂(N) = N/2 − N(N−1)/12 = N(7−N)/12. -/
theorem hit_empty_product (N : ℚ) : N/2 - N*(N-1)/12 = N*(7-N)/12 := by ring

/-- q₂'s contact pattern {0, 3, 4} (the active LP vertex). -/
theorem q2_contact : q2 0 = gQ 0 ∧ q2 3 = gQ 3 ∧ q2 4 = gQ 4 := by
  refine ⟨?_, ?_, ?_⟩ <;> norm_num [gQ, q2]

/-- q₃'s contact pattern {0, 1, 3, 7}. -/
theorem q3_contact : q3 0 = gQ 0 ∧ q3 1 = gQ 1 ∧ q3 3 = gQ 3 ∧ q3 7 = gQ 7 := by
  refine ⟨?_, ?_, ?_, ?_⟩ <;> norm_num [gQ, q3]

/-- **THM-711's threshold**: 12·(1 − capRat 10) = 432/91. -/
theorem k9_threshold : 12 * (1 - capRat 10) = 432/91 := by norm_num [capRat]

/-- **THM-712's requirement constant**: 1 − capRat 9 = 2025/4004. -/
theorem k8_threshold : 1 - capRat 9 = 2025/4004 := by norm_num [capRat]

end MajorantCerts
end LRC14
end LonelyRunner
