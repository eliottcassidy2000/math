/-
  TournamentH7.LRCEigenTransfer — THM-710 formalized: the factorial moments are
  EIGENVECTORS of the far-element occupancy transfer (mac-mini cont.38).

  The operator (kps THM-700/701): T p j = ((7−j)/7)·p j + ((j+1)/7)·p (j+1)
  (with p 8 := 0).  The falling-factorial moments m_r = Σ_j j^(r) p j satisfy
  m_r(T p) = ((7−r)/7)·m_r(p) EXACTLY — proved here for r = 1, 2, 3 (the ladder's
  degrees) by full expansion over the eight cells and `ring`.  The rung facts
  (THM-710's propagation): (6/7)·cap_{k+1} + 1/7 ≤ cap_{k+2}, exact rationals.
-/
import Mathlib
import TournamentH7.LRCFourteenSkeleton

namespace LonelyRunner
namespace LRC14
namespace EigenTransfer

/-- The occupancy distribution as an 8-vector (cells N = 0..7). -/
def T (p : Fin 8 → ℝ) : Fin 8 → ℝ :=
  fun j => ((7 - (j : ℝ)) / 7) * p j +
    (((j : ℝ) + 1) / 7) * (if h : (j : ℕ) + 1 < 8 then p ⟨(j : ℕ) + 1, h⟩ else 0)

/-- First factorial moment. -/
def m1 (p : Fin 8 → ℝ) : ℝ := ∑ j : Fin 8, (j : ℝ) * p j
/-- Second factorial moment. -/
def m2 (p : Fin 8 → ℝ) : ℝ := ∑ j : Fin 8, (j : ℝ) * ((j : ℝ) - 1) * p j
/-- Third factorial moment. -/
def m3 (p : Fin 8 → ℝ) : ℝ :=
  ∑ j : Fin 8, (j : ℝ) * ((j : ℝ) - 1) * ((j : ℝ) - 2) * p j

/-- **THM-710, r = 1**: `m1` is an eigenvector of the transfer, eigenvalue `6/7`. -/
theorem m1_eigen (p : Fin 8 → ℝ) : m1 (T p) = (6 / 7) * m1 p := by
  simp only [m1, T, Fin.sum_univ_succ]
  norm_num
  ring

/-- **THM-710, r = 2**: eigenvalue `5/7`. -/
theorem m2_eigen (p : Fin 8 → ℝ) : m2 (T p) = (5 / 7) * m2 p := by
  simp only [m2, T, Fin.sum_univ_succ]
  norm_num
  ring

/-- **THM-710, r = 3**: eigenvalue `4/7`. -/
theorem m3_eigen (p : Fin 8 → ℝ) : m3 (T p) = (4 / 7) * m3 p := by
  simp only [m3, T, Fin.sum_univ_succ]
  norm_num
  ring

/-- **The majorant-value transfer**: V = 1 − m1/2 + m2/12 obeys
`V(Tp) = (6/7)·V(p) + 1/7 − m2(p)/84`. -/
theorem V_transfer (p : Fin 8 → ℝ) :
    1 - m1 (T p) / 2 + m2 (T p) / 12
      = (6 / 7) * (1 - m1 p / 2 + m2 p / 12) + 1 / 7 - m2 p / 84 := by
  rw [m1_eigen, m2_eigen]
  ring

/-- **THM-710's rung propagation, exact**: `(6/7)·capRat (k+1) + 1/7 ≤ capRat (k+2)`
for every binding rung `k = 8..11`. -/
theorem rung_propagation : ∀ k ∈ [8, 9, 10, 11],
    (6 / 7) * (capRat (k + 1) : ℝ) + 1 / 7 ≤ (capRat (k + 2) : ℝ) := by
  intro k hk
  fin_cases hk <;> norm_num [capRat]

end EigenTransfer
end LRC14
end LonelyRunner
