/-
  TournamentH7.LRCResidueLiar — THE RESIDUE-LIAR FAMILY, CLOSED BY FORMULA
  (kind-pasteur-2026-07-04-S3, HYP-4078).

  The census-hard single-swap coverer family {1,…,11,13,12K} (12 | X) — the shape that
  mimics the arithmetic progression's residues but stays strictly loose — is lonely at the
  explicit rational time t = (5K+2)/(12K+5), with margin M = K/(12K+5) > 1/14 for every
  K ≥ 3.  This is the LRC(14) analogue of the tournament H=7 gap (THM-029): ONE
  combinatorial obstruction, closed by an explicit forcing (here a residue table), giving
  an INFINITE family of lonely certificates with NO native_decide.

  Proof: at t = (5K+2)/(12K+5) each runner v sits at r_v/(12K+5) where r_v = v·(5K+2) mod
  (12K+5) is a LINEAR polynomial in K, and K ≤ r_v ≤ (12K+5) − K for every runner (checked
  termwise, binding at v=5 with r=K and v=12K with r=11K+5), so the lattice distance is ≥ K,
  and 14·K ≥ 12K+5 ⟺ K ≥ 3.

  Kernel-pure (propext, Classical.choice, Quot.sound); no sorry, no native_decide.
-/
import Mathlib
import TournamentH7.LonelyRunner
import TournamentH7.LRCSpread13

namespace LonelyRunner
namespace LRC14

/-- **LATTICE-DISTANCE FLOOR**: if `N = qq·Q + r` with `κ ≤ r ≤ Q − κ` (`0 < Q`, `0 ≤ κ`),
then `N` is integer-distance `≥ κ` from every multiple of `Q`. -/
lemma lattice_dist_ge {Q kappa N : ℤ} (qq r : ℤ) (hQ : 0 < Q)
    (hN : N = qq * Q + r) (hlo : kappa ≤ r) (hhi : r ≤ Q - kappa) (m : ℤ) :
    kappa ≤ |N - m * Q| := by
  have hrw : N - m * Q = (qq - m) * Q + r := by rw [hN]; ring
  rw [le_abs]
  rcases lt_trichotomy (qq - m) 0 with h | h | h
  · -- qq − m ≤ −1 : N − mQ ≤ −κ, so κ ≤ −(N − mQ)
    right
    have hle : qq - m ≤ -1 := by omega
    nlinarith [hrw, hhi, mul_le_mul_of_nonneg_right hle hQ.le]
  · -- qq = m : N − mQ = r ≥ κ
    left
    rw [hrw, h, zero_mul, zero_add]; exact hlo
  · -- qq − m ≥ 1 : N − mQ ≥ Q + κ ≥ κ
    left
    have hge : 1 ≤ qq - m := by omega
    nlinarith [hrw, hlo, hQ, mul_le_mul_of_nonneg_right hge hQ.le]

/-- The residue-liar / single-swap coverer family `{1,…,11,13,12K}` (13 speeds). -/
def residueLiar (K : ℤ) : Fin 13 → ℤ :=
  ![1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 12 * K]

/-- **THE RESIDUE-LIAR FAMILY IS LONELY** at `t = (5K+2)/(12K+5)` for every `K ≥ 3`:
an infinite family of LRC(14) lonely certificates in closed form (margin `K/(12K+5) > 1/14`),
kernel-pure, no `native_decide`. -/
theorem residueLiar_lonely (K : ℤ) (hK : 3 ≤ K) :
    Lonely 14 (residueLiar K) (((5 * K + 2 : ℤ) : ℝ) / ((12 * K + 5 : ℤ) : ℝ)) := by
  apply lonely14_of_ratio (residueLiar K) (5 * K + 2) (12 * K + 5) (by omega)
  intro i m
  -- reusable closer: value `val`, quotient `qq`, residue `r`
  have key : ∀ val qq r : ℤ, val * (5 * K + 2) = qq * (12 * K + 5) + r →
      K ≤ r → r ≤ (12 * K + 5) - K →
      (12 * K + 5) ≤ 14 * |val * (5 * K + 2) - m * (12 * K + 5)| := by
    intro val qq r hEq hlo hhi
    have hm : K ≤ |val * (5 * K + 2) - m * (12 * K + 5)| :=
      lattice_dist_ge qq r (by omega) hEq hlo hhi m
    nlinarith [hm, hK]
  fin_cases i <;> simp only [residueLiar]
  · exact key 1 0 (5 * K + 2) (by ring) (by omega) (by omega)
  · exact key 2 0 (10 * K + 4) (by ring) (by omega) (by omega)
  · exact key 3 1 (3 * K + 1) (by ring) (by omega) (by omega)
  · exact key 4 1 (8 * K + 3) (by ring) (by omega) (by omega)
  · exact key 5 2 K (by ring) (by omega) (by omega)
  · exact key 6 2 (6 * K + 2) (by ring) (by omega) (by omega)
  · exact key 7 2 (11 * K + 4) (by ring) (by omega) (by omega)
  · exact key 8 3 (4 * K + 1) (by ring) (by omega) (by omega)
  · exact key 9 3 (9 * K + 3) (by ring) (by omega) (by omega)
  · exact key 10 4 (2 * K) (by ring) (by omega) (by omega)
  · exact key 11 4 (7 * K + 2) (by ring) (by omega) (by omega)
  · exact key 13 5 (5 * K + 1) (by ring) (by omega) (by omega)
  · exact key (12 * K) (5 * K - 1) (11 * K + 5) (by ring) (by omega) (by omega)

/-- Concrete instance: the first COVERING residue-liar `{1,…,11,13,84}` (K=7, `84 = 6·14`,
`12K+5 = 89 = F₁₁`) is lonely at `37/89`, margin `7/89 > 1/14`. -/
theorem residueLiar84_lonely :
    Lonely 14 (residueLiar 7) ((37 : ℝ) / 89) := by
  have h := residueLiar_lonely 7 (by norm_num)
  norm_num at h
  exact h

-- Kernel-purity: expect only [propext, Classical.choice, Quot.sound] — NO native_decide/ofReduceBool.
#print axioms residueLiar_lonely
#print axioms residueLiar84_lonely

end LRC14
end LonelyRunner
