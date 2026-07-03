/-
  TournamentH7.LRCFarCutSplit — THE FAR-COUNT-7 DISPATCH (klein-2026-07-02-S115, HYP-4020).

  Integrates the fleet's endgame into ONE reduction, splitting the remaining LRC(14) content at
  the union-bound wall `j = 7 = 1/(2·(1/14))` (OPEN-Q-108): each danger arc has measure `1/7`, so
  up to six far arcs leave room (union bound / simulpeel / block-6 work) but seven can tile the
  circle — only the JOINT rate (mac-mini JointRateCore) breaks past.

  This dispatch reduces `LRC14Statement` to exactly four legs:
    (0) the LRC(≤13) citation node (owner-sanctioned);
    (1) hwindow  — bounded families `|v| ≤ 22`  (the machine-checked window census, kps/klein);
    (2) hle6     — COMPRESSED families (no ratio-13 dominant) with ≤ 6 far entries
                   (mac-mini `lonely_of_simul_peel` / kps block-6 territory — union bound holds);
    (3) hge7     — COMPRESSED families with ≥ 7 far entries
                   (the JointRateCore target — the last genuinely-open crux).
  The DOMINANT case (some entry ≥ 13× every other) is discharged unconditionally by klein-S114's
  `top_ratio_lonely_13'`.  Pure logic glue over the proved closers; sorry-free.
-/
import TournamentH7.LRCTopRatioPeel13
import TournamentH7.LRC14CertRoute

namespace LonelyRunner
namespace LRC14

/-- The number of "far" entries (beyond the window cut `22`). -/
def farCount22 (v : Fin 13 → ℤ) : ℕ :=
  (Finset.univ.filter (fun i => 22 < |v i|)).card

/-- **The far-count-7 dispatch.**  LRC(14) from the citation node + the window census + the two
compressed residuals split at the union-bound wall (`≤ 6` far vs `≥ 7` far). -/
theorem lrc14_of_farcut_split (cite : LRCUpTo13)
    (hwindow : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → (∀ i, |v i| ≤ 22) →
      ∃ t : ℝ, Lonely 14 v t)
    (hle6 : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → CoveringFamily v →
      (∀ i0 : Fin 13, ∃ i, i ≠ i0 ∧ |v i0| < 13 * |v i|) →
      farCount22 v ≤ 6 →
      ∃ t : ℝ, Lonely 14 v t)
    (hge7 : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → CoveringFamily v →
      (∀ i0 : Fin 13, ∃ i, i ≠ i0 ∧ |v i0| < 13 * |v i|) →
      7 ≤ farCount22 v →
      ∃ t : ℝ, Lonely 14 v t) :
    LRC14Statement := by
  apply lrc14_of_covering_lonely
  intro v hv hcov
  by_cases hbnd : ∀ i, |v i| ≤ 22
  · exact hwindow v hv hbnd
  · by_cases hpeel : ∃ i0 : Fin 13, ∀ i, i ≠ i0 → 13 * |v i| ≤ |v i0|
    · obtain ⟨i0, hp⟩ := hpeel
      exact top_ratio_lonely_13' cite v hv i0 hp
    · push Not at hpeel
      have hcomp : ∀ i0 : Fin 13, ∃ i, i ≠ i0 ∧ |v i0| < 13 * |v i| := by
        intro i0; obtain ⟨i, hi, hlt⟩ := hpeel i0; exact ⟨i, hi, by omega⟩
      by_cases hc : farCount22 v ≤ 6
      · exact hle6 v hv hcov hcomp hc
      · exact hge7 v hv hcov hcomp (by omega)

#print axioms lrc14_of_farcut_split

end LRC14
end LonelyRunner
