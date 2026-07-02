/-
klein-2026-07-02-S112 (HYP-4017) — THE MINIMAL FINAL SURFACE: LRC(14) from
{LRC(≤13) citation node} + {covering far families are lonely}.

Sharpens kps-S17's `lrc14_of_peel20`: the peel-transport hypothesis is replaced by the
honest open-math name. A far family that is NOT covering never reaches the remaining
hypothesis — it is lonely at t = 1/q by the denominator sieve (`lrc14_of_covering_lonely`,
THM-523 route); a family with no far entry is lonely by the closed window leg
(`hwindow20_closed`, kps-S17 band data + citation). What remains is EXACTLY the covering
floor: covering families with an entry beyond the censused window. The comparison theorem
below discharges kps's peel20 hypothesis from this surface, so the two endgames are
interderivable — but `CoveringFarLonely 20` is the smaller, transport-free assumption.
-/
import TournamentH7.LRCWindowData
import TournamentH7.LRC14CertRoute

namespace LonelyRunner.LRC14

open LonelyRunner WindowData

/-- **The covering-far surface at cut W**: every covering family with an entry beyond
the window is lonely. This is the single remaining analytic hypothesis of LRC(14). -/
def CoveringFarLonely (W : ℤ) : Prop :=
  ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → CoveringFamily v → 0 < farCountW W v →
    ∃ t : ℝ, Lonely 14 v t

/-- **LRC(14) from the citation node + the covering-far surface at W = 20.**
Strictly leaner than `lrc14_of_peel20`: no loneliness-transport (peel) hypothesis. -/
theorem lrc14_of_covering_far (cite : LRCUpTo13) (hcovfar : CoveringFarLonely 20) :
    LRC14Statement := by
  apply lrc14_of_covering_lonely
  intro v hv hcov
  by_cases hfar : 0 < farCountW 20 v
  · exact hcovfar v hv hcov hfar
  · have h0 : farCountW 20 v = 0 := by omega
    exact hwindow20_closed cite v hv ((farCountW_eq_zero_iff 20 v).mp h0)

/-- **Surface comparison**: the covering-far surface discharges kps-S17's `peel20`
hypothesis outright — for a far family, loneliness holds directly (covering → the
surface; non-covering → the sieve), and the all-ones family witnesses the far-count
descent with a constant transport. Hence `lrc14_of_peel20`'s hypothesis pair
{cite, peel20} is implied by {cite, CoveringFarLonely 20}. -/
theorem peel20_of_covering_far (hcovfar : CoveringFarLonely 20) :
    ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → 0 < farCountW 20 v →
      ∃ v' : Fin 13 → ℤ, (∀ i, v' i ≠ 0) ∧ farCountW 20 v' < farCountW 20 v ∧
        ((∃ t : ℝ, Lonely 14 v' t) → ∃ t : ℝ, Lonely 14 v t) := by
  intro v hv hfar
  have hlon : ∃ t : ℝ, Lonely 14 v t := by
    by_cases hc : CoveringFamily v
    · exact hcovfar v hv hc hfar
    · unfold CoveringFamily at hc
      push Not at hc
      obtain ⟨q, hq2, hq14, hdiv⟩ := hc
      exact ⟨(1 : ℝ) / q, sieve_one_div 14 q v hq14 (by omega) hdiv⟩
  refine ⟨fun _ => 1, fun _ => one_ne_zero, ?_, fun _ => hlon⟩
  have hone : farCountW 20 (fun _ => (1 : ℤ)) = 0 := by simp [farCountW]
  omega

end LonelyRunner.LRC14
