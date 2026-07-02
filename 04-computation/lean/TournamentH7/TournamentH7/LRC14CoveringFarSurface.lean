/-
klein-2026-07-02-S112 (HYP-4017) — THE MINIMAL FINAL SURFACE: LRC(14) from
{LRC(≤13) citation node} + {covering far families are lonely}.

Sharpens kps-S17's `lrc14_of_peel20`: the peel-transport hypothesis is replaced by the
honest open-math name. A far family that is NOT covering never reaches the remaining
hypothesis — it is lonely at t = 1/q by the denominator sieve (`lrc14_of_covering_lonely`,
THM-523 route); a family with no far entry is lonely by the closed window leg
(`hwindow20_closed`, kps-S17 band data + citation). What remains is EXACTLY the covering
floor: covering families with an entry beyond the censused window. The comparison theorem
discharges kps's peel20 hypothesis from this surface, so the two endgames are
interderivable — but `CoveringFarLonely` is the smaller, transport-free assumption.

The instances section records the certified sub-universe of the surface: the consecutive
blocks (an INFINITE covering-far family, via the AP tail certificate) and the ladder's
deep-well family ![1,…,12,182]. These are the shapes the middle-band program generalizes.
-/
import TournamentH7.LRCWindowData
import TournamentH7.LRC14CertRoute
import TournamentH7.LadderPackData

namespace LonelyRunner.LRC14

open LonelyRunner WindowData

/-- **The covering-far surface at cut W**: every covering family with an entry beyond
the window is lonely. This is the single remaining analytic hypothesis of LRC(14). -/
def CoveringFarLonely (W : ℤ) : Prop :=
  ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → CoveringFamily v → 0 < farCountW W v →
    ∃ t : ℝ, Lonely 14 v t

/-- **The W-generic endgame**: LRC(14) from any closed bounded window at cut W together
with the covering-far surface at the same cut. Raising W with new band data strictly
shrinks the remaining hypothesis. -/
theorem lrc14_of_covering_far_of_window (W : ℤ)
    (hwindow : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → (∀ i, |v i| ≤ W) →
      ∃ t : ℝ, Lonely 14 v t)
    (hcovfar : CoveringFarLonely W) : LRC14Statement := by
  apply lrc14_of_covering_lonely
  intro v hv hcov
  by_cases hfar : 0 < farCountW W v
  · exact hcovfar v hv hcov hfar
  · have h0 : farCountW W v = 0 := by omega
    exact hwindow v hv ((farCountW_eq_zero_iff W v).mp h0)

/-- **LRC(14) from the citation node + the covering-far surface at W = 20.**
Strictly leaner than `lrc14_of_peel20`: no loneliness-transport (peel) hypothesis. -/
theorem lrc14_of_covering_far (cite : LRCUpTo13) (hcovfar : CoveringFarLonely 20) :
    LRC14Statement :=
  lrc14_of_covering_far_of_window 20 (hwindow20_closed cite) hcovfar

/-- **Surface comparison**: the covering-far surface discharges kps-S17's `peel20`
hypothesis outright — for a far family, loneliness holds directly (covering → the
surface; non-covering → the sieve), and the all-ones family witnesses the far-count
descent with a constant transport. Hence the endgame pair {cite, peel20} is implied
by {cite, CoveringFarLonely 20}. -/
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

/-! ## Certified instances of the surface

What is already PROVED inside `CoveringFarLonely 20`'s domain. -/

/-- **The consecutive blocks are certified covering-far instances**: for every V ≥ 21
with V % 14 ≠ 13, the family (V, V−1, …, V−12) is covering, has a far entry, and is
machine-checked lonely (AP tail certificate). An infinite certified slice of the
remaining surface. -/
theorem coveringFar_blockFamily (V : ℤ) (hV : 21 ≤ V) (hmod : V % 14 ≠ 13) :
    CoveringFamily (blockFamily V) ∧ 0 < farCountW 20 (blockFamily V) ∧
      ∃ t : ℝ, Lonely 14 (blockFamily V) t := by
  obtain ⟨hcov, hlon⟩ := covering_block_lonely V (by omega) hmod
  refine ⟨hcov, ?_, hlon⟩
  have hb0 : blockFamily V 0 = V := by
    unfold blockFamily
    norm_num
  unfold farCountW
  rw [Finset.card_pos]
  refine ⟨0, Finset.mem_filter.mpr ⟨Finset.mem_univ _, ?_⟩⟩
  rw [hb0, abs_of_pos (by omega : (0 : ℤ) < V)]
  omega

/-- **The ladder's deep well is a certified covering-far instance**: ![1,…,12,182] is
covering (each q ≤ 12 divides itself; 13 and 14 divide 182), has the far entry 182,
and is lonely via the fuel-checked ladder certificate. -/
theorem coveringFar_deepWell :
    CoveringFamily LadderPackData.deepWellV ∧
      0 < farCountW 20 LadderPackData.deepWellV ∧
      ∃ t : ℝ, Lonely 14 LadderPackData.deepWellV t := by
  refine ⟨?_, ?_, LadderPackData.deepWell_lonely⟩
  · intro q hq2 hq14
    interval_cases q
    · exact ⟨1, by decide⟩
    · exact ⟨2, by decide⟩
    · exact ⟨3, by decide⟩
    · exact ⟨4, by decide⟩
    · exact ⟨5, by decide⟩
    · exact ⟨6, by decide⟩
    · exact ⟨7, by decide⟩
    · exact ⟨8, by decide⟩
    · exact ⟨9, by decide⟩
    · exact ⟨10, by decide⟩
    · exact ⟨11, by decide⟩
    · exact ⟨12, by decide⟩
    · exact ⟨12, by decide⟩
  · unfold farCountW
    rw [Finset.card_pos]
    exact ⟨12, Finset.mem_filter.mpr ⟨Finset.mem_univ _, by decide⟩⟩

end LonelyRunner.LRC14
