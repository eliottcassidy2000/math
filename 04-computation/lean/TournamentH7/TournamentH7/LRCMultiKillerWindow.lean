/-
  TournamentH7.LRCMultiKillerWindow — THE MULTI-KILLER WINDOW + THE CITATION MARGIN
  (klein-2026-07-05-S134, HYP-4099).

  The composition that powers the PEEL TOWER: a base with margin β at one point holds a
  Lipschitz window of half-width `(β − 1/14)/B`; up to SIX peeled tops are then handled
  SIMULTANEOUSLY inside that window by the Hunter block step with ZERO pair credits —
  each top's teeth eat at most `1/7` of the window plus a `3/(7·|top|)` boundary fee
  (`teeth_mass`), and `≤ 6` tops leave a positive remainder whenever

      Σ_{tops} 3/(7·|v j|)  <  2·((β − 1/14)/B)·((7 − #tops)/7).

  With citation floors `β = 1/(14 − #tops)` this yields the tower thresholds
  `(k+1)/(13−k)·B`, whose product down the tower telescopes to `C(13,6) = 1716` — an
  explicit a-priori spread bound for the tower's finite leaf.  `cite_margin` extracts
  the β-margin point for any index-subset of ≤ 12 runners from the citation node.

  Kernel-pure; consumes kps's `hunter_block_step`/`teeth_mass` and the S132 transport.
-/
import TournamentH7.LRCRealRegions
import TournamentH7.LRCOneWindowPeel

namespace LonelyRunner
namespace MultiKillerWindow

open LRC14 RealRegion

/-- **The citation margin**: any subset of ≤ 12 of the runners has a common point with
margin `1/(card + 1)`, straight from the LRC(≤13) citation node. -/
theorem cite_margin (cite : LRCUpTo13) (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0)
    (T : Finset (Fin 13)) (hT : T.card ≤ 12) :
    ∃ t : ℝ, ∀ i ∈ T, ∀ m : ℤ, (1 : ℝ) / (T.card + 1) ≤ |(v i : ℝ) * t - m| := by
  let e : Fin T.card ≃ T := T.equivFin.symm
  set w : Fin T.card → ℤ := fun m => v (e m : Fin 13) with hw
  have hwne : ∀ m, w m ≠ 0 := fun m => hv _
  obtain ⟨t, hL⟩ := cite T.card hT w hwne
  refine ⟨t, fun i hi m => ?_⟩
  set m0 : Fin T.card := e.symm ⟨i, hi⟩ with hm0
  have hwm0 : w m0 = v i := by
    rw [hw, hm0]
    simp only
    rw [Equiv.apply_symm_apply]
  have hstep := hL m0 m
  rw [hwm0] at hstep
  have hcast : ((T.card + 1 : ℕ) : ℝ) = (T.card : ℝ) + 1 := by push_cast; ring
  rw [hcast] at hstep
  have : ((T.card : ℕ) : ℝ) + 1 = (T.card : ℝ) + 1 := by norm_num
  linarith [hstep]

/-- **THE MULTI-KILLER WINDOW**: base margin β at one point + up to six tops clearing
the summed boundary-fee criterion ⟹ the whole 13-family is lonely.  The tops are
handled simultaneously by the Hunter block step (zero credits) inside the base's
Lipschitz window: each eats ≤ 1/7 of the window plus its boundary fee. -/
theorem lonely_of_window_multi (v : Fin 13 → ℤ) (tops : Finset (Fin 13))
    (hcard : tops.card ≤ 6) (tstar β B : ℝ)
    (hβ : 1/14 < β)
    (hbase : ∀ i, i ∉ tops → ∀ m : ℤ, β ≤ |(v i : ℝ) * tstar - m|)
    (hB : ∀ i, i ∉ tops → |(v i : ℝ)| ≤ B) (hBpos : 0 < B)
    (hv : ∀ i, v i ≠ 0)
    (hcrit : ((tops.toList.map fun j => 3 / (7 * |(v j : ℝ)|)).sum)
      < 2 * ((β - 1/14) / B) * ((7 - (tops.card : ℝ)) / 7)) :
    ∃ t : ℝ, Lonely 14 v t := by
  set δ : ℝ := (β - 1/14) / B with hδ
  have hδpos : 0 < δ := div_pos (by linarith) hBpos
  set a : ℝ := tstar - δ with ha
  set b : ℝ := tstar + δ with hb
  have hab : a ≤ b := by rw [ha, hb]; linarith
  have hLab : b - a = 2 * δ := by rw [ha, hb]; ring
  set ws : List ℤ := tops.toList.map (fun j => |v j|) with hws
  have hwpos : ∀ w ∈ ws, 0 < w := by
    intro w hw
    rw [hws, List.mem_map] at hw
    obtain ⟨j, -, rfl⟩ := hw
    exact abs_pos.mpr (hv j)
  have hlen : ws.length = tops.card := by
    rw [hws, List.length_map, Finset.length_toList]
  -- the fee sum over the list equals the Finset criterion sum
  have hfees : ((ws.map fun (w : ℤ) => 3 / (7 * (w : ℝ))).sum)
      = ((tops.toList.map fun j => 3 / (7 * |(v j : ℝ)|)).sum) := by
    rw [hws, List.map_map]
    congr 1
    apply List.map_congr_left
    intro j _
    simp only [Function.comp_apply]
    rw [Int.cast_abs]
  -- the zero-credit ledger
  have hledger : 0 < (b - a)
      - ((ws.map fun (w : ℤ) => rlength (rinter [(a, b)] (teeth w a b))).sum)
      + pairCredits [(a, b)] (ws.map fun (w : ℤ) => teeth w a b) := by
    have hcred := pairCredits_nonneg (ws.map fun (w : ℤ) => teeth w a b) [(a, b)]
    have hsingles : ((ws.map fun (w : ℤ) => rlength (rinter [(a, b)] (teeth w a b))).sum)
        ≤ ((ws.map fun (w : ℤ) => (b - a) / 7 + 3 / (7 * (w : ℝ))).sum) := by
      apply List.sum_le_sum
      intro w hw
      rw [rlength_inter_window_clipsum]
      exact teeth_mass (hwpos w hw) a b hab
    have hsplit : ((ws.map fun (w : ℤ) => (b - a) / 7 + 3 / (7 * (w : ℝ))).sum)
        = (ws.length : ℝ) * ((b - a) / 7)
          + ((ws.map fun (w : ℤ) => 3 / (7 * (w : ℝ))).sum) := by
      induction ws with
      | nil => simp
      | cons x t ih =>
          simp only [List.map_cons, List.sum_cons, List.length_cons, ih]
          push_cast
          ring
    rw [hsplit, hfees, hlen] at hsingles
    have hcardR : ((tops.card : ℕ) : ℝ) ≤ 6 := by exact_mod_cast hcard
    have harith : (tops.card : ℝ) * ((b - a) / 7)
        + ((tops.toList.map fun j => 3 / (7 * |(v j : ℝ)|)).sum) < b - a := by
      rw [hLab]
      have h2 : 2 * ((β - 1/14) / B) * ((7 - (tops.card : ℝ)) / 7)
          = 2 * δ - (tops.card : ℝ) * (2 * δ / 7) := by
        rw [hδ]; ring
      rw [h2] at hcrit
      linarith
    linarith
  obtain ⟨t, hta, htb, hgood⟩ := hunter_block_step ws hwpos a b hab hledger
  have htw : |t - tstar| ≤ δ := by
    rw [abs_le]
    constructor
    · rw [ha] at hta; linarith
    · rw [hb] at htb; linarith
  refine ⟨t, fun i m => ?_⟩
  show (1 : ℝ) / (14 : ℕ) ≤ |(v i : ℝ) * t - (m : ℝ)|
  push_cast
  by_cases hi : i ∈ tops
  · have hwmem : |v i| ∈ ws := by
      rw [hws, List.mem_map]
      exact ⟨i, Finset.mem_toList.mpr hi, rfl⟩
    rcases abs_cases (v i) with ⟨heq, -⟩ | ⟨heq, -⟩
    · have h14 := hgood |v i| hwmem m
      rw [heq] at h14
      exact_mod_cast h14
    · have h14 := hgood |v i| hwmem (-m)
      rw [heq] at h14
      have hcast : ((-(v i) : ℤ) : ℝ) * t - ((-m : ℤ) : ℝ) = -((v i : ℝ) * t - (m : ℝ)) := by
        push_cast; ring
      rw [hcast, abs_neg] at h14
      exact h14
  · have hbi := hbase i hi m
    have hBi := hB i hi
    have key : |(v i : ℝ) * tstar - m| ≤ |(v i : ℝ) * t - m| + |(v i : ℝ)| * |tstar - t| := by
      calc |(v i : ℝ) * tstar - m|
          = |((v i : ℝ) * t - m) + (v i : ℝ) * (tstar - t)| := by congr 1; ring
        _ ≤ |(v i : ℝ) * t - m| + |(v i : ℝ) * (tstar - t)| := abs_add_le _ _
        _ = |(v i : ℝ) * t - m| + |(v i : ℝ)| * |tstar - t| := by rw [abs_mul]
    have habs' : |(v i : ℝ)| * |tstar - t| ≤ B * δ := by
      apply mul_le_mul hBi ?_ (abs_nonneg _) (le_of_lt hBpos)
      rwa [abs_sub_comm]
    have hBδ : B * δ = β - 1/14 := by
      rw [hδ]
      field_simp
    linarith

#print axioms cite_margin
#print axioms lonely_of_window_multi

end MultiKillerWindow
end LonelyRunner
