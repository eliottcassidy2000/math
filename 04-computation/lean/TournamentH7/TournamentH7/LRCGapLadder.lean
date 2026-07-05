/-
  TournamentH7.LRCGapLadder — THE l ≥ 2 PEEL LADDER AT THE GAP THRESHOLD
  (kind-pasteur-2026-07-05-S7, HYP-4115; renumbered from 4113, opus-S80 first-committed).

  Extends the S6 top-compression (`gap_compressed_24`, the `l = 1` rung) to the
  full order-statistic ladder, by running klein's radius-parametric window stack
  (S136, `LRCTeethR`) at `ρ = 2/25` instead of `1/13`:

  * `gap_tower_step` — the ρ = 2/25 mirror of klein's `tower_step_12`: remove
    `l ≥ 1` tops, cite the `12 − l` sub-base at `1/(13−l)` (which exceeds `2/25`
    for EVERY `l ≥ 1`), and if the tops' teeth-mass fees clear the window, the
    WHOLE family has a `2/25`-margin point — the dichotomy's loose branch.

  * `gap_ladder_rung` — the contrapositive, integer form: a gap violator (no
    `2/25`-point) satisfies, for EVERY subset `S` with `1 ≤ |S| = l ≤ 6` and
    every bound `B` on the complement:

        ∃ j ∈ S,  (2l−1)·(25−4l)·|vⱼ| ≤ 150·l·(13−l)·B.

    Taking `S` = the top-`l` runners: `w₍ₗ₎ ≤ C_l·w₍ₗ₊₁₎` with
    `C_l = 150·l·(13−l)/((2l−1)(25−4l))` ≈ 64.7, 69.2, 85.7, 133.3, 572.7 for
    `l = 2..6` (the `l = 1` rung has the sharper one-tooth constant 24 from S6).
    The fee budget dies at `l = 7` (`4l ≥ 25`): six rungs is the whole ladder.

  A gap violator's order statistics are thus chained top-down: no runner can be
  more than a bounded factor above the next `l` runners, for every `l ≤ 6`,
  while the S2 spread gate forces the whole set to be wide.  Squeezed at every
  scale.  Mirror-don't-reinvent: the machinery is klein's `margin_of_window_of_fees`
  + `teethR_mass` + `rlength_inter_window_clipsum`; this file only re-aims it.

  Kernel-pure; no `sorry`, no `native_decide`.
-/
import Mathlib
import TournamentH7.LRCTeethR

namespace LonelyRunner
namespace GapLadder

open TeethR RealRegion

private lemma sum_map_const {α : Type*} (L : List α) (c : ℝ) :
    (L.map fun _ => c).sum = L.length * c := by
  induction L with
  | nil => simp
  | cons x xs ih =>
    simp only [List.map_cons, List.sum_cons, List.length_cons, ih]
    push_cast
    ring

/-- **THE GAP-LEVEL TOWER STEP** (ρ = 2/25 mirror of klein's `tower_step_12`):
a 12-family whose `l ≥ 1` tops clear the fee criterion has a common `2/25`-margin
point — the sub-base of `12 − l` runners is cited at `1/(13−l) > 2/25`, and the
tops ride the ρ = 2/25 window.  The conclusion is exactly the dichotomy's loose
branch. -/
theorem gap_tower_step (cite : LRCUpTo13) (v : Fin 12 → ℤ) (hv : ∀ i, v i ≠ 0)
    (tops : List (Fin 12)) (fee : Fin 12 → ℝ) (B : ℝ) (hBpos : 0 < B)
    (hl1 : tops ≠ [])
    (hB : ∀ i, i ∉ tops → |(v i : ℝ)| ≤ B)
    (hfee : ∀ (tstar : ℝ), ∀ j ∈ tops,
      rlength (rinter
        [(tstar - (1/(13 - (tops.toFinset.card : ℝ)) - 2/25) / B,
          tstar + (1/(13 - (tops.toFinset.card : ℝ)) - 2/25) / B)]
        (teethR (2/25) |v j|
          (tstar - (1/(13 - (tops.toFinset.card : ℝ)) - 2/25) / B)
          (tstar + (1/(13 - (tops.toFinset.card : ℝ)) - 2/25) / B))) ≤ fee j)
    (hcrit : (tops.map fee).sum
      < 2 * ((1/(13 - (tops.toFinset.card : ℝ)) - 2/25) / B)) :
    ∃ t : ℝ, ∀ i, ∀ m : ℤ, (2 : ℝ) / 25 ≤ |(v i : ℝ) * t - m| := by
  classical
  set T : Finset (Fin 12) := Finset.univ \ tops.toFinset with hT
  have hTcard : T.card = 12 - tops.toFinset.card := by
    rw [hT, Finset.card_univ_diff, Fintype.card_fin]
  have htopspos : 0 < tops.toFinset.card := by
    rcases tops with _ | ⟨j, rest⟩
    · exact absurd rfl hl1
    · apply Finset.card_pos.mpr
      exact ⟨j, by simp⟩
  have htops12 : tops.toFinset.card ≤ 12 := by
    calc tops.toFinset.card ≤ Finset.univ.card := Finset.card_le_card (Finset.subset_univ _)
      _ = 12 := by rw [Finset.card_univ, Fintype.card_fin]
  have hTle : T.card ≤ 12 := by omega
  obtain ⟨tstar, hmargin⟩ := cite_margin_gen cite v T hTle (fun i _ => hv i)
  set β : ℝ := 1 / (13 - (tops.toFinset.card : ℝ)) with hβdef
  have hcard13 : (tops.toFinset.card : ℝ) ≤ 12 := by exact_mod_cast htops12
  have hcard1 : (1 : ℝ) ≤ (tops.toFinset.card : ℝ) := by exact_mod_cast htopspos
  have hβgt : (2 : ℝ) / 25 < β := by
    rw [hβdef, div_lt_div_iff₀ (by norm_num) (by linarith)]
    linarith
  have hbase : ∀ i, i ∉ tops → ∀ m : ℤ, β ≤ |(v i : ℝ) * tstar - m| := by
    intro i hi m
    have hiT : i ∈ T := by
      rw [hT, Finset.mem_sdiff]
      exact ⟨Finset.mem_univ i, fun hc => hi (List.mem_toFinset.mp hc)⟩
    have hstep := hmargin i hiT m
    have hcast : (T.card : ℝ) + 1 = 13 - (tops.toFinset.card : ℝ) := by
      rw [hTcard]
      push_cast [Nat.cast_sub htops12]
      ring
    rw [hβdef, ← hcast]
    exact hstep
  exact margin_of_window_of_fees v tops fee (2/25) tstar β B
    (by norm_num) (by norm_num) hβgt hbase hB hBpos (fun j _ => hv j) (hfee tstar) hcrit

/-- **THE LADDER RUNG (integer form).**  A gap violator (no `2/25`-margin point)
satisfies, for every subset `S` with `1 ≤ |S| ≤ 6` and every positive bound `B`
on the complement: some `j ∈ S` has

    `(2·|S|−1)·(25−4·|S|)·|vⱼ| ≤ 150·|S|·(13−|S|)·B`.

Taking `S` = the top-`l` runners bounds each order statistic by `C_l` times the
next: `C_2 ≈ 64.7, C_3 ≈ 69.2, C_4 ≈ 85.7, C_5 ≈ 133.3, C_6 ≈ 572.7`. -/
theorem gap_ladder_rung (cite : LRCUpTo13) (v : Fin 12 → ℤ) (hv : ∀ i, v i ≠ 0)
    (hnl : ∀ t : ℝ, ∃ i, ∃ m : ℤ, |(v i : ℝ) * t - m| < 2 / 25)
    (S : Finset (Fin 12)) (hS0 : S.Nonempty) (hS6 : S.card ≤ 6)
    (B : ℤ) (hB0 : 0 < B) (hB : ∀ i, i ∉ S → |v i| ≤ B) :
    ∃ j ∈ S, (2 * (S.card : ℤ) - 1) * (25 - 4 * S.card) * |v j|
      ≤ 150 * S.card * (13 - (S.card : ℤ)) * B := by
  classical
  by_contra hcon
  push_neg at hcon
  -- notation
  set l : ℕ := S.card with hldef
  have hl1 : 1 ≤ l := Finset.card_pos.mpr hS0
  have hlR1 : (1 : ℝ) ≤ (l : ℝ) := by exact_mod_cast hl1
  have hlR6 : (l : ℝ) ≤ 6 := by exact_mod_cast hS6
  have hBR : (0 : ℝ) < (B : ℝ) := by exact_mod_cast hB0
  -- positive ladder factors (over ℝ)
  have hf1 : (0 : ℝ) < 2 * l - 1 := by linarith
  have hf2 : (0 : ℝ) < 25 - 4 * l := by linarith
  have hf3 : (0 : ℝ) < 13 - (l : ℝ) := by linarith
  have hf4 : (0 : ℝ) < 150 * l * (13 - (l : ℝ)) * B := by positivity
  -- the tops list
  set tops : List (Fin 12) := S.toList with htopsdef
  have htopsFinset : tops.toFinset = S := by
    rw [htopsdef]; exact Finset.toList_toFinset S
  have htopslen : tops.length = l := by
    rw [htopsdef, Finset.length_toList, hldef]
  have htopsne : tops ≠ [] := by
    intro h
    have hempty : S = ∅ := by
      rw [← htopsFinset, h]
      rfl
    exact Finset.nonempty_iff_ne_empty.mp hS0 hempty
  have hmemtops : ∀ j, j ∈ tops ↔ j ∈ S := by
    intro j
    rw [htopsdef, Finset.mem_toList]
  -- the window half-width δ and the fee function
  set δ : ℝ := (1 / (13 - (l : ℝ)) - 2 / 25) / (B : ℝ) with hδdef
  have hδnum : (0 : ℝ) < 1 / (13 - (l : ℝ)) - 2 / 25 := by
    rw [sub_pos, div_lt_div_iff₀ (by norm_num) hf3]
    linarith
  have hδ0 : 0 < δ := div_pos hδnum hBR
  set fee : Fin 12 → ℝ := fun j => 2 * (2/25) * (2 * δ) + 6 * (2/25) / (|v j| : ℝ)
    with hfeedef
  -- per-top speed positivity and the reciprocal bound from the contradiction
  have hvpos : ∀ j ∈ S, (0 : ℝ) < (|v j| : ℝ) := by
    intro j hj
    exact_mod_cast abs_pos.mpr (hv j)
  have hrecip : ∀ j ∈ S, 6 * (2/25) / (|v j| : ℝ)
      < (12/25) * ((2 * l - 1) * (25 - 4 * l)) / (150 * l * (13 - (l : ℝ)) * B) := by
    intro j hj
    have hcj := hcon j hj
    -- integer strict: 150 l (13−l) B < (2l−1)(25−4l)|v j|
    have hcjR : 150 * l * (13 - (l : ℝ)) * B < (2 * l - 1) * (25 - 4 * l) * (|v j| : ℝ) := by
      have h1 : (150 * (S.card : ℤ) * (13 - (S.card : ℤ)) * B : ℤ)
          < (2 * (S.card : ℤ) - 1) * (25 - 4 * S.card) * |v j| := hcj
      calc (150 : ℝ) * l * (13 - (l : ℝ)) * B
          = ((150 * (S.card : ℤ) * (13 - (S.card : ℤ)) * B : ℤ) : ℝ) := by
            push_cast [hldef]; ring
        _ < (((2 * (S.card : ℤ) - 1) * (25 - 4 * S.card) * |v j| : ℤ) : ℝ) := by
            exact_mod_cast h1
        _ = (2 * l - 1) * (25 - 4 * l) * (|v j| : ℝ) := by
            push_cast [hldef]; ring
    rw [div_lt_div_iff₀ (hvpos j hj) hf4]
    have hexpand : (12/25) * ((2 * l - 1) * (25 - 4 * l)) * (|v j| : ℝ)
        = (12/25) * ((2 * l - 1) * (25 - 4 * l) * (|v j| : ℝ)) := by ring
    rw [hexpand]
    have h2 : 6 * (2/25) * (150 * l * (13 - (l : ℝ)) * B)
        = (12/25) * (150 * l * (13 - (l : ℝ)) * B) := by ring
    rw [h2]
    apply mul_lt_mul_of_pos_left hcjR
    norm_num
  -- the fee criterion: per-element STRICT (from the integer bump) + the exact
  -- identity l·K = 2δ (the fee budget saturates exactly at the reciprocal bound)
  set K : ℝ := 2 * (2/25) * (2 * δ)
    + (12/25) * ((2 * l - 1) * (25 - 4 * l)) / (150 * l * (13 - (l : ℝ)) * B)
    with hKdef
  have hδval : δ = (2 * l - 1) / (25 * (13 - (l : ℝ)) * B) := by
    have h13 : (13 - (l : ℝ)) ≠ 0 := ne_of_gt hf3
    have hBne : (B : ℝ) ≠ 0 := ne_of_gt hBR
    rw [hδdef]
    field_simp
    ring
  have hKeq : (l : ℝ) * K = 2 * δ := by
    have h13 : (13 - (l : ℝ)) ≠ 0 := ne_of_gt hf3
    have hBne : (B : ℝ) ≠ 0 := ne_of_gt hBR
    have hlne : (l : ℝ) ≠ 0 := by linarith
    rw [hKdef, hδval]
    field_simp
    ring
  have hcrit : (tops.map fee).sum < 2 * δ := by
    obtain ⟨j₀, hj₀⟩ := List.exists_mem_of_ne_nil tops htopsne
    have hstrict : (tops.map fee).sum < (tops.map fun _ => K).sum := by
      apply List.sum_lt_sum
      · intro j hj
        have hr := hrecip j ((hmemtops j).mp hj)
        simp only [hfeedef, hKdef]
        linarith
      · refine ⟨j₀, hj₀, ?_⟩
        have hr := hrecip j₀ ((hmemtops j₀).mp hj₀)
        simp only [hfeedef, hKdef]
        linarith
    calc (tops.map fee).sum < (tops.map fun _ => K).sum := hstrict
      _ = l * K := by rw [sum_map_const, htopslen]
      _ = 2 * δ := hKeq
  -- assemble: the tower step gives a 2/25 point, contradicting hnl
  have hBbound : ∀ i, i ∉ tops → |(v i : ℝ)| ≤ (B : ℝ) := by
    intro i hi
    have : i ∉ S := fun hc => hi ((hmemtops i).mpr hc)
    have := hB i this
    rw [← Int.cast_abs]
    exact_mod_cast this
  have hfee : ∀ (tstar : ℝ), ∀ j ∈ tops,
      rlength (rinter
        [(tstar - (1/(13 - (tops.toFinset.card : ℝ)) - 2/25) / (B : ℝ),
          tstar + (1/(13 - (tops.toFinset.card : ℝ)) - 2/25) / (B : ℝ))]
        (teethR (2/25) |v j|
          (tstar - (1/(13 - (tops.toFinset.card : ℝ)) - 2/25) / (B : ℝ))
          (tstar + (1/(13 - (tops.toFinset.card : ℝ)) - 2/25) / (B : ℝ)))) ≤ fee j := by
    intro tstar j hj
    have hcard : (tops.toFinset.card : ℝ) = (l : ℝ) := by
      rw [htopsFinset, hldef]
    rw [hcard]
    rw [rlength_inter_window_clipsum]
    have hwj : 0 < |v j| := abs_pos.mpr (hv j)
    have hmass := teethR_mass (by norm_num : (0:ℝ) ≤ 2/25) hwj
      (tstar - δ) (tstar + δ) (by linarith)
    have hwidth : (tstar + δ) - (tstar - δ) = 2 * δ := by ring
    rw [hwidth, Int.cast_abs] at hmass
    rw [hfeedef]
    simp only
    -- match the window expressions: (1/(13−l) − 2/25)/B = δ
    have hδmatch : (1/(13 - (l : ℝ)) - 2/25) / (B : ℝ) = δ := by rw [hδdef]
    rw [hδmatch]
    exact hmass
  have hcrit' : (tops.map fee).sum
      < 2 * ((1/(13 - (tops.toFinset.card : ℝ)) - 2/25) / (B : ℝ)) := by
    have hcard : (tops.toFinset.card : ℝ) = (l : ℝ) := by
      rw [htopsFinset, hldef]
    rw [hcard]
    have hδmatch : (1/(13 - (l : ℝ)) - 2/25) / (B : ℝ) = δ := by rw [hδdef]
    rw [hδmatch]
    exact hcrit
  obtain ⟨t, hgood⟩ := gap_tower_step cite v hv tops fee (B : ℝ) hBR htopsne
    hBbound hfee hcrit'
  obtain ⟨i, m, hlt⟩ := hnl t
  exact absurd (hgood i m) (not_le.mpr hlt)

#print axioms gap_tower_step
#print axioms gap_ladder_rung

end GapLadder
end LonelyRunner
