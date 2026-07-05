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

/-! ## The far-top sharpening (S134 named-next; HYP-4104)

A window SHORTER than one inter-tooth gap `6/(7w)` meets at most ONE tooth of runner
`w`, so its clipped mass is `≤ 1/(7w)` — no density term, no boundary fee.  Far tops
therefore compose WITHOUT the 7-wall: ANY number of them fit one window. -/

/-- A list whose elements pairwise cannot both be positive sums to at most the cap. -/
theorem sum_le_single_bound (l : List ℝ) (c : ℝ) (hc : 0 ≤ c)
    (hnn : ∀ x ∈ l, 0 ≤ x) (hcap : ∀ x ∈ l, x ≤ c)
    (hpair : l.Pairwise (fun x y => x = 0 ∨ y = 0)) :
    l.sum ≤ c := by
  induction l with
  | nil => simpa using hc
  | cons x t ih =>
      rw [List.pairwise_cons] at hpair
      rcases eq_or_ne x 0 with rfl | hx0
      · rw [List.sum_cons, zero_add]
        exact ih (fun y hy => hnn y (List.mem_cons_of_mem _ hy))
          (fun y hy => hcap y (List.mem_cons_of_mem _ hy)) hpair.2
      · have hzero : ∀ y ∈ t, y = 0 := by
          intro y hy
          rcases hpair.1 y hy with h | h
          · exact absurd h hx0
          · exact h
        have hts : t.sum = 0 := List.sum_eq_zero hzero
        rw [List.sum_cons, hts, add_zero]
        exact hcap x (List.mem_cons_self ..)

/-- Positive clip forces the tooth's edges to straddle the window. -/
theorem tooth_clip_pos_edges {w : ℤ} (hw : 0 < w) {a b : ℝ} {n : ℤ}
    (hn : 0 < clipLen (tooth w n) a b) :
    ((n : ℝ) - 1/14) / w < b ∧ a < ((n : ℝ) + 1/14) / w := by
  unfold clipLen tooth at hn
  simp only at hn
  have hn' : 0 < min (((n : ℝ) + 1/14) / w) b - max (((n : ℝ) - 1/14) / w) a := by
    by_contra hle
    push Not at hle
    rw [max_eq_left hle] at hn
    exact lt_irrefl 0 hn
  constructor
  · calc ((n : ℝ) - 1/14) / w ≤ max (((n : ℝ) - 1/14) / w) a := le_max_left _ _
      _ < min (((n : ℝ) + 1/14) / w) b := by linarith
      _ ≤ b := min_le_right _ _
  · calc a ≤ max (((n : ℝ) - 1/14) / w) a := le_max_right _ _
      _ < min (((n : ℝ) + 1/14) / w) b := by linarith
      _ ≤ ((n : ℝ) + 1/14) / w := min_le_left _ _

/-- Two different teeth cannot both clip positively into a window shorter than the
inter-tooth gap `6/(7w)`. -/
theorem tooth_clip_disjoint {w : ℤ} (hw : 0 < w) {a b : ℝ}
    (hfar : 7 * (b - a) * (w : ℝ) ≤ 6) {m m' : ℤ} (hlt : m < m')
    (hpos : 0 < clipLen (tooth w m) a b) (hpos' : 0 < clipLen (tooth w m') a b) :
    False := by
  have hwR : (0 : ℝ) < (w : ℝ) := by exact_mod_cast hw
  obtain ⟨-, hra⟩ := tooth_clip_pos_edges hw hpos
  obtain ⟨hlb, -⟩ := tooth_clip_pos_edges hw hpos'
  rw [lt_div_iff₀ hwR] at hra
  rw [div_lt_iff₀ hwR] at hlb
  have hm1 : (m : ℝ) + 1 ≤ (m' : ℝ) := by exact_mod_cast hlt
  have hexp : 7 * (b - a) * (w : ℝ) = 7 * (b * w) - 7 * (a * w) := by ring
  rw [hexp] at hfar
  linarith

/-- **`teeth_mass_far`**: a sub-gap window carries at most one tooth of runner `w` —
clipped mass `≤ 1/(7w)`, with NO density term and NO boundary fee. -/
theorem teeth_mass_far {w : ℤ} (hw : 0 < w) (a b : ℝ)
    (hfar : 7 * (b - a) * (w : ℝ) ≤ 6) :
    ((teeth w a b).map fun p => clipLen p a b).sum ≤ 1 / (7 * (w : ℝ)) := by
  have hwR : (0 : ℝ) < (w : ℝ) := by exact_mod_cast hw
  apply sum_le_single_bound
  · positivity
  · intro x hx
    rw [List.mem_map] at hx
    obtain ⟨p, -, rfl⟩ := hx
    exact clipLen_nonneg p a b
  · intro x hx
    rw [List.mem_map] at hx
    obtain ⟨p, hp, rfl⟩ := hx
    unfold teeth at hp
    rw [List.mem_map] at hp
    obtain ⟨i, -, rfl⟩ := hp
    exact clipLen_tooth_le hw _ a b
  · unfold teeth
    rw [List.map_map]
    refine List.Pairwise.map _ (fun i j hij => ?_) List.pairwise_lt_range
    by_contra hcon
    push Not at hcon
    obtain ⟨h1, h2⟩ := hcon
    have hp1 : 0 < clipLen (tooth w (⌈(w : ℝ) * a⌉ - 1 + (i : ℤ))) a b :=
      lt_of_le_of_ne (clipLen_nonneg _ _ _) (Ne.symm h1)
    have hp2 : 0 < clipLen (tooth w (⌈(w : ℝ) * a⌉ - 1 + (j : ℤ))) a b :=
      lt_of_le_of_ne (clipLen_nonneg _ _ _) (Ne.symm h2)
    exact tooth_clip_disjoint hw hfar (by omega) hp1 hp2

/-- **The window with abstract fees**: any per-top mass bounds summing below the window
length close the family; the S134 criterion and the far-top sharpening are instances. -/
theorem lonely_of_window_of_fees (v : Fin 13 → ℤ) (tops : List (Fin 13))
    (fee : Fin 13 → ℝ) (tstar β B : ℝ)
    (hβ : 1/14 < β)
    (hbase : ∀ i, i ∉ tops → ∀ m : ℤ, β ≤ |(v i : ℝ) * tstar - m|)
    (hB : ∀ i, i ∉ tops → |(v i : ℝ)| ≤ B) (hBpos : 0 < B)
    (hv : ∀ i, v i ≠ 0)
    (hfee : ∀ j ∈ tops,
      rlength (rinter [(tstar - (β - 1/14) / B, tstar + (β - 1/14) / B)]
        (teeth |v j| (tstar - (β - 1/14) / B) (tstar + (β - 1/14) / B))) ≤ fee j)
    (hcrit : (tops.map fee).sum < 2 * ((β - 1/14) / B)) :
    ∃ t : ℝ, Lonely 14 v t := by
  set δ : ℝ := (β - 1/14) / B with hδ
  have hδpos : 0 < δ := div_pos (by linarith) hBpos
  set a : ℝ := tstar - δ with ha
  set b : ℝ := tstar + δ with hb
  have hab : a ≤ b := by rw [ha, hb]; linarith
  have hLab : b - a = 2 * δ := by rw [ha, hb]; ring
  set ws : List ℤ := tops.map (fun j => |v j|) with hws
  have hwpos : ∀ w ∈ ws, 0 < w := by
    intro w hw
    rw [hws, List.mem_map] at hw
    obtain ⟨j, -, rfl⟩ := hw
    exact abs_pos.mpr (hv j)
  have hledger : 0 < (b - a)
      - ((ws.map fun (w : ℤ) => rlength (rinter [(a, b)] (teeth w a b))).sum)
      + pairCredits [(a, b)] (ws.map fun (w : ℤ) => teeth w a b) := by
    have hcred := pairCredits_nonneg (ws.map fun (w : ℤ) => teeth w a b) [(a, b)]
    have hsingles : ((ws.map fun (w : ℤ) => rlength (rinter [(a, b)] (teeth w a b))).sum)
        ≤ (tops.map fee).sum := by
      rw [hws, List.map_map]
      apply List.sum_le_sum
      intro j hj
      simp only [Function.comp_apply]
      exact hfee j hj
    have : (tops.map fee).sum < b - a := by rw [hLab]; exact hcrit
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
      exact ⟨i, hi, rfl⟩
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

/-- **THE SPARSE-TOP MULTI-KILLER WINDOW** (the ~3× fee sharpening, and more): tops in
the SUB-GAP regime (`7·2δ·|v j| ≤ 6`, i.e. the window is shorter than one inter-tooth
gap — the mid-scale killers where teeth_mass's +3 fee was dominant) carry NO density
term — ANY number of them closes with just `Σ 1/(7|v j|) < 2δ`.  The 7-wall vanishes
in this regime; a single such top recovers EXACTLY the S132 sharp threshold
`B/(14(β−1/14))` through the ledger. -/
theorem lonely_of_window_multi_far (v : Fin 13 → ℤ) (tops : List (Fin 13))
    (tstar β B : ℝ)
    (hβ : 1/14 < β)
    (hbase : ∀ i, i ∉ tops → ∀ m : ℤ, β ≤ |(v i : ℝ) * tstar - m|)
    (hB : ∀ i, i ∉ tops → |(v i : ℝ)| ≤ B) (hBpos : 0 < B)
    (hv : ∀ i, v i ≠ 0)
    (hfar : ∀ j ∈ tops, 7 * (2 * ((β - 1/14) / B)) * |(v j : ℝ)| ≤ 6)
    (hcrit : (tops.map fun j => 1 / (7 * |(v j : ℝ)|)).sum < 2 * ((β - 1/14) / B)) :
    ∃ t : ℝ, Lonely 14 v t := by
  apply lonely_of_window_of_fees v tops (fun j => 1 / (7 * |(v j : ℝ)|)) tstar β B
    hβ hbase hB hBpos hv ?_ hcrit
  intro j hj
  have hja : (0 : ℤ) < |v j| := abs_pos.mpr (hv j)
  rw [rlength_inter_window_clipsum]
  have hlen : (tstar + (β - 1/14) / B) - (tstar - (β - 1/14) / B)
      = 2 * ((β - 1/14) / B) := by ring
  have habs : ((|v j| : ℤ) : ℝ) = |(v j : ℝ)| := by
    push_cast
    rfl
  have hfar' : 7 * ((tstar + (β - 1/14) / B) - (tstar - (β - 1/14) / B)) * ((|v j| : ℤ) : ℝ) ≤ 6 := by
    rw [hlen, habs]
    exact hfar j hj
  have := teeth_mass_far hja (tstar - (β - 1/14) / B) (tstar + (β - 1/14) / B) hfar'
  rw [habs] at this
  exact this

#print axioms cite_margin
#print axioms lonely_of_window_multi
#print axioms teeth_mass_far
#print axioms lonely_of_window_of_fees
#print axioms lonely_of_window_multi_far

end MultiKillerWindow
end LonelyRunner
