/-
  TournamentH7.LRCTeethR — THE RADIUS-PARAMETRIC WINDOW STACK (klein-2026-07-05-S136,
  HYP-4107; the Fin-12/ρ=1/13 transcription requested by mac-mini-S52 and kps-S3).

  Everything the window/ledger machinery needs, with the danger radius ρ a PARAMETER
  (0 < ρ < 1/2) instead of the hardcoded 1/14: teeth, mass bounds (density AND sparse),
  the Hunter block step, and the ι-generic multi-top window whose conclusion is a
  MARGIN statement — `Lonely 14` is the ρ = 1/14, ι = Fin 13 instance, and the n = 13
  rigidity level (ρ = 1/13, β = 1/12, Fin 12) is another: mac-mini's fee table
  `T_l = 156(13−l)/(13−2l)` lives here, with the `2l < 13` wall as the ρ-generic
  `(1 − 2lρ) > 0` remainder.

  Sources mirrored: LRCBlockSix (tooth/teeth/mass, 1/14 → ρ), LRCRealRegions
  (hunter_block_step; the region calculus and hunter_ledger are radius-free and
  consumed as-is), LRCMultiKillerWindow (the sparse bound and the fee composition).
-/
import TournamentH7.LRCRealRegions
import TournamentH7.LRCMultiKillerWindow

namespace LonelyRunner
namespace TeethR

open LRC14 RealRegion

/-- One radius-ρ tooth of runner `w` at integer `m`. -/
noncomputable def toothR (ρ : ℝ) (w : ℤ) (m : ℤ) : ℝ × ℝ :=
  (((m : ℝ) - ρ) / w, ((m : ℝ) + ρ) / w)

/-- The radius-ρ teeth of runner `w` meeting the window `[a, b]`. -/
noncomputable def teethR (ρ : ℝ) (w : ℤ) (a b : ℝ) : List (ℝ × ℝ) :=
  (List.range (⌊(w : ℝ) * b⌋ + 1 - (⌈(w : ℝ) * a⌉ - 1) + 1).toNat).map fun (i : ℕ) =>
    toothR ρ w (⌈(w : ℝ) * a⌉ - 1 + (i : ℤ))

theorem teethR_live {ρ : ℝ} (hρ : 0 ≤ ρ) {w : ℤ} (hw : 0 < w) (a b : ℝ) :
    ∀ r ∈ teethR ρ w a b, r.1 ≤ r.2 := by
  intro r hr
  rw [teethR, List.mem_map] at hr
  obtain ⟨i, -, rfl⟩ := hr
  have hwR : (0 : ℝ) < (w : ℝ) := by exact_mod_cast hw
  unfold toothR
  simp only
  gcongr
  linarith

theorem teethR_sortedSep {ρ : ℝ} (hρ : ρ ≤ 1/2) {w : ℤ} (hw : 0 < w) (a b : ℝ) :
    SortedSep (teethR ρ w a b) := by
  have hwR : (0 : ℝ) < (w : ℝ) := by exact_mod_cast hw
  unfold teethR
  apply sortedSep_map_range
  intro i j hij _
  unfold toothR
  simp only
  have hij' : (i : ℝ) + 1 ≤ (j : ℝ) := by exact_mod_cast hij
  gcongr
  push_cast
  linarith

/-- Each clipped radius-ρ tooth carries at most `2ρ/w`. -/
theorem clipLen_toothR_le {ρ : ℝ} (hρ : 0 ≤ ρ) {w : ℤ} (hw : 0 < w) (m : ℤ) (a b : ℝ) :
    clipLen (toothR ρ w m) a b ≤ 2 * ρ / (w : ℝ) := by
  have hwR : (0 : ℝ) < (w : ℝ) := by exact_mod_cast hw
  unfold clipLen toothR
  simp only
  rcases le_total (min (((m : ℝ) + ρ) / w) b - max (((m : ℝ) - ρ) / w) a) 0 with h | h
  · rw [max_eq_left h]
    positivity
  · rw [max_eq_right h]
    have h1 : min (((m : ℝ) + ρ) / w) b ≤ ((m : ℝ) + ρ) / w := min_le_left _ _
    have h2 : ((m : ℝ) - ρ) / w ≤ max (((m : ℝ) - ρ) / w) a := le_max_left _ _
    have h3 : ((m : ℝ) + ρ) / w - ((m : ℝ) - ρ) / w = 2 * ρ / (w : ℝ) := by
      field_simp
      ring
    linarith

/-- **Density mass bound**: the clipped radius-ρ teeth of `w` carry at most
`2ρ(b−a) + 6ρ/w` over the window. -/
theorem teethR_mass {ρ : ℝ} (hρ : 0 ≤ ρ) {w : ℤ} (hw : 0 < w) (a b : ℝ) (hab : a ≤ b) :
    ((teethR ρ w a b).map fun p => clipLen p a b).sum
      ≤ 2 * ρ * (b - a) + 6 * ρ / (w : ℝ) := by
  have hwR : (0 : ℝ) < (w : ℝ) := by exact_mod_cast hw
  have hper : ∀ x ∈ (teethR ρ w a b).map fun p => clipLen p a b,
      x ≤ 2 * ρ / (w : ℝ) := by
    intro x hx
    rw [List.mem_map] at hx
    obtain ⟨p, hp, rfl⟩ := hx
    rw [teethR, List.mem_map] at hp
    obtain ⟨i, -, rfl⟩ := hp
    exact clipLen_toothR_le hρ hw _ a b
  have hlen : ((teethR ρ w a b).map fun p => clipLen p a b).length
      = (⌊(w : ℝ) * b⌋ + 1 - (⌈(w : ℝ) * a⌉ - 1) + 1).toNat := by
    rw [List.length_map, teethR, List.length_map, List.length_range]
  have hsum : ((teethR ρ w a b).map fun p => clipLen p a b).sum
      ≤ (((⌊(w : ℝ) * b⌋ + 1 - (⌈(w : ℝ) * a⌉ - 1) + 1).toNat : ℝ)) * (2 * ρ / (w : ℝ)) := by
    rw [← hlen]
    exact List.sum_le_card_nsmul _ _ hper |>.trans (by rw [nsmul_eq_mul])
  have hcount : (((⌊(w : ℝ) * b⌋ + 1 - (⌈(w : ℝ) * a⌉ - 1) + 1).toNat : ℝ))
      ≤ (w : ℝ) * (b - a) + 3 := by
    have h1 : ((⌊(w : ℝ) * b⌋ : ℤ) : ℝ) ≤ (w : ℝ) * b := Int.floor_le _
    have h2 : (w : ℝ) * a ≤ ((⌈(w : ℝ) * a⌉ : ℤ) : ℝ) := Int.le_ceil _
    rcases le_or_gt (⌊(w : ℝ) * b⌋ + 1 - (⌈(w : ℝ) * a⌉ - 1) + 1) 0 with h | h
    · rw [Int.toNat_of_nonpos h]
      norm_num
      nlinarith [mul_nonneg hwR.le (sub_nonneg.mpr hab)]
    · have hchain : (((⌊(w : ℝ) * b⌋ + 1 - (⌈(w : ℝ) * a⌉ - 1) + 1).toNat : ℕ) : ℝ)
          = ((⌊(w : ℝ) * b⌋ + 1 - (⌈(w : ℝ) * a⌉ - 1) + 1 : ℤ) : ℝ) := by
        rw [← Int.cast_natCast, Int.toNat_of_nonneg (le_of_lt h)]
      rw [hchain]
      push_cast
      have hexp : (w : ℝ) * (b - a) = (w : ℝ) * b - (w : ℝ) * a := by ring
      linarith
  calc ((teethR ρ w a b).map fun p => clipLen p a b).sum
      ≤ _ * (2 * ρ / (w : ℝ)) := hsum
    _ ≤ ((w : ℝ) * (b - a) + 3) * (2 * ρ / (w : ℝ)) := by
        apply mul_le_mul_of_nonneg_right hcount
        positivity
    _ = 2 * ρ * (b - a) + 6 * ρ / (w : ℝ) := by
        field_simp
        ring

/-! ## The sparse bound at radius ρ (mirror of S135's teeth_mass_far) -/

/-- Positive clip forces the ρ-tooth's edges to straddle the window. -/
theorem toothR_clip_pos_edges {ρ : ℝ} {w : ℤ} (hw : 0 < w) {a b : ℝ} {n : ℤ}
    (hn : 0 < clipLen (toothR ρ w n) a b) :
    ((n : ℝ) - ρ) / w < b ∧ a < ((n : ℝ) + ρ) / w := by
  unfold clipLen toothR at hn
  simp only at hn
  have hn' : 0 < min (((n : ℝ) + ρ) / w) b - max (((n : ℝ) - ρ) / w) a := by
    by_contra hle
    push Not at hle
    rw [max_eq_left hle] at hn
    exact lt_irrefl 0 hn
  constructor
  · calc ((n : ℝ) - ρ) / w ≤ max (((n : ℝ) - ρ) / w) a := le_max_left _ _
      _ < min (((n : ℝ) + ρ) / w) b := by linarith
      _ ≤ b := min_le_right _ _
  · calc a ≤ max (((n : ℝ) - ρ) / w) a := le_max_right _ _
      _ < min (((n : ℝ) + ρ) / w) b := by linarith
      _ ≤ ((n : ℝ) + ρ) / w := min_le_left _ _

/-- Two different ρ-teeth cannot both clip positively into a window shorter than the
inter-tooth gap `(1 − 2ρ)/w`. -/
theorem toothR_clip_disjoint {ρ : ℝ} (hρ : ρ ≤ 1/2) {w : ℤ} (hw : 0 < w) {a b : ℝ}
    (hfar : (b - a) * (w : ℝ) ≤ 1 - 2 * ρ) {m m' : ℤ} (hlt : m < m')
    (hpos : 0 < clipLen (toothR ρ w m) a b) (hpos' : 0 < clipLen (toothR ρ w m') a b) :
    False := by
  have hwR : (0 : ℝ) < (w : ℝ) := by exact_mod_cast hw
  obtain ⟨-, hra⟩ := toothR_clip_pos_edges hw hpos
  obtain ⟨hlb, -⟩ := toothR_clip_pos_edges hw hpos'
  rw [lt_div_iff₀ hwR] at hra
  rw [div_lt_iff₀ hwR] at hlb
  have hm1 : (m : ℝ) + 1 ≤ (m' : ℝ) := by exact_mod_cast hlt
  have hexp : (b - a) * (w : ℝ) = b * w - a * w := by ring
  rw [hexp] at hfar
  linarith

/-- **Sparse mass bound at radius ρ**: a sub-gap window carries at most one ρ-tooth. -/
theorem teethR_mass_far {ρ : ℝ} (hρ : 0 ≤ ρ) (hρ2 : ρ ≤ 1/2) {w : ℤ} (hw : 0 < w)
    (a b : ℝ) (hfar : (b - a) * (w : ℝ) ≤ 1 - 2 * ρ) :
    ((teethR ρ w a b).map fun p => clipLen p a b).sum ≤ 2 * ρ / (w : ℝ) := by
  have hwR : (0 : ℝ) < (w : ℝ) := by exact_mod_cast hw
  apply MultiKillerWindow.sum_le_single_bound
  · positivity
  · intro x hx
    rw [List.mem_map] at hx
    obtain ⟨p, -, rfl⟩ := hx
    exact clipLen_nonneg p a b
  · intro x hx
    rw [List.mem_map] at hx
    obtain ⟨p, hp, rfl⟩ := hx
    rw [teethR, List.mem_map] at hp
    obtain ⟨i, -, rfl⟩ := hp
    exact clipLen_toothR_le hρ hw _ a b
  · unfold teethR
    rw [List.map_map]
    refine List.Pairwise.map _ (fun i j hij => ?_) List.pairwise_lt_range
    by_contra hcon
    push Not at hcon
    obtain ⟨h1, h2⟩ := hcon
    have hp1 : 0 < clipLen (toothR ρ w (⌈(w : ℝ) * a⌉ - 1 + (i : ℤ))) a b :=
      lt_of_le_of_ne (clipLen_nonneg _ _ _) (Ne.symm h1)
    have hp2 : 0 < clipLen (toothR ρ w (⌈(w : ℝ) * a⌉ - 1 + (j : ℤ))) a b :=
      lt_of_le_of_ne (clipLen_nonneg _ _ _) (Ne.symm h2)
    exact toothR_clip_disjoint hρ2 hw hfar (by omega) hp1 hp2

/-! ## The avoidance bridge at radius ρ (mirror of good_of_avoid_teeth) -/

/-- A point of the window avoiding every listed ρ-tooth is ρ-far from every integer
mark of the runner. -/
theorem good_of_avoid_teethR {ρ : ℝ} (hρ : 0 < ρ) (hρ2 : ρ ≤ 1/2) {w : ℤ} (hw : 0 < w)
    {a b t : ℝ} (hta : a ≤ t) (htb : t ≤ b)
    (havoid : ∀ p ∈ teethR ρ w a b, t ≤ p.1 ∨ p.2 ≤ t) :
    ∀ m : ℤ, ρ ≤ |(w : ℝ) * t - m| := by
  intro m
  have hwR : (0 : ℝ) < (w : ℝ) := by exact_mod_cast hw
  set mlo : ℤ := ⌈(w : ℝ) * a⌉ - 1 with hmlo
  set mhi : ℤ := ⌊(w : ℝ) * b⌋ + 1 with hmhi
  have hmloR : (mlo : ℝ) < (w : ℝ) * a := by
    rw [hmlo]
    push_cast
    linarith [Int.ceil_lt_add_one ((w : ℝ) * a)]
  have hmhiR : (w : ℝ) * b < (mhi : ℝ) := by
    rw [hmhi]
    push_cast
    linarith [Int.lt_floor_add_one ((w : ℝ) * b)]
  rcases lt_or_ge m mlo with hm | hm
  · have h2 : (w : ℝ) * a ≤ (w : ℝ) * t := mul_le_mul_of_nonneg_left hta hwR.le
    have h3 : (m : ℝ) ≤ (mlo : ℝ) - 1 := by
      have hmz : m ≤ mlo - 1 := by omega
      calc (m : ℝ) ≤ ((mlo - 1 : ℤ) : ℝ) := by exact_mod_cast hmz
        _ = (mlo : ℝ) - 1 := by push_cast; ring
    rw [abs_of_nonneg (by linarith)]
    linarith
  · rcases le_or_gt m mhi with hm2 | hm2
    · have hmem : toothR ρ w m ∈ teethR ρ w a b := by
        unfold teethR
        rw [List.mem_map]
        refine ⟨(m - mlo).toNat, ?_, ?_⟩
        · rw [List.mem_range]
          rw [hmlo] at hm ⊢
          rw [hmhi] at hm2
          omega
        · congr 1
          rw [hmlo] at hm ⊢
          omega
      rcases havoid _ hmem with hleft | hright
      · unfold toothR at hleft
        simp only at hleft
        rw [le_div_iff₀ hwR] at hleft
        rw [abs_of_nonpos (by nlinarith)]
        nlinarith
      · unfold toothR at hright
        simp only at hright
        rw [div_le_iff₀ hwR] at hright
        rw [abs_of_nonneg (by nlinarith)]
        nlinarith
    · have h2 : (w : ℝ) * t ≤ (w : ℝ) * b := mul_le_mul_of_nonneg_left htb hwR.le
      have h3 : (mhi : ℝ) + 1 ≤ (m : ℝ) := by
        have hmz : mhi + 1 ≤ m := by omega
        calc (mhi : ℝ) + 1 = ((mhi + 1 : ℤ) : ℝ) := by push_cast; ring
          _ ≤ (m : ℝ) := by exact_mod_cast hmz
      rw [abs_of_nonpos (by linarith)]
      linarith

/-! ## The Hunter block step at radius ρ -/

/-- **The ρ-parametric Hunter block step**: positive zero-or-more-credit ledger over the
window ⟹ a common ρ-good point for the whole block.  The region calculus and the
Hunter ledger are radius-free; only the teeth and the avoidance bridge carry ρ. -/
theorem hunter_block_stepR {ρ : ℝ} (hρ : 0 < ρ) (hρ2 : ρ ≤ 1/2)
    (ws : List ℤ) (hpos : ∀ w ∈ ws, 0 < w) (a b : ℝ) (hab : a ≤ b)
    (hledger : 0 < (b - a)
        - ((ws.map fun (w : ℤ) => rlength (rinter [(a, b)] (teethR ρ w a b))).sum)
        + pairCredits [(a, b)] (ws.map fun (w : ℤ) => teethR ρ w a b)) :
    ∃ t : ℝ, a ≤ t ∧ t ≤ b ∧
      ∀ w ∈ ws, ∀ m : ℤ, ρ ≤ |(w : ℝ) * t - m| := by
  set Ds : List RRegion := ws.map fun (w : ℤ) => teethR ρ w a b with hDs
  have hlive : ∀ D ∈ Ds, ∀ r ∈ D, r.1 ≤ r.2 := by
    intro D hD
    rw [hDs, List.mem_map] at hD
    obtain ⟨w, hw, rfl⟩ := hD
    exact teethR_live (le_of_lt hρ) (hpos w hw) a b
  have hsep : ∀ D ∈ Ds, SortedSep D := by
    intro D hD
    rw [hDs, List.mem_map] at hD
    obtain ⟨w, hw, rfl⟩ := hD
    exact teethR_sortedSep hρ2 (hpos w hw) a b
  have hIlen : rlength [(a, b)] = b - a := by
    unfold rlength
    simp only [List.map_cons, List.map_nil, List.sum_cons, List.sum_nil, add_zero]
    exact max_eq_right (by linarith)
  have hled := hunter_ledger Ds hlive hsep [(a, b)]
  have hsum_eq : (Ds.map fun D => rlength (rinter [(a, b)] D)).sum
      = (ws.map fun (w : ℤ) => rlength (rinter [(a, b)] (teethR ρ w a b))).sum := by
    rw [hDs, List.map_map]
    rfl
  have hpos' : 0 < rlength (Ds.foldl rdiff [(a, b)]) := by
    rw [hIlen, hsum_eq] at hled
    linarith
  obtain ⟨p', hp'mem, hp'lt⟩ := exists_pos_interval_of_rlength_pos hpos'
  obtain ⟨⟨p, hpI, hpt1, hpt2⟩, havoid⟩ :=
    rdiff_chain_point_good Ds [(a, b)] p' hp'mem p'.1 (le_refl _) (le_of_lt hp'lt)
  have hpab : p = (a, b) := by
    rcases List.mem_cons.mp hpI with rfl | hf
    · rfl
    · exact absurd hf List.not_mem_nil
  subst hpab
  simp only at hpt1 hpt2
  refine ⟨p'.1, hpt1, hpt2, ?_⟩
  intro w hw m
  apply good_of_avoid_teethR hρ hρ2 (hpos w hw) hpt1 hpt2
  intro r hr
  exact havoid (teethR ρ w a b) (by rw [hDs, List.mem_map]; exact ⟨w, hw, rfl⟩) r hr

/-! ## The ι-generic margin window -/

/-- **THE ρ-PARAMETRIC MULTI-TOP WINDOW** (ι-generic, margin-form conclusion): a base
with margin β > ρ at one point plus tops whose per-top mass fees sum below the window
give a common ρ-margin point for the WHOLE family.  `Lonely 14` is the
`ι = Fin 13, ρ = 1/14` instance; the n = 13 rigidity level is `ι = Fin 12, ρ = 1/13,
β = 1/12`.  Fees are supplied by `teethR_mass` (density) or `teethR_mass_far`
(sparse) per top. -/
theorem margin_of_window_of_fees {ι : Type*} (v : ι → ℤ)
    (tops : List ι) (fee : ι → ℝ) (ρ tstar β B : ℝ)
    (hρ : 0 < ρ) (hρ2 : ρ ≤ 1/2) (hβ : ρ < β)
    (hbase : ∀ i, i ∉ tops → ∀ m : ℤ, β ≤ |(v i : ℝ) * tstar - m|)
    (hB : ∀ i, i ∉ tops → |(v i : ℝ)| ≤ B) (hBpos : 0 < B)
    (hv : ∀ j ∈ tops, v j ≠ 0)
    (hfee : ∀ j ∈ tops,
      rlength (rinter [(tstar - (β - ρ) / B, tstar + (β - ρ) / B)]
        (teethR ρ |v j| (tstar - (β - ρ) / B) (tstar + (β - ρ) / B))) ≤ fee j)
    (hcrit : (tops.map fee).sum < 2 * ((β - ρ) / B)) :
    ∃ t : ℝ, ∀ i, ∀ m : ℤ, ρ ≤ |(v i : ℝ) * t - m| := by
  set δ : ℝ := (β - ρ) / B with hδ
  have hδpos : 0 < δ := div_pos (by linarith) hBpos
  set a : ℝ := tstar - δ with ha
  set b : ℝ := tstar + δ with hb
  have hab : a ≤ b := by rw [ha, hb]; linarith
  have hLab : b - a = 2 * δ := by rw [ha, hb]; ring
  set ws : List ℤ := tops.map (fun j => |v j|) with hws
  have hwpos : ∀ w ∈ ws, 0 < w := by
    intro w hw
    rw [hws, List.mem_map] at hw
    obtain ⟨j, hj, rfl⟩ := hw
    exact abs_pos.mpr (hv j hj)
  have hledger : 0 < (b - a)
      - ((ws.map fun (w : ℤ) => rlength (rinter [(a, b)] (teethR ρ w a b))).sum)
      + pairCredits [(a, b)] (ws.map fun (w : ℤ) => teethR ρ w a b) := by
    have hcred := pairCredits_nonneg (ws.map fun (w : ℤ) => teethR ρ w a b) [(a, b)]
    have hsingles : ((ws.map fun (w : ℤ) => rlength (rinter [(a, b)] (teethR ρ w a b))).sum)
        ≤ (tops.map fee).sum := by
      rw [hws, List.map_map]
      apply List.sum_le_sum
      intro j hj
      simp only [Function.comp_apply]
      exact hfee j hj
    have hlt : (tops.map fee).sum < b - a := by rw [hLab]; exact hcrit
    linarith
  obtain ⟨t, hta, htb, hgood⟩ := hunter_block_stepR hρ hρ2 ws hwpos a b hab hledger
  have htw : |t - tstar| ≤ δ := by
    rw [abs_le]
    constructor
    · rw [ha] at hta; linarith
    · rw [hb] at htb; linarith
  refine ⟨t, fun i m => ?_⟩
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
    have hBδ : B * δ = β - ρ := by
      rw [hδ]
      field_simp
    linarith

#print axioms teethR_mass
/-! ## The citation margin, ι-generic, and the level-12 tower step -/

/-- **The citation margin, ι-generic**: any ≤ 12 runners have a common point with
margin `1/(card + 1)`, from the LRC(≤13) citation node. -/
theorem cite_margin_gen (cite : LRCUpTo13) {ι : Type*} [DecidableEq ι] (v : ι → ℤ)
    (T : Finset ι) (hT : T.card ≤ 12) (hv : ∀ i ∈ T, v i ≠ 0) :
    ∃ t : ℝ, ∀ i ∈ T, ∀ m : ℤ, (1 : ℝ) / (T.card + 1) ≤ |(v i : ℝ) * t - m| := by
  let e : Fin T.card ≃ T := T.equivFin.symm
  set w : Fin T.card → ℤ := fun m => v (e m : ι) with hw
  have hwne : ∀ m, w m ≠ 0 := fun m => hv _ (e m).2
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
  linarith [hstep]

/-- **THE LEVEL-12 TOWER STEP** (the transcription mac-mini-S52/kps-S3 requested):
a 12-family whose `l ≥ 1` tops clear the fee criterion has a common `1/13`-margin
point — the sub-base of `12 − l` runners is cited at `1/(13−l) > 1/13`, and the tops
ride the ρ = 1/13 window.  Fees per top from `teethR_mass`/`teethR_mass_far`;
mac-mini's `T_l` thresholds are instances of `hcrit`. -/
theorem tower_step_12 (cite : LRCUpTo13) (v : Fin 12 → ℤ) (hv : ∀ i, v i ≠ 0)
    (tops : List (Fin 12)) (fee : Fin 12 → ℝ) (B : ℝ) (hBpos : 0 < B)
    (hl1 : tops ≠ [])
    (hB : ∀ i, i ∉ tops → |(v i : ℝ)| ≤ B)
    (hfee : ∀ (tstar : ℝ), ∀ j ∈ tops,
      rlength (rinter
        [(tstar - (1/(13 - (tops.toFinset.card : ℝ)) - 1/13) / B,
          tstar + (1/(13 - (tops.toFinset.card : ℝ)) - 1/13) / B)]
        (teethR (1/13) |v j|
          (tstar - (1/(13 - (tops.toFinset.card : ℝ)) - 1/13) / B)
          (tstar + (1/(13 - (tops.toFinset.card : ℝ)) - 1/13) / B))) ≤ fee j)
    (hcrit : (tops.map fee).sum
      < 2 * ((1/(13 - (tops.toFinset.card : ℝ)) - 1/13) / B)) :
    ∃ t : ℝ, ∀ i, ∀ m : ℤ, (1 : ℝ) / 13 ≤ |(v i : ℝ) * t - m| := by
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
  have hβgt : (1 : ℝ) / 13 < β := by
    rw [hβdef]
    apply div_lt_div_of_pos_left (by norm_num) (by linarith) (by linarith)
  have hbase : ∀ i, i ∉ tops → ∀ m : ℤ, β ≤ |(v i : ℝ) * tstar - m| := by
    intro i hi m
    have hiT : i ∈ T := by
      rw [hT, Finset.mem_sdiff]
      exact ⟨Finset.mem_univ i, fun hc => hi (List.mem_toFinset.mp hc)⟩
    have := hmargin i hiT m
    have hcast : (T.card : ℝ) + 1 = 13 - (tops.toFinset.card : ℝ) := by
      rw [hTcard]
      push_cast [Nat.cast_sub htops12]  -- guard; adjusted below
      ring
    rw [hβdef, ← hcast]
    exact this
  exact margin_of_window_of_fees v tops fee (1/13) tstar β B
    (by norm_num) (by norm_num) hβgt hbase hB hBpos (fun j _ => hv j) (hfee tstar) hcrit

/-! ## The 1/14 bridge: the legacy stack is definitionally the ρ = 1/14 instance -/

theorem toothR_fourteen (w m : ℤ) : toothR (1/14) w m = tooth w m := rfl

theorem teethR_fourteen (w : ℤ) (a b : ℝ) : teethR (1/14) w a b = teeth w a b := rfl

/-- `Lonely 14` from the parametric window: the ρ = 1/14, `Fin 13` instance of
`margin_of_window_of_fees` (subsumes the S134/S135 window lemmas). -/
theorem lonely14_of_margin_window (v : Fin 13 → ℤ) (tops : List (Fin 13))
    (fee : Fin 13 → ℝ) (tstar β B : ℝ)
    (hβ : 1/14 < β)
    (hbase : ∀ i, i ∉ tops → ∀ m : ℤ, β ≤ |(v i : ℝ) * tstar - m|)
    (hB : ∀ i, i ∉ tops → |(v i : ℝ)| ≤ B) (hBpos : 0 < B)
    (hv : ∀ j ∈ tops, v j ≠ 0)
    (hfee : ∀ j ∈ tops,
      rlength (rinter [(tstar - (β - 1/14) / B, tstar + (β - 1/14) / B)]
        (teethR (1/14) |v j| (tstar - (β - 1/14) / B) (tstar + (β - 1/14) / B))) ≤ fee j)
    (hcrit : (tops.map fee).sum < 2 * ((β - 1/14) / B)) :
    ∃ t : ℝ, Lonely 14 v t := by
  obtain ⟨t, ht⟩ := margin_of_window_of_fees v tops fee (1/14) tstar β B
    (by norm_num) (by norm_num) hβ hbase hB hBpos hv hfee hcrit
  refine ⟨t, fun i m => ?_⟩
  show (1 : ℝ) / (14 : ℕ) ≤ _
  push_cast
  exact ht i m

#print axioms teethR_mass
