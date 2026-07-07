/-
  TournamentH7.LRCLadderD1  (mac-mini-2026-07-06-S33, HYP-4632/opus-S123)

  THE d=1 LADDER BOUND, formalized — one of the two remaining pieces of (C).

  opus-S123: (C) at N=12 reduces to the d=1 and d=2 defect strata (d>=3 GREEN via
  kps's mod-25 floor).  The d=1 stratum is the canonical family

        V(x) = {1, 2, …, 11} ∪ {x}          (11-term AP + one outlier)

  and the ladder bound is: for x ≠ 12,  M(V(x)) ≥ 2/25  (so V(x) is not in the
  open gap (1/13, 2/25); x = 12 completes the AP, giving M = 1/13).

  TWO-WITNESS PROOF (verified, lrc_d1_ladder_witnesses_macmini_S33):
    · 12 ∤ x        : witness t = 1/12,  every speed clears 1/12 ≥ 2/25.
    · x = 12k, k≥2  : witness t = k/(12k+1),  every speed clears k/(12k+1) ≥ 2/25.
  Since every positive x ≠ 12 falls in one case, M(V(x)) ≥ 2/25 throughout.

  This file formalizes both witnesses (kernel-pure), giving
  `d1_reach_ge : x ≠ 12 → 2/25 ≤ sSup (margin (V x) '' [0,1])` — the d=1 half of
  the AP-rigidity closure `M < 2/25 ⟹ AP`.
-/
import TournamentH7.LRCWitnessAttainment

open TournamentH7.LRCWitness

namespace TournamentH7.LRCLadderD1

/-- The canonical d=1 family `{1,…,11, x}` as `Fin 12 → ℤ`. -/
def V (x : ℤ) : Fin 12 → ℤ := fun i =>
  if h : i.val = 11 then x else (i.val : ℤ) + 1

lemma V_last (x : ℤ) : V x 11 = x := by simp [V]

/-- `distZ y ≤ 1/2`. -/
lemma distZ_le_half (y : ℝ) : distZ y ≤ 1 / 2 := by
  have hmem : ((round y : ℤ) : ℝ) ∈ Set.range ((↑) : ℤ → ℝ) := ⟨round y, rfl⟩
  calc distZ y ≤ dist y ((round y : ℤ) : ℝ) := Metric.infDist_le_dist_of_mem hmem
    _ = |y - (round y : ℝ)| := Real.dist_eq _ _
    _ ≤ 1 / 2 := abs_sub_round y

lemma margin_le_half (v : Fin 12 → ℤ) (t : ℝ) : margin v t ≤ 1 / 2 :=
  le_trans (Finset.inf'_le _ (Finset.mem_univ (0 : Fin 12))) (distZ_le_half _)

/-- **Rational → real bridge** at denominator `q`, numerator threshold `c`.
If `c ≤ |v·b − m·q|` for all integers `m`, then `c/q ≤ |v·(b/q) − m|`. -/
lemma bridge (q b c v : ℤ) (hq : 0 < q) (hb : ∀ m : ℤ, c ≤ |v * b - m * q|) (m : ℤ) :
    (c : ℝ) / q ≤ |(v : ℝ) * ((b : ℝ) / q) - (m : ℝ)| := by
  have key : (c : ℝ) ≤ |(v : ℝ) * b - (m : ℝ) * q| := by
    have h : (c : ℝ) ≤ ((|v * b - m * q| : ℤ) : ℝ) := by exact_mod_cast hb m
    rw [Int.cast_abs] at h; push_cast at h; exact h
  have hqR : (0 : ℝ) < q := by exact_mod_cast hq
  have hfac : (v : ℝ) * ((b : ℝ) / q) - (m : ℝ) = ((v : ℝ) * b - (m : ℝ) * q) / q := by
    field_simp
  rw [hfac, abs_div, abs_of_pos hqR, le_div_iff₀ hqR, div_mul_cancel₀]
  · exact key
  · exact ne_of_gt hqR

/-! ### Case 1 — the generic witness `t = 1/12` (12 ∤ x) -/

/-- If `12 ∤ w` then `1 ≤ |w - m·12|` for every `m` (nonzero integer has abs ≥ 1). -/
lemma one_le_abs_of_not_dvd {w : ℤ} (hw : ¬ (12 ∣ w)) (m : ℤ) : (1 : ℤ) ≤ |w - m * 12| := by
  have hne : w - m * 12 ≠ 0 := by
    intro h; apply hw; exact ⟨m, by omega⟩
  exact Int.one_le_abs hne

/-- **Generic d=1 witness.**  If `12 ∤ x`, every speed of `V x` is `≥ 1/12` from ℤ at
`t = 1/12`; so `margin (V x) (1/12) ≥ 1/12`. -/
theorem d1_generic_margin (x : ℤ) (hx : ¬ (12 ∣ x)) :
    (1 : ℝ) / 12 ≤ margin (V x) (1 / 12) := by
  rw [le_margin_iff]
  intro i m
  -- every speed w = V x i satisfies 12 ∤ w, so the bridge (q=12,b=1,c=1) applies
  have hbridge : ∀ w : ℤ, ¬ (12 ∣ w) → (1 : ℝ) / 12 ≤ |(w : ℝ) * (1 / 12 : ℝ) - (m : ℝ)| := by
    intro w hw
    have := bridge 12 1 1 w (by norm_num) (fun m => by
      have := one_le_abs_of_not_dvd hw m; simpa using this) m
    simpa using this
  have hval : ¬ (12 ∣ V x i) := by
    have hi12 : i.val < 12 := i.isLt
    unfold V
    split
    · exact hx
    · rename_i h
      -- V x i = i.val + 1 ∈ {1,…,11}, not a multiple of 12
      rintro ⟨c, hc⟩; omega
  exact hbridge (V x i) hval

/-- Generic reach bound: `12 ∤ x ⟹ reach(V x) ≥ 1/12 > 2/25`. -/
theorem d1_generic_reach (x : ℤ) (hx : ¬ (12 ∣ x)) :
    (2 : ℝ) / 25 < sSup (margin (V x) '' Set.Icc (0 : ℝ) 1) := by
  have hmem : margin (V x) (1 / 12) ∈ margin (V x) '' Set.Icc (0 : ℝ) 1 :=
    ⟨1 / 12, Set.mem_Icc.mpr ⟨by norm_num, by norm_num⟩, rfl⟩
  have hbdd : BddAbove (margin (V x) '' Set.Icc (0 : ℝ) 1) :=
    ⟨1 / 2, by rintro _ ⟨t, _, rfl⟩; exact margin_le_half (V x) t⟩
  have hle : margin (V x) (1 / 12) ≤ sSup (margin (V x) '' Set.Icc (0 : ℝ) 1) :=
    le_csSup hbdd hmem
  have := d1_generic_margin x hx
  have h2 : (2 : ℝ) / 25 < (1 : ℝ) / 12 := by norm_num
  linarith

/-! ### Case 2 — the resonant witness `t = k/(12k+1)` (x = 12k, k ≥ 2)

The integer core: for `1 ≤ i ≤ 11` and `k ≥ 1`, `k ≤ |i·k − m·(12k+1)|` for all `m`
(sign split on `m`); and the outlier `12k` gives `12k·k ≡ −k (mod 12k+1)`. -/

/-- Integer clearance for the small speeds: `1 ≤ i ≤ 11, 1 ≤ k ⟹ k ≤ |i·k − m·(12k+1)|`. -/
lemma small_speed_clear (i k m : ℤ) (hi1 : 1 ≤ i) (hi11 : i ≤ 11) (hk : 1 ≤ k) :
    k ≤ |i * k - m * (12 * k + 1)| := by
  by_cases hm : m ≤ 0
  · -- m ≤ 0 : i*k - m*(12k+1) = i*k + (-m)(12k+1) ≥ i*k ≥ k, positive
    have h1 : (0 : ℤ) ≤ (i - 1) * k := mul_nonneg (by linarith) (by linarith)
    have h2 : (0 : ℤ) ≤ (-m) * (12 * k + 1) := mul_nonneg (by linarith) (by linarith)
    rw [abs_of_nonneg (by nlinarith)]; nlinarith
  · -- m ≥ 1 : i*k - m*(12k+1) ≤ 11k - (12k+1) < 0, abs = m*(12k+1) - i*k ≥ k+1
    push_neg at hm
    have h1 : (0 : ℤ) ≤ (11 - i) * k := mul_nonneg (by linarith) (by linarith)
    have h2 : (0 : ℤ) ≤ (m - 1) * (12 * k + 1) := mul_nonneg (by linarith) (by linarith)
    rw [abs_of_neg (by nlinarith)]; nlinarith

/-- The outlier `12k` clears: `k ≤ |12k·k − m·(12k+1)|` (since `12k² ≡ −k`, min at `m=k`). -/
lemma outlier_clear (k m : ℤ) (hk : 1 ≤ k) :
    k ≤ |12 * k * k - m * (12 * k + 1)| := by
  -- 12k²·1 − m(12k+1) = (k−m)(12k+1) − k, so A + k = (k−m)(12k+1)
  rcases lt_trichotomy m k with h | h | h
  · -- m < k : (k−m) ≥ 1 ⟹ A ≥ 11k+1 > 0
    have h1 : (12 * k + 1) ≤ (k - m) * (12 * k + 1) :=
      le_mul_of_one_le_left (by linarith) (by linarith)
    rw [abs_of_nonneg (by nlinarith)]; nlinarith
  · -- m = k : A = −k
    rw [h, show 12 * k * k - k * (12 * k + 1) = -k by ring, abs_neg, abs_of_nonneg (by linarith)]
  · -- m > k : (m−k) ≥ 1 ⟹ A ≤ −13k−1 < 0
    have h1 : (12 * k + 1) ≤ (m - k) * (12 * k + 1) :=
      le_mul_of_one_le_left (by linarith) (by linarith)
    rw [abs_of_neg (by nlinarith)]; nlinarith

/-- **Resonant d=1 witness.**  For `k ≥ 2`, every speed of `V (12k)` is `≥ 2/25` from ℤ
at `t = k/(12k+1)`; so `margin (V (12k)) (k/(12k+1)) ≥ 2/25`. -/
theorem d1_resonant_margin (k : ℤ) (hk : 2 ≤ k) :
    (2 : ℝ) / 25 ≤ margin (V (12 * k)) ((k : ℝ) / (12 * k + 1)) := by
  have hk1 : (1 : ℤ) ≤ k := by omega
  have hqpos : (0 : ℤ) < 12 * k + 1 := by omega
  -- first: k/(12k+1) is the clearance; and 2/25 ≤ k/(12k+1) for k ≥ 2
  have hclear : (k : ℝ) / (12 * k + 1) ≤ margin (V (12 * k)) ((k : ℝ) / (12 * k + 1)) := by
    rw [le_margin_iff]
    intro i m
    have hi12 : i.val < 12 := i.isLt
    by_cases hlast : i.val = 11
    · -- outlier speed = 12k
      have hv : V (12 * k) i = 12 * k := by simp [V, hlast]
      rw [hv]
      have := bridge (12 * k + 1) k k (12 * k) hqpos
        (fun m => outlier_clear k m hk1) m
      simpa using this
    · -- small speed = i.val + 1 ∈ {1,…,11}
      have hv : V (12 * k) i = (i.val : ℤ) + 1 := by simp [V, hlast]
      have hi11 : (i.val : ℤ) + 1 ≤ 11 := by
        have h10 : i.val ≤ 10 := by omega
        have : (i.val : ℤ) ≤ 10 := by exact_mod_cast h10
        linarith
      have hi1 : (1 : ℤ) ≤ (i.val : ℤ) + 1 := by
        have : (0 : ℤ) ≤ (i.val : ℤ) := Int.natCast_nonneg _
        linarith
      rw [hv]
      have := bridge (12 * k + 1) k k ((i.val : ℤ) + 1) hqpos
        (fun m => small_speed_clear ((i.val : ℤ) + 1) k m hi1 hi11 hk1) m
      simpa using this
  have hkR : (2 : ℝ) ≤ (k : ℝ) := by exact_mod_cast hk
  have hqR : (0 : ℝ) < 12 * (k : ℝ) + 1 := by linarith
  have hge : (2 : ℝ) / 25 ≤ (k : ℝ) / (12 * k + 1) := by
    rw [le_div_iff₀ hqR]
    have : (2 : ℝ) / 25 * (12 * k + 1) ≤ k := by nlinarith
    linarith
  linarith

/-- Resonant reach bound: `k ≥ 2 ⟹ reach(V (12k)) ≥ 2/25`. -/
theorem d1_resonant_reach (k : ℤ) (hk : 2 ≤ k) :
    (2 : ℝ) / 25 ≤ sSup (margin (V (12 * k)) '' Set.Icc (0 : ℝ) 1) := by
  have hqpos : (0 : ℝ) < 12 * (k : ℝ) + 1 := by
    have : (2 : ℝ) ≤ (k : ℝ) := by exact_mod_cast hk
    linarith
  have ht01 : (k : ℝ) / (12 * k + 1) ∈ Set.Icc (0 : ℝ) 1 := by
    have hkpos : (0 : ℝ) ≤ (k : ℝ) := by positivity
    refine Set.mem_Icc.mpr ⟨by positivity, ?_⟩
    rw [div_le_one hqpos]; linarith
  have hmem : margin (V (12 * k)) ((k : ℝ) / (12 * k + 1)) ∈
      margin (V (12 * k)) '' Set.Icc (0 : ℝ) 1 := ⟨_, ht01, rfl⟩
  have hbdd : BddAbove (margin (V (12 * k)) '' Set.Icc (0 : ℝ) 1) :=
    ⟨1 / 2, by rintro _ ⟨t, _, rfl⟩; exact margin_le_half (V (12 * k)) t⟩
  exact le_trans (d1_resonant_margin k hk) (le_csSup hbdd hmem)

/-- **The d=1 ladder bound (headline).**  For every positive `x ≠ 12`, the canonical
d=1 family `V x = {1,…,11, x}` has reach `≥ 2/25` — so it is NOT in the open gap
`(1/13, 2/25)`.  (`x = 12` completes the AP, giving reach `1/13`, the gap's lower
edge.)  Combined with kps's mod-25 floor (`d ≥ 3`) and the d=2 bound, this is a
piece of the AP-rigidity closure `M < 2/25 ⟹ AP` = (C). -/
theorem d1_reach_ge (x : ℤ) (hx0 : 0 < x) (hx12 : x ≠ 12) :
    (2 : ℝ) / 25 ≤ sSup (margin (V x) '' Set.Icc (0 : ℝ) 1) := by
  by_cases h : (12 : ℤ) ∣ x
  · obtain ⟨k, rfl⟩ := h
    have hk2 : 2 ≤ k := by
      rcases lt_trichotomy k 1 with hk | hk | hk
      · omega
      · exact absurd (by rw [hk]; ring) hx12
      · omega
    exact d1_resonant_reach k hk2
  · exact le_of_lt (d1_generic_reach x h)

#print axioms d1_generic_reach
#print axioms d1_resonant_reach
#print axioms d1_reach_ge

end TournamentH7.LRCLadderD1
