/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: opus (LRC multi-agent project, 2026-07-05-S76)
-/
import TournamentH7.LRC14CertRoute

/-!
# The one-window peel at band 1/13 (hdich lift-rigidity, window leg — HYP-4099)

klein's `lonely_of_window_margin` (S132) peels killers at band `1/14`.  The `hdich` lift
rigidity (opus-S74 trichotomy, leg 2) needs the SAME lemma one level down: 12-runner families
at band `1/13`, base margin `β = 1/12` from the LRC(12) citation.  This file transcribes the
window peel to `n = 13`:

  * `good_point13_in_long_interval` — an interval longer than one bad arc `2/(13u)` contains a
    `1/13`-good point (constructive: the left endpoint or `(m₀ + 1/13)/u`);
  * `lonely13_of_window_margin` — base margin `β > 1/13` at one point beats any killer above
    `B/(13(β − 1/13))`;
  * `lifted_lonely13_of_margin` — the rigidity consumer at `β = 1/12`, `B = 12`: a lifted
    runner above `144` is peeled outright.  Sieve-surviving lifts `r + 13k ≤ 144` (k ≤ 11)
    are exactly mac-mini's finite CRT sweep (S51).

Kernel-pure; no native_decide.  With residue pinning (S75) and this leg, `hdich` =
sieve (`sieve_one_div`) + this window + the finite sweep.
-/

namespace LonelyRunner
namespace WindowMargin13

open LRC14

/-- **A long interval contains a 1/13-good point**: if `ℓ > 2/(13u)` (one bad arc), some
`t ∈ [c, c+ℓ]` has `|u·t − m| ≥ 1/13` for every integer `m`. -/
theorem good_point13_in_long_interval (u c ℓ : ℝ) (hu : 0 < u) (hℓ : 2 / (13 * u) < ℓ) :
    ∃ t : ℝ, c ≤ t ∧ t ≤ c + ℓ ∧ ∀ m : ℤ, 1/13 ≤ |u * t - m| := by
  have hℓpos : 0 < ℓ := lt_trans (by positivity) hℓ
  by_cases h0 : ∀ m : ℤ, 1/13 ≤ |u * c - m|
  · exact ⟨c, le_refl c, by linarith, h0⟩
  · push Not at h0
    obtain ⟨m₀, hm₀⟩ := h0
    have h1 : u * c - m₀ < 1/13 := (abs_lt.mp hm₀).2
    have h2 : -(1/13) < u * c - m₀ := (abs_lt.mp hm₀).1
    have h3 : 2 / 13 < ℓ * u := by
      rw [div_lt_iff₀ (by positivity : (0:ℝ) < 13 * u)] at hℓ
      linarith
    refine ⟨((m₀ : ℝ) + 1/13) / u, ?_, ?_, ?_⟩
    · rw [le_div_iff₀ hu]
      nlinarith
    · rw [div_le_iff₀ hu]
      nlinarith
    · intro m
      have hval : u * (((m₀ : ℝ) + 1/13) / u) = (m₀ : ℝ) + 1/13 := by
        field_simp
      rw [hval]
      rcases eq_or_ne m m₀ with rfl | hne
      · rw [show (m : ℝ) + 1/13 - m = 1/13 by ring]
        norm_num
      · have h4 : (1 : ℝ) ≤ |(m₀ : ℝ) - m| := by
          have hZ : (1 : ℤ) ≤ |m₀ - m| :=
            Int.one_le_abs (sub_ne_zero.mpr (fun h => hne h.symm))
          calc (1 : ℝ) = ((1 : ℤ) : ℝ) := by norm_num
            _ ≤ ((|m₀ - m| : ℤ) : ℝ) := by exact_mod_cast hZ
            _ = |((m₀ - m : ℤ) : ℝ)| := by rw [Int.cast_abs]
            _ = |(m₀ : ℝ) - m| := by push_cast; ring_nf
        have h5 : |(m₀ : ℝ) - m| ≤ |(m₀ : ℝ) + 1/13 - m| + 1/13 := by
          calc |(m₀ : ℝ) - m| = |((m₀ : ℝ) + 1/13 - m) + (-(1/13))| := by ring_nf
            _ ≤ |(m₀ : ℝ) + 1/13 - m| + |(-(1/13) : ℝ)| := abs_add_le _ _
            _ = |(m₀ : ℝ) + 1/13 - m| + 1/13 := by norm_num
        linarith

/-- **The one-window peel at band 1/13**: 11 base runners with margin `β > 1/13` at a single
point `t*` beat any killer above the threshold `B/(13(β − 1/13))`. -/
theorem lonely13_of_window_margin (v : Fin 12 → ℤ) (istar : Fin 12) (tstar β B : ℝ)
    (hβ : 1/13 < β)
    (hbase : ∀ i, i ≠ istar → ∀ m : ℤ, β ≤ |(v i : ℝ) * tstar - m|)
    (hB : ∀ i, i ≠ istar → |(v i : ℝ)| ≤ B) (hBpos : 0 < B)
    (hkne : v istar ≠ 0)
    (hkill : B < 13 * (β - 1/13) * |(v istar : ℝ)|) :
    ∃ t : ℝ, Lonely 13 v t := by
  set u : ℝ := |(v istar : ℝ)| with hu
  have hupos : 0 < u := abs_pos.mpr (by exact_mod_cast hkne)
  set δ : ℝ := (β - 1/13) / B with hδ
  have hδpos : 0 < δ := div_pos (by linarith) hBpos
  have hlen : 2 / (13 * u) < 2 * δ := by
    rw [div_lt_iff₀ (by positivity : (0:ℝ) < 13 * u), hδ]
    rw [show 2 * ((β - 1/13) / B) * (13 * u) = 2 * (13 * (β - 1/13) * u) / B by ring]
    rw [lt_div_iff₀ hBpos]
    nlinarith
  obtain ⟨t, ht1, ht2, htk⟩ :=
    good_point13_in_long_interval u (tstar - δ) (2 * δ) hupos hlen
  have htw : |t - tstar| ≤ δ := by
    rw [abs_le]
    constructor <;> linarith
  refine ⟨t, fun i m => ?_⟩
  show (1 : ℝ) / (13 : ℕ) ≤ |(v i : ℝ) * t - m|
  push_cast
  by_cases hi : i = istar
  · subst hi
    rcases abs_cases ((v i : ℝ)) with ⟨habs, -⟩ | ⟨habs, -⟩
    · have := htk m
      rw [hu, habs] at this
      linarith
    · have := htk (-m)
      rw [hu, habs] at this
      have heq : |-(v i : ℝ) * t - (-m : ℤ)| = |(v i : ℝ) * t - m| := by
        push_cast
        rw [show -(v i : ℝ) * t - (-(m : ℝ)) = -((v i : ℝ) * t - m) by ring, abs_neg]
      rw [heq] at this
      linarith
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
    have hBδ : B * δ = β - 1/13 := by
      rw [hδ]
      field_simp
    linarith

/-- **The rigidity consumer** (hdich lift-rigidity, window leg at the citation margin):
if the 11 non-lifted runners have margin `1/12` at some point and speeds `≤ 12`, any lifted
runner above `144` is peeled — sieve-surviving lifts `r + 13k ≤ 144` are the finite CRT
sweep. -/
theorem lifted_lonely13_of_margin (v : Fin 12 → ℤ) (istar : Fin 12) (tstar : ℝ)
    (hbase : ∀ i, i ≠ istar → ∀ m : ℤ, (1:ℝ)/12 ≤ |(v i : ℝ) * tstar - m|)
    (hB : ∀ i, i ≠ istar → |(v i : ℝ)| ≤ 12)
    (hkne : v istar ≠ 0)
    (hkill : (144 : ℝ) < |(v istar : ℝ)|) :
    ∃ t : ℝ, Lonely 13 v t := by
  apply lonely13_of_window_margin v istar tstar (1/12) 12 (by norm_num) hbase hB
    (by norm_num) hkne
  have : (13 : ℝ) * (1/12 - 1/13) = 1/12 := by norm_num
  rw [this]
  linarith

end WindowMargin13
end LonelyRunner
