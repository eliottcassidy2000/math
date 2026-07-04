/-
  TournamentH7.LRCLargeTightener — THE LARGE-TIGHTENER DISCREPANCY CORE (THM-615 Lemma 3)
  (klein-2026-07-04-S125, HYP-4079).

  opus-S65 PROVED THM-615 Lemma 3 (the loose end of the m=2,f=2 confinement): if a tightener is
  large enough (`max(w₁,w₂) > u_max/(6(M(U)−1/12))`), then `M(2U ∪ {w₁,w₂}) ≥ 1/12`.  Mechanism:
  a large tightener's orbit `{w₁ t}` is dense enough on the good band `I₀ = {g_E ≥ 1/12}` (length
  `≥ (M(U)−1/12)/u_max`, from the margin, THM-613 = klein LRCMarginMeasure) to hit a MODERATE value
  (`reach ∈ [1/12, 5/12]`), which is NON-extremity, so `Ψ ≥ 1/12` (folding, klein LRCFolding).

  This file formalizes that DISCREPANCY CORE, sorry-free, on `LRCFolding.reach`:
   * `exists_reach_eq_sixth` : for `w ≥ 1`, EVERY length-`1/w` interval contains a `t` with
       `reach (w·t) = 1/6`  (`w·t` sweeps a full period ⇒ hits the midpoint value `1/6`).
   * `large_tightener_folded_ge` : if the even-part band `gv ≥ 1/12` holds on a length-`1/w₁`
       interval, then some `t` there has `min(gv(t), Ψ(t)) ≥ 1/12` — the folded value clears `1/12`.

  Combined with `LRCMarginMeasure.lonely_of_margin` (the band `I₀` has length `≥ (M(U)−1/12)/u_max`),
  `1/w₁ ≤ |I₀|` (i.e. `w₁` large) gives Lemma 3.  This is the LOOSE-END half of the confinement.

  SCOPE (HONEST): the large-tightener case only.  The small-tightener × near-AP residual (opus's
  argmax barrier) is unproved and LRC(14)-equivalent — NOT here.
-/
import TournamentH7.LRCFolding

namespace LonelyRunner.LargeTightener

open LonelyRunner.Folding

/-- **The discrepancy hit.**  For `w ≥ 1` and any `c`, the length-`1/w` interval `[c, c+1/w]`
  contains a point `t` with `reach (w·t) = 1/6` (a moderate value in `[1/12, 5/12]`).  Because
  `w·t` ranges over an interval of length `1` there, so it meets the point `k₀ + 1/6`. -/
theorem exists_reach_eq_sixth (w : ℤ) (hw : 1 ≤ w) (c : ℝ) :
    ∃ t ∈ Set.Icc c (c + 1 / (w : ℝ)), reach ((w : ℝ) * t) = 1 / 6 := by
  have hwR : (0 : ℝ) < (w : ℝ) := by exact_mod_cast hw
  set k₀ : ℤ := ⌈(w : ℝ) * c - 1 / 6⌉ with hk
  have hceil_lo : (w : ℝ) * c - 1 / 6 ≤ (k₀ : ℝ) := Int.le_ceil _
  have hceil_hi : (k₀ : ℝ) < (w : ℝ) * c - 1 / 6 + 1 := Int.ceil_lt_add_one _
  refine ⟨((k₀ : ℝ) + 1 / 6) / (w : ℝ), ⟨?_, ?_⟩, ?_⟩
  · -- c ≤ (k₀ + 1/6)/w
    rw [le_div_iff₀ hwR]; nlinarith [hceil_lo]
  · -- (k₀ + 1/6)/w ≤ c + 1/w
    rw [div_le_iff₀ hwR]
    have hexp : (c + 1 / (w : ℝ)) * (w : ℝ) = c * (w : ℝ) + 1 := by field_simp
    rw [hexp]; nlinarith [hceil_hi]
  · -- reach (w · (k₀+1/6)/w) = reach (k₀ + 1/6) = reach (1/6) = 1/6
    rw [mul_div_cancel₀ _ (ne_of_gt hwR),
      show ((k₀ : ℝ) + 1 / 6) = 1 / 6 + (k₀ : ℤ) by push_cast; ring, reach_add_int]
    unfold reach
    rw [Int.fract_eq_self.mpr ⟨by norm_num, by norm_num⟩]
    norm_num

/-- **The large-tightener fold clears 1/12 (THM-615 Lemma 3, core).**  If the even-part view `gv`
  is `≥ 1/12` throughout a length-`1/w₁` interval `[c, c+1/w₁]` (`w₁ ≥ 1`), then for ANY second
  tightener `w₂`, some `t` in that interval has the folded value
  `min( gv(t), Ψ(reach(w₁t), reach(w₂t)) ) ≥ 1/12`  (`Ψ = max(min a b, ½ − max a b)`).
  So a tightener whose comb fits inside the good band cannot force extremity everywhere. -/
theorem large_tightener_folded_ge (w₁ w₂ : ℤ) (hw₁ : 1 ≤ w₁) (c : ℝ)
    (gv : ℝ → ℝ) (hgv : ∀ t ∈ Set.Icc c (c + 1 / (w₁ : ℝ)), (1 : ℝ) / 12 ≤ gv t) :
    ∃ t ∈ Set.Icc c (c + 1 / (w₁ : ℝ)),
      (1 : ℝ) / 12 ≤ min (gv t)
        (max (min (reach ((w₁ : ℝ) * t)) (reach ((w₂ : ℝ) * t)))
             (1 / 2 - max (reach ((w₁ : ℝ) * t)) (reach ((w₂ : ℝ) * t)))) := by
  obtain ⟨t, ht, hr⟩ := exists_reach_eq_sixth w₁ hw₁ c
  refine ⟨t, ht, le_min (hgv t ht) ?_⟩
  -- Ψ ≥ 1/12: reach(w₁ t) = 1/6 ∈ [1/12, 5/12] is not extremity.
  apply psi_ge_of_not_extremity
  rintro (⟨hlt, _⟩ | ⟨_, hgt⟩)
  · rw [hr] at hlt; norm_num at hlt
  · rw [hr] at hgt; norm_num at hgt

#print axioms exists_reach_eq_sixth
#print axioms large_tightener_folded_ge

end LonelyRunner.LargeTightener
