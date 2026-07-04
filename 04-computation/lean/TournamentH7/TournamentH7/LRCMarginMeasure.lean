/-
  TournamentH7.LRCMarginMeasure — THE MARGIN → MEASURE BRIDGE, FORMALIZED (THM-613)
  (klein-2026-07-03-S120, HYP-4068).

  opus-S60 identified "the one problem: measure floor = margin rigidity".  THM-613 (paper) is
  the quantitative half: `F(t) = min_v ||v t||` is `v_max`-Lipschitz and peaks at `M(S)`, so the
  super-level set `{t : F(t) >= b}` (= the loneliness set at threshold `b`) contains an interval
  of half-width `(b - 1/n)/v_max` around the optimal time, whence

     meas{ t : Lonely n v t }  >=  2 (b - 1/n) / v_max        (b = M(S) the sharpest).

  This file FORMALIZES that bridge, sorry-free, GENERAL, and instantiates it at the covering-min
  extremizer (`deepWell14`, `M = 14/183`): a concrete positive measure floor
  `meas(lonely deepWell14) >= 13/233142`.

  SCOPE (honest): THM-613 is the bridge that CONVERTS the covering-min lower bound `M >= 14/183`
  (the rigidity's OUTPUT) into a measure floor.  It does NOT prove `M >= 14/183` — the covering-min
  lower-bound RIGIDITY (tight `M=1/n` => tight locus {AP,GW}; the confinement descent, THM-612)
  is LRC(14)-equivalent and remains OPEN.  This is the measure-floor side of that "one problem",
  made rigorous.
-/
import TournamentH7.LonelyRunner
import TournamentH7.LRCDeepWellWitness
import TournamentH7.LRCDeepWellLonely

namespace LonelyRunner.Margin

open MeasureTheory Set

/-- **Margin ⇒ lonely on an interval (Lipschitz kernel of THM-613).**  If at `tstar` every runner
  is `≥ b` from the integers (loneliness at threshold `b`), all speeds have `|v i| ≤ V`, and
  `1/n ≤ b`, then every `t` within `(b - 1/n)/V` of `tstar` is `n`-lonely.  Pure reverse-triangle
  inequality: `|v_i t − m| ≥ |v_i tstar − m| − |v_i|·|t − tstar| ≥ b − (b − 1/n) = 1/n`. -/
theorem lonely_of_margin {ι : Type*} (n : ℕ) (v : ι → ℤ) (tstar t b V : ℝ)
    (hb : ∀ i, ∀ m : ℤ, b ≤ |(v i : ℝ) * tstar - m|)
    (hV : ∀ i, |(v i : ℝ)| ≤ V)
    (hn : (1 : ℝ) / n ≤ b)
    (hVpos : 0 < V)
    (hd : |t - tstar| ≤ (b - 1 / n) / V) :
    Lonely n v t := by
  intro i m
  have hbi : b ≤ |(v i : ℝ) * tstar - m| := hb i m
  have hstep : |(v i : ℝ) * tstar - m| - |(v i : ℝ) * t - m| ≤ |(v i : ℝ) * (tstar - t)| := by
    calc |(v i : ℝ) * tstar - m| - |(v i : ℝ) * t - m|
        ≤ |((v i : ℝ) * tstar - m) - ((v i : ℝ) * t - m)| := abs_sub_abs_le_abs_sub _ _
      _ = |(v i : ℝ) * (tstar - t)| := by congr 1; ring
  have hbound : |(v i : ℝ) * (tstar - t)| ≤ V * |t - tstar| := by
    rw [abs_mul, abs_sub_comm tstar t]
    exact mul_le_mul_of_nonneg_right (hV i) (abs_nonneg _)
  have hVd : V * |t - tstar| ≤ b - 1 / n := by
    calc V * |t - tstar| ≤ V * ((b - 1 / n) / V) := mul_le_mul_of_nonneg_left hd (le_of_lt hVpos)
      _ = b - 1 / n := by field_simp
  linarith [hbi, hstep, hbound, hVd]

/-- **The margin → measure floor (THM-613), formalized.**  The `n`-loneliness set of `v` has
  Lebesgue measure `≥ 2 (b − 1/n)/V` for any threshold `1/n ≤ b` witnessed at `tstar` and any
  speed bound `V`.  Taking `b = M(S)` gives the sharpest floor `2(M(S) − 1/n)/v_max`. -/
theorem margin_measure_floor {ι : Type*} (n : ℕ) (v : ι → ℤ) (tstar b V : ℝ)
    (hb : ∀ i, ∀ m : ℤ, b ≤ |(v i : ℝ) * tstar - m|)
    (hV : ∀ i, |(v i : ℝ)| ≤ V)
    (hn : (1 : ℝ) / n ≤ b)
    (hVpos : 0 < V) :
    ENNReal.ofReal (2 * ((b - 1 / n) / V)) ≤ volume {t : ℝ | Lonely n v t} := by
  set w : ℝ := (b - 1 / n) / V with hw
  have hsub : Set.Icc (tstar - w) (tstar + w) ⊆ {t : ℝ | Lonely n v t} := by
    intro t ht
    have hdist : |t - tstar| ≤ w := by rw [abs_le]; constructor <;> [linarith [ht.1]; linarith [ht.2]]
    exact lonely_of_margin n v tstar t b V hb hV hn hVpos hdist
  calc ENNReal.ofReal (2 * w)
      = volume (Set.Icc (tstar - w) (tstar + w)) := by rw [Real.volume_Icc]; congr 1; ring
    _ ≤ volume {t : ℝ | Lonely n v t} := measure_mono hsub

/-! ### Instantiation at the covering-min extremizer `deepWell14 = {1,…,12,182}` -/

open DeepWell

/-- The deep well is `≥ 14/183` from the integers at `t* = 14/183` (the sharp margin, `> 1/14`).
  Real form at threshold `14/183`, uniform over runners — feeds `margin_measure_floor`. -/
theorem deepWell14_margin (i : Fin 13) (m : ℤ) :
    (14 : ℝ) / 183 ≤ |(deepWell14 i : ℝ) * ((14 : ℝ) / 183) - m| := by
  -- integer bound: 14 ≤ |v*14 - m*183| for v = deepWell14 i, then divide by 183.
  have hbridge : ∀ (v : ℤ), (∀ k : ℤ, (14 : ℤ) ≤ |v * 14 - k * 183|) →
      (14 : ℝ) / 183 ≤ |(v : ℝ) * ((14 : ℝ) / 183) - m| := by
    intro v hb
    have key : (14 : ℝ) ≤ |(v : ℝ) * 14 - m * 183| := by
      have h : (14 : ℝ) ≤ ((|v * 14 - m * 183| : ℤ) : ℝ) := by exact_mod_cast hb m
      rw [Int.cast_abs] at h; push_cast at h; exact h
    have hfac : (v : ℝ) * ((14 : ℝ) / 183) - m = ((v : ℝ) * 14 - m * 183) / 183 := by
      rw [eq_div_iff (by norm_num : (183 : ℝ) ≠ 0)]; ring
    rw [hfac, abs_div, show |(183 : ℝ)| = 183 by norm_num,
      le_div_iff₀ (by norm_num : (0 : ℝ) < 183)]
    rw [show (14 : ℝ) / 183 * 183 = 14 by norm_num]; exact key
  by_cases h : i.val = 12
  · have hv : deepWell14 i = 182 := by simp [deepWell14, h]
    rw [hv]
    refine hbridge 182 (fun k => ?_)
    have hd := defect_runner_lonely (n := 14) (by norm_num) k
    rw [phi6_14] at hd
    have he : (14 : ℤ) * (14 - 1) * 14 = 182 * 14 := by norm_num
    rw [he] at hd; exact hd
  · have hi : i.val ≤ 11 := by omega
    have hv : deepWell14 i = (i.val : ℤ) + 1 := by simp [deepWell14, h]
    rw [hv]
    refine hbridge ((i.val : ℤ) + 1) (fun k => ?_)
    have ha := ap_runner_lonely (n := 14) (j := (i.val : ℤ) + 1) (by norm_num) (by omega) (by omega) k
    rw [phi6_14] at ha; exact ha

/-- Every deep-well speed has `|v i| ≤ 182`. -/
theorem deepWell14_speed_le (i : Fin 13) : |(deepWell14 i : ℝ)| ≤ 182 := by
  by_cases h : i.val = 12
  · have hv : deepWell14 i = 182 := by simp [deepWell14, h]
    rw [hv]; norm_num
  · have hi : i.val ≤ 11 := by omega
    have hv : deepWell14 i = (i.val : ℤ) + 1 := by simp [deepWell14, h]
    rw [hv, abs_of_nonneg (by exact_mod_cast (show (0 : ℤ) ≤ (i.val : ℤ) + 1 by omega))]
    have hle : ((i.val : ℤ) : ℝ) ≤ 11 := by exact_mod_cast hi
    push_cast at hle ⊢; linarith

/-- **Concrete measure floor for the covering-min extremizer.**  The loneliness set of the deep
  well `{1,…,12,182}` (the unique covering-min family, `M = 14/183 > 1/14`) has positive Lebesgue
  measure `≥ 13/233142`  (`= 2·(14/183 − 1/14)/182 = (13/1281)/182`) — a rigorous instance of
  THM-613's floor, given the (open) covering-min value.  -/
theorem deepWell14_measure_floor :
    ENNReal.ofReal (13 / 233142) ≤ volume {t : ℝ | Lonely 14 deepWell14 t} := by
  have h := margin_measure_floor 14 deepWell14 ((14 : ℝ) / 183) ((14 : ℝ) / 183) 182
    deepWell14_margin deepWell14_speed_le (by norm_num) (by norm_num)
  refine le_trans ?_ h
  apply ENNReal.ofReal_le_ofReal
  norm_num

#print axioms lonely_of_margin
#print axioms margin_measure_floor
#print axioms deepWell14_measure_floor

end LonelyRunner.Margin
