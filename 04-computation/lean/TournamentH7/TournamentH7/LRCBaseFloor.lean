/-
  TournamentH7.LRCBaseFloor — STEP 1 of the far-element peel: the base good-region floor
  (mac-mini-2026-07-03-S25, THM-609).

  opus-S49's far-peel of `CoveringFarLonely 22` reduced to ONE genuine lemma (step 1): the
  ≤12-speed base has a POSITIVE-length good region.  kps-S30 reduced that to producing ONE
  rational strictly-good point (`goodRegion2_length_pos_of_strict`).  This file supplies it —
  the real→rational density bridge from the LRC(≤13) citation.

  From the citation the base is `Lonely 13` at some real `t₀`: `dist(s·t₀) ≥ 1/13 = 1/14 + 1/182`
  (opus `slack_of_lonely13`).  A rational `x` within `1/(364·V)` of `t₀` (`V = max|base|`) has, by
  the reverse triangle bound, `dist(s·x) ≥ 1/14 + 1/364 > 1/14` — STRICTLY good.  Reducing `x` to
  `[0,1)` (loneliness is 1-periodic) and feeding it to `goodRegion2_length_pos_of_strict` closes
  step 1.  This is THM-609 in its member form.

  Kernel-pure (`propext, Classical.choice, Quot.sound`); no `sorry`, no `native_decide`.
-/
import Mathlib
import TournamentH7.LonelyRunner
import TournamentH7.LRCFarPeelGood
import TournamentH7.LRCScaleSeparation

namespace LonelyRunner
namespace RatIntervals

/-- **The real→rational strict-good bridge**: given a real `t₀` at which every base speed keeps
integer-distance `≥ 1/14 + 1/182` (the LRC(≤13) slack), there is a rational `x ∈ [0,1)` at which
every base speed keeps STRICT integer-distance `> 1/14`. -/
theorem exists_strict_good_rat (base : List ℤ) (V : ℤ) (hVpos : 0 < V)
    (hV : ∀ s ∈ base, |s| ≤ V) (t0 : ℝ)
    (hslack : ∀ s ∈ base, ∀ m : ℤ, (1 : ℝ) / 14 + 1 / 182 ≤ |(s : ℝ) * t0 - (m : ℝ)|) :
    ∃ x : ℚ, 0 ≤ x ∧ x < 1 ∧ ∀ s ∈ base, ∀ m : ℤ, (1 : ℚ) / 14 < |(s : ℚ) * x - (m : ℚ)| := by
  have hVR : (0 : ℝ) < (V : ℝ) := by exact_mod_cast hVpos
  set ε : ℝ := 1 / (364 * (V : ℝ)) with hεdef
  have hεpos : 0 < ε := by rw [hεdef]; positivity
  obtain ⟨q, hq⟩ := exists_rat_near t0 hεpos
  -- strict-good for `q` itself, at every integer
  have key : ∀ s ∈ base, ∀ m : ℤ, (1 : ℚ) / 14 < |(s : ℚ) * (q : ℚ) - (m : ℚ)| := by
    intro s hs m
    have hsV : |(s : ℝ)| ≤ (V : ℝ) := by
      have := hV s hs; rw [← Int.cast_abs]; exact_mod_cast this
    have hreal : (1 : ℝ) / 14 < |(s : ℝ) * (q : ℝ) - (m : ℝ)| := by
      have h1 : (1 : ℝ) / 14 + 1 / 182 ≤ |(s : ℝ) * t0 - (m : ℝ)| := hslack s hs m
      have h2 : |(s : ℝ) * t0 - (m : ℝ)|
          ≤ |(s : ℝ) * (q : ℝ) - (m : ℝ)| + |(s : ℝ) * t0 - (s : ℝ) * (q : ℝ)| := by
        have he : (s : ℝ) * t0 - (m : ℝ)
            = ((s : ℝ) * (q : ℝ) - (m : ℝ)) + ((s : ℝ) * t0 - (s : ℝ) * (q : ℝ)) := by ring
        rw [he]; exact abs_add_le _ _
      have hq' : |t0 - (q : ℝ)| < ε := hq
      have h3 : |(s : ℝ) * t0 - (s : ℝ) * (q : ℝ)| < 1 / 364 := by
        have heq : |(s : ℝ) * t0 - (s : ℝ) * (q : ℝ)| = |(s : ℝ)| * |t0 - (q : ℝ)| := by
          rw [← abs_mul]; congr 1; ring
        have hVε : (V : ℝ) * ε = 1 / 364 := by rw [hεdef]; field_simp
        rw [heq, ← hVε]
        calc |(s : ℝ)| * |t0 - (q : ℝ)|
            ≤ (V : ℝ) * |t0 - (q : ℝ)| := by gcongr
          _ < (V : ℝ) * ε := mul_lt_mul_of_pos_left hq' hVR
      linarith
    -- transfer the strict bound ℝ → ℚ
    have hcast : ((1 : ℚ) / 14 : ℝ) < ((|(s : ℚ) * (q : ℚ) - (m : ℚ)| : ℚ) : ℝ) := by
      push_cast [Rat.cast_abs]; exact hreal
    exact_mod_cast hcast
  -- move `q` into `[0,1)` by its fractional part (distance is 1-periodic)
  refine ⟨Int.fract q, Int.fract_nonneg q, Int.fract_lt_one q, ?_⟩
  intro s hs m
  have hfr : Int.fract q = q - (⌊q⌋ : ℚ) := by
    have := Int.floor_add_fract q; linarith
  have hstep : (s : ℚ) * Int.fract q - (m : ℚ)
      = (s : ℚ) * (q : ℚ) - ((s * ⌊q⌋ + m : ℤ) : ℚ) := by
    rw [hfr]; push_cast; ring
  rw [hstep]
  exact key s hs (s * ⌊q⌋ + m)

/-- **STEP 1 (the base good-region floor), from the citation.** A base of `≤ 12` positive speeds has a
positive-length good region at band `1/14`, by the `LRC(≤13)` citation.  This is THM-609; it discharges the
`hbig`-free step-1 obligation of the far-element peel (`goodRegion2_length_pos_of_strict`). -/
theorem base_goodRegion_floor (cite : LonelyRunner.LRCUpTo13) (base : List ℤ)
    (hpos : ∀ s ∈ base, 0 < s) (hlen : base.length ≤ 12) :
    0 < length (goodRegion2 base (1 / 14)) := by
  -- a uniform bound `V` on `|s|`
  obtain ⟨V, hVpos, hV⟩ : ∃ V : ℤ, 0 < V ∧ ∀ s ∈ base, |s| ≤ V := by
    have hnn : ∀ x ∈ base.map fun s => |s|, (0 : ℤ) ≤ x := by
      intro x hx
      simp only [List.mem_map] at hx
      obtain ⟨a, _, rfl⟩ := hx
      exact abs_nonneg a
    refine ⟨1 + (base.map fun s => |s|).sum, by linarith [List.sum_nonneg hnn], ?_⟩
    intro s hs
    have hle : |s| ≤ (base.map fun s => |s|).sum :=
      List.single_le_sum hnn _ (List.mem_map_of_mem hs)
    linarith
  -- the citation point `t₀` for the ≤12-speed base
  obtain ⟨t0, hL⟩ := cite base.length hlen (fun i => base.get i)
    (fun i => (hpos _ (base.get_mem i)).ne')
  -- LRC(≤13) slack: `1/14 + 1/182 ≤ dist` for every base speed
  have hslack : ∀ s ∈ base, ∀ m : ℤ, (1 : ℝ) / 14 + 1 / 182 ≤ |(s : ℝ) * t0 - (m : ℝ)| := by
    intro s hs m
    obtain ⟨i, rfl⟩ := List.mem_iff_get.mp hs
    have hi := hL i m
    have hkey : (1 : ℝ) / 14 + 1 / 182 ≤ 1 / ((base.length + 1 : ℕ) : ℝ) := by
      have h13 : ((base.length + 1 : ℕ) : ℝ) ≤ 13 := by
        have : base.length + 1 ≤ 13 := by omega
        exact_mod_cast this
      have hpos1 : (0 : ℝ) < ((base.length + 1 : ℕ) : ℝ) := by positivity
      have := one_div_le_one_div_of_le hpos1 h13
      norm_num at this ⊢; linarith
    calc (1 : ℝ) / 14 + 1 / 182 ≤ 1 / ((base.length + 1 : ℕ) : ℝ) := hkey
      _ ≤ |(base.get i : ℝ) * t0 - (m : ℝ)| := hi
  -- the strictly-good rational point, then kps's positivity lemma
  obtain ⟨x, hx0, hx1, hgood⟩ := exists_strict_good_rat base V hVpos hV t0 hslack
  exact goodRegion2_length_pos_of_strict hpos hx0 hx1 hgood

#print axioms exists_strict_good_rat
#print axioms base_goodRegion_floor

end RatIntervals
end LonelyRunner
