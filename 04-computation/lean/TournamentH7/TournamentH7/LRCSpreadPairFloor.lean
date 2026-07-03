/-
  TournamentH7.LRCSpreadPairFloor — THE SPREAD PAIR-FLOOR, Stage 1 (klein-2026-07-02-S118,
  HYP-4023).  The first Hunter pair credit, lower-bounded on the PRISTINE window for a
  spread consecutive pair — the quantitative input kps-S24's clustered/spread repair
  routes to the pair lane.

  THE MECHANISM (Python-pinned, lrc14_spread_pair_floor_pin_klein_S118.py: identity
  PASS 200/200 random cases, adversarial floor ≥ 0.92·(L/49) at D·L ≥ 1):
   * TWO-TOOTH TRAPEZOID IDENTITY (this stage): the overlap of `tooth w₁ n` and
     `tooth w₂ m` is EXACTLY `max 0 (min (2h/w₂) ((S − |r|)/(w₁w₂)))` with
     `r = n·w₂ − m·w₁` the lattice residue and `S = h(w₁+w₂)` the support half-width.
     The proof is case-free: min-sub-max turns the clip length into a four-way min.
   * THE RESIDUE WALK (stage 3): `r` advances by `−D (mod w₂)`, `D = w₂ − w₁`, per
     w₂-tooth; `D·L ≥ 2` forces wraps, each landing in the half-support and staying
     `~S/(2D)` consecutive teeth, each worth `≥ h/w₂` — the floor `≥ L/392`.
  The clustered complement (`D·L < 2`) is kps-S24's C0-combo lane; the wide complement
  (`14D > w₁`) is mac-mini's dense-events regime.
-/
import TournamentH7.LRCRealRegions

namespace LonelyRunner
namespace SpreadPairFloor

open LRC14 RealRegion

/-! ## Stage 1: the two-tooth trapezoid identity -/

/-- Four-way min: `min a b − max c d` distributes into the four differences. -/
theorem min_sub_max (a b c d : ℝ) :
    min a b - max c d = min (min (a - c) (a - d)) (min (b - c) (b - d)) := by
  simp only [min_def, max_def]
  split_ifs <;> linarith

/-- Balanced min: `min (S + r) (S − r) = S − |r|`. -/
theorem min_add_sub_abs (S r : ℝ) : min (S + r) (S - r) = S - |r| := by
  rcases le_total r 0 with hr | hr
  · rw [abs_of_nonpos hr, min_eq_left (by linarith)]; ring
  · rw [abs_of_nonneg hr, min_eq_right (by linarith)]

/-- The trapezoid: the overlap profile of two teeth at lattice residue `r`.
Plateau `2·(1/14)/w₂` (the narrower tooth), support `|r| < (1/14)(w₁+w₂)`. -/
noncomputable def trap (w₁ w₂ r : ℝ) : ℝ :=
  max 0 (min (2 * (1/14) / w₂) (((1/14) * (w₁ + w₂) - |r|) / (w₁ * w₂)))

theorem trap_nonneg (w₁ w₂ r : ℝ) : 0 ≤ trap w₁ w₂ r := le_max_left _ _

/-- **The two-tooth trapezoid identity**: the clip length of `tooth w₁ n` against the
window given by `tooth w₂ m` is exactly the trapezoid at the residue `n·w₂ − m·w₁`. -/
theorem clipLen_tooth_tooth (w₁ w₂ : ℤ) (hw₁ : 0 < w₁) (hw₁₂ : w₁ ≤ w₂) (n m : ℤ) :
    clipLen (tooth w₁ n) ((tooth w₂ m).1) ((tooth w₂ m).2)
      = trap (w₁ : ℝ) (w₂ : ℝ) (((n * w₂ - m * w₁ : ℤ) : ℝ)) := by
  have hw₁R : (0 : ℝ) < (w₁ : ℝ) := by exact_mod_cast hw₁
  have hw₂R : (0 : ℝ) < (w₂ : ℝ) := lt_of_lt_of_le hw₁R (by exact_mod_cast hw₁₂)
  have hw₁₂R : (w₁ : ℝ) ≤ (w₂ : ℝ) := by exact_mod_cast hw₁₂
  have hmul : (0 : ℝ) < (w₁ : ℝ) * w₂ := mul_pos hw₁R hw₂R
  unfold clipLen tooth trap
  simp only
  rw [min_sub_max]
  push_cast
  congr 1
  have e₁ : ((n : ℝ) + 1/14) / w₁ - ((n : ℝ) - 1/14) / w₁ = 2 * (1/14) / (w₁ : ℝ) := by
    field_simp; ring
  have e₂ : ((m : ℝ) + 1/14) / w₂ - ((m : ℝ) - 1/14) / w₂ = 2 * (1/14) / (w₂ : ℝ) := by
    field_simp; ring
  have e₃ : ((n : ℝ) + 1/14) / w₁ - ((m : ℝ) - 1/14) / w₂
      = ((1/14) * ((w₁ : ℝ) + w₂) + ((n : ℝ) * w₂ - (m : ℝ) * w₁)) / ((w₁ : ℝ) * w₂) := by
    field_simp; ring
  have e₄ : ((m : ℝ) + 1/14) / w₂ - ((n : ℝ) - 1/14) / w₁
      = ((1/14) * ((w₁ : ℝ) + w₂) - ((n : ℝ) * w₂ - (m : ℝ) * w₁)) / ((w₁ : ℝ) * w₂) := by
    field_simp; ring
  rw [e₁, e₂, e₃, e₄]
  set P : ℝ := (w₁ : ℝ) * w₂ with hP
  set S : ℝ := (1/14) * ((w₁ : ℝ) + w₂) with hS
  set r : ℝ := (n : ℝ) * w₂ - (m : ℝ) * w₁ with hr
  -- RHS's (S − |r|)/P is the min of the two leg values
  have habs : (S - |r|) / P = min ((S + r) / P) ((S - r) / P) := by
    rw [← min_add_sub_abs S r]
    rcases le_total (S + r) (S - r) with hxy | hxy
    · rw [min_eq_left hxy, min_eq_left (by apply div_le_div_of_nonneg_right hxy hmul)]
    · rw [min_eq_right hxy, min_eq_right (by apply div_le_div_of_nonneg_right hxy hmul)]
  have hplateau : 2 * (1/14) / (w₂ : ℝ) ≤ 2 * (1/14) / (w₁ : ℝ) :=
    div_le_div_of_nonneg_left (by norm_num) hw₁R hw₁₂R
  rw [habs]
  -- min (min (2h/w₁) ((S+r)/P)) (min ((S−r)/P) (2h/w₂)) = min (2h/w₂) (min ((S+r)/P) ((S−r)/P))
  apply le_antisymm
  · apply le_min
    · exact (min_le_right _ _).trans (min_le_right _ _)
    · apply le_min
      · exact (min_le_left _ _).trans (min_le_right _ _)
      · exact (min_le_right _ _).trans (min_le_left _ _)
  · apply le_min
    · apply le_min
      · exact (min_le_left _ _).trans hplateau
      · exact (min_le_right _ _).trans (min_le_left _ _)
    · apply le_min
      · exact (min_le_right _ _).trans (min_le_right _ _)
      · exact min_le_left _ _

end SpreadPairFloor
end LonelyRunner
