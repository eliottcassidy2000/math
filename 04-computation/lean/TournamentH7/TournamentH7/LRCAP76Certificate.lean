/-
  TournamentH7.LRCAP76Certificate — THE TIGHT AP₇₆ DENSITY-FLOOR CERTIFICATE.
  opus-2026-07-07-S145 (HYP-5197), executing boxeph-S2's handoff with their bridge:
  the 24 roof-superlevel intervals of the q ≤ 6 Farey-76 cells, each in Good(1/7)(AP₇₆)
  by one `good_of_roof_gt` call, pairwise disjoint (sorted), summing EXACTLY to
  2314528732/40290957525 — the published exact value of μ_{1/7}(AP₇₆).
  This discharges `TailDiameter.AP76Certificate`, making the k=13 diameter ≤ 75 leg
  (`hlarge_floor_of_diam_le`) unconditional.
  Cell data verified exactly: 04-computation/gen_ap76_lean_opus_S145.py.
  Kernel-pure: no sorry, no native_decide.
-/
import Mathlib
import TournamentH7.LRCFareyRoofBridge
import TournamentH7.LRCTailDiameter

namespace LonelyRunner
namespace AP76Certificate

open TailDiameter FareyRoofBridge MeasureTheory
open scoped ENNReal

/-- Interval 1: Farey-76 cell `(0/1, 1/76)`, roof-superlevel `(0, 2 / 175)`. -/
theorem I1_good :
    Set.Ioo (0 : ℝ) (2 / 175) ⊆ Good (1 / 7) (Finset.Icc (1 : ℤ) 76) ∩ Set.Icc (0 : ℝ) 1 := by
  intro x hx; obtain ⟨hxl, hxu⟩ := hx
  refine ⟨good_of_roof_gt (p := (0 : ℤ)) (q := (1 : ℤ)) (p' := (1 : ℤ)) (q' := (76 : ℤ))
    (k := (76 : ℤ)) (x := x) (θ := (1 / 7 : ℝ))
    (by norm_num) (by norm_num) (by norm_num) (by norm_num)
    (by push_cast; nlinarith [hxl]) (by push_cast; nlinarith [hxu]) (by push_cast; nlinarith [hxu]),
    le_of_lt hxl, by nlinarith [hxu]⟩

/-- Interval 2: Farey-76 cell `(12/73, 1/6)`, roof-superlevel `(78 / 469, 1 / 6)`. -/
theorem I2_good :
    Set.Ioo (78 / 469 : ℝ) (1 / 6) ⊆ Good (1 / 7) (Finset.Icc (1 : ℤ) 76) ∩ Set.Icc (0 : ℝ) 1 := by
  intro x hx; obtain ⟨hxl, hxu⟩ := hx
  refine ⟨good_of_roof_gt (p := (12 : ℤ)) (q := (73 : ℤ)) (p' := (1 : ℤ)) (q' := (6 : ℤ))
    (k := (76 : ℤ)) (x := x) (θ := (1 / 7 : ℝ))
    (by norm_num) (by norm_num) (by norm_num) (by norm_num)
    (by push_cast; nlinarith [hxl]) (by push_cast; nlinarith [hxu]) (by push_cast; nlinarith [hxl]),
    by nlinarith [hxl], by nlinarith [hxu]⟩

/-- Interval 3: Farey-76 cell `(1/6, 12/71)`, roof-superlevel `(1 / 6, 76 / 455)`. -/
theorem I3_good :
    Set.Ioo (1 / 6 : ℝ) (76 / 455) ⊆ Good (1 / 7) (Finset.Icc (1 : ℤ) 76) ∩ Set.Icc (0 : ℝ) 1 := by
  intro x hx; obtain ⟨hxl, hxu⟩ := hx
  refine ⟨good_of_roof_gt (p := (1 : ℤ)) (q := (6 : ℤ)) (p' := (12 : ℤ)) (q' := (71 : ℤ))
    (k := (76 : ℤ)) (x := x) (θ := (1 / 7 : ℝ))
    (by norm_num) (by norm_num) (by norm_num) (by norm_num)
    (by push_cast; nlinarith [hxl]) (by push_cast; nlinarith [hxu]) (by push_cast; nlinarith [hxu]),
    by nlinarith [hxl], by nlinarith [hxu]⟩

/-- Interval 4: Farey-76 cell `(15/76, 1/5)`, roof-superlevel `(99 / 497, 1 / 5)`. -/
theorem I4_good :
    Set.Ioo (99 / 497 : ℝ) (1 / 5) ⊆ Good (1 / 7) (Finset.Icc (1 : ℤ) 76) ∩ Set.Icc (0 : ℝ) 1 := by
  intro x hx; obtain ⟨hxl, hxu⟩ := hx
  refine ⟨good_of_roof_gt (p := (15 : ℤ)) (q := (76 : ℤ)) (p' := (1 : ℤ)) (q' := (5 : ℤ))
    (k := (76 : ℤ)) (x := x) (θ := (1 / 7 : ℝ))
    (by norm_num) (by norm_num) (by norm_num) (by norm_num)
    (by push_cast; nlinarith [hxl]) (by push_cast; nlinarith [hxu]) (by push_cast; nlinarith [hxl]),
    by nlinarith [hxl], by nlinarith [hxu]⟩

/-- Interval 5: Farey-76 cell `(1/5, 15/74)`, roof-superlevel `(1 / 5, 97 / 483)`. -/
theorem I5_good :
    Set.Ioo (1 / 5 : ℝ) (97 / 483) ⊆ Good (1 / 7) (Finset.Icc (1 : ℤ) 76) ∩ Set.Icc (0 : ℝ) 1 := by
  intro x hx; obtain ⟨hxl, hxu⟩ := hx
  refine ⟨good_of_roof_gt (p := (1 : ℤ)) (q := (5 : ℤ)) (p' := (15 : ℤ)) (q' := (74 : ℤ))
    (k := (76 : ℤ)) (x := x) (θ := (1 / 7 : ℝ))
    (by norm_num) (by norm_num) (by norm_num) (by norm_num)
    (by push_cast; nlinarith [hxl]) (by push_cast; nlinarith [hxu]) (by push_cast; nlinarith [hxu]),
    by nlinarith [hxl], by nlinarith [hxu]⟩

/-- Interval 6: Farey-76 cell `(18/73, 1/4)`, roof-superlevel `(40 / 161, 1 / 4)`. -/
theorem I6_good :
    Set.Ioo (40 / 161 : ℝ) (1 / 4) ⊆ Good (1 / 7) (Finset.Icc (1 : ℤ) 76) ∩ Set.Icc (0 : ℝ) 1 := by
  intro x hx; obtain ⟨hxl, hxu⟩ := hx
  refine ⟨good_of_roof_gt (p := (18 : ℤ)) (q := (73 : ℤ)) (p' := (1 : ℤ)) (q' := (4 : ℤ))
    (k := (76 : ℤ)) (x := x) (θ := (1 / 7 : ℝ))
    (by norm_num) (by norm_num) (by norm_num) (by norm_num)
    (by push_cast; nlinarith [hxl]) (by push_cast; nlinarith [hxu]) (by push_cast; nlinarith [hxl]),
    by nlinarith [hxl], by nlinarith [hxu]⟩

/-- Interval 7: Farey-76 cell `(1/4, 19/75)`, roof-superlevel `(1 / 4, 125 / 497)`. -/
theorem I7_good :
    Set.Ioo (1 / 4 : ℝ) (125 / 497) ⊆ Good (1 / 7) (Finset.Icc (1 : ℤ) 76) ∩ Set.Icc (0 : ℝ) 1 := by
  intro x hx; obtain ⟨hxl, hxu⟩ := hx
  refine ⟨good_of_roof_gt (p := (1 : ℤ)) (q := (4 : ℤ)) (p' := (19 : ℤ)) (q' := (75 : ℤ))
    (k := (76 : ℤ)) (x := x) (θ := (1 / 7 : ℝ))
    (by norm_num) (by norm_num) (by norm_num) (by norm_num)
    (by push_cast; nlinarith [hxl]) (by push_cast; nlinarith [hxu]) (by push_cast; nlinarith [hxu]),
    by nlinarith [hxl], by nlinarith [hxu]⟩

/-- Interval 8: Farey-76 cell `(25/76, 1/3)`, roof-superlevel `(169 / 511, 1 / 3)`. -/
theorem I8_good :
    Set.Ioo (169 / 511 : ℝ) (1 / 3) ⊆ Good (1 / 7) (Finset.Icc (1 : ℤ) 76) ∩ Set.Icc (0 : ℝ) 1 := by
  intro x hx; obtain ⟨hxl, hxu⟩ := hx
  refine ⟨good_of_roof_gt (p := (25 : ℤ)) (q := (76 : ℤ)) (p' := (1 : ℤ)) (q' := (3 : ℤ))
    (k := (76 : ℤ)) (x := x) (θ := (1 / 7 : ℝ))
    (by norm_num) (by norm_num) (by norm_num) (by norm_num)
    (by push_cast; nlinarith [hxl]) (by push_cast; nlinarith [hxu]) (by push_cast; nlinarith [hxl]),
    by nlinarith [hxl], by nlinarith [hxu]⟩

/-- Interval 9: Farey-76 cell `(1/3, 25/74)`, roof-superlevel `(1 / 3, 167 / 497)`. -/
theorem I9_good :
    Set.Ioo (1 / 3 : ℝ) (167 / 497) ⊆ Good (1 / 7) (Finset.Icc (1 : ℤ) 76) ∩ Set.Icc (0 : ℝ) 1 := by
  intro x hx; obtain ⟨hxl, hxu⟩ := hx
  refine ⟨good_of_roof_gt (p := (1 : ℤ)) (q := (3 : ℤ)) (p' := (25 : ℤ)) (q' := (74 : ℤ))
    (k := (76 : ℤ)) (x := x) (θ := (1 / 7 : ℝ))
    (by norm_num) (by norm_num) (by norm_num) (by norm_num)
    (by push_cast; nlinarith [hxl]) (by push_cast; nlinarith [hxu]) (by push_cast; nlinarith [hxu]),
    by nlinarith [hxl], by nlinarith [hxu]⟩

/-- Interval 10: Farey-76 cell `(29/73, 2/5)`, roof-superlevel `(95 / 238, 2 / 5)`. -/
theorem I10_good :
    Set.Ioo (95 / 238 : ℝ) (2 / 5) ⊆ Good (1 / 7) (Finset.Icc (1 : ℤ) 76) ∩ Set.Icc (0 : ℝ) 1 := by
  intro x hx; obtain ⟨hxl, hxu⟩ := hx
  refine ⟨good_of_roof_gt (p := (29 : ℤ)) (q := (73 : ℤ)) (p' := (2 : ℤ)) (q' := (5 : ℤ))
    (k := (76 : ℤ)) (x := x) (θ := (1 / 7 : ℝ))
    (by norm_num) (by norm_num) (by norm_num) (by norm_num)
    (by push_cast; nlinarith [hxl]) (by push_cast; nlinarith [hxu]) (by push_cast; nlinarith [hxl]),
    by nlinarith [hxl], by nlinarith [hxu]⟩

/-- Interval 11: Farey-76 cell `(2/5, 29/72)`, roof-superlevel `(2 / 5, 188 / 469)`. -/
theorem I11_good :
    Set.Ioo (2 / 5 : ℝ) (188 / 469) ⊆ Good (1 / 7) (Finset.Icc (1 : ℤ) 76) ∩ Set.Icc (0 : ℝ) 1 := by
  intro x hx; obtain ⟨hxl, hxu⟩ := hx
  refine ⟨good_of_roof_gt (p := (2 : ℤ)) (q := (5 : ℤ)) (p' := (29 : ℤ)) (q' := (72 : ℤ))
    (k := (76 : ℤ)) (x := x) (θ := (1 / 7 : ℝ))
    (by norm_num) (by norm_num) (by norm_num) (by norm_num)
    (by push_cast; nlinarith [hxl]) (by push_cast; nlinarith [hxu]) (by push_cast; nlinarith [hxu]),
    by nlinarith [hxl], by nlinarith [hxu]⟩

/-- Interval 12: Farey-76 cell `(37/75, 1/2)`, roof-superlevel `(253 / 511, 1 / 2)`. -/
theorem I12_good :
    Set.Ioo (253 / 511 : ℝ) (1 / 2) ⊆ Good (1 / 7) (Finset.Icc (1 : ℤ) 76) ∩ Set.Icc (0 : ℝ) 1 := by
  intro x hx; obtain ⟨hxl, hxu⟩ := hx
  refine ⟨good_of_roof_gt (p := (37 : ℤ)) (q := (75 : ℤ)) (p' := (1 : ℤ)) (q' := (2 : ℤ))
    (k := (76 : ℤ)) (x := x) (θ := (1 / 7 : ℝ))
    (by norm_num) (by norm_num) (by norm_num) (by norm_num)
    (by push_cast; nlinarith [hxl]) (by push_cast; nlinarith [hxu]) (by push_cast; nlinarith [hxl]),
    by nlinarith [hxl], by nlinarith [hxu]⟩

/-- Interval 13: Farey-76 cell `(1/2, 38/75)`, roof-superlevel `(1 / 2, 258 / 511)`. -/
theorem I13_good :
    Set.Ioo (1 / 2 : ℝ) (258 / 511) ⊆ Good (1 / 7) (Finset.Icc (1 : ℤ) 76) ∩ Set.Icc (0 : ℝ) 1 := by
  intro x hx; obtain ⟨hxl, hxu⟩ := hx
  refine ⟨good_of_roof_gt (p := (1 : ℤ)) (q := (2 : ℤ)) (p' := (38 : ℤ)) (q' := (75 : ℤ))
    (k := (76 : ℤ)) (x := x) (θ := (1 / 7 : ℝ))
    (by norm_num) (by norm_num) (by norm_num) (by norm_num)
    (by push_cast; nlinarith [hxl]) (by push_cast; nlinarith [hxu]) (by push_cast; nlinarith [hxu]),
    by nlinarith [hxl], by nlinarith [hxu]⟩

/-- Interval 14: Farey-76 cell `(43/72, 3/5)`, roof-superlevel `(281 / 469, 3 / 5)`. -/
theorem I14_good :
    Set.Ioo (281 / 469 : ℝ) (3 / 5) ⊆ Good (1 / 7) (Finset.Icc (1 : ℤ) 76) ∩ Set.Icc (0 : ℝ) 1 := by
  intro x hx; obtain ⟨hxl, hxu⟩ := hx
  refine ⟨good_of_roof_gt (p := (43 : ℤ)) (q := (72 : ℤ)) (p' := (3 : ℤ)) (q' := (5 : ℤ))
    (k := (76 : ℤ)) (x := x) (θ := (1 / 7 : ℝ))
    (by norm_num) (by norm_num) (by norm_num) (by norm_num)
    (by push_cast; nlinarith [hxl]) (by push_cast; nlinarith [hxu]) (by push_cast; nlinarith [hxl]),
    by nlinarith [hxl], by nlinarith [hxu]⟩

/-- Interval 15: Farey-76 cell `(3/5, 44/73)`, roof-superlevel `(3 / 5, 143 / 238)`. -/
theorem I15_good :
    Set.Ioo (3 / 5 : ℝ) (143 / 238) ⊆ Good (1 / 7) (Finset.Icc (1 : ℤ) 76) ∩ Set.Icc (0 : ℝ) 1 := by
  intro x hx; obtain ⟨hxl, hxu⟩ := hx
  refine ⟨good_of_roof_gt (p := (3 : ℤ)) (q := (5 : ℤ)) (p' := (44 : ℤ)) (q' := (73 : ℤ))
    (k := (76 : ℤ)) (x := x) (θ := (1 / 7 : ℝ))
    (by norm_num) (by norm_num) (by norm_num) (by norm_num)
    (by push_cast; nlinarith [hxl]) (by push_cast; nlinarith [hxu]) (by push_cast; nlinarith [hxu]),
    by nlinarith [hxl], by nlinarith [hxu]⟩

/-- Interval 16: Farey-76 cell `(49/74, 2/3)`, roof-superlevel `(330 / 497, 2 / 3)`. -/
theorem I16_good :
    Set.Ioo (330 / 497 : ℝ) (2 / 3) ⊆ Good (1 / 7) (Finset.Icc (1 : ℤ) 76) ∩ Set.Icc (0 : ℝ) 1 := by
  intro x hx; obtain ⟨hxl, hxu⟩ := hx
  refine ⟨good_of_roof_gt (p := (49 : ℤ)) (q := (74 : ℤ)) (p' := (2 : ℤ)) (q' := (3 : ℤ))
    (k := (76 : ℤ)) (x := x) (θ := (1 / 7 : ℝ))
    (by norm_num) (by norm_num) (by norm_num) (by norm_num)
    (by push_cast; nlinarith [hxl]) (by push_cast; nlinarith [hxu]) (by push_cast; nlinarith [hxl]),
    by nlinarith [hxl], by nlinarith [hxu]⟩

/-- Interval 17: Farey-76 cell `(2/3, 51/76)`, roof-superlevel `(2 / 3, 342 / 511)`. -/
theorem I17_good :
    Set.Ioo (2 / 3 : ℝ) (342 / 511) ⊆ Good (1 / 7) (Finset.Icc (1 : ℤ) 76) ∩ Set.Icc (0 : ℝ) 1 := by
  intro x hx; obtain ⟨hxl, hxu⟩ := hx
  refine ⟨good_of_roof_gt (p := (2 : ℤ)) (q := (3 : ℤ)) (p' := (51 : ℤ)) (q' := (76 : ℤ))
    (k := (76 : ℤ)) (x := x) (θ := (1 / 7 : ℝ))
    (by norm_num) (by norm_num) (by norm_num) (by norm_num)
    (by push_cast; nlinarith [hxl]) (by push_cast; nlinarith [hxu]) (by push_cast; nlinarith [hxu]),
    by nlinarith [hxl], by nlinarith [hxu]⟩

/-- Interval 18: Farey-76 cell `(56/75, 3/4)`, roof-superlevel `(372 / 497, 3 / 4)`. -/
theorem I18_good :
    Set.Ioo (372 / 497 : ℝ) (3 / 4) ⊆ Good (1 / 7) (Finset.Icc (1 : ℤ) 76) ∩ Set.Icc (0 : ℝ) 1 := by
  intro x hx; obtain ⟨hxl, hxu⟩ := hx
  refine ⟨good_of_roof_gt (p := (56 : ℤ)) (q := (75 : ℤ)) (p' := (3 : ℤ)) (q' := (4 : ℤ))
    (k := (76 : ℤ)) (x := x) (θ := (1 / 7 : ℝ))
    (by norm_num) (by norm_num) (by norm_num) (by norm_num)
    (by push_cast; nlinarith [hxl]) (by push_cast; nlinarith [hxu]) (by push_cast; nlinarith [hxl]),
    by nlinarith [hxl], by nlinarith [hxu]⟩

/-- Interval 19: Farey-76 cell `(3/4, 55/73)`, roof-superlevel `(3 / 4, 121 / 161)`. -/
theorem I19_good :
    Set.Ioo (3 / 4 : ℝ) (121 / 161) ⊆ Good (1 / 7) (Finset.Icc (1 : ℤ) 76) ∩ Set.Icc (0 : ℝ) 1 := by
  intro x hx; obtain ⟨hxl, hxu⟩ := hx
  refine ⟨good_of_roof_gt (p := (3 : ℤ)) (q := (4 : ℤ)) (p' := (55 : ℤ)) (q' := (73 : ℤ))
    (k := (76 : ℤ)) (x := x) (θ := (1 / 7 : ℝ))
    (by norm_num) (by norm_num) (by norm_num) (by norm_num)
    (by push_cast; nlinarith [hxl]) (by push_cast; nlinarith [hxu]) (by push_cast; nlinarith [hxu]),
    by nlinarith [hxl], by nlinarith [hxu]⟩

/-- Interval 20: Farey-76 cell `(59/74, 4/5)`, roof-superlevel `(386 / 483, 4 / 5)`. -/
theorem I20_good :
    Set.Ioo (386 / 483 : ℝ) (4 / 5) ⊆ Good (1 / 7) (Finset.Icc (1 : ℤ) 76) ∩ Set.Icc (0 : ℝ) 1 := by
  intro x hx; obtain ⟨hxl, hxu⟩ := hx
  refine ⟨good_of_roof_gt (p := (59 : ℤ)) (q := (74 : ℤ)) (p' := (4 : ℤ)) (q' := (5 : ℤ))
    (k := (76 : ℤ)) (x := x) (θ := (1 / 7 : ℝ))
    (by norm_num) (by norm_num) (by norm_num) (by norm_num)
    (by push_cast; nlinarith [hxl]) (by push_cast; nlinarith [hxu]) (by push_cast; nlinarith [hxl]),
    by nlinarith [hxl], by nlinarith [hxu]⟩

/-- Interval 21: Farey-76 cell `(4/5, 61/76)`, roof-superlevel `(4 / 5, 398 / 497)`. -/
theorem I21_good :
    Set.Ioo (4 / 5 : ℝ) (398 / 497) ⊆ Good (1 / 7) (Finset.Icc (1 : ℤ) 76) ∩ Set.Icc (0 : ℝ) 1 := by
  intro x hx; obtain ⟨hxl, hxu⟩ := hx
  refine ⟨good_of_roof_gt (p := (4 : ℤ)) (q := (5 : ℤ)) (p' := (61 : ℤ)) (q' := (76 : ℤ))
    (k := (76 : ℤ)) (x := x) (θ := (1 / 7 : ℝ))
    (by norm_num) (by norm_num) (by norm_num) (by norm_num)
    (by push_cast; nlinarith [hxl]) (by push_cast; nlinarith [hxu]) (by push_cast; nlinarith [hxu]),
    by nlinarith [hxl], by nlinarith [hxu]⟩

/-- Interval 22: Farey-76 cell `(59/71, 5/6)`, roof-superlevel `(379 / 455, 5 / 6)`. -/
theorem I22_good :
    Set.Ioo (379 / 455 : ℝ) (5 / 6) ⊆ Good (1 / 7) (Finset.Icc (1 : ℤ) 76) ∩ Set.Icc (0 : ℝ) 1 := by
  intro x hx; obtain ⟨hxl, hxu⟩ := hx
  refine ⟨good_of_roof_gt (p := (59 : ℤ)) (q := (71 : ℤ)) (p' := (5 : ℤ)) (q' := (6 : ℤ))
    (k := (76 : ℤ)) (x := x) (θ := (1 / 7 : ℝ))
    (by norm_num) (by norm_num) (by norm_num) (by norm_num)
    (by push_cast; nlinarith [hxl]) (by push_cast; nlinarith [hxu]) (by push_cast; nlinarith [hxl]),
    by nlinarith [hxl], by nlinarith [hxu]⟩

/-- Interval 23: Farey-76 cell `(5/6, 61/73)`, roof-superlevel `(5 / 6, 391 / 469)`. -/
theorem I23_good :
    Set.Ioo (5 / 6 : ℝ) (391 / 469) ⊆ Good (1 / 7) (Finset.Icc (1 : ℤ) 76) ∩ Set.Icc (0 : ℝ) 1 := by
  intro x hx; obtain ⟨hxl, hxu⟩ := hx
  refine ⟨good_of_roof_gt (p := (5 : ℤ)) (q := (6 : ℤ)) (p' := (61 : ℤ)) (q' := (73 : ℤ))
    (k := (76 : ℤ)) (x := x) (θ := (1 / 7 : ℝ))
    (by norm_num) (by norm_num) (by norm_num) (by norm_num)
    (by push_cast; nlinarith [hxl]) (by push_cast; nlinarith [hxu]) (by push_cast; nlinarith [hxu]),
    by nlinarith [hxl], by nlinarith [hxu]⟩

/-- Interval 24: Farey-76 cell `(75/76, 1/1)`, roof-superlevel `(173 / 175, 1)`. -/
theorem I24_good :
    Set.Ioo (173 / 175 : ℝ) (1) ⊆ Good (1 / 7) (Finset.Icc (1 : ℤ) 76) ∩ Set.Icc (0 : ℝ) 1 := by
  intro x hx; obtain ⟨hxl, hxu⟩ := hx
  refine ⟨good_of_roof_gt (p := (75 : ℤ)) (q := (76 : ℤ)) (p' := (1 : ℤ)) (q' := (1 : ℤ))
    (k := (76 : ℤ)) (x := x) (θ := (1 / 7 : ℝ))
    (by norm_num) (by norm_num) (by norm_num) (by norm_num)
    (by push_cast; nlinarith [hxl]) (by push_cast; nlinarith [hxu]) (by push_cast; nlinarith [hxl]),
    by nlinarith [hxl], le_of_lt hxu⟩

/-- disjointness of two `Ioo` when the first ends at or before the second starts. -/
private theorem ioo_disj {a b c e : ℝ} (h : b ≤ c) :
    Disjoint (Set.Ioo a b) (Set.Ioo c e) := by
  apply Set.disjoint_left.mpr
  intro x hx1 hx2
  exact absurd (lt_of_lt_of_le hx1.2 (le_trans h (le_of_lt hx2.1))) (lt_irrefl x)

/-- **THE AP₇₆ CERTIFICATE on `{1,…,76}`.**  The 24 disjoint roof intervals sum to
the exact value. -/
theorem ap76_sum :
    ENNReal.ofReal ((2314528732 : ℝ) / 40290957525) ≤
      muGood (1 / 7) (Finset.Icc (1 : ℤ) 76) := by
  set I1 := Set.Ioo (0 : ℝ) (2 / 175) with hI1
  set I2 := Set.Ioo (78 / 469 : ℝ) (1 / 6) with hI2
  set I3 := Set.Ioo (1 / 6 : ℝ) (76 / 455) with hI3
  set I4 := Set.Ioo (99 / 497 : ℝ) (1 / 5) with hI4
  set I5 := Set.Ioo (1 / 5 : ℝ) (97 / 483) with hI5
  set I6 := Set.Ioo (40 / 161 : ℝ) (1 / 4) with hI6
  set I7 := Set.Ioo (1 / 4 : ℝ) (125 / 497) with hI7
  set I8 := Set.Ioo (169 / 511 : ℝ) (1 / 3) with hI8
  set I9 := Set.Ioo (1 / 3 : ℝ) (167 / 497) with hI9
  set I10 := Set.Ioo (95 / 238 : ℝ) (2 / 5) with hI10
  set I11 := Set.Ioo (2 / 5 : ℝ) (188 / 469) with hI11
  set I12 := Set.Ioo (253 / 511 : ℝ) (1 / 2) with hI12
  set I13 := Set.Ioo (1 / 2 : ℝ) (258 / 511) with hI13
  set I14 := Set.Ioo (281 / 469 : ℝ) (3 / 5) with hI14
  set I15 := Set.Ioo (3 / 5 : ℝ) (143 / 238) with hI15
  set I16 := Set.Ioo (330 / 497 : ℝ) (2 / 3) with hI16
  set I17 := Set.Ioo (2 / 3 : ℝ) (342 / 511) with hI17
  set I18 := Set.Ioo (372 / 497 : ℝ) (3 / 4) with hI18
  set I19 := Set.Ioo (3 / 4 : ℝ) (121 / 161) with hI19
  set I20 := Set.Ioo (386 / 483 : ℝ) (4 / 5) with hI20
  set I21 := Set.Ioo (4 / 5 : ℝ) (398 / 497) with hI21
  set I22 := Set.Ioo (379 / 455 : ℝ) (5 / 6) with hI22
  set I23 := Set.Ioo (5 / 6 : ℝ) (391 / 469) with hI23
  set I24 := Set.Ioo (173 / 175 : ℝ) (1) with hI24
  -- pairwise disjointness (intervals are sorted)
  have d1_2 : Disjoint I1 I2 := ioo_disj (by norm_num)
  have d1_3 : Disjoint I1 I3 := ioo_disj (by norm_num)
  have d2_3 : Disjoint I2 I3 := ioo_disj (by norm_num)
  have d1_4 : Disjoint I1 I4 := ioo_disj (by norm_num)
  have d2_4 : Disjoint I2 I4 := ioo_disj (by norm_num)
  have d3_4 : Disjoint I3 I4 := ioo_disj (by norm_num)
  have d1_5 : Disjoint I1 I5 := ioo_disj (by norm_num)
  have d2_5 : Disjoint I2 I5 := ioo_disj (by norm_num)
  have d3_5 : Disjoint I3 I5 := ioo_disj (by norm_num)
  have d4_5 : Disjoint I4 I5 := ioo_disj (by norm_num)
  have d1_6 : Disjoint I1 I6 := ioo_disj (by norm_num)
  have d2_6 : Disjoint I2 I6 := ioo_disj (by norm_num)
  have d3_6 : Disjoint I3 I6 := ioo_disj (by norm_num)
  have d4_6 : Disjoint I4 I6 := ioo_disj (by norm_num)
  have d5_6 : Disjoint I5 I6 := ioo_disj (by norm_num)
  have d1_7 : Disjoint I1 I7 := ioo_disj (by norm_num)
  have d2_7 : Disjoint I2 I7 := ioo_disj (by norm_num)
  have d3_7 : Disjoint I3 I7 := ioo_disj (by norm_num)
  have d4_7 : Disjoint I4 I7 := ioo_disj (by norm_num)
  have d5_7 : Disjoint I5 I7 := ioo_disj (by norm_num)
  have d6_7 : Disjoint I6 I7 := ioo_disj (by norm_num)
  have d1_8 : Disjoint I1 I8 := ioo_disj (by norm_num)
  have d2_8 : Disjoint I2 I8 := ioo_disj (by norm_num)
  have d3_8 : Disjoint I3 I8 := ioo_disj (by norm_num)
  have d4_8 : Disjoint I4 I8 := ioo_disj (by norm_num)
  have d5_8 : Disjoint I5 I8 := ioo_disj (by norm_num)
  have d6_8 : Disjoint I6 I8 := ioo_disj (by norm_num)
  have d7_8 : Disjoint I7 I8 := ioo_disj (by norm_num)
  have d1_9 : Disjoint I1 I9 := ioo_disj (by norm_num)
  have d2_9 : Disjoint I2 I9 := ioo_disj (by norm_num)
  have d3_9 : Disjoint I3 I9 := ioo_disj (by norm_num)
  have d4_9 : Disjoint I4 I9 := ioo_disj (by norm_num)
  have d5_9 : Disjoint I5 I9 := ioo_disj (by norm_num)
  have d6_9 : Disjoint I6 I9 := ioo_disj (by norm_num)
  have d7_9 : Disjoint I7 I9 := ioo_disj (by norm_num)
  have d8_9 : Disjoint I8 I9 := ioo_disj (by norm_num)
  have d1_10 : Disjoint I1 I10 := ioo_disj (by norm_num)
  have d2_10 : Disjoint I2 I10 := ioo_disj (by norm_num)
  have d3_10 : Disjoint I3 I10 := ioo_disj (by norm_num)
  have d4_10 : Disjoint I4 I10 := ioo_disj (by norm_num)
  have d5_10 : Disjoint I5 I10 := ioo_disj (by norm_num)
  have d6_10 : Disjoint I6 I10 := ioo_disj (by norm_num)
  have d7_10 : Disjoint I7 I10 := ioo_disj (by norm_num)
  have d8_10 : Disjoint I8 I10 := ioo_disj (by norm_num)
  have d9_10 : Disjoint I9 I10 := ioo_disj (by norm_num)
  have d1_11 : Disjoint I1 I11 := ioo_disj (by norm_num)
  have d2_11 : Disjoint I2 I11 := ioo_disj (by norm_num)
  have d3_11 : Disjoint I3 I11 := ioo_disj (by norm_num)
  have d4_11 : Disjoint I4 I11 := ioo_disj (by norm_num)
  have d5_11 : Disjoint I5 I11 := ioo_disj (by norm_num)
  have d6_11 : Disjoint I6 I11 := ioo_disj (by norm_num)
  have d7_11 : Disjoint I7 I11 := ioo_disj (by norm_num)
  have d8_11 : Disjoint I8 I11 := ioo_disj (by norm_num)
  have d9_11 : Disjoint I9 I11 := ioo_disj (by norm_num)
  have d10_11 : Disjoint I10 I11 := ioo_disj (by norm_num)
  have d1_12 : Disjoint I1 I12 := ioo_disj (by norm_num)
  have d2_12 : Disjoint I2 I12 := ioo_disj (by norm_num)
  have d3_12 : Disjoint I3 I12 := ioo_disj (by norm_num)
  have d4_12 : Disjoint I4 I12 := ioo_disj (by norm_num)
  have d5_12 : Disjoint I5 I12 := ioo_disj (by norm_num)
  have d6_12 : Disjoint I6 I12 := ioo_disj (by norm_num)
  have d7_12 : Disjoint I7 I12 := ioo_disj (by norm_num)
  have d8_12 : Disjoint I8 I12 := ioo_disj (by norm_num)
  have d9_12 : Disjoint I9 I12 := ioo_disj (by norm_num)
  have d10_12 : Disjoint I10 I12 := ioo_disj (by norm_num)
  have d11_12 : Disjoint I11 I12 := ioo_disj (by norm_num)
  have d1_13 : Disjoint I1 I13 := ioo_disj (by norm_num)
  have d2_13 : Disjoint I2 I13 := ioo_disj (by norm_num)
  have d3_13 : Disjoint I3 I13 := ioo_disj (by norm_num)
  have d4_13 : Disjoint I4 I13 := ioo_disj (by norm_num)
  have d5_13 : Disjoint I5 I13 := ioo_disj (by norm_num)
  have d6_13 : Disjoint I6 I13 := ioo_disj (by norm_num)
  have d7_13 : Disjoint I7 I13 := ioo_disj (by norm_num)
  have d8_13 : Disjoint I8 I13 := ioo_disj (by norm_num)
  have d9_13 : Disjoint I9 I13 := ioo_disj (by norm_num)
  have d10_13 : Disjoint I10 I13 := ioo_disj (by norm_num)
  have d11_13 : Disjoint I11 I13 := ioo_disj (by norm_num)
  have d12_13 : Disjoint I12 I13 := ioo_disj (by norm_num)
  have d1_14 : Disjoint I1 I14 := ioo_disj (by norm_num)
  have d2_14 : Disjoint I2 I14 := ioo_disj (by norm_num)
  have d3_14 : Disjoint I3 I14 := ioo_disj (by norm_num)
  have d4_14 : Disjoint I4 I14 := ioo_disj (by norm_num)
  have d5_14 : Disjoint I5 I14 := ioo_disj (by norm_num)
  have d6_14 : Disjoint I6 I14 := ioo_disj (by norm_num)
  have d7_14 : Disjoint I7 I14 := ioo_disj (by norm_num)
  have d8_14 : Disjoint I8 I14 := ioo_disj (by norm_num)
  have d9_14 : Disjoint I9 I14 := ioo_disj (by norm_num)
  have d10_14 : Disjoint I10 I14 := ioo_disj (by norm_num)
  have d11_14 : Disjoint I11 I14 := ioo_disj (by norm_num)
  have d12_14 : Disjoint I12 I14 := ioo_disj (by norm_num)
  have d13_14 : Disjoint I13 I14 := ioo_disj (by norm_num)
  have d1_15 : Disjoint I1 I15 := ioo_disj (by norm_num)
  have d2_15 : Disjoint I2 I15 := ioo_disj (by norm_num)
  have d3_15 : Disjoint I3 I15 := ioo_disj (by norm_num)
  have d4_15 : Disjoint I4 I15 := ioo_disj (by norm_num)
  have d5_15 : Disjoint I5 I15 := ioo_disj (by norm_num)
  have d6_15 : Disjoint I6 I15 := ioo_disj (by norm_num)
  have d7_15 : Disjoint I7 I15 := ioo_disj (by norm_num)
  have d8_15 : Disjoint I8 I15 := ioo_disj (by norm_num)
  have d9_15 : Disjoint I9 I15 := ioo_disj (by norm_num)
  have d10_15 : Disjoint I10 I15 := ioo_disj (by norm_num)
  have d11_15 : Disjoint I11 I15 := ioo_disj (by norm_num)
  have d12_15 : Disjoint I12 I15 := ioo_disj (by norm_num)
  have d13_15 : Disjoint I13 I15 := ioo_disj (by norm_num)
  have d14_15 : Disjoint I14 I15 := ioo_disj (by norm_num)
  have d1_16 : Disjoint I1 I16 := ioo_disj (by norm_num)
  have d2_16 : Disjoint I2 I16 := ioo_disj (by norm_num)
  have d3_16 : Disjoint I3 I16 := ioo_disj (by norm_num)
  have d4_16 : Disjoint I4 I16 := ioo_disj (by norm_num)
  have d5_16 : Disjoint I5 I16 := ioo_disj (by norm_num)
  have d6_16 : Disjoint I6 I16 := ioo_disj (by norm_num)
  have d7_16 : Disjoint I7 I16 := ioo_disj (by norm_num)
  have d8_16 : Disjoint I8 I16 := ioo_disj (by norm_num)
  have d9_16 : Disjoint I9 I16 := ioo_disj (by norm_num)
  have d10_16 : Disjoint I10 I16 := ioo_disj (by norm_num)
  have d11_16 : Disjoint I11 I16 := ioo_disj (by norm_num)
  have d12_16 : Disjoint I12 I16 := ioo_disj (by norm_num)
  have d13_16 : Disjoint I13 I16 := ioo_disj (by norm_num)
  have d14_16 : Disjoint I14 I16 := ioo_disj (by norm_num)
  have d15_16 : Disjoint I15 I16 := ioo_disj (by norm_num)
  have d1_17 : Disjoint I1 I17 := ioo_disj (by norm_num)
  have d2_17 : Disjoint I2 I17 := ioo_disj (by norm_num)
  have d3_17 : Disjoint I3 I17 := ioo_disj (by norm_num)
  have d4_17 : Disjoint I4 I17 := ioo_disj (by norm_num)
  have d5_17 : Disjoint I5 I17 := ioo_disj (by norm_num)
  have d6_17 : Disjoint I6 I17 := ioo_disj (by norm_num)
  have d7_17 : Disjoint I7 I17 := ioo_disj (by norm_num)
  have d8_17 : Disjoint I8 I17 := ioo_disj (by norm_num)
  have d9_17 : Disjoint I9 I17 := ioo_disj (by norm_num)
  have d10_17 : Disjoint I10 I17 := ioo_disj (by norm_num)
  have d11_17 : Disjoint I11 I17 := ioo_disj (by norm_num)
  have d12_17 : Disjoint I12 I17 := ioo_disj (by norm_num)
  have d13_17 : Disjoint I13 I17 := ioo_disj (by norm_num)
  have d14_17 : Disjoint I14 I17 := ioo_disj (by norm_num)
  have d15_17 : Disjoint I15 I17 := ioo_disj (by norm_num)
  have d16_17 : Disjoint I16 I17 := ioo_disj (by norm_num)
  have d1_18 : Disjoint I1 I18 := ioo_disj (by norm_num)
  have d2_18 : Disjoint I2 I18 := ioo_disj (by norm_num)
  have d3_18 : Disjoint I3 I18 := ioo_disj (by norm_num)
  have d4_18 : Disjoint I4 I18 := ioo_disj (by norm_num)
  have d5_18 : Disjoint I5 I18 := ioo_disj (by norm_num)
  have d6_18 : Disjoint I6 I18 := ioo_disj (by norm_num)
  have d7_18 : Disjoint I7 I18 := ioo_disj (by norm_num)
  have d8_18 : Disjoint I8 I18 := ioo_disj (by norm_num)
  have d9_18 : Disjoint I9 I18 := ioo_disj (by norm_num)
  have d10_18 : Disjoint I10 I18 := ioo_disj (by norm_num)
  have d11_18 : Disjoint I11 I18 := ioo_disj (by norm_num)
  have d12_18 : Disjoint I12 I18 := ioo_disj (by norm_num)
  have d13_18 : Disjoint I13 I18 := ioo_disj (by norm_num)
  have d14_18 : Disjoint I14 I18 := ioo_disj (by norm_num)
  have d15_18 : Disjoint I15 I18 := ioo_disj (by norm_num)
  have d16_18 : Disjoint I16 I18 := ioo_disj (by norm_num)
  have d17_18 : Disjoint I17 I18 := ioo_disj (by norm_num)
  have d1_19 : Disjoint I1 I19 := ioo_disj (by norm_num)
  have d2_19 : Disjoint I2 I19 := ioo_disj (by norm_num)
  have d3_19 : Disjoint I3 I19 := ioo_disj (by norm_num)
  have d4_19 : Disjoint I4 I19 := ioo_disj (by norm_num)
  have d5_19 : Disjoint I5 I19 := ioo_disj (by norm_num)
  have d6_19 : Disjoint I6 I19 := ioo_disj (by norm_num)
  have d7_19 : Disjoint I7 I19 := ioo_disj (by norm_num)
  have d8_19 : Disjoint I8 I19 := ioo_disj (by norm_num)
  have d9_19 : Disjoint I9 I19 := ioo_disj (by norm_num)
  have d10_19 : Disjoint I10 I19 := ioo_disj (by norm_num)
  have d11_19 : Disjoint I11 I19 := ioo_disj (by norm_num)
  have d12_19 : Disjoint I12 I19 := ioo_disj (by norm_num)
  have d13_19 : Disjoint I13 I19 := ioo_disj (by norm_num)
  have d14_19 : Disjoint I14 I19 := ioo_disj (by norm_num)
  have d15_19 : Disjoint I15 I19 := ioo_disj (by norm_num)
  have d16_19 : Disjoint I16 I19 := ioo_disj (by norm_num)
  have d17_19 : Disjoint I17 I19 := ioo_disj (by norm_num)
  have d18_19 : Disjoint I18 I19 := ioo_disj (by norm_num)
  have d1_20 : Disjoint I1 I20 := ioo_disj (by norm_num)
  have d2_20 : Disjoint I2 I20 := ioo_disj (by norm_num)
  have d3_20 : Disjoint I3 I20 := ioo_disj (by norm_num)
  have d4_20 : Disjoint I4 I20 := ioo_disj (by norm_num)
  have d5_20 : Disjoint I5 I20 := ioo_disj (by norm_num)
  have d6_20 : Disjoint I6 I20 := ioo_disj (by norm_num)
  have d7_20 : Disjoint I7 I20 := ioo_disj (by norm_num)
  have d8_20 : Disjoint I8 I20 := ioo_disj (by norm_num)
  have d9_20 : Disjoint I9 I20 := ioo_disj (by norm_num)
  have d10_20 : Disjoint I10 I20 := ioo_disj (by norm_num)
  have d11_20 : Disjoint I11 I20 := ioo_disj (by norm_num)
  have d12_20 : Disjoint I12 I20 := ioo_disj (by norm_num)
  have d13_20 : Disjoint I13 I20 := ioo_disj (by norm_num)
  have d14_20 : Disjoint I14 I20 := ioo_disj (by norm_num)
  have d15_20 : Disjoint I15 I20 := ioo_disj (by norm_num)
  have d16_20 : Disjoint I16 I20 := ioo_disj (by norm_num)
  have d17_20 : Disjoint I17 I20 := ioo_disj (by norm_num)
  have d18_20 : Disjoint I18 I20 := ioo_disj (by norm_num)
  have d19_20 : Disjoint I19 I20 := ioo_disj (by norm_num)
  have d1_21 : Disjoint I1 I21 := ioo_disj (by norm_num)
  have d2_21 : Disjoint I2 I21 := ioo_disj (by norm_num)
  have d3_21 : Disjoint I3 I21 := ioo_disj (by norm_num)
  have d4_21 : Disjoint I4 I21 := ioo_disj (by norm_num)
  have d5_21 : Disjoint I5 I21 := ioo_disj (by norm_num)
  have d6_21 : Disjoint I6 I21 := ioo_disj (by norm_num)
  have d7_21 : Disjoint I7 I21 := ioo_disj (by norm_num)
  have d8_21 : Disjoint I8 I21 := ioo_disj (by norm_num)
  have d9_21 : Disjoint I9 I21 := ioo_disj (by norm_num)
  have d10_21 : Disjoint I10 I21 := ioo_disj (by norm_num)
  have d11_21 : Disjoint I11 I21 := ioo_disj (by norm_num)
  have d12_21 : Disjoint I12 I21 := ioo_disj (by norm_num)
  have d13_21 : Disjoint I13 I21 := ioo_disj (by norm_num)
  have d14_21 : Disjoint I14 I21 := ioo_disj (by norm_num)
  have d15_21 : Disjoint I15 I21 := ioo_disj (by norm_num)
  have d16_21 : Disjoint I16 I21 := ioo_disj (by norm_num)
  have d17_21 : Disjoint I17 I21 := ioo_disj (by norm_num)
  have d18_21 : Disjoint I18 I21 := ioo_disj (by norm_num)
  have d19_21 : Disjoint I19 I21 := ioo_disj (by norm_num)
  have d20_21 : Disjoint I20 I21 := ioo_disj (by norm_num)
  have d1_22 : Disjoint I1 I22 := ioo_disj (by norm_num)
  have d2_22 : Disjoint I2 I22 := ioo_disj (by norm_num)
  have d3_22 : Disjoint I3 I22 := ioo_disj (by norm_num)
  have d4_22 : Disjoint I4 I22 := ioo_disj (by norm_num)
  have d5_22 : Disjoint I5 I22 := ioo_disj (by norm_num)
  have d6_22 : Disjoint I6 I22 := ioo_disj (by norm_num)
  have d7_22 : Disjoint I7 I22 := ioo_disj (by norm_num)
  have d8_22 : Disjoint I8 I22 := ioo_disj (by norm_num)
  have d9_22 : Disjoint I9 I22 := ioo_disj (by norm_num)
  have d10_22 : Disjoint I10 I22 := ioo_disj (by norm_num)
  have d11_22 : Disjoint I11 I22 := ioo_disj (by norm_num)
  have d12_22 : Disjoint I12 I22 := ioo_disj (by norm_num)
  have d13_22 : Disjoint I13 I22 := ioo_disj (by norm_num)
  have d14_22 : Disjoint I14 I22 := ioo_disj (by norm_num)
  have d15_22 : Disjoint I15 I22 := ioo_disj (by norm_num)
  have d16_22 : Disjoint I16 I22 := ioo_disj (by norm_num)
  have d17_22 : Disjoint I17 I22 := ioo_disj (by norm_num)
  have d18_22 : Disjoint I18 I22 := ioo_disj (by norm_num)
  have d19_22 : Disjoint I19 I22 := ioo_disj (by norm_num)
  have d20_22 : Disjoint I20 I22 := ioo_disj (by norm_num)
  have d21_22 : Disjoint I21 I22 := ioo_disj (by norm_num)
  have d1_23 : Disjoint I1 I23 := ioo_disj (by norm_num)
  have d2_23 : Disjoint I2 I23 := ioo_disj (by norm_num)
  have d3_23 : Disjoint I3 I23 := ioo_disj (by norm_num)
  have d4_23 : Disjoint I4 I23 := ioo_disj (by norm_num)
  have d5_23 : Disjoint I5 I23 := ioo_disj (by norm_num)
  have d6_23 : Disjoint I6 I23 := ioo_disj (by norm_num)
  have d7_23 : Disjoint I7 I23 := ioo_disj (by norm_num)
  have d8_23 : Disjoint I8 I23 := ioo_disj (by norm_num)
  have d9_23 : Disjoint I9 I23 := ioo_disj (by norm_num)
  have d10_23 : Disjoint I10 I23 := ioo_disj (by norm_num)
  have d11_23 : Disjoint I11 I23 := ioo_disj (by norm_num)
  have d12_23 : Disjoint I12 I23 := ioo_disj (by norm_num)
  have d13_23 : Disjoint I13 I23 := ioo_disj (by norm_num)
  have d14_23 : Disjoint I14 I23 := ioo_disj (by norm_num)
  have d15_23 : Disjoint I15 I23 := ioo_disj (by norm_num)
  have d16_23 : Disjoint I16 I23 := ioo_disj (by norm_num)
  have d17_23 : Disjoint I17 I23 := ioo_disj (by norm_num)
  have d18_23 : Disjoint I18 I23 := ioo_disj (by norm_num)
  have d19_23 : Disjoint I19 I23 := ioo_disj (by norm_num)
  have d20_23 : Disjoint I20 I23 := ioo_disj (by norm_num)
  have d21_23 : Disjoint I21 I23 := ioo_disj (by norm_num)
  have d22_23 : Disjoint I22 I23 := ioo_disj (by norm_num)
  have d1_24 : Disjoint I1 I24 := ioo_disj (by norm_num)
  have d2_24 : Disjoint I2 I24 := ioo_disj (by norm_num)
  have d3_24 : Disjoint I3 I24 := ioo_disj (by norm_num)
  have d4_24 : Disjoint I4 I24 := ioo_disj (by norm_num)
  have d5_24 : Disjoint I5 I24 := ioo_disj (by norm_num)
  have d6_24 : Disjoint I6 I24 := ioo_disj (by norm_num)
  have d7_24 : Disjoint I7 I24 := ioo_disj (by norm_num)
  have d8_24 : Disjoint I8 I24 := ioo_disj (by norm_num)
  have d9_24 : Disjoint I9 I24 := ioo_disj (by norm_num)
  have d10_24 : Disjoint I10 I24 := ioo_disj (by norm_num)
  have d11_24 : Disjoint I11 I24 := ioo_disj (by norm_num)
  have d12_24 : Disjoint I12 I24 := ioo_disj (by norm_num)
  have d13_24 : Disjoint I13 I24 := ioo_disj (by norm_num)
  have d14_24 : Disjoint I14 I24 := ioo_disj (by norm_num)
  have d15_24 : Disjoint I15 I24 := ioo_disj (by norm_num)
  have d16_24 : Disjoint I16 I24 := ioo_disj (by norm_num)
  have d17_24 : Disjoint I17 I24 := ioo_disj (by norm_num)
  have d18_24 : Disjoint I18 I24 := ioo_disj (by norm_num)
  have d19_24 : Disjoint I19 I24 := ioo_disj (by norm_num)
  have d20_24 : Disjoint I20 I24 := ioo_disj (by norm_num)
  have d21_24 : Disjoint I21 I24 := ioo_disj (by norm_num)
  have d22_24 : Disjoint I22 I24 := ioo_disj (by norm_num)
  have d23_24 : Disjoint I23 I24 := ioo_disj (by norm_num)
  -- union ⊆ Good ∩ [0,1]
  have hu : I1 ∪ I2 ∪ I3 ∪ I4 ∪ I5 ∪ I6 ∪ I7 ∪ I8 ∪ I9 ∪ I10 ∪ I11 ∪ I12 ∪ I13 ∪ I14 ∪ I15 ∪ I16 ∪ I17 ∪ I18 ∪ I19 ∪ I20 ∪ I21 ∪ I22 ∪ I23 ∪ I24 ⊆ Good (1 / 7) (Finset.Icc (1 : ℤ) 76) ∩ Set.Icc (0 : ℝ) 1 := by
    exact Set.union_subset (Set.union_subset (Set.union_subset (Set.union_subset (Set.union_subset (Set.union_subset (Set.union_subset (Set.union_subset (Set.union_subset (Set.union_subset (Set.union_subset (Set.union_subset (Set.union_subset (Set.union_subset (Set.union_subset (Set.union_subset (Set.union_subset (Set.union_subset (Set.union_subset (Set.union_subset (Set.union_subset (Set.union_subset (Set.union_subset I1_good I2_good) I3_good) I4_good) I5_good) I6_good) I7_good) I8_good) I9_good) I10_good) I11_good) I12_good) I13_good) I14_good) I15_good) I16_good) I17_good) I18_good) I19_good) I20_good) I21_good) I22_good) I23_good) I24_good
  -- measure of the disjoint union = sum of lengths (peel from the right)
  have hvol : volume (I1 ∪ I2 ∪ I3 ∪ I4 ∪ I5 ∪ I6 ∪ I7 ∪ I8 ∪ I9 ∪ I10 ∪ I11 ∪ I12 ∪ I13 ∪ I14 ∪ I15 ∪ I16 ∪ I17 ∪ I18 ∪ I19 ∪ I20 ∪ I21 ∪ I22 ∪ I23 ∪ I24) = volume I1 + volume I2 + volume I3 + volume I4 + volume I5 + volume I6 + volume I7 + volume I8 + volume I9 + volume I10 + volume I11 + volume I12 + volume I13 + volume I14 + volume I15 + volume I16 + volume I17 + volume I18 + volume I19 + volume I20 + volume I21 + volume I22 + volume I23 + volume I24 := by
    rw [measure_union ((((((((((((((((((((((d1_24.union_left d2_24).union_left d3_24).union_left d4_24).union_left d5_24).union_left d6_24).union_left d7_24).union_left d8_24).union_left d9_24).union_left d10_24).union_left d11_24).union_left d12_24).union_left d13_24).union_left d14_24).union_left d15_24).union_left d16_24).union_left d17_24).union_left d18_24).union_left d19_24).union_left d20_24).union_left d21_24).union_left d22_24).union_left d23_24) measurableSet_Ioo,
        measure_union (((((((((((((((((((((d1_23.union_left d2_23).union_left d3_23).union_left d4_23).union_left d5_23).union_left d6_23).union_left d7_23).union_left d8_23).union_left d9_23).union_left d10_23).union_left d11_23).union_left d12_23).union_left d13_23).union_left d14_23).union_left d15_23).union_left d16_23).union_left d17_23).union_left d18_23).union_left d19_23).union_left d20_23).union_left d21_23).union_left d22_23) measurableSet_Ioo,
        measure_union ((((((((((((((((((((d1_22.union_left d2_22).union_left d3_22).union_left d4_22).union_left d5_22).union_left d6_22).union_left d7_22).union_left d8_22).union_left d9_22).union_left d10_22).union_left d11_22).union_left d12_22).union_left d13_22).union_left d14_22).union_left d15_22).union_left d16_22).union_left d17_22).union_left d18_22).union_left d19_22).union_left d20_22).union_left d21_22) measurableSet_Ioo,
        measure_union (((((((((((((((((((d1_21.union_left d2_21).union_left d3_21).union_left d4_21).union_left d5_21).union_left d6_21).union_left d7_21).union_left d8_21).union_left d9_21).union_left d10_21).union_left d11_21).union_left d12_21).union_left d13_21).union_left d14_21).union_left d15_21).union_left d16_21).union_left d17_21).union_left d18_21).union_left d19_21).union_left d20_21) measurableSet_Ioo,
        measure_union ((((((((((((((((((d1_20.union_left d2_20).union_left d3_20).union_left d4_20).union_left d5_20).union_left d6_20).union_left d7_20).union_left d8_20).union_left d9_20).union_left d10_20).union_left d11_20).union_left d12_20).union_left d13_20).union_left d14_20).union_left d15_20).union_left d16_20).union_left d17_20).union_left d18_20).union_left d19_20) measurableSet_Ioo,
        measure_union (((((((((((((((((d1_19.union_left d2_19).union_left d3_19).union_left d4_19).union_left d5_19).union_left d6_19).union_left d7_19).union_left d8_19).union_left d9_19).union_left d10_19).union_left d11_19).union_left d12_19).union_left d13_19).union_left d14_19).union_left d15_19).union_left d16_19).union_left d17_19).union_left d18_19) measurableSet_Ioo,
        measure_union ((((((((((((((((d1_18.union_left d2_18).union_left d3_18).union_left d4_18).union_left d5_18).union_left d6_18).union_left d7_18).union_left d8_18).union_left d9_18).union_left d10_18).union_left d11_18).union_left d12_18).union_left d13_18).union_left d14_18).union_left d15_18).union_left d16_18).union_left d17_18) measurableSet_Ioo,
        measure_union (((((((((((((((d1_17.union_left d2_17).union_left d3_17).union_left d4_17).union_left d5_17).union_left d6_17).union_left d7_17).union_left d8_17).union_left d9_17).union_left d10_17).union_left d11_17).union_left d12_17).union_left d13_17).union_left d14_17).union_left d15_17).union_left d16_17) measurableSet_Ioo,
        measure_union ((((((((((((((d1_16.union_left d2_16).union_left d3_16).union_left d4_16).union_left d5_16).union_left d6_16).union_left d7_16).union_left d8_16).union_left d9_16).union_left d10_16).union_left d11_16).union_left d12_16).union_left d13_16).union_left d14_16).union_left d15_16) measurableSet_Ioo,
        measure_union (((((((((((((d1_15.union_left d2_15).union_left d3_15).union_left d4_15).union_left d5_15).union_left d6_15).union_left d7_15).union_left d8_15).union_left d9_15).union_left d10_15).union_left d11_15).union_left d12_15).union_left d13_15).union_left d14_15) measurableSet_Ioo,
        measure_union ((((((((((((d1_14.union_left d2_14).union_left d3_14).union_left d4_14).union_left d5_14).union_left d6_14).union_left d7_14).union_left d8_14).union_left d9_14).union_left d10_14).union_left d11_14).union_left d12_14).union_left d13_14) measurableSet_Ioo,
        measure_union (((((((((((d1_13.union_left d2_13).union_left d3_13).union_left d4_13).union_left d5_13).union_left d6_13).union_left d7_13).union_left d8_13).union_left d9_13).union_left d10_13).union_left d11_13).union_left d12_13) measurableSet_Ioo,
        measure_union ((((((((((d1_12.union_left d2_12).union_left d3_12).union_left d4_12).union_left d5_12).union_left d6_12).union_left d7_12).union_left d8_12).union_left d9_12).union_left d10_12).union_left d11_12) measurableSet_Ioo,
        measure_union (((((((((d1_11.union_left d2_11).union_left d3_11).union_left d4_11).union_left d5_11).union_left d6_11).union_left d7_11).union_left d8_11).union_left d9_11).union_left d10_11) measurableSet_Ioo,
        measure_union ((((((((d1_10.union_left d2_10).union_left d3_10).union_left d4_10).union_left d5_10).union_left d6_10).union_left d7_10).union_left d8_10).union_left d9_10) measurableSet_Ioo,
        measure_union (((((((d1_9.union_left d2_9).union_left d3_9).union_left d4_9).union_left d5_9).union_left d6_9).union_left d7_9).union_left d8_9) measurableSet_Ioo,
        measure_union ((((((d1_8.union_left d2_8).union_left d3_8).union_left d4_8).union_left d5_8).union_left d6_8).union_left d7_8) measurableSet_Ioo,
        measure_union (((((d1_7.union_left d2_7).union_left d3_7).union_left d4_7).union_left d5_7).union_left d6_7) measurableSet_Ioo,
        measure_union ((((d1_6.union_left d2_6).union_left d3_6).union_left d4_6).union_left d5_6) measurableSet_Ioo,
        measure_union (((d1_5.union_left d2_5).union_left d3_5).union_left d4_5) measurableSet_Ioo,
        measure_union ((d1_4.union_left d2_4).union_left d3_4) measurableSet_Ioo,
        measure_union (d1_3.union_left d2_3) measurableSet_Ioo,
        measure_union d1_2 measurableSet_Ioo]
  have hbound : volume I1 + volume I2 + volume I3 + volume I4 + volume I5 + volume I6 + volume I7 + volume I8 + volume I9 + volume I10 + volume I11 + volume I12 + volume I13 + volume I14 + volume I15 + volume I16 + volume I17 + volume I18 + volume I19 + volume I20 + volume I21 + volume I22 + volume I23 + volume I24
      ≤ muGood (1 / 7) (Finset.Icc (1 : ℤ) 76) := by
    rw [← hvol]; exact measure_mono hu
  refine le_trans ?_ hbound
  simp only [hI1, hI2, hI3, hI4, hI5, hI6, hI7, hI8, hI9, hI10, hI11, hI12, hI13, hI14, hI15, hI16, hI17, hI18, hI19, hI20, hI21, hI22, hI23, hI24, Real.volume_Ioo]
  rw [← ENNReal.ofReal_add (by norm_num) (by norm_num),
      ← ENNReal.ofReal_add (by norm_num) (by norm_num),
      ← ENNReal.ofReal_add (by norm_num) (by norm_num),
      ← ENNReal.ofReal_add (by norm_num) (by norm_num),
      ← ENNReal.ofReal_add (by norm_num) (by norm_num),
      ← ENNReal.ofReal_add (by norm_num) (by norm_num),
      ← ENNReal.ofReal_add (by norm_num) (by norm_num),
      ← ENNReal.ofReal_add (by norm_num) (by norm_num),
      ← ENNReal.ofReal_add (by norm_num) (by norm_num),
      ← ENNReal.ofReal_add (by norm_num) (by norm_num),
      ← ENNReal.ofReal_add (by norm_num) (by norm_num),
      ← ENNReal.ofReal_add (by norm_num) (by norm_num),
      ← ENNReal.ofReal_add (by norm_num) (by norm_num),
      ← ENNReal.ofReal_add (by norm_num) (by norm_num),
      ← ENNReal.ofReal_add (by norm_num) (by norm_num),
      ← ENNReal.ofReal_add (by norm_num) (by norm_num),
      ← ENNReal.ofReal_add (by norm_num) (by norm_num),
      ← ENNReal.ofReal_add (by norm_num) (by norm_num),
      ← ENNReal.ofReal_add (by norm_num) (by norm_num),
      ← ENNReal.ofReal_add (by norm_num) (by norm_num),
      ← ENNReal.ofReal_add (by norm_num) (by norm_num),
      ← ENNReal.ofReal_add (by norm_num) (by norm_num),
      ← ENNReal.ofReal_add (by norm_num) (by norm_num)]
  exact ENNReal.ofReal_le_ofReal (by norm_num)

/-- **`AP76Certificate` HOLDS** — the Prop `LRCTailDiameter` consumes, on `{0,…,75}`
via translation.  With `hlarge_floor_of_diam_le`, the k=13 diameter ≤ 75 leg is now
unconditional. -/
theorem ap76_certificate : TailDiameter.AP76Certificate := by
  unfold TailDiameter.AP76Certificate
  have hE : (Finset.Icc (1 : ℤ) 76).image (fun e => e - 1) = Finset.Icc (0 : ℤ) 75 := by
    ext n; simp only [Finset.mem_image, Finset.mem_Icc]
    constructor
    · rintro ⟨a, ⟨h1, h2⟩, rfl⟩; omega
    · intro h; exact ⟨n + 1, ⟨by omega, by omega⟩, by omega⟩
  have htr : muGood (1 / 7) (Finset.Icc (0 : ℤ) 75)
      = muGood (1 / 7) (Finset.Icc (1 : ℤ) 76) := by
    have := muGood_translate (1 / 7) (Finset.Icc (1 : ℤ) 76) 1
    rw [hE] at this; exact this
  rw [htr]; exact ap76_sum

/-- The k=13 `hlarge` floor, UNCONDITIONAL: every integer family inside a translated
window of diameter ≤ 75 has `μ_{1/7} ≥ m_P`. -/
theorem hlarge_floor_diam75_unconditional
    (E : Finset ℤ) (m D : ℤ) (hD : D ≤ 75)
    (hE : ∀ e ∈ E, e - m ∈ Finset.Icc (0 : ℤ) D) :
    ENNReal.ofReal ((14249 : ℝ) / 252252) ≤ muGood (1 / 7) E :=
  TailDiameter.hlarge_floor_of_diam_le ap76_certificate E m D hD hE

end AP76Certificate
end LonelyRunner
