/-
klein-2026-07-02-S96 (HYP-4004) — the L_y census decide-table (DAG-ledger N1, [TABLE]→[LEAN]).

Provenance: full bounded-spread window censuses of the THM-534 moment-LP functional L_y
(float scan + exact top-8 confirmation, exact rational arithmetic at the maxima):
ly_full_window_census_klein.py (k = 8, 9) and ly_windows_k10_13_klein.py (k = 10..13).
ZERO shapes over cap at every row. This module records the certified maxima and verifies
the cap inequalities max L_y ≤ cap_k in Lean (norm_num — kernel-checked rational
arithmetic). The census enumerations themselves are [TABLE] nodes (Python, exact); their
Lean re-enumeration is the remaining step for a full [LEAN] flip and follows the
NestDecidable Helly/decide pattern.
-/
import Mathlib

namespace LRCLyDecideTable

/-- k=8: max L_y = 2633/7350 (at the dilated AP {0,2,...,14}) ≤ cap_8 = 2243/5880. -/
theorem row8  : (2633/7350 : ℚ) ≤ 2243/5880 := by norm_num
/-- k=9: max L_y = 26083/52920 (consec) ≤ cap_9 = 1979/4004 — the tight rung. -/
theorem row9  : (26083/52920 : ℚ) ≤ 1979/4004 := by norm_num
/-- k=10: max L_y = 45253/79380 (consec) ≤ cap_10 = 55/91. -/
theorem row10 : (45253/79380 : ℚ) ≤ 55/91 := by norm_num
/-- k=11: max L_y = 12073/17640 ≤ cap_11 = 66/91. -/
theorem row11 : (12073/17640 : ℚ) ≤ 66/91 := by norm_num
/-- k=12: max L_y = 688/945 ≤ cap_12 = 78/91. -/
theorem row12 : (688/945 : ℚ) ≤ 78/91 := by norm_num
/-- k=13: max L_y = 1037591/1375920 ≤ cap_13 = 1. -/
theorem row13 : (1037591/1375920 : ℚ) ≤ 1 := by norm_num

/-- The tight rung's exact margin: cap_9 − max_9 = 10441/7567560 ≈ 0.00138. -/
theorem tight_rung_margin : (1979/4004 : ℚ) - 26083/52920 = 10441/7567560 := by norm_num

end LRCLyDecideTable
