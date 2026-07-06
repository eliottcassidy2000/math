/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: opus (LRC multi-agent project, 2026-07-06-S95)
-/
import TournamentH7.LRCClearCert

/-!
# The |B| = 5 production table, Lean sample (HYP-4236)

The generator (`cert_table_gen_opus_S95.py`) emitted zone-clear component certificates
for ALL 792 five-element bases — 792/792 toothMiss-verified, zero failures
(`cert_table_B5_opus_S95.tsv`).  The components CLUSTER: `(3/106, 4/106)` serves the six
smallest-component bases and `(29/94, 30/94)` the next four, so TWO in-kernel certificate
lemmas cover the ten worst bases of the whole table.  This file is that sample — the
production-consumption shape for klein's assembly plug: look up the base's row, apply the
interval's `clear` lemma, feed `strict_lonely_of_clear_component`.
-/

namespace LonelyRunner
namespace ClearRowsB5

open ClearCert

/-- Certificates on `(3/106, 4/106)` for every runner appearing in its six bases. -/
theorem cert_106 (w : ℤ)
    (hw : w = 3 ∨ w = 5 ∨ w = 6 ∨ w = 8 ∨ w = 9 ∨ w = 10 ∨ w = 11 ∨ w = 12) :
    ∀ m : ℤ, w * 3 - 106 ≤ 106 * m → 106 * m ≤ w * (3 + 1) + 106 →
      toothMiss w 3 1 106 m := by
  rcases hw with rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl <;>
  · intro m h1 h2
    have hlo : -1 ≤ m := by omega
    have hhi : m ≤ 1 := by omega
    interval_cases m <;> decide

/-- Certificates on `(29/94, 30/94)` for every runner appearing in its four bases. -/
theorem cert_94 (w : ℤ)
    (hw : w = 1 ∨ w = 2 ∨ w = 4 ∨ w = 6 ∨ w = 8 ∨ w = 9 ∨ w = 10 ∨ w = 11) :
    ∀ m : ℤ, w * 29 - 94 ≤ 94 * m → 94 * m ≤ w * (29 + 1) + 94 →
      toothMiss w 29 1 94 m := by
  rcases hw with rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl <;>
  · intro m h1 h2
    have hlo : -1 ≤ m := by omega
    have hhi : m ≤ 5 := by omega
    interval_cases m <;> decide

/-- The `(3/106, 4/106)` component is strictly clear for all its table runners —
covers the six worst |B| = 5 bases at once. -/
theorem clear_106 (w : ℤ)
    (hw : w = 3 ∨ w = 5 ∨ w = 6 ∨ w = 8 ∨ w = 9 ∨ w = 10 ∨ w = 11 ∨ w = 12) :
    ∀ t : ℝ, (3 : ℝ) / 106 < t → t < ((3 : ℝ) + 1) / 106 →
      ∀ m : ℤ, (1 : ℝ) / 13 < |(w : ℝ) * t - m| := by
  intro t ht1 ht2 m
  exact clear_of_cert w 3 1 106
    (by rcases hw with rfl|rfl|rfl|rfl|rfl|rfl|rfl|rfl <;> norm_num)
    (by norm_num) (by norm_num) (cert_106 w hw) t (by push_cast at ht1 ⊢; linarith)
    (by push_cast at ht2 ⊢; linarith) m

/-- The `(29/94, 30/94)` component, likewise — covers its four bases. -/
theorem clear_94 (w : ℤ)
    (hw : w = 1 ∨ w = 2 ∨ w = 4 ∨ w = 6 ∨ w = 8 ∨ w = 9 ∨ w = 10 ∨ w = 11) :
    ∀ t : ℝ, (29 : ℝ) / 94 < t → t < ((29 : ℝ) + 1) / 94 →
      ∀ m : ℤ, (1 : ℝ) / 13 < |(w : ℝ) * t - m| := by
  intro t ht1 ht2 m
  exact clear_of_cert w 29 1 94
    (by rcases hw with rfl|rfl|rfl|rfl|rfl|rfl|rfl|rfl <;> norm_num)
    (by norm_num) (by norm_num) (cert_94 w hw) t (by push_cast at ht1 ⊢; linarith)
    (by push_cast at ht2 ⊢; linarith) m

#print axioms clear_106
#print axioms clear_94

end ClearRowsB5
end LonelyRunner
