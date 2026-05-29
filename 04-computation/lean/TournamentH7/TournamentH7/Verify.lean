/-
  TournamentH7.Verify — Axiom dependency audit

  Inspect after `lake build` — the `#print axioms` output should list
  exactly the axioms documented in OCF.lean, Forbidden.lean, H21.lean,
  H63.lean, and Redei.lean (plus Lean's foundational `propext`,
  `Classical.choice`, `Quot.sound`).
-/

import TournamentH7.HSpectrum

open Tournament

/-! ### Audit each individual theorem -/

theorem H_ne_seven_audit {n : ℕ} (T : Tournament n) : H T ≠ 7 := H_ne_seven T
#print axioms H_ne_seven_audit

theorem H_ne_twentyone_audit {n : ℕ} (T : Tournament n) : H T ≠ 21 := H_ne_twentyone T
#print axioms H_ne_twentyone_audit

theorem H_ne_sixtythree_audit {n : ℕ} (T : Tournament n) : H T ≠ 63 := H_ne_sixtythree T
#print axioms H_ne_sixtythree_audit

theorem H_pos_audit {n : ℕ} (hn : 1 ≤ n) (T : Tournament n) : H T ≠ 0 := H_pos hn T
#print axioms H_pos_audit

theorem forbidden_trio_audit {n : ℕ} (T : Tournament n) :
    H T ≠ 7 ∧ H T ≠ 21 ∧ H T ≠ 63 := H_not_in_forbidden_trio T
#print axioms forbidden_trio_audit
