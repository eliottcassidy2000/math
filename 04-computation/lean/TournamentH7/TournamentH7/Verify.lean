/-
  TournamentH7.Verify — Axiom dependency audit

  Inspect after `lake build` — the `#print axioms` output should list
  exactly the axioms documented in OCF.lean, Forbidden.lean, H21.lean,
  H63.lean, and Redei.lean (plus Lean's foundational `propext`,
  `Classical.choice`, `Quot.sound`).
-/

import TournamentH7.HSpectrum
import TournamentH7.Tilings
import TournamentH7.GridReflection

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

/-! ### Project-novel results — audit -/

/-- THM-280 corollary: grid-symmetric ↔ self-complementary via vertex
    reversal. -/
theorem gridSym_iff_audit {n : ℕ} (T : Tournament n)
    (hbp : HasBasePath T) (hn : 1 ≤ n) :
    IsGridSymmetric T ↔
    (∀ i j : Fin n,
       T.arc i j = (op T).arc (vertexReversal n i) (vertexReversal n j)) :=
  gridSym_iff_sc_via_reversal T hbp hn
#print axioms gridSym_iff_audit

/-- Score-formula corollary: regular tournaments are not self-flip
    (project-novel, oracle-2026-05-11-S1).  Used to prove
    Paley(p) ∉ SF for p ≡ 3 (mod 4). -/
theorem regular_not_SF_audit {n : ℕ} (T : Tournament n)
    (hbp : HasBasePath T) (hn : 3 ≤ n) (hreg : IsRegular T)
    (v0 : Fin n) (hv0 : v0.val = 0) :
    (tilde T).outDegree v0 ≠ T.outDegree v0 :=
  regular_not_SF T hbp hn hreg v0 hv0
#print axioms regular_not_SF_audit
