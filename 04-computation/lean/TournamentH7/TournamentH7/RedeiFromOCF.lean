/-
  TournamentH7.RedeiFromOCF — Rédei's H ≥ 1 derived from OCF

  Project canon Rédei's theorem has two parts:
    (R1) Every tournament has at least one Hamilton path: H(T) ≥ 1.
    (R2) The number of Hamilton paths is always odd.

  Currently both are axiomatised in `Redei.lean`.  But (R1) actually
  follows from OCF (THM-326): H(T) = 1 + 2 α_1(T) + 4 α_2(T) + …
  Since each α_k(T) ≥ 0, this is at least 1.

  This module formalises this observation.
-/

import TournamentH7.OCF
import TournamentH7.SCC

namespace Tournament

variable {n : ℕ}

/-! ### Rédei's existence (H ≥ 1) — derived from OCF, no axiom needed -/

/-- **Theorem (Rédei's existence, oracle-S10).**

    For every tournament T, `H(T) ≥ 1`. Derived from OCF
    (without needing the separate `redei_existence` axiom). -/
theorem H_ge_one (T : Tournament n) : 1 ≤ H T := by
  have hocf := ocf T
  -- hocf : H T = 1 + 2 α_1 + 4 α_2 + 8 α_3 + 16 α_4
  -- Each α_k ≥ 0 (natural numbers), so RHS ≥ 1.
  omega

/-! ### Corollaries -/

/-- H(T) ≠ 0 (Rédei's existence). -/
theorem H_ne_zero (T : Tournament n) : H T ≠ 0 := by
  have := H_ge_one T
  omega

end Tournament
