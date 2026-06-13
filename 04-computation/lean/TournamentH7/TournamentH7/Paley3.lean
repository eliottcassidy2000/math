/-
  TournamentH7.Paley3 — The Paley tournament at p = 3 IS the threeCycle.

  The Paley tournament Paley(3) on ZMod 3 has arc i → j iff (j - i) is a
  non-zero quadratic residue mod 3.  At p = 3: the squares are {0, 1}.
  Non-zero squares: {1}.  So arc i → j iff (j - i) ≡ 1 (mod 3).
   • 0 → 1: difference 1. ✓
   • 1 → 2: difference 1. ✓
   • 2 → 0: difference -2 ≡ 1. ✓
   So Paley(3) = the 3-cycle 0 → 1 → 2 → 0.

  This module records: threeCycle ≅ Paley(3) (axiomatically), and
  derives consequences via the threeCycle's known properties.
-/

import TournamentH7.IsoCharacterizations
import TournamentH7.PaleyAxiomatic
import TournamentH7.SmallTournaments

namespace Tournament

/-! ### Axiomatic: Paley(3) ≅ threeCycle -/

/-- **Axiom (canonical).** Paley(3) is the threeCycle tournament.
    Verified by hand: Paley arc (i → j iff (j-i) is non-zero quadratic
    residue mod 3 = 1) coincides with the 3-cycle. -/
axiom paley_3_iso_threeCycle :
    ∃ P : PaleyType 3, P.T ≅ threeCycle

/-! ### Consequences -/

/-- The threeCycle is regular (already proved in SmallTournaments). -/
example : IsRegular threeCycle := threeCycle_isRegular

/-- The threeCycle has H = 3 (matches the project canon H(Paley(3)) = 3). -/
example : H threeCycle = 3 := H_threeCycle_eq_three

/-- Any Paley(3)-type tournament has H = 3. -/
theorem paley_3_H_eq_three (P : PaleyType 3) (h : P.T ≅ threeCycle) :
    H P.T = 3 := by
  rw [H_iso_invariant P.T threeCycle h]
  exact H_threeCycle_eq_three

end Tournament
