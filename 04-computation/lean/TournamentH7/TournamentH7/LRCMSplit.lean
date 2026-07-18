/-
  TournamentH7.LRCMSplit — the M-split, and the full reduction of LRC(14) to the
  inverse theorem.  boxeph-2026-07-18-S108.

  The maximizer `M(V)` never needs to be defined as a supremum: for the LRC thresholds,
  `M(V) < 1/n  ⟺  ¬ ∃ t, Lonely n V t`.  So the split "M ≥ 1/14 (immediately lonely) vs
  M < 1/13 (the crux)" is just `Classical` case analysis on `∃ t, Lonely 14 V t`, using
  the band-shrink monotonicity `Lonely 13 ⟹ Lonely 14`.

  `M_split` (PROVED): to prove a 13-family is `1/14`-lonely it suffices to handle the
  `M < 1/13` sub-case (`¬∃ Lonely 13 ⟹ ∃ Lonely 14`).  The `M ≥ 1/14` families are
  immediate; families with no `1/14`-lonely time have `M < 1/14 < 1/13`.

  `crux_of_dominance` (PROVED): the crux follows from the inverse theorem in dominance
  form (`M<1/13 ⟹ ρ≥13`) plus `ap_core_bridge` (S105).

  `lonely14_of_dominance` (PROVED): the COMPLETE reduction — LRC(14) for `v` from
  LRC(≤13) [cited] and the dominance inverse theorem `INV` [open].  Combined with the
  sieve dispatch (S106) and the density discharge (S107), the elementary skeleton of
  LRC(14) is now assembled in the kernel down to the single open `INV`.
-/
import Mathlib
import TournamentH7.LonelyRunner
import TournamentH7.LRCSieveDispatch
import TournamentH7.LRCAPCoreBridge

namespace LonelyRunner

/-- **The M-split (PROVED).**  To prove a 13-family is `1/14`-lonely it suffices to
handle the `M < 1/13` sub-case.  If some time is already `1/14`-lonely we are done;
otherwise no time is `1/14`-lonely, hence (a fortiori, `1/14 < 1/13`) none is
`1/13`-lonely, and the crux applies.  This is the split
"`M ≥ 1/14` (immediate) vs `M < 1/13` (crux)". -/
theorem M_split (v : Fin 13 → ℤ)
    (crux : (¬ ∃ t, Lonely 13 v t) → ∃ t, Lonely 14 v t) :
    ∃ t, Lonely 14 v t := by
  by_cases h : ∃ t, Lonely 14 v t
  · exact h
  · refine crux ?_
    rintro ⟨t, hlt⟩
    exact h ⟨t, lonely14_of_lonely_le (by norm_num) (by norm_num) hlt⟩

/-- **The crux from dominance (PROVED).**  If `M < 1/13` (`¬∃ Lonely 13`) forces the
`ρ ≥ 13` dominance (the inverse theorem `INV`), then `ap_core_bridge` (S105) makes the
family `1/14`-lonely. -/
theorem crux_of_dominance (cite : LRCUpTo13) (v : Fin 13 → ℤ) (hpos : ∀ i, 0 < v i)
    (inv_dom : (¬ ∃ t, Lonely 13 v t) →
      ∃ vstar : Fin 13, ∀ i, i ≠ vstar → 13 * v i ≤ v vstar)
    (hnl : ¬ ∃ t, Lonely 13 v t) : ∃ t, Lonely 14 v t := by
  obtain ⟨vstar, hdom⟩ := inv_dom hnl
  exact ap_core_bridge cite v hpos vstar hdom

/-- **The complete reduction of LRC(14) to INV (PROVED).**  Every 13-family of positive
speeds is `1/14`-lonely, given the LRC(≤13) citation and the inverse theorem in dominance
form `INV : M<1/13 ⟹ ρ≥13`.  The M-split reduces to the crux; the crux to
`INV + ap_core_bridge`.  (Non-covering families are absorbed automatically: `¬∃Lonely13`
already entails covering, so no separate sieve case is needed here.)

So `LRC(14) ⟸ LRC(≤13)[cited] + INV[open]`, kernel-checked. -/
theorem lonely14_of_dominance (cite : LRCUpTo13) (v : Fin 13 → ℤ) (hpos : ∀ i, 0 < v i)
    (inv_dom : (¬ ∃ t, Lonely 13 v t) →
      ∃ vstar : Fin 13, ∀ i, i ≠ vstar → 13 * v i ≤ v vstar) :
    ∃ t, Lonely 14 v t :=
  M_split v (crux_of_dominance cite v hpos inv_dom)

/-- **LRC(14), reduced to LRC(≤13) + INV (PROVED, universally quantified).**  The whole
conjecture follows from the citation and the inverse theorem, stated for all 13-families
at once. -/
theorem LRC14_of_INV (cite : LRCUpTo13)
    (INV : ∀ v : Fin 13 → ℤ, (∀ i, 0 < v i) → (¬ ∃ t, Lonely 13 v t) →
      ∃ vstar : Fin 13, ∀ i, i ≠ vstar → 13 * v i ≤ v vstar) :
    ∀ v : Fin 13 → ℤ, (∀ i, 0 < v i) → ∃ t, Lonely 14 v t :=
  fun v hpos => lonely14_of_dominance cite v hpos (INV v hpos)

end LonelyRunner

#print axioms LonelyRunner.M_split
#print axioms LonelyRunner.crux_of_dominance
#print axioms LonelyRunner.lonely14_of_dominance
#print axioms LonelyRunner.LRC14_of_INV
