/-
  TournamentH7.LRCMSplit — the M-split, and the full reduction of LRC(14) to the
  covering inverse theorem.  boxeph-2026-07-18-S108; corrected codex-2026-07-18.

  The maximizer `M(V)` never needs to be defined as a supremum: for the LRC thresholds,
  `M(V) < 1/n  ⟺  ¬ ∃ t, Lonely n V t`.  So the split "M ≥ 1/14 (immediately lonely) vs
  M < 1/13 (the crux)" is just `Classical` case analysis on `∃ t, Lonely 14 V t`, using
  the band-shrink monotonicity `Lonely 13 ⟹ Lonely 14`.

  `M_split` (PROVED): to prove a 13-family is `1/14`-lonely it suffices to handle the
  covering `M < 1/13` sub-case.  The `M ≥ 1/14` families are immediate; a family
  with no `1/14`-lonely time has `M < 1/14 < 1/13` and, by the divisor sieve, covers
  every modulus `2..14`.

  `crux_of_dominance` (PROVED): the crux follows from the covering inverse theorem in
  dominance form (`Covering(2..14) ∧ M<1/13 ⟹ ρ≥13`) plus `ap_core_bridge` (S105).

  `lonely14_of_dominance` (PROVED): the COMPLETE reduction — LRC(14) for `v` from
  LRC(≤13) [cited] and the covering dominance inverse theorem `INVcov` [open].  The
  covering premise is essential: `¬∃ Lonely 13` alone forces only divisibility coverage
  through 13, whereas the inverse theorem needs coverage through 14.
-/
import Mathlib
import TournamentH7.LonelyRunner
import TournamentH7.LRCSieveDispatch
import TournamentH7.LRCAPCoreBridge

namespace LonelyRunner

/-- **The covering M-split (PROVED).**  To prove a 13-family is `1/14`-lonely it
suffices to handle the covering `M < 1/13` sub-case.  If some time is already
`1/14`-lonely we are done.  Otherwise `counterexample_needs_all_divisors` gives
coverage of every modulus `2..14`, and band monotonicity gives no `1/13`-lonely
time.  This is the honest split "`M ≥ 1/14` (immediate) vs
`Covering(2..14) ∧ M < 1/13` (crux)". -/
theorem M_split (v : Fin 13 → ℤ)
    (crux : Covering v → (¬ ∃ t, Lonely 13 v t) → ∃ t, Lonely 14 v t) :
    ∃ t, Lonely 14 v t := by
  by_cases h : ∃ t, Lonely 14 v t
  · exact h
  · have hcover : Covering v := by
      intro n hn2 hn14
      exact counterexample_needs_all_divisors 14 v (fun t ht => h ⟨t, ht⟩) n hn2 hn14
    refine crux hcover ?_
    rintro ⟨t, hlt⟩
    exact h ⟨t, lonely14_of_lonely_le (by norm_num) (by norm_num) hlt⟩

/-- **The covering inverse theorem (OPEN, dominance form).**  Every positive covering
13-family with no `1/13`-lonely time has a speed dominating every other speed
13-fold.  Coverage means every modulus `2..14` divides some speed; it cannot be
dropped from this statement. -/
def INVcov : Prop :=
  ∀ v : Fin 13 → ℤ, (∀ i, 0 < v i) → Covering v → (¬ ∃ t, Lonely 13 v t) →
    ∃ vstar : Fin 13, ∀ i, i ≠ vstar → 13 * v i ≤ v vstar

/-- **The covering crux from dominance (PROVED).**  If
`Covering(2..14) ∧ M < 1/13` forces `ρ ≥ 13` (the inverse theorem `INVcov`), then
`ap_core_bridge` (S105) makes the family `1/14`-lonely. -/
theorem crux_of_dominance (cite : LRCUpTo13) (v : Fin 13 → ℤ) (hpos : ∀ i, 0 < v i)
    (hcover : Covering v)
    (inv_dom : Covering v → (¬ ∃ t, Lonely 13 v t) →
      ∃ vstar : Fin 13, ∀ i, i ≠ vstar → 13 * v i ≤ v vstar)
    (hnl : ¬ ∃ t, Lonely 13 v t) : ∃ t, Lonely 14 v t := by
  obtain ⟨vstar, hdom⟩ := inv_dom hcover hnl
  exact ap_core_bridge cite v hpos vstar hdom

/-- **The complete reduction of LRC(14) to `INVcov` (PROVED).**  Every 13-family of
positive speeds is `1/14`-lonely, given the LRC(≤13) citation and the covering inverse
theorem in dominance form
`INVcov : Covering(2..14) ∧ M<1/13 ⟹ ρ≥13`.  In the no-`Lonely14` branch the divisor
sieve supplies the required `Covering(2..14)` premise; band monotonicity supplies
`M<1/13`; then `ap_core_bridge` closes the branch.

So `LRC(14) ⟸ LRC(≤13)[cited] + INVcov[open]`, kernel-checked. -/
theorem lonely14_of_dominance (cite : LRCUpTo13) (v : Fin 13 → ℤ) (hpos : ∀ i, 0 < v i)
    (inv_dom : Covering v → (¬ ∃ t, Lonely 13 v t) →
      ∃ vstar : Fin 13, ∀ i, i ≠ vstar → 13 * v i ≤ v vstar) :
    ∃ t, Lonely 14 v t :=
  M_split v (fun hcover => crux_of_dominance cite v hpos hcover inv_dom)

/-- **LRC(14), reduced to LRC(≤13) + `INVcov` (PROVED, universally quantified).**
The whole conjecture follows from the citation and the covering inverse theorem. -/
theorem LRC14_of_INVcov (cite : LRCUpTo13) (inv : INVcov) :
    ∀ v : Fin 13 → ℤ, (∀ i, 0 < v i) → ∃ t, Lonely 14 v t :=
  fun v hpos => lonely14_of_dominance cite v hpos (inv v hpos)

end LonelyRunner

#print axioms LonelyRunner.M_split
#print axioms LonelyRunner.crux_of_dominance
#print axioms LonelyRunner.lonely14_of_dominance
#print axioms LonelyRunner.LRC14_of_INVcov
