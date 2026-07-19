/-
  TournamentH7.LRCMSplit — the covering M-split and two honest inverse interfaces
  for LRC(14).  boxeph-2026-07-18-S108; corrected codex-2026-07-18.

  The maximizer `M(V)` never needs to be defined as a supremum: for the LRC thresholds,
  `M(V) < 1/n  ⟺  ¬ ∃ t, Lonely n V t`.  Thus the split "M ≥ 1/14 (immediately
  lonely) vs M < 1/13 (the crux)" is just case analysis on `∃ t, Lonely 14 V t`,
  followed by band-shrink monotonicity `Lonely 13 ⟹ Lonely 14`.

  `M_split` (PROVED) makes the structural domain explicit.  A family with no
  `1/14`-lonely time covers every modulus `2..14` by the divisor sieve and has no
  `1/13`-lonely time by band monotonicity.

  `INVcov` is a historical sufficient premise, now REFUTED by the doubled AP
  `v_i=2(i+1)` (THM-1153):

      positive + Covering(2..14) + no Lonely13  ⟹  13-fold dominance.

  Together with the cited LRC(≤13) bridge it conditionally implies LRC(14), but the
  premise is false: the doubled AP is Covering, has no Lonely13 time, and has no
  13-dominant speed.  The consumer theorems are retained as kernel-valid historical
  implications, not as a live reduction.

  `ResidualINV` is the exact counterexample-structural target:

      positive + Covering(2..14) + no Lonely14  ⟹  13-fold dominance.

  It is an honest interface, but not by itself a mathematical reduction: with the
  cited AP-core bridge it is logically equivalent to the working LRC(14) statement.
  The omitted-covering global premise `(¬ ∃ Lonely 13) ⟹ dominance` is false;
  `{1,...,13}` is a counterexample to that premise.
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
time. -/
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

/-- **The historical covering inverse premise (REFUTED, dominance form).**
Every positive covering 13-family with no `1/13`-lonely time was conjectured to
have a speed dominating every other speed 13-fold.  THM-1153 refutes this exact
Prop with `2*{1,...,13}`.  It remains named only so the already proved conditional
implications have a stable interface. -/
def INVcov : Prop :=
  ∀ v : Fin 13 → ℤ, (∀ i, 0 < v i) → Covering v → (¬ ∃ t, Lonely 13 v t) →
    ∃ vstar : Fin 13, ∀ i, i ≠ vstar → 13 * v i ≤ v vstar

/-- **The covering crux from dominance (PROVED).**  If
`Covering(2..14) ∧ M < 1/13` forces `ρ ≥ 13`, then `ap_core_bridge` makes the
family `1/14`-lonely. -/
theorem crux_of_dominance (cite : LRCUpTo13) (v : Fin 13 → ℤ) (hpos : ∀ i, 0 < v i)
    (hcover : Covering v)
    (inv_dom : Covering v → (¬ ∃ t, Lonely 13 v t) →
      ∃ vstar : Fin 13, ∀ i, i ≠ vstar → 13 * v i ≤ v vstar)
    (hnl : ¬ ∃ t, Lonely 13 v t) : ∃ t, Lonely 14 v t := by
  obtain ⟨vstar, hdom⟩ := inv_dom hcover hnl
  exact ap_core_bridge cite v hpos vstar hdom

/-- **Conditional implication from the refuted `INVcov` premise (PROVED).**  In the
no-`Lonely14` branch the divisor sieve supplies `Covering(2..14)` and band
monotonicity supplies `M < 1/13`; `ap_core_bridge` then closes the branch. -/
theorem lonely14_of_dominance (cite : LRCUpTo13) (v : Fin 13 → ℤ) (hpos : ∀ i, 0 < v i)
    (inv_dom : Covering v → (¬ ∃ t, Lonely 13 v t) →
      ∃ vstar : Fin 13, ∀ i, i ≠ vstar → 13 * v i ≤ v vstar) :
    ∃ t, Lonely 14 v t :=
  M_split v (fun hcover => crux_of_dominance cite v hpos hcover inv_dom)

/-- **Historical conditional capstone (PROVED, but `INVcov` is refuted).** -/
theorem LRC14_of_INVcov (cite : LRCUpTo13) (inv : INVcov) :
    ∀ v : Fin 13 → ℤ, (∀ i, 0 < v i) → ∃ t, Lonely 14 v t :=
  fun v hpos => lonely14_of_dominance cite v hpos (inv v hpos)

/-- **The exact residual inverse target.**  Dominance is requested only for an actual
counterexample candidate: a positive covering family having no `1/14`-lonely time.
The `Covering` argument is forced by the final argument via the divisor sieve, but
retaining it records the intended structural domain. -/
def ResidualINV : Prop :=
  ∀ v : Fin 13 → ℤ, (∀ i, 0 < v i) → Covering v →
    (¬ ∃ t, Lonely 14 v t) →
    ∃ vstar : Fin 13, ∀ i, i ≠ vstar → 13 * v i ≤ v vstar

/-- The refuted but formally stronger premise `INVcov` supplies the exact residual
target.  This implication remains logically valid and is not a live proof route. -/
theorem residualINV_of_INVcov (inv : INVcov) : ResidualINV := by
  intro v hpos hcover h14
  have h13 : ¬ ∃ t, Lonely 13 v t := by
    rintro ⟨t, ht⟩
    exact h14 ⟨t, lonely14_of_lonely_le (by norm_num) (by norm_num) ht⟩
  exact inv v hpos hcover h13

/-- **The exact hard-strip disjunction (PROVED from `ResidualINV`).**  Every positive
family is either already `1/14`-lonely or has a `13`-dominant speed.  In the second
branch Covering is derived from the absence of a lonely time by the sieve. -/
theorem lonely14_or_dominance_of_residual_INV (inv : ResidualINV)
    (v : Fin 13 → ℤ) (hpos : ∀ i, 0 < v i) :
    (∃ t, Lonely 14 v t) ∨
      ∃ vstar : Fin 13, ∀ i, i ≠ vstar → 13 * v i ≤ v vstar := by
  by_cases h14 : ∃ t, Lonely 14 v t
  · exact Or.inl h14
  · right
    have hcover : Covering v := by
      by_contra hnc
      exact h14 (sieve_dispatch v hnc)
    exact inv v hpos hcover h14

/-- **LRC(14) from the exact residual inverse (PROVED).**  The already-lonely branch
is immediate; a dominant speed in the residual branch is dispatched by
`ap_core_bridge`. -/
theorem lonely14_of_residual_INV (cite : LRCUpTo13) (inv : ResidualINV)
    (v : Fin 13 → ℤ) (hpos : ∀ i, 0 < v i) : ∃ t, Lonely 14 v t := by
  rcases lonely14_or_dominance_of_residual_INV inv v hpos with h14 | ⟨vstar, hdom⟩
  · exact h14
  · exact ap_core_bridge cite v hpos vstar hdom

/-- **LRC(14), reduced to the residual inverse (PROVED, universally quantified).** -/
theorem LRC14_of_residual_INV (cite : LRCUpTo13) (inv : ResidualINV) :
    ∀ v : Fin 13 → ℤ, (∀ i, 0 < v i) → ∃ t, Lonely 14 v t :=
  fun v hpos => lonely14_of_residual_INV cite inv v hpos

/-- Any proof of the working LRC(14) statement makes `ResidualINV` vacuous. -/
theorem residualINV_of_LRC14
    (hLRC14 : ∀ v : Fin 13 → ℤ, (∀ i, 0 < v i) → ∃ t, Lonely 14 v t) :
    ResidualINV := by
  intro v hpos _ h14
  exact False.elim (h14 (hLRC14 v hpos))

/-- With the cited AP-core bridge, `ResidualINV` is logically equivalent to the
working LRC(14) statement.  This records why the residual interface is exact rather
than a noncircular reduction. -/
theorem residualINV_iff_LRC14 (cite : LRCUpTo13) :
    ResidualINV ↔
      ∀ v : Fin 13 → ℤ, (∀ i, 0 < v i) → ∃ t, Lonely 14 v t :=
  ⟨LRC14_of_residual_INV cite, residualINV_of_LRC14⟩

end LonelyRunner

#print axioms LonelyRunner.M_split
#print axioms LonelyRunner.crux_of_dominance
#print axioms LonelyRunner.lonely14_of_dominance
#print axioms LonelyRunner.LRC14_of_INVcov
#print axioms LonelyRunner.residualINV_of_INVcov
#print axioms LonelyRunner.lonely14_or_dominance_of_residual_INV
#print axioms LonelyRunner.lonely14_of_residual_INV
#print axioms LonelyRunner.LRC14_of_residual_INV
#print axioms LonelyRunner.residualINV_of_LRC14
#print axioms LonelyRunner.residualINV_iff_LRC14
