/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: opus (LRC multi-agent project, 2026-07-02-S43)
-/
import TournamentH7.LRCFourteenSkeleton
import TournamentH7.LRC14Assembly

/-!
# The completeness surface: `hfloor` reduced to census + peel by well-founded assembly

**The honest gap this file pins (owner-flagged):** the corpus's censuses certify the families
they CONTAIN; no theorem yet proves those families EXHAUST the covering universe.  The
`hfloor` parameter of `lrc14_endgame` silently includes that exhaustion.  This file makes the
completeness leg PRECISE and discharges its glue: the well-founded assembly
`hfloor_of_census_and_peel` proves that `hfloor` follows from exactly two named legs —

* **census** (`farCount v = 0`, the bounded-spread branch): every ingredient is already
  machine-checked (klein's `LyWindowEnum` k = 8..13 rows; the N1 leg (a) census) — what
  remains is the format wiring, a finite matter;
* **peel** (`farCount v > 0`): one application of klein's far-element rate lemma
  (HYP-4001(b), proved on paper, python-verified exact) + the damped comparison margins
  (leg (c), decide-shaped, computed) + the w-band sweeps (leg (d), claimed HYP-4005) —
  the ONE genuinely unformalized mathematical statement in the chain.

The induction itself — the only piece that was neither stated nor proved anywhere — is pure
logic and is proved here sorry-free.  After this file, the completeness question is not a
vague doubt about `hfloor`: it is the conjunction of two finite-shaped parameters whose
provenance is named.  (The deeper measure-route exhaustiveness — THM-602's trichotomy and
THM-595's case tree — remains the PAPER architecture behind the same split; the witness route
needs only this counted form.)
-/

namespace LonelyRunner.LRC14.Completeness

open LonelyRunner.LRC14

/-- **The completeness assembly.**  If (census) every nonzero 13-tuple with no far element
satisfies the witness floor, and (peel) every tuple with a far element reduces to one with
strictly fewer far elements while transporting the floor back, then the floor holds for ALL
nonzero 13-tuples — i.e. `lrc14_endgame`'s `hfloor` parameter is discharged.  `farCount` is
the abstract far-element counter (instantiated by the spread cut of the N1 restructure);
the induction is strong induction on it. -/
theorem hfloor_of_census_and_peel
    (farCount : (Fin 13 → ℤ) → ℕ)
    (census : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → farCount v = 0 →
      (witnessMP : ℝ) ≤ witnessG2 (shapeOf v))
    (peel : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → 0 < farCount v →
      ∃ v' : Fin 13 → ℤ, (∀ i, v' i ≠ 0) ∧ farCount v' < farCount v ∧
        ((witnessMP : ℝ) ≤ witnessG2 (shapeOf v') →
         (witnessMP : ℝ) ≤ witnessG2 (shapeOf v))) :
    ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) →
      (witnessMP : ℝ) ≤ witnessG2 (shapeOf v) := by
  suffices H : ∀ n : ℕ, ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → farCount v = n →
      (witnessMP : ℝ) ≤ witnessG2 (shapeOf v) by
    exact fun v hv => H (farCount v) v hv rfl
  intro n
  induction n using Nat.strong_induction_on with
  | _ n ih =>
    intro v hv hn
    rcases Nat.eq_zero_or_pos n with rfl | hpos
    · exact census v hv hn
    · obtain ⟨v', hv', hlt, himp⟩ := peel v hv (hn ▸ hpos)
      exact himp (ih (farCount v') (by omega) v' hv' rfl)

/-- **The endgame through the completeness surface**: LRC(14) from census + peel + `hpartA`.
The two-parameter surface of `lrc14_endgame` is now (census, peel, hpartA), each finite-shaped
with named provenance. -/
theorem lrc14_of_census_peel_partA
    (farCount : (Fin 13 → ℤ) → ℕ)
    (census : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → farCount v = 0 →
      (witnessMP : ℝ) ≤ witnessG2 (shapeOf v))
    (peel : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → 0 < farCount v →
      ∃ v' : Fin 13 → ℤ, (∀ i, v' i ≠ 0) ∧ farCount v' < farCount v ∧
        ((witnessMP : ℝ) ≤ witnessG2 (shapeOf v') →
         (witnessMP : ℝ) ≤ witnessG2 (shapeOf v)))
    (hpartA : ∀ v : Fin 13 → ℤ, 0 < witnessG2 (shapeOf v) →
      (1 : ℝ) / 14 ≤ Mreach v) :
    LRC14Statement :=
  TournamentH7.LRC14Assembly.lrc14_endgame
    (hfloor_of_census_and_peel farCount census peel) hpartA

end LonelyRunner.LRC14.Completeness
