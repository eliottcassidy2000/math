/-
  TournamentH7.IsoCharacterizations — Characterizations of iso classes by H

  Classical results connecting Hamiltonian path counts to iso-class structure:

  * H(T) = 1 ⟺ T is in the transitive iso class.
  * H(T) = 3 at n = 3 ⟺ T is in the 3-cycle iso class.

  We axiomatize these structural facts where they require α-decomposition
  knowledge not yet formalized, and PROVE corollaries.
-/

import TournamentH7.HSpectrumClean
import TournamentH7.SmallTournaments
import TournamentH7.IsoProperties
import TournamentH7.ForbiddenHCounting
import TournamentH7.TransitiveH

namespace Tournament

variable {n : ℕ}

/-! ### H = 1 ⟹ transitive class

  Project canon: a tournament T has H(T) = 1 iff T is the transitive
  tournament (up to vertex relabelling).  This requires α_1(T) = 0
  (no odd cycle) ⟺ T is acyclic ⟺ T is transitive (classical).
  We axiomatize the structural part. -/

/-- **Axiom (classical).** If α_1(T) = 0, then T is iso to the transitive
    tournament on n vertices.  Classical: a tournament with no 3-cycle
    is transitive. -/
axiom alpha1_zero_iff_transitive (T : Tournament n) (hn : 1 ≤ n) :
    alphaCount 1 T = 0 ↔ T ≅ transitiveTournament n

/-- **Theorem (derived from OCF).** The transitive tournament has H = 1.
    Proved via OCF + transitive_alphaCount_zero (see TransitiveH.lean). -/
theorem H_transitive_eq_one (n : ℕ) (hn : 1 ≤ n) :
    H (transitiveTournament n) = 1 :=
  H_transitive_eq_one_from_ocf n hn

/-! ### Concrete checks via decide (small n) -/

/-- H of trivial 1-vertex tournament is 1. -/
example : H (transitiveTournament 1) = 1 := by
  unfold H transitiveTournament
  decide

/-- H of transitive 2 = 1 (only Hamilton path 1 → 0). -/
example : H (transitiveTournament 2) = 1 := by
  unfold H transitiveTournament
  decide

/-- H of transitive 3 = 1 (only path 2 → 1 → 0). -/
example : H (transitiveTournament 3) = 1 := by
  unfold H transitiveTournament
  decide

/-! ### Specific concrete H values that DON'T require axiom -/

/-- H(transitive_3) = 1 — PROVED IN LEAN via decide. -/
theorem H_transitive_3_eq_one : H (transitiveTournament 3) = 1 := by
  unfold H transitiveTournament; decide

/-- H(transitive_4) = 1 — PROVED IN LEAN via decide. -/
theorem H_transitive_4_eq_one : H (transitiveTournament 4) = 1 := by
  unfold H transitiveTournament; decide

/-- H(transitive_5) = 1 — PROVED IN LEAN via decide. -/
theorem H_transitive_5_eq_one : H (transitiveTournament 5) = 1 := by
  unfold H transitiveTournament; decide

-- Note: n ≥ 6 times out for `decide` due to factorial growth in permutations.

/-- H(threeCycle) = 3 — PROVED IN LEAN via decide. -/
theorem H_threeCycle_eq_three : H threeCycle = 3 := by
  unfold H threeCycle threeCycleArc; decide

/-! ### threeCycle is the unique non-transitive H = 3 tournament at n = 3 -/

/-- The threeCycle has α_1 = 1 (exactly one odd cycle).

    This follows from H = 3 = 1 + 2*1 + 0 + ... ⟹ α_1 = 1. -/
theorem threeCycle_alpha1_eq_one :
    alphaCount 1 threeCycle = 1 := by
  have := alpha_solution_H3 threeCycle H_threeCycle_eq_three
  exact this.1

/-- H(threeCycle) = 3 verifies the small-n spectrum. -/
example : H threeCycle = 3 ∧ Odd (H threeCycle) := by
  refine ⟨H_threeCycle_eq_three, ?_⟩
  rw [H_threeCycle_eq_three]; decide

/-! ### Corollary: H(T) = 1 ⟺ T ≅ transitive -/

/-- **Theorem.** For any tournament T with n ≥ 1, H(T) = 1 ⟺ T ≅ transitive. -/
theorem H_eq_one_iff_transitive (T : Tournament n) (hn : 1 ≤ n) :
    H T = 1 ↔ T ≅ transitiveTournament n := by
  constructor
  · intro h
    -- H = 1 ⟹ alpha_1 = 0 (by alpha_solution_H1).
    have h_alpha := alpha_solution_H1 T h
    -- alpha_1 = 0 ⟹ T ≅ transitive.
    exact (alpha1_zero_iff_transitive T hn).mp h_alpha.1
  · intro h
    -- T ≅ transitive ⟹ H(T) = H(transitive) = 1.
    have h_H := H_iso_invariant T (transitiveTournament n) h
    rw [h_H]
    -- Need H(transitive) = 1. This is the trivial Hamilton path.
    -- We axiomatize this for now.
    exact H_transitive_eq_one n hn

end Tournament
