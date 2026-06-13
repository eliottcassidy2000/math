/-
  TournamentH7.SmallTournaments — concrete tournament constructors

  Provides explicit tournament constructions and elementary lemmas.
  The transitive tournament is built generically; the 3-cycle is built
  by truth table.

  This module shows the formalisation is *computational* — the
  axioms (irreflexive, total, asymmetric) are checked at definition time.
-/

import TournamentH7.Basic
import TournamentH7.SCC
import TournamentH7.Tilings
import Mathlib.Tactic.FinCases

namespace Tournament

/-! ### The transitive tournament -/

/-- `transitiveTournament n` — the totally ordered tournament where
    higher-indexed vertices beat lower-indexed ones. -/
def transitiveTournament (n : ℕ) : Tournament n where
  arc i j := decide (j.val < i.val)
  irrefl i := by simp
  total i j hij := by
    have hij_val : i.val ≠ j.val := fun H => hij (Fin.ext H)
    by_cases h : j.val < i.val
    · left; simp [h]
    · right
      have : i.val < j.val := by omega
      simp [this]
  asym i j := by
    intro ⟨h1, h2⟩
    simp at h1 h2; omega

/-- The transitive tournament has the base path. -/
lemma transitive_hasBasePath (n : ℕ) :
    HasBasePath (transitiveTournament n) := by
  intro i hi
  show decide (i.val < i.val + 1) = true
  simp

/-! ### The 3-cycle by truth table -/

/-- Arc function of the 3-cycle: 0→1, 1→2, 2→0. -/
def threeCycleArc : Fin 3 → Fin 3 → Bool := fun i j =>
  (i.val = 0 ∧ j.val = 1) ∨ (i.val = 1 ∧ j.val = 2) ∨ (i.val = 2 ∧ j.val = 0)

/-- The 3-cycle tournament: 0 → 1 → 2 → 0. -/
def threeCycle : Tournament 3 where
  arc := threeCycleArc
  irrefl i := by
    unfold threeCycleArc
    fin_cases i <;> decide
  total i j hij := by
    unfold threeCycleArc
    fin_cases i <;> fin_cases j <;> simp_all (config := { decide := true })
  asym i j h := by
    have ⟨h1, h2⟩ := h
    unfold threeCycleArc at h1 h2
    fin_cases i <;> fin_cases j <;> simp_all (config := { decide := true })

/-- The 3-cycle is `IsRegular` (every score equals 1 = (3-1)/2). -/
lemma threeCycle_isRegular : IsRegular threeCycle := by
  intro v
  show 2 * threeCycle.outDegree v = 3 - 1
  unfold Tournament.outDegree
  show 2 * (Finset.univ.filter (fun w => threeCycle.arc v w = true)).card = 2
  fin_cases v
  all_goals
    show 2 * (Finset.univ.filter _).card = 2
    decide

/-! ### Score sequence of the transitive tournament -/

/-- The score *sequence* of `transitiveTournament n` is `(0, 1, 2, …, n-1)`.
    In particular, no vertex has score (n-1)/2 unless that's an integer
    AND only one vertex hits it — so the transitive tournament is *never*
    regular for n ≥ 2. -/
lemma transitive_not_regular (n : ℕ) (hn : 2 ≤ n) :
    ¬ IsRegular (transitiveTournament n) := by
  intro hreg
  -- vertex 0 has out-degree 0; vertex 1 has out-degree 1.  In a regular
  -- tournament, both would equal (n − 1)/2.  Contradiction unless n ≤ 1.
  have h0 := hreg ⟨0, by omega⟩
  have h1 := hreg ⟨1, by omega⟩
  -- 2 * (outDegree 0) = n - 1 and 2 * (outDegree 1) = n - 1.
  -- We use: outDegree v ≤ v.val (since arcs only go down).
  have hb0 : (transitiveTournament n).outDegree ⟨0, by omega⟩ = 0 := by
    have hn' : 0 < n := by omega
    show (Finset.univ.filter (fun w : Fin n =>
      (transitiveTournament n).arc ⟨0, hn'⟩ w = true)).card = 0
    rw [Finset.card_eq_zero]
    ext w
    simp only [Finset.mem_filter, Finset.mem_univ, true_and,
      Finset.notMem_empty, iff_false]
    show ¬ decide (w.val < (⟨0, hn'⟩ : Fin n).val) = true
    simp
  rw [hb0] at h0
  -- h0 : 2 * 0 = n - 1 ⟹ n = 1.
  omega

/-! ### Note on threeCycle and the base path

    The 3-cycle `threeCycle` is `0 → 1 → 2 → 0`, NOT the base-path
    direction `2 → 1 → 0`.  So `threeCycle` does *not* satisfy
    `HasBasePath` (which requires `(i+1) → i`).  The opposite tournament
    `op threeCycle` is `0 → 2 → 1 → 0`, which also doesn't satisfy
    HasBasePath (in fact one needs the LINEAR direction).  -/

/-! ### Examples of axiom-free results from SmallTournaments -/

/-- Transitive tournament `Trans_2` (n = 2) has scores (0, 1). -/
example : (transitiveTournament 2).outDegree ⟨0, by omega⟩ = 0 := by
  show (Finset.univ.filter
    (fun w : Fin 2 => (transitiveTournament 2).arc ⟨0, by omega⟩ w = true)).card = 0
  decide

example : (transitiveTournament 2).outDegree ⟨1, by omega⟩ = 1 := by
  show (Finset.univ.filter
    (fun w : Fin 2 => (transitiveTournament 2).arc ⟨1, by omega⟩ w = true)).card = 1
  decide

/-! ### Named tiny score facts -/

/-- The source of `Trans_3` has score 0. -/
lemma transitive_three_score_zero :
    (transitiveTournament 3).outDegree ⟨0, by omega⟩ = 0 := by
  show (Finset.univ.filter
    (fun w : Fin 3 => (transitiveTournament 3).arc ⟨0, by omega⟩ w = true)).card = 0
  decide

/-- The middle vertex of `Trans_3` has score 1. -/
lemma transitive_three_score_one :
    (transitiveTournament 3).outDegree ⟨1, by omega⟩ = 1 := by
  show (Finset.univ.filter
    (fun w : Fin 3 => (transitiveTournament 3).arc ⟨1, by omega⟩ w = true)).card = 1
  decide

/-- The sink of `Trans_3` in the base-path orientation has score 2. -/
lemma transitive_three_score_two :
    (transitiveTournament 3).outDegree ⟨2, by omega⟩ = 2 := by
  show (Finset.univ.filter
    (fun w : Fin 3 => (transitiveTournament 3).arc ⟨2, by omega⟩ w = true)).card = 2
  decide

/-- The transitive tournament on three vertices has score vector `(0, 1, 2)`. -/
lemma transitive_three_score_vector :
    (transitiveTournament 3).outDegree ⟨0, by omega⟩ = 0 ∧
    (transitiveTournament 3).outDegree ⟨1, by omega⟩ = 1 ∧
    (transitiveTournament 3).outDegree ⟨2, by omega⟩ = 2 := by
  exact ⟨transitive_three_score_zero,
    transitive_three_score_one,
    transitive_three_score_two⟩

end Tournament
