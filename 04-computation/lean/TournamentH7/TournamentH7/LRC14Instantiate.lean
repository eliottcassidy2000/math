/-
klein-2026-07-02-S109 (HYP-4014) — the concrete instantiation of the strong-induction
wrapper: `farCount` = the spread-20 cut (matching the censused window packs, max ≤ 20 =
6084 kernel-gated rows), so the endgame's remaining obligations become exactly two
CONCRETE statements: `census20` (bounded tuples clear the G2 floor — the ledger/table
leg) and `peel20` (a far element peels with floor transport — the rate-lemma leg,
F3-sharp's Lean face). This file adds NO axioms and NO sorries: it fixes the vocabulary
and re-exports the wrapper at the concrete cut.
-/
import Mathlib
import TournamentH7.LRC14CompletenessSurface

namespace LonelyRunner.LRC14

open TournamentH7 LonelyRunner.LRC14.Completeness

/-- The concrete far-element counter: entries beyond the censused window `20`. -/
def farCount20 (v : Fin 13 → ℤ) : ℕ :=
  (Finset.univ.filter fun i => 20 < |v i|).card

/-- `farCount20 v = 0` iff every entry is within the censused window. -/
theorem farCount20_eq_zero_iff (v : Fin 13 → ℤ) :
    farCount20 v = 0 ↔ ∀ i, |v i| ≤ 20 := by
  unfold farCount20
  rw [Finset.card_eq_zero, Finset.filter_eq_empty_iff]
  constructor
  · intro h i
    have := h (Finset.mem_univ i)
    omega
  · intro h i _
    have := h i
    omega

/-- **The concrete endgame surface**: LRC(14) from exactly
(1) `census20` — every nonzero 13-tuple inside the window clears the G2 floor
    (the table/ledger leg over the 6084 censused rows and their shape classes);
(2) `peel20` — any tuple with an entry beyond 20 peels one far element with floor
    transport (the rate-lemma leg; F3-sharp);
(3) `hpartA` — positive G2 implies reach (kps's dissolution).
No other obligations remain. -/
theorem lrc14_concrete
    (census20 : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → (∀ i, |v i| ≤ 20) →
      (witnessMP : ℝ) ≤ witnessG2 (shapeOf v))
    (peel20 : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → 0 < farCount20 v →
      ∃ v' : Fin 13 → ℤ, (∀ i, v' i ≠ 0) ∧ farCount20 v' < farCount20 v ∧
        ((witnessMP : ℝ) ≤ witnessG2 (shapeOf v') →
         (witnessMP : ℝ) ≤ witnessG2 (shapeOf v)))
    (hpartA : ∀ v : Fin 13 → ℤ, 0 < witnessG2 (shapeOf v) →
      (1 : ℝ) / 14 ≤ Mreach v) :
    LRC14Statement :=
  lrc14_of_census_peel_partA farCount20
    (fun v hv h0 => census20 v hv ((farCount20_eq_zero_iff v).mp h0))
    peel20 hpartA

end LonelyRunner.LRC14
