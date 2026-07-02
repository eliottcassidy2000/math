/-
klein-2026-07-02-S111 (HYP-4016) — the composed window leg at W = 20: the REPEAT case is
formally eliminated through the wiring via the LRC(≤13) citation node, so `hwindow` at
the censused cut reduces to exactly the DISTINCT census (= the 6084 pack rows' quantified
form). Surface shrink: hwindow20 ⟸ {hdistinct20} + LRCUpTo13.
-/
import Mathlib
import TournamentH7.LRC14WindowWiring
import TournamentH7.LRC13Citation

namespace LonelyRunner.LRC14

open LonelyRunner

/-- A 13-tuple is strictly monotone (hence has distinct entries). -/
def StrictMono13 (u : Fin 13 → ℤ) : Prop := StrictMono u

/-- **The composed window leg at the censused cut**: bounded-window loneliness from
(1) the DISTINCT census — every strictly monotone positive primitive covering tuple
    ≤ 20 is lonely (the 6084 kernel-gated pack rows, quantified), and
(2) the LRC(≤13) citation node — which absorbs every monotone tuple with a REPEAT
    (≤ 12 distinct speeds ⟹ lonely at 1/13 ≥ 1/14, `lonely14_of_repeat`). -/
theorem hwindow20
    (cite : LRCUpTo13)
    (hdistinct20 : ∀ u : Fin 13 → ℤ, (∀ i, 0 < u i) → StrictMono u →
      (∀ i, u i ≤ 20) → CoveringFamily u → tupleGcd u = 1 → ∃ t : ℝ, Lonely 14 u t) :
    ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → (∀ i, |v i| ≤ 20) → ∃ t : ℝ, Lonely 14 v t := by
  apply hwindow_of_normalized_census 20
  intro u hpos hmono hle hcov hprim
  by_cases hrep : ∃ i j : Fin 13, i ≠ j ∧ u i = u j
  · obtain ⟨i, j, hij, hval⟩ := hrep
    exact lonely14_of_repeat cite u (fun k => (hpos k).ne') hij hval
  · push_neg at hrep
    have hsm : StrictMono u := by
      intro a b hab
      rcases lt_or_eq_of_le (hmono (le_of_lt hab)) with h | h
      · exact h
      · exact absurd h (hrep a b (ne_of_lt hab))
    exact hdistinct20 u hpos hsm hle hcov hprim

end LonelyRunner.LRC14
