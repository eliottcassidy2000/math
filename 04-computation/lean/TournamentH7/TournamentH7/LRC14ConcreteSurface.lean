/-
  TournamentH7.LRC14ConcreteSurface — the spread-20 endgame in DISCHARGEABLE vocabulary
  (kind-pasteur-2026-07-02-S16, HYP-3973).  Single-writer satellite.

  klein-S109's `lrc14_concrete` fixes the induction counter (`farCount20`) but phrases its
  `census20`/`peel20`/`hpartA` obligations over the skeleton's OPAQUE `witnessG2`/`shapeOf`.
  Opaque constants have no provable properties, so those parameters cannot be discharged by
  the pack rows (which prove `Lonely` facts) — the same trap flagged in the S9 frontier note.

  THIS FILE: the SAME spread-20 surface with every obligation in final `Lonely` vocabulary,
  so each is genuinely dischargeable by existing machinery:

    * `census20` — exactly the S15 census class at `W = 20` (monotone, positive, primitive,
      covering, bounded): the 6084 kernel-gated pack rows + the repeat sweep;
    * `peel20`   — far-element peel with loneliness transport: the S14 peel-gate leg
      (`damped_peel`/`goodRegion2_pos_of_peel` + threshold data).

  `lrc14_concrete_surface` : census20 → peel20 → LRC14Statement, kernel-pure glue on top of
  the S10 strong-induction shell and the S15 window wiring.  No opaque symbols anywhere.
-/
import TournamentH7.LRC14WindowWiring

namespace LonelyRunner
namespace LRC14

/-- The concrete far-element counter at the censused window `W = 20`. -/
def farCountW (W : ℤ) (v : Fin 13 → ℤ) : ℕ :=
  (Finset.univ.filter fun i => W < |v i|).card

theorem farCountW_eq_zero_iff (W : ℤ) (v : Fin 13 → ℤ) :
    farCountW W v = 0 ↔ ∀ i, |v i| ≤ W := by
  unfold farCountW
  rw [Finset.card_eq_zero, Finset.filter_eq_empty_iff]
  constructor
  · intro h i
    have := h (Finset.mem_univ i)
    omega
  · intro h i _
    have := h i
    omega

/-- **The spread-`W` endgame surface, fully dischargeable**: LRC(14) from
(1) the normalized bounded census at `W` (the packs' exact class), and
(2) the far-element peel with loneliness transport (the peel-gate leg).
Both obligations are in final `Lonely` vocabulary — no opaque symbols. -/
theorem lrc14_of_census_peel_concrete (W : ℤ)
    (census : ∀ u : Fin 13 → ℤ, (∀ i, 0 < u i) → Monotone u → (∀ i, u i ≤ W) →
      CoveringFamily u → tupleGcd u = 1 → ∃ t : ℝ, Lonely 14 u t)
    (peel : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → 0 < farCountW W v →
      ∃ v' : Fin 13 → ℤ, (∀ i, v' i ≠ 0) ∧ farCountW W v' < farCountW W v ∧
        ((∃ t : ℝ, Lonely 14 v' t) → ∃ t : ℝ, Lonely 14 v t)) :
    LRC14Statement := by
  apply lrc14_of_lonely_census_peel (farCountW W)
  · -- base: inside the window, the S15 wiring reduces to the census class
    intro v hv h0
    exact hwindow_of_normalized_census W census v hv ((farCountW_eq_zero_iff W v).mp h0)
  · exact peel

/-- The spread-20 instance: the counter and window of the current pack corpus
(6084 kernel-gated rows at `max ≤ 20`). -/
theorem lrc14_of_census20_peel20
    (census20 : ∀ u : Fin 13 → ℤ, (∀ i, 0 < u i) → Monotone u → (∀ i, u i ≤ 20) →
      CoveringFamily u → tupleGcd u = 1 → ∃ t : ℝ, Lonely 14 u t)
    (peel20 : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → 0 < farCountW 20 v →
      ∃ v' : Fin 13 → ℤ, (∀ i, v' i ≠ 0) ∧ farCountW 20 v' < farCountW 20 v ∧
        ((∃ t : ℝ, Lonely 14 v' t) → ∃ t : ℝ, Lonely 14 v t)) :
    LRC14Statement :=
  lrc14_of_census_peel_concrete 20 census20 peel20

end LRC14
end LonelyRunner
