/-
  TournamentH7.LRCFinsetBridge — connect `LRC14_of_INV`'s Lonely-form to the
  `Finset ℕ` / `t ∈ [0,1]` target `LRC14.LRC14` of `LRC14Ledger`.
  boxeph-2026-07-18-S109.

  `LRC14_of_INV` (S108) concludes `∀ v : Fin 13 → ℤ, (0<v) → ∃ t, Lonely 14 v t`.
  The ledger target quantifies over a `Finset ℕ` of card 13, asks the lonely time in
  `[0,1]`, and uses `∀ w ∈ W`.  Two elementary transfers:
    · enumerate `W` as `v : Fin 13 → ℤ` via `W.equivFinOfCardEq`;
    · loneliness is shift-invariant, so replace `t` by `Int.fract t ∈ [0,1)`.

  `LRC14_finset_of_INV` (PROVED, kernel-pure): `LRC(≤13) + INV ⟹ LRC14.LRC14`.
-/
import Mathlib
import TournamentH7.LonelyRunner
import TournamentH7.LRCMSplit
import TournamentH7.LRC14Ledger

open Set

namespace LonelyRunner

/-- **Shift invariance.** A `1/n`-lonely time reduces to `[0,1)`: `Lonely n v t` gives
`Lonely n v (Int.fract t)`, since `v_i·fract t` differs from `v_i·t` by the integer
`v_i·⌊t⌋`, absorbed by the integer `m`. -/
theorem lonely_fract {ι : Type*} {n : ℕ} (v : ι → ℤ) (t : ℝ)
    (h : Lonely n v t) : Lonely n v (Int.fract t) := by
  intro i m
  have hrw : (v i : ℝ) * Int.fract t - m
      = (v i : ℝ) * t - ((m + v i * ⌊t⌋ : ℤ) : ℝ) := by
    rw [Int.fract]; push_cast; ring
  rw [hrw]; exact h i (m + v i * ⌊t⌋)

/-- **The Finset target from INV (PROVED).**  `LRC14_of_INV`'s `Fin 13 → ℤ` conclusion
transfers to the `Finset ℕ`, `t ∈ [0,1]` statement `LRC14.LRC14`: enumerate `W`, apply
`LRC14_of_INV`, and reduce the lonely time to `[0,1)` by `lonely_fract`. -/
theorem LRC14_finset_of_INV (cite : LRCUpTo13)
    (INV : ∀ v : Fin 13 → ℤ, (∀ i, 0 < v i) → (¬ ∃ t, Lonely 13 v t) →
      ∃ vstar : Fin 13, ∀ i, i ≠ vstar → 13 * v i ≤ v vstar) :
    LRC14.LRC14 := by
  intro W hWpos hWcard
  -- enumerate `W` as `v : Fin 13 → ℤ`
  set e : Fin 13 ≃ ↥W := (W.equivFinOfCardEq hWcard).symm with he
  set v : Fin 13 → ℤ := fun i => ((e i : ℕ) : ℤ) with hv
  have hvpos : ∀ i, 0 < v i := by
    intro i
    have hpos := hWpos (e i) (e i).2
    simp only [hv]; exact_mod_cast hpos
  -- LRC(14) in the Lonely form
  obtain ⟨t, ht⟩ := LRC14_of_INV cite INV v hvpos
  -- reduce the lonely time to [0,1)
  refine ⟨Int.fract t, mem_Icc.mpr ⟨Int.fract_nonneg t, le_of_lt (Int.fract_lt_one t)⟩, ?_⟩
  intro w hw a
  set i : Fin 13 := e.symm ⟨w, hw⟩ with hi
  have hei : (e i : ℕ) = w := by simp only [hi, Equiv.apply_symm_apply]
  have hvi : (v i : ℝ) = (w : ℝ) := by
    simp only [hv, hei]; push_cast; ring
  have hl := lonely_fract v t ht i a
  rwa [hvi] at hl

end LonelyRunner

#print axioms LonelyRunner.lonely_fract
#print axioms LonelyRunner.LRC14_finset_of_INV
