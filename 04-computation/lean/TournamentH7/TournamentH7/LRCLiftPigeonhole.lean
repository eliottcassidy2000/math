/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: opus (LRC multi-agent project, 2026-07-05-S80)
-/
import Mathlib.Data.Finset.Card
import Mathlib.Data.Nat.GCD.Basic
import Mathlib.Tactic

/-!
# The l ≥ 7 structural lemmas (HYP-4107 in Lean; feeds `tower_step_12`)

The two 3-line facts behind the l ≥ 7 lift-rigidity closure (opus-S78):

  * `exists_high_of_card_ge_seven` — any ≥ 7 lifted positions include one of the six
    unique-multiple indices `{6, …, 11}` (0-indexed residues 7..12): pigeonhole against the
    ≤ 5 unlifted slots.
  * `dvd_lift_height` — a lifted unique-multiple residue `r ∈ [7, 12]` keeps modulus `r`
    covered only if `r ∣ k`; hence its lifted value is `≥ r + 13r = 14r ≥ 98`.

Wiring note (discovered attempting the `tower_step_12` instantiation, for the fee-table
owners): at `l = 7` the window half-width is `δ = (1/6 − 1/13)/12 = 7/936`, and the GENERAL
fee `6ρ/w` (`teethR_mass`) for a mid-scale top `w ∈ [29, 97]` (too fast for sub-gap
`13·2δ·w ≤ 11`, too slow for a thin tooth to matter) already exceeds the whole budget
`2δ = 7/468`.  The S78 mathematics is safe (numeric floor `3/19`), but the `T_l` fee table
must handle the mid-scale stratum separately — e.g. per-top windows, or the citation margin
taken at `1/(13−l')` for the sub-gap tops only.  Handed to mac-mini/klein with this file.
-/

namespace LonelyRunner
namespace LiftPigeonhole

/-- **Pigeonhole**: a set of ≥ 7 positions in `Fin 12` meets the six high indices
`{6, …, 11}` (the unique-multiple residues 7..12). -/
theorem exists_high_of_card_ge_seven (T : Finset (Fin 12)) (hT : 7 ≤ T.card) :
    ∃ i ∈ T, 6 ≤ (i : ℕ) := by
  by_contra hcon
  push_neg at hcon
  have hsub : T ⊆ Finset.univ.filter (fun i : Fin 12 => (i : ℕ) < 6) := by
    intro i hi
    have := hcon i hi
    exact Finset.mem_filter.mpr ⟨Finset.mem_univ i, by omega⟩
  have hcard : (Finset.univ.filter (fun i : Fin 12 => (i : ℕ) < 6)).card = 6 := by decide
  have := Finset.card_le_card hsub
  omega

/-- **Height forcing**: for a unique-multiple residue `r ∈ [7, 12]`, the lifted value
`r + 13k` is a multiple of `r` only if `r ∣ k` — so a sieve-surviving lift has
`k ≥ r`, hence value `≥ 14r ≥ 98`. -/
theorem dvd_lift_height {r k : ℕ} (hr7 : 7 ≤ r) (hr12 : r ≤ 12)
    (hdvd : r ∣ (r + 13 * k)) : r ∣ k := by
  have h13 : r ∣ 13 * k := (Nat.dvd_add_right (dvd_refl r)).mp hdvd
  have hco : Nat.Coprime r 13 := by
    interval_cases r <;> decide
  exact (Nat.Coprime.dvd_of_dvd_mul_left hco h13)

/-- The forced height, spelled out: a sieve-surviving lifted unique-multiple residue sits at
value at least `14r ≥ 98`. -/
theorem lift_value_ge {r k : ℕ} (hr7 : 7 ≤ r) (hr12 : r ≤ 12) (hk : 0 < k)
    (hdvd : r ∣ (r + 13 * k)) : 14 * r ≤ r + 13 * k := by
  have hrk := dvd_lift_height hr7 hr12 hdvd
  have hkr : r ≤ k := Nat.le_of_dvd hk hrk
  nlinarith

end LiftPigeonhole
end LonelyRunner
