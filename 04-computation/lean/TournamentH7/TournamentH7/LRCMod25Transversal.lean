/-
  TournamentH7.LRCMod25Transversal — the mod-25 clearing rotation EXISTS when the
  family MISSES a ±-pair (mac-mini-2026-07-06-S33, HYP-4642).

  Supplies the EXISTENCE half that kps-S41's `mod25_covering_floor`
  (`LRCMod25Floor`) assumed as a hypothesis: it takes a rotation `c` off `{0,±1}`
  mod 25 and concludes `M ≥ 2/25`, but did not produce `c`.  Here `c` is produced
  explicitly.

  The mod-25 crux dichotomy (mac-mini S32 / concurrent pair-blocking, kps-S42):
  a `(ℤ/25)*` rotation clearing every speed off `{0,±1}` EXISTS  ⟺  the family is
  NOT a full transversal mod 25 — i.e. its unit-speeds MISS one of the ten antipodal
  pairs `{a, -a}`.  This file proves the `⟸` direction with the closed-form witness

        miss the pair `{a, -a}`  (no speed `≡ 0`)   ⟹   `c := a⁻¹ mod 25` clears,

  because `v·c ≡ 0,1,-1 (mod 25) ⟺ v ≡ 0, a, -a (mod 25)`, all excluded.  Composed
  with `loose_of_mod25_covering` this PROVES the non-transversal (miss-a-pair) branch
  of the LRC(14) crux (G): a 12-family that is not a full transversal mod 25 has
  `M ≥ 2/25` (above the first gap `(1/13, 2/25)`).  The remaining crux is then the
  full-transversal residual (pair-blocker ⟹ dilated AP).

  Kernel-pure; no `sorry`, no `native_decide`.
-/
import Mathlib
import TournamentH7.LRCMod25Floor

namespace LonelyRunner
namespace Mod25Transversal

open Mod25Floor

/-- **The clearing rotation, from a missed pair.**
If `c` is a mod-25 inverse of `a` (`a·c ≡ 1`), no speed is `≡ 0 (mod 25)`, and the
speeds miss the antipodal pair `{a, -a}` mod 25, then `c` rotates every speed into the
safe band `[2, 23]` mod 25.  This is exactly the hypothesis of
`mod25_covering_floor`. -/
theorem covering_of_misses_pair {ι : Type*} (v : ι → ℤ) (a c : ℤ)
    (hinv : (a : ZMod 25) * (c : ZMod 25) = 1)
    (hnz : ∀ i, (v i : ZMod 25) ≠ 0)
    (hmiss : ∀ i, (v i : ZMod 25) ≠ (a : ZMod 25) ∧ (v i : ZMod 25) ≠ -(a : ZMod 25)) :
    ∀ i, 2 ≤ (v i * c) % 25 ∧ (v i * c) % 25 ≤ 23 := by
  intro i
  -- Step 1: in `ZMod 25`, `(v i)·c ∉ {0, 1, -1}`, since `v i ∉ {0, a, -a}` and `a·c = 1`.
  have h0 : (v i : ZMod 25) * (c : ZMod 25) ≠ 0 := by
    intro h; exact hnz i (by linear_combination (a : ZMod 25) * h - (v i : ZMod 25) * hinv)
  have h1 : (v i : ZMod 25) * (c : ZMod 25) ≠ 1 := by
    intro h; exact (hmiss i).1 (by linear_combination (a : ZMod 25) * h - (v i : ZMod 25) * hinv)
  have h2 : (v i : ZMod 25) * (c : ZMod 25) ≠ -1 := by
    intro h; exact (hmiss i).2 (by linear_combination (a : ZMod 25) * h - (v i : ZMod 25) * hinv)
  -- Step 2: transfer to non-divisibility of the integer `v i * c`.
  have d0 : ¬ (((25 : ℕ) : ℤ) ∣ v i * c) := by
    rw [← ZMod.intCast_zmod_eq_zero_iff_dvd]; push_cast; exact h0
  have d1 : ¬ (((25 : ℕ) : ℤ) ∣ (v i * c - 1)) := by
    rw [← ZMod.intCast_zmod_eq_zero_iff_dvd]; push_cast; rw [sub_eq_zero]; exact h1
  have d2 : ¬ (((25 : ℕ) : ℤ) ∣ (v i * c + 1)) := by
    rw [← ZMod.intCast_zmod_eq_zero_iff_dvd]; push_cast; rw [add_eq_zero_iff_eq_neg]; exact h2
  -- Step 3: `omega` reads `∤` and `%` and closes the band.
  simp only [Nat.cast_ofNat] at d0 d1 d2
  omega

/-- **The non-transversal (miss-a-pair) branch of the crux, PROVED.**
A 12-speed family that misses an antipodal unit-pair `{a, -a}` mod 25 (and has no
multiple of 25) is LOOSE: there is a common `2/25`-margin point `t = c/25`, so
`M ≥ 2/25` — above the first gap. -/
theorem loose_of_misses_pair {ι : Type*} (v : ι → ℤ) (a c : ℤ)
    (hinv : (a : ZMod 25) * (c : ZMod 25) = 1)
    (hnz : ∀ i, (v i : ZMod 25) ≠ 0)
    (hmiss : ∀ i, (v i : ZMod 25) ≠ (a : ZMod 25) ∧ (v i : ZMod 25) ≠ -(a : ZMod 25)) :
    ∃ t : ℝ, ∀ i, ∀ m : ℤ, (2 : ℝ) / 25 ≤ |(v i : ℝ) * t - m| :=
  loose_of_mod25_covering v c (covering_of_misses_pair v a c hinv hnz hmiss)

#print axioms covering_of_misses_pair
#print axioms loose_of_misses_pair

end Mod25Transversal
end LonelyRunner
