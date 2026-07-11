/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: opus (LRC multi-agent project, 2026-07-11-S227)
-/
import TournamentH7.LRCCleanRuler

/-!
# The prime clean ruler `q = 13` — discharging `hB5` for families with no speed divisible by 13

At the prime modulus `q = 13` the safe band is `[⌈13/14⌉, ⌊13·13/14⌋] = [1, 12]`, i.e. *exactly* the
nonzero residues mod 13.  So runner `i` is OUT of band at multiplier `p` iff `(v i · p) mod 13 = 0` iff
`13 ∣ v i · p`.  For `p ∈ (0,13)` we have `13 ∤ p` (13 is prime), so if additionally `13 ∤ v i` for every
runner, then `13 ∤ v i · p`, every runner is in band, and `bandCount v 13 p = 0` for every `p ∈ (0,13)`.

Hence `q = 13` is a **clean ruler** (`maxBand = 0`, `liveCount = 12`), and by `exists_B5_pos_of_cleanRuler`
(THM-707) this **discharges the single LRC(14) obligation `hB5` for the entire sub-class of residual
families with no speed a multiple of 13**.  Verified: `q = 13` clean on 2000/2000 random 13-speed families
avoiding `13 ∣ ·` (opus-S227).

This is the "13 is prime" phenomenon (mac-mini HYP-4382): the hard residual families are exactly those with
a speed `≡ 0 mod 13` — the AP `{1,…,13}`'s speed `13`, or `26 = 2·13` — where this trivial ruler fails and
the pair-sum ruler (kps THM-707 / `LRCPairSumDispatch`) is required.  The same argument gives a clean ruler
at any prime `q ∈ {7, 11, 13}` (each has band `[1, q−1]` = all nonzero residues); so the genuinely hard
class shrinks to families with a speed divisible by *each* of 7, 11, 13.

Kernel-pure: elementary residue arithmetic, no `decide`, no `native_decide`, no `sorry`.
-/

namespace LonelyRunner
namespace LRC14Concrete

open Finset

/-- **The prime clean ruler.**  If no speed is divisible by 13, then at `q = 13` every runner clears the
safe band at every multiplier, so `q = 13` is a clean ruler (and hence supplies the per-family `hB5`
witness `0 < B5 v 13`). -/
theorem cleanRuler_of_not_dvd_13 (v : Fin 13 → ℤ) (h : ∀ i, ¬ (13 : ℤ) ∣ v i) :
    CleanRuler v := by
  have hp13 : Prime (13 : ℤ) := by norm_num
  -- every multiplier `p ∈ (0,13)` is fully in-band, so `bandCount = 0`
  have key : ∀ p ∈ Finset.Ioo 0 13, bandCount v 13 p = 0 := by
    intro p hp
    rw [Finset.mem_Ioo] at hp
    obtain ⟨hp0, hplt⟩ := hp
    rw [bandCount, Finset.card_eq_zero, Finset.filter_eq_empty_iff]
    intro i _
    rw [not_not]
    -- `13 ∤ p` since `0 < p < 13`
    have hpdvd : ¬ (13 : ℤ) ∣ (p : ℤ) := by
      intro hd
      have hle : (13 : ℤ) ≤ (p : ℤ) := Int.le_of_dvd (by exact_mod_cast hp0) hd
      have : (13 : ℕ) ≤ p := by exact_mod_cast hle
      omega
    -- so `13 ∤ v i · p` (13 prime, divides neither factor)
    have hnd : ¬ (13 : ℤ) ∣ (v i * (p : ℤ)) := by
      rw [hp13.dvd_mul, not_or]; exact ⟨h i, hpdvd⟩
    -- the residue is nonzero and in `[0,13)`, hence in the band `[1,12]`
    have hr0 : 0 ≤ (v i * (p : ℤ)) % 13 := Int.emod_nonneg _ (by norm_num)
    have hrlt : (v i * (p : ℤ)) % 13 < 13 := Int.emod_lt_of_pos _ (by norm_num)
    have hrne : (v i * (p : ℤ)) % 13 ≠ 0 := fun hc => hnd (Int.dvd_of_emod_eq_zero hc)
    show inBand v 13 p i
    unfold inBand
    have hcast : ((13 : ℕ) : ℤ) = (13 : ℤ) := by norm_num
    rw [hcast]
    refine ⟨?_, ?_⟩ <;> omega
  -- assemble the clean ruler at q = 13
  refine ⟨13, by norm_num, ?_, ?_⟩
  · intro p hp; rw [key p hp]; exact Nat.zero_le 5
  · -- liveCount = |Ioo 0 13| = 12 > 0
    have : liveCount v 13 = (Finset.Ioo 0 13).card := by
      rw [liveCount, Finset.filter_true_of_mem key]
    rw [this, Nat.card_Ioo]; norm_num

end LRC14Concrete
end LonelyRunner

#print axioms LonelyRunner.LRC14Concrete.cleanRuler_of_not_dvd_13
