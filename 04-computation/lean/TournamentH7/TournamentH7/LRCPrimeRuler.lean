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

/-- **The prime clean ruler (general).**  For any prime `q ≤ 14`, the safe band at `q` is `[1, q−1]` = the
nonzero residues mod `q`.  If no speed is divisible by `q`, then every runner clears the band at every
multiplier `p ∈ (0,q)` (since `q ∤ p` and `q ∤ v i`, so `q ∤ v i · p`), giving `bandCount ≡ 0`: `q` is a
clean ruler, and hence supplies the per-family `hB5` witness `0 < B5 v q`. -/
theorem cleanRuler_of_prime_not_dvd (v : Fin 13 → ℤ) (q : ℕ)
    (hq : q.Prime) (hq14 : q ≤ 14) (h : ∀ i, ¬ (q : ℤ) ∣ v i) :
    CleanRuler v := by
  have hq2 : 2 ≤ q := hq.two_le
  have hqZ : Prime (q : ℤ) := Nat.prime_iff_prime_int.mp hq
  have hqZpos : (0 : ℤ) < (q : ℤ) := by exact_mod_cast (by omega : 0 < q)
  have hqZne : (q : ℤ) ≠ 0 := ne_of_gt hqZpos
  have hq14Z : (q : ℤ) ≤ 14 := by exact_mod_cast hq14
  -- every multiplier `p ∈ (0,q)` is fully in-band, so `bandCount = 0`
  have key : ∀ p ∈ Finset.Ioo 0 q, bandCount v q p = 0 := by
    intro p hp
    rw [Finset.mem_Ioo] at hp
    obtain ⟨hp0, hplt⟩ := hp
    rw [bandCount, Finset.card_eq_zero, Finset.filter_eq_empty_iff]
    intro i _
    rw [not_not]
    -- `q ∤ p` since `0 < p < q`
    have hpdvd : ¬ (q : ℤ) ∣ (p : ℤ) := by
      intro hd
      have hle : (q : ℤ) ≤ (p : ℤ) := Int.le_of_dvd (by exact_mod_cast hp0) hd
      have : q ≤ p := by exact_mod_cast hle
      omega
    -- so `q ∤ v i · p` (q prime, divides neither factor)
    have hnd : ¬ (q : ℤ) ∣ (v i * (p : ℤ)) := by
      rw [hqZ.dvd_mul, not_or]; exact ⟨h i, hpdvd⟩
    -- the residue is nonzero and in `[0,q)`, hence in the band `[1, q−1]`
    have hr0 : 0 ≤ (v i * (p : ℤ)) % (q : ℤ) := Int.emod_nonneg _ hqZne
    have hrlt : (v i * (p : ℤ)) % (q : ℤ) < (q : ℤ) := Int.emod_lt_of_pos _ hqZpos
    have hrne : (v i * (p : ℤ)) % (q : ℤ) ≠ 0 := fun hc => hnd (Int.dvd_of_emod_eq_zero hc)
    show inBand v q p i
    unfold inBand
    refine ⟨?_, ?_⟩ <;> omega
  -- assemble the clean ruler at `q`
  refine ⟨q, by omega, ?_, ?_⟩
  · intro p hp; rw [key p hp]; exact Nat.zero_le 5
  · -- liveCount = |Ioo 0 q| = q − 1 > 0
    have hlc : liveCount v q = (Finset.Ioo 0 q).card := by
      rw [liveCount, Finset.filter_true_of_mem key]
    rw [hlc, Nat.card_Ioo]; omega

/-- The `q = 13` instance (opus-S227): no speed divisible by 13 ⟹ clean ruler. -/
theorem cleanRuler_of_not_dvd_13 (v : Fin 13 → ℤ) (h : ∀ i, ¬ (13 : ℤ) ∣ v i) :
    CleanRuler v :=
  cleanRuler_of_prime_not_dvd v 13 (by norm_num) (by norm_num) h

/-- The `q = 11` instance. -/
theorem cleanRuler_of_not_dvd_11 (v : Fin 13 → ℤ) (h : ∀ i, ¬ (11 : ℤ) ∣ v i) :
    CleanRuler v :=
  cleanRuler_of_prime_not_dvd v 11 (by norm_num) (by norm_num) h

/-- The `q = 7` instance. -/
theorem cleanRuler_of_not_dvd_7 (v : Fin 13 → ℤ) (h : ∀ i, ¬ (7 : ℤ) ∣ v i) :
    CleanRuler v :=
  cleanRuler_of_prime_not_dvd v 7 (by norm_num) (by norm_num) h

/-- **The composite (opus-S228).**  If the family avoids `0 mod q` for *some* prime `q ∈ {7, 11, 13}`, it has
a clean ruler.  So the only residuals NOT discharged by a prime clean ruler are those with a speed divisible
by *each* of 7, 11, 13 — the structured hard core handled by the pair-sum ruler (kps THM-707). -/
theorem cleanRuler_of_avoids_some_prime (v : Fin 13 → ℤ)
    (h : (∀ i, ¬ (7 : ℤ) ∣ v i) ∨ (∀ i, ¬ (11 : ℤ) ∣ v i) ∨ (∀ i, ¬ (13 : ℤ) ∣ v i)) :
    CleanRuler v := by
  rcases h with h | h | h
  · exact cleanRuler_of_not_dvd_7 v h
  · exact cleanRuler_of_not_dvd_11 v h
  · exact cleanRuler_of_not_dvd_13 v h

end LRC14Concrete
end LonelyRunner

#print axioms LonelyRunner.LRC14Concrete.cleanRuler_of_prime_not_dvd
#print axioms LonelyRunner.LRC14Concrete.cleanRuler_of_not_dvd_13
#print axioms LonelyRunner.LRC14Concrete.cleanRuler_of_avoids_some_prime
