/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: kind-pasteur (LRC multi-agent project, 2026-07-11-S127 cont.31)
-/
import TournamentH7.LRCCleanRuler

/-!
# The prime clean ruler — every prime `P ≤ 13` discharges `hB5` for families avoiding `P`

opus-S227's THM-709 showed the modulus `q = 13` discharges `hB5` for residual families with no speed
divisible by `13` (via a margin/pinning argument). This file generalizes it to **every prime `P ≤ 13`**, as a
two-step corollary of the clean-ruler lemma `b5_pos_of_clean` (`LRCCleanRuler`):

> **`b5_pos_of_prime_ndvd`**: for a prime `P ≤ 13`, if `P ∤ vᵢ` for all `i`, then `q = P` is a *perfectly*
> clean ruler — every multiplier is fully safe (`bandCount = 0`), so `liveCount = P − 1` and `B5 v P = P − 1 > 0`.

The mechanism is transparent: at modulus `P ≤ 13` the safe band `[P/14, 13P/14]` contains **all** nonzero
residues `{1,…,P−1}`, so a runner is unsafe only at a residue `0`, i.e. only if `P ∣ vᵢ·p'`. With `P` prime,
`P ∤ vᵢ`, and `0 < p' < P`, that never happens. So `q = 13` (opus) is the top instance of a six-member family
`P ∈ {2,3,5,7,11,13}`, and a residual family is discharged whenever *some* prime `≤ 13` divides none of its
speeds.
-/

namespace LonelyRunner
namespace LRC14Concrete

open Finset

/-- At a prime modulus `P ≤ 13` with no speed divisible by `P`, every nonzero multiplier `p' ∈ (0,P)` is
fully safe: `bandCount v P p' = 0`.  (The `1/7` safe band at modulus `P ≤ 13` covers all nonzero residues.) -/
lemma bandCount_eq_zero_of_prime_ndvd (v : Fin 13 → ℤ) (P : ℕ)
    (hP : P.Prime) (hP13 : P ≤ 13) (hnd : ∀ i, ¬ (P : ℤ) ∣ v i)
    (p' : ℕ) (hp' : p' ∈ Finset.Ioo 0 P) :
    bandCount v P p' = 0 := by
  rw [bandCount, Finset.card_eq_zero, Finset.filter_eq_empty_iff]
  intro i _
  rw [not_not]
  rw [Finset.mem_Ioo] at hp'
  obtain ⟨hp'0, hp'P⟩ := hp'
  -- `P` does not divide `v i * p'`
  have hPz : Prime (P : ℤ) := Nat.prime_iff_prime_int.mp hP
  have hnp' : ¬ (P : ℤ) ∣ (p' : ℤ) := by
    intro hdvd
    rw [Int.natCast_dvd_natCast] at hdvd
    exact absurd (Nat.le_of_dvd hp'0 hdvd) (by omega)
  have hnmul : ¬ (P : ℤ) ∣ (v i * (p' : ℤ)) := by
    rw [hPz.dvd_mul]; push_neg; exact ⟨hnd i, hnp'⟩
  -- so the residue `r = (v i * p') % P` lies in `[1, P-1]`
  set r : ℤ := (v i * (p' : ℤ)) % (P : ℤ) with hr
  have hPpos : (0 : ℤ) < P := by exact_mod_cast hP.pos
  have hr0 : 0 ≤ r := Int.emod_nonneg _ (by positivity)
  have hrP : r < P := Int.emod_lt_of_pos _ hPpos
  have hrne : r ≠ 0 := by
    rw [hr, Ne, ← Int.dvd_iff_emod_eq_zero]; exact hnmul
  have hr1 : 1 ≤ r := by omega
  have hP13' : (P : ℤ) ≤ 13 := by exact_mod_cast hP13
  -- the band inequalities `P ≤ 14 r ≤ 13 P` follow by omega
  refine ⟨?_, ?_⟩
  · show (P : ℤ) ≤ 14 * r; omega
  · show 14 * r ≤ 13 * (P : ℤ); omega

/-- **The prime clean ruler (generalizing opus-S227's THM-709 `q = 13`).**  For a prime `P ≤ 13` with no
speed divisible by `P`, `q = P` is a clean ruler and `0 < B5 v P` (in fact `B5 v P = liveCount v P = P − 1`). -/
theorem b5_pos_of_prime_ndvd (v : Fin 13 → ℤ) (P : ℕ)
    (hP : P.Prime) (hP13 : P ≤ 13) (hnd : ∀ i, ¬ (P : ℤ) ∣ v i) :
    0 < B5 v P := by
  have hclean : ∀ p' ∈ Finset.Ioo 0 P, bandCount v P p' ≤ 5 := fun p' hp' => by
    rw [bandCount_eq_zero_of_prime_ndvd v P hP hP13 hnd p' hp']; exact Nat.zero_le 5
  have hlive : 0 < liveCount v P := by
    rw [liveCount, Finset.card_pos]
    refine ⟨1, ?_⟩
    rw [Finset.mem_filter, Finset.mem_Ioo]
    have h1P : 1 < P := hP.one_lt
    exact ⟨⟨Nat.one_pos, h1P⟩, bandCount_eq_zero_of_prime_ndvd v P hP hP13 hnd 1
      (Finset.mem_Ioo.mpr ⟨Nat.one_pos, h1P⟩)⟩
  exact b5_pos_of_clean v P hclean hlive

/-- **The prime clean ruler supplies the `hB5` witness** for the sub-class avoiding `P`. -/
theorem exists_B5_pos_of_prime_ndvd (v : Fin 13 → ℤ) (P : ℕ)
    (hP : P.Prime) (hP13 : P ≤ 13) (hnd : ∀ i, ¬ (P : ℤ) ∣ v i) :
    ∃ q : ℕ, 0 < q ∧ 0 < B5 v q :=
  ⟨P, hP.pos, b5_pos_of_prime_ndvd v P hP hP13 hnd⟩

end LRC14Concrete
end LonelyRunner

#print axioms LonelyRunner.LRC14Concrete.b5_pos_of_prime_ndvd
