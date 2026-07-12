/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: kind-pasteur (LRC multi-agent project, 2026-07-11-S127 cont.38)
-/
import TournamentH7.LRCCleanRuler

/-!
# The general clean divisibility ruler for `q ≤ 14` (generalizing THM-712 to composites)

For `q ≤ 14` the safe band `[q/14, 13q/14]` contains **every** nonzero residue mod `q`, so a runner is unsafe
at multiplier `p` exactly when `q ∣ vᵢ·p`.  Writing `d = q / gcd(q,p)` (a divisor of `q` with `d ≥ 2` when
`0 < p < q`), Euclid gives `q ∣ vᵢ·p ⟹ d ∣ vᵢ`, so `bandCount(v,q,p) ≤ #{i : d ∣ vᵢ}`.  Hence:

> **`b5_pos_of_div_clean`** — if `q ∤ vᵢ` for all `i` (LIVE: `p = 1` covers) and every divisor `d ≥ 2` of `q`
> divides at most `5` of the speeds (CLEAN: `maxBand ≤ 5`), then `0 < B5 v q`.

This extends the prime clean ruler (`b5_pos_of_prime_ndvd`, THM-712) from primes to the composite moduli
`q ∈ {8,9,10,12,14}` — the tier-1 "rigorous divisibility" part of the bounded-window covering (HYP-6035).
-/

namespace LonelyRunner
namespace LRC14Concrete

open Finset

/-- For `2 ≤ q ≤ 14`, a runner is unsafe (`¬ inBand`) exactly when `q ∣ vᵢ·p`. -/
lemma not_inBand_iff_dvd (v : Fin 13 → ℤ) (q p : ℕ) (hq2 : 2 ≤ q) (hq14 : q ≤ 14) (i : Fin 13) :
    ¬ inBand v q p i ↔ (q : ℤ) ∣ v i * (p : ℤ) := by
  unfold inBand
  have hq0 : (0 : ℤ) < q := by exact_mod_cast (by omega : 0 < q)
  set r : ℤ := (v i * (p : ℤ)) % (q : ℤ) with hr
  have hr0 : 0 ≤ r := Int.emod_nonneg _ (by positivity)
  have hrq : r < q := Int.emod_lt_of_pos _ hq0
  have hq14' : (q : ℤ) ≤ 14 := by exact_mod_cast hq14
  have hq2' : (2 : ℤ) ≤ q := by exact_mod_cast hq2
  have hdvd : (q : ℤ) ∣ v i * (p : ℤ) ↔ r = 0 := by rw [hr, Int.dvd_iff_emod_eq_zero]
  rw [hdvd]
  constructor
  · intro h
    by_contra hne
    exact h ⟨by omega, by omega⟩
  · intro h ⟨h1, _⟩
    omega

/-- For `2 ≤ q ≤ 14`, `bandCount` counts the runners with `q ∣ vᵢ·p`. -/
lemma bandCount_eq_dvd (v : Fin 13 → ℤ) (q p : ℕ) (hq2 : 2 ≤ q) (hq14 : q ≤ 14) :
    bandCount v q p = (univ.filter (fun i => (q : ℤ) ∣ v i * (p : ℤ))).card := by
  unfold bandCount
  apply congrArg Finset.card
  apply Finset.filter_congr
  intro i _
  exact not_inBand_iff_dvd v q p hq2 hq14 i

/-- The divisibility bound on `bandCount`: `#{i : q ∣ vᵢ·p} ≤ #{i : d ∣ vᵢ}` for `d = q / gcd(q,p) ≥ 2`. -/
lemma bandCount_le_of_div (v : Fin 13 → ℤ) (q p : ℕ) (hq2 : 2 ≤ q) (hq14 : q ≤ 14)
    (hp : p ∈ Finset.Ioo 0 q)
    (hclean : ∀ d, d ∣ q → 2 ≤ d → (univ.filter (fun i => (d : ℤ) ∣ v i)).card ≤ 5) :
    bandCount v q p ≤ 5 := by
  rw [bandCount_eq_dvd v q p hq2 hq14]
  rw [Finset.mem_Ioo] at hp
  obtain ⟨hp0, hpq⟩ := hp
  set g := Nat.gcd q p with hg
  have hgpos : 0 < g := Nat.gcd_pos_of_pos_left p (by omega)
  have hgq : g ∣ q := Nat.gcd_dvd_left q p
  have hgp : g ∣ p := Nat.gcd_dvd_right q p
  obtain ⟨d, hd⟩ := hgq          -- q = g * d
  obtain ⟨p', hp'⟩ := hgp         -- p = g * p'
  have hd_dvd_q : d ∣ q := ⟨g, by rw [hd]; ring⟩
  have hg_le_p : g ≤ p := Nat.le_of_dvd (by omega) ⟨p', hp'⟩
  have hd2 : 2 ≤ d := by
    rcases Nat.lt_or_ge d 2 with h | h
    · interval_cases d <;> omega
    · exact h
  have hcop : Nat.Coprime d p' := by
    have := Nat.coprime_div_gcd_div_gcd (m := q) (n := p) hgpos
    rwa [← hg, hd, hp', Nat.mul_div_cancel_left _ hgpos, Nat.mul_div_cancel_left _ hgpos] at this
  have hsub : (univ.filter (fun i => (q : ℤ) ∣ v i * (p : ℤ)))
      ⊆ (univ.filter (fun i => (d : ℤ) ∣ v i)) := by
    intro i hi
    rw [Finset.mem_filter] at hi ⊢
    refine ⟨hi.1, ?_⟩
    have hqp : (q : ℤ) ∣ v i * (p : ℤ) := hi.2
    have : (g : ℤ) * (d : ℤ) ∣ (g : ℤ) * (v i * (p' : ℤ)) := by
      have e1 : ((q : ℕ) : ℤ) = (g : ℤ) * (d : ℤ) := by exact_mod_cast hd
      have e2 : ((p : ℕ) : ℤ) = (g : ℤ) * (p' : ℤ) := by exact_mod_cast hp'
      rw [e1] at hqp; rw [e2] at hqp
      calc (g : ℤ) * (d : ℤ) ∣ v i * ((g : ℤ) * (p' : ℤ)) := hqp
        _ = (g : ℤ) * (v i * (p' : ℤ)) := by ring
    have hg0 : (g : ℤ) ≠ 0 := by exact_mod_cast (by omega : g ≠ 0)
    have hdd : (d : ℤ) ∣ v i * (p' : ℤ) := (mul_dvd_mul_iff_left hg0).mp this
    have hcopZ : IsCoprime (d : ℤ) (p' : ℤ) := Nat.isCoprime_iff_coprime.mpr hcop
    exact hcopZ.dvd_of_dvd_mul_right hdd
  calc (univ.filter (fun i => (q : ℤ) ∣ v i * (p : ℤ))).card
      ≤ (univ.filter (fun i => (d : ℤ) ∣ v i)).card := Finset.card_le_card hsub
    _ ≤ 5 := hclean d hd_dvd_q hd2

/-- **The composite clean divisibility ruler.**  For `2 ≤ q ≤ 14`, if no speed is divisible by `q` and every
divisor `d ≥ 2` of `q` divides at most `5` speeds, then `0 < B5 v q`. -/
theorem b5_pos_of_div_clean (v : Fin 13 → ℤ) (q : ℕ) (hq2 : 2 ≤ q) (hq14 : q ≤ 14)
    (hlive : ∀ i, ¬ (q : ℤ) ∣ v i)
    (hclean : ∀ d, d ∣ q → 2 ≤ d → (univ.filter (fun i => (d : ℤ) ∣ v i)).card ≤ 5) :
    0 < B5 v q := by
  apply b5_pos_of_clean v q
  · intro p hp
    exact bandCount_le_of_div v q p hq2 hq14 hp hclean
  · rw [liveCount, Finset.card_pos]
    refine ⟨1, ?_⟩
    rw [Finset.mem_filter, Finset.mem_Ioo]
    refine ⟨⟨Nat.one_pos, by omega⟩, ?_⟩
    rw [bandCount_eq_dvd v q 1 hq2 hq14, Finset.card_eq_zero, Finset.filter_eq_empty_iff]
    intro i _
    simpa using hlive i

end LRC14Concrete
end LonelyRunner

#print axioms LonelyRunner.LRC14Concrete.b5_pos_of_div_clean
