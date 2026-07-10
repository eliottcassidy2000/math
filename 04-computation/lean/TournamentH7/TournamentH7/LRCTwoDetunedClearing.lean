/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: opus (LRC multi-agent project, 2026-07-10-S211)
-/
import Mathlib
import TournamentH7.LonelyRunner
import TournamentH7.LRC13Citation
import TournamentH7.LRCDetunedD2
import TournamentH7.LRCIntervalCount

/-!
# THM-678 `d = 2` — the counting `TwoDetunedClearing`, PROVED (opus-S211)

Discharges `DetunedD2.TwoDetunedClearing`: for distinct-support detuned `δ₁, δ₂` (neither divisible by
`g ≥ 1`, `(q₁,q₂) ≠ (2,2)`), every real `u` admits a branch `c ∈ [0,g)` clearing both detuned phases.

Union bound over `c ∈ [0, g)`:
* `LRCIntervalCount.bad_count_le` — the `c` failing coordinate `j` number `≤ dⱼ·(⌊qⱼ/7⌋+1)`
  (the de-circled `ψ`-injection count);
* `LRCIntervalCount.sum_lt` — those two bounds sum to `< g` exactly when `(q₁,q₂) ≠ (2,2)`;
* so the two bad sets cannot cover `[0,g)`, and any uncovered `c` clears both (a non-`< 1/14` phase is
  `≥ 1/14` from every integer).

Composing with `DetunedD2.lonely14_of_two_detuned` (opus-S210) this makes the generic `d = 2` detuned
dispatch unconditional from the LRC(≤13) citation — kernel-pure, `[propext, Classical.choice, Quot.sound]`.
-/

namespace LonelyRunner
namespace DetunedD2

open scoped Classical

/-- **THM-678 `d = 2`, the counting — PROVED.** The two-detuned clearing obligation holds. -/
theorem twoDetunedClearing : TwoDetunedClearing := by
  intro δ₁ δ₂ g hg hδ₁ hδ₂ hne u
  have hg1 : (1 : ℤ) ≤ g := hg
  -- the two bad-branch sets and their counts
  set bad₁ := (Finset.Ico (0 : ℤ) g).filter
      (fun (c : ℤ) => ∃ n : ℤ, |(δ₁ : ℝ) * ((u + (c : ℝ)) / (g : ℝ)) - n| < 1 / 14) with hbad₁
  set bad₂ := (Finset.Ico (0 : ℤ) g).filter
      (fun (c : ℤ) => ∃ n : ℤ, |(δ₂ : ℝ) * ((u + (c : ℝ)) / (g : ℝ)) - n| < 1 / 14) with hbad₂
  have hb₁ : bad₁.card ≤ (Int.gcd δ₁ g) * ((g / (Int.gcd δ₁ g : ℤ)).toNat / 7 + 1) :=
    LRCIntervalCount.bad_count_le δ₁ g hg1 u
  have hb₂ : bad₂.card ≤ (Int.gcd δ₂ g) * ((g / (Int.gcd δ₂ g : ℤ)).toNat / 7 + 1) :=
    LRCIntervalCount.bad_count_le δ₂ g hg1 u
  -- naturals `dⱼ, Qⱼ` with `g.toNat = dⱼ·Qⱼ` and `Qⱼ ≥ 2`
  set d₁ := Int.gcd δ₁ g with hd₁
  set d₂ := Int.gcd δ₂ g with hd₂
  set Q₁ := (g / (Int.gcd δ₁ g : ℤ)).toNat with hQ₁def
  set Q₂ := (g / (Int.gcd δ₂ g : ℤ)).toNat with hQ₂def
  have hd₁pos : 0 < d₁ := by rw [hd₁, Int.gcd_pos_iff]; right; omega
  have hd₂pos : 0 < d₂ := by rw [hd₂, Int.gcd_pos_iff]; right; omega
  -- `g.toNat = dⱼ·Qⱼ`
  have key : ∀ (δ : ℤ), ¬ (g ∣ δ) →
      g.toNat = (Int.gcd δ g) * (g / (Int.gcd δ g : ℤ)).toNat ∧
      2 ≤ (g / (Int.gcd δ g : ℤ)).toNat ∧
      ((g / (Int.gcd δ g : ℤ)).toNat = 2 ↔ g / (Int.gcd δ g : ℤ) = 2) := by
    intro δ hδ
    have hdvd : ((Int.gcd δ g : ℤ)) ∣ g := Int.gcd_dvd_right δ g
    have hdz : (0 : ℤ) < (Int.gcd δ g : ℤ) := by
      have : 0 < Int.gcd δ g := by rw [Int.gcd_pos_iff]; right; omega
      exact_mod_cast this
    have hqnn : (0 : ℤ) ≤ g / (Int.gcd δ g : ℤ) := Int.ediv_nonneg (le_of_lt hg) (le_of_lt hdz)
    have heq : g = (Int.gcd δ g : ℤ) * (g / (Int.gcd δ g : ℤ)) :=
      (Int.mul_ediv_cancel' hdvd).symm
    refine ⟨?_, ?_, ?_⟩
    · -- `g.toNat = gcd · Q`
      have hcast : (g.toNat : ℤ) = ((Int.gcd δ g) * (g / (Int.gcd δ g : ℤ)).toNat : ℕ) := by
        rw [Int.toNat_of_nonneg (le_of_lt hg)]
        push_cast
        rw [Int.toNat_of_nonneg hqnn]
        exact heq
      exact_mod_cast hcast
    · -- `2 ≤ Q`
      rcases Nat.lt_or_ge (g / (Int.gcd δ g : ℤ)).toNat 2 with h | h
      · exfalso
        interval_cases hQ : (g / (Int.gcd δ g : ℤ)).toNat
        · -- Q = 0 ⟹ g ≤ 0
          have : g / (Int.gcd δ g : ℤ) = 0 := by
            have := Int.toNat_of_nonneg hqnn; omega
          rw [this, mul_zero] at heq; omega
        · -- Q = 1 ⟹ g = gcd ⟹ g ∣ δ
          have hq1 : g / (Int.gcd δ g : ℤ) = 1 := by
            have := Int.toNat_of_nonneg hqnn; omega
          rw [hq1, mul_one] at heq
          exact hδ (heq ▸ Int.gcd_dvd_left δ g)
      · exact h
    · -- `Q = 2 ↔ g/gcd = 2`
      constructor
      · intro h; omega
      · intro h; rw [h]; rfl
  obtain ⟨hgQ₁, hQ₁2, hiff₁⟩ := key δ₁ hδ₁
  obtain ⟨hgQ₂, hQ₂2, hiff₂⟩ := key δ₂ hδ₂
  -- `(q₁,q₂) ≠ (2,2)` transported to naturals
  have hne_nat : ¬ (Q₁ = 2 ∧ Q₂ = 2) := by
    intro ⟨e1, e2⟩; exact hne ⟨hiff₁.mp e1, hiff₂.mp e2⟩
  -- the sum bound
  have hsum : d₁ * (Q₁ / 7 + 1) + d₂ * (Q₂ / 7 + 1) < g.toNat :=
    LRCIntervalCount.sum_lt d₁ Q₁ d₂ Q₂ g.toNat hd₁pos hd₂pos hQ₁2 hQ₂2 hne_nat hgQ₁ hgQ₂
  -- the two bad sets cannot cover `[0, g)`
  have hcard_lt : (bad₁ ∪ bad₂).card < (Finset.Ico (0 : ℤ) g).card := by
    have hIco : (Finset.Ico (0 : ℤ) g).card = g.toNat := by
      rw [Int.card_Ico]; congr 1; omega
    calc (bad₁ ∪ bad₂).card ≤ bad₁.card + bad₂.card := Finset.card_union_le _ _
      _ ≤ d₁ * (Q₁ / 7 + 1) + d₂ * (Q₂ / 7 + 1) := Nat.add_le_add hb₁ hb₂
      _ < g.toNat := hsum
      _ = (Finset.Ico (0 : ℤ) g).card := hIco.symm
  -- an uncovered branch
  have hnotsub : ¬ (Finset.Ico (0 : ℤ) g ⊆ bad₁ ∪ bad₂) := by
    intro h; exact absurd (Finset.card_le_card h) (by omega)
  rw [Finset.not_subset] at hnotsub
  obtain ⟨c, hcIco, hcnot⟩ := hnotsub
  rw [Finset.mem_union, not_or] at hcnot
  obtain ⟨hc1, hc2⟩ := hcnot
  -- at `c` both detuned phases clear `1/14`
  refine ⟨c, ?_, ?_⟩
  · intro n
    refine not_lt.mp (fun h => hc1 ?_)
    rw [hbad₁, Finset.mem_filter]
    exact ⟨hcIco, n, h⟩
  · intro n
    refine not_lt.mp (fun h => hc2 ?_)
    rw [hbad₂, Finset.mem_filter]
    exact ⟨hcIco, n, h⟩

/-- **THM-678 `d = 2` generic, unconditional.** A family `v = g·H ∪ {δ₁, δ₂}` (`g ≥ 2` dividing all but
`i₁, i₂`, `(q₁,q₂) ≠ (2,2)`) is lonely, from the LRC(≤13) citation alone — the counting is now discharged. -/
theorem lonely14_of_two_detuned' (cite : LRCUpTo13)
    (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0) (g : ℤ) (hg : 2 ≤ g)
    (i₁ i₂ : Fin 13) (hne : i₁ ≠ i₂)
    (hdvd : ∀ j, j ≠ i₁ → j ≠ i₂ → g ∣ v j) (hδ1 : ¬ g ∣ v i₁) (hδ2 : ¬ g ∣ v i₂)
    (hq : ¬ (g / (Int.gcd (v i₁) g : ℤ) = 2 ∧ g / (Int.gcd (v i₂) g : ℤ) = 2)) :
    ∃ t : ℝ, Lonely 14 v t :=
  lonely14_of_two_detuned cite twoDetunedClearing v hv g hg i₁ i₂ hne hdvd hδ1 hδ2 hq

end DetunedD2
end LonelyRunner

