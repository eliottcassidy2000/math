/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: opus (LRC multi-agent project, 2026-07-10-S212)
-/
import Mathlib
import TournamentH7.LonelyRunner
import TournamentH7.LRC13Citation
import TournamentH7.LRCDetunedD3

/-!
# THM-678 `d = 3` — the clean coarse corollary (opus-S212)

kps-S127's `LRCDetunedD3` proves the `d = 3` dispatch under the raw counting hypothesis
`Σⱼ badCount δⱼ g < g.toNat`. This file adds the clean **all-fine corollary**: when every detuning
denominator `qⱼ = g/gcd(δⱼ,g) ≥ 8`, the counting hypothesis is automatic (`Nⱼ/qⱼ ≤ 15/56`, so
`Σ ≤ 3·15/56 = 45/56 < 1`) — dispatching without any per-instance sum check. This is the common case
(dissociated detuning, `δⱼ` fairly coprime to `g`); the surviving exceptions are all small-`q`
(`kps`'s `lrc14_three_detuned_exceptional`), the `d = 3` analogue of `d = 2`'s `(2,2)` residual.

Kernel-pure: no `sorry`, no `native_decide`. Axioms: `[propext, Classical.choice, Quot.sound]`.
-/

namespace LonelyRunner
namespace DetunedD3

open scoped Classical

/-- `g.toNat` factors as `gcd(δ,g)·(g/gcd(δ,g))` (both nonneg). -/
private theorem gcd_toNat_mul (δ g : ℤ) (hg : 0 < g) :
    g.toNat = (Int.gcd δ g) * (g / (Int.gcd δ g : ℤ)).toNat := by
  have hdvd : ((Int.gcd δ g : ℤ)) ∣ g := Int.gcd_dvd_right δ g
  have hdz : (0 : ℤ) < (Int.gcd δ g : ℤ) := by
    have : 0 < Int.gcd δ g := by rw [Int.gcd_pos_iff]; right; omega
    exact_mod_cast this
  have hqnn : (0 : ℤ) ≤ g / (Int.gcd δ g : ℤ) := Int.ediv_nonneg (le_of_lt hg) (le_of_lt hdz)
  have heq : g = (Int.gcd δ g : ℤ) * (g / (Int.gcd δ g : ℤ)) := (Int.mul_ediv_cancel' hdvd).symm
  have hcast : (g.toNat : ℤ) = ((Int.gcd δ g) * (g / (Int.gcd δ g : ℤ)).toNat : ℕ) := by
    rw [Int.toNat_of_nonneg (le_of_lt hg)]; push_cast; rw [Int.toNat_of_nonneg hqnn]; exact heq
  exact_mod_cast hcast

/-- The per-coordinate "coarse" bound: `56·badCount δ g ≤ 15·g` when `q = g/gcd(δ,g) ≥ 8`
(i.e. `N/q ≤ 15/56`). -/
private theorem coarse_badCount_bound (δ g : ℤ) (hg : 0 < g)
    (hq : 8 ≤ g / (Int.gcd δ g : ℤ)) : 56 * badCount δ g ≤ 15 * g.toNat := by
  have hQ8 : 8 ≤ (g / (Int.gcd δ g : ℤ)).toNat := by omega
  have hG : g.toNat = (Int.gcd δ g) * (g / (Int.gcd δ g : ℤ)).toNat := gcd_toNat_mul δ g hg
  have h56 : 56 * ((g / (Int.gcd δ g : ℤ)).toNat / 7 + 1)
      ≤ 15 * (g / (Int.gcd δ g : ℤ)).toNat := by omega
  calc 56 * badCount δ g
      = (Int.gcd δ g) * (56 * ((g / (Int.gcd δ g : ℤ)).toNat / 7 + 1)) := by rw [badCount]; ring
    _ ≤ (Int.gcd δ g) * (15 * (g / (Int.gcd δ g : ℤ)).toNat) := Nat.mul_le_mul_left _ h56
    _ = 15 * ((Int.gcd δ g) * (g / (Int.gcd δ g : ℤ)).toNat) := by ring
    _ = 15 * g.toNat := by rw [hG]

/-- **THM-678 `d = 3`, the clean coarse corollary.** If all three detuning denominators
`qᵢ = g/gcd(δᵢ,g) ≥ 8` (the all-fine case, `Σ Nᵢ/qᵢ ≤ 3·15/56 = 45/56 < 1`), the family
`v = g·H ∪ {δ₁,δ₂,δ₃}` is lonely — from the LRC(≤13) citation, no per-instance sum check. -/
theorem lonely14_of_three_detuned_coarse (cite : LRCUpTo13)
    (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0) (g : ℤ) (hg : 2 ≤ g)
    (i₁ i₂ i₃ : Fin 13) (h12 : i₁ ≠ i₂) (h13 : i₁ ≠ i₃) (h23 : i₂ ≠ i₃)
    (hdvd : ∀ j, j ≠ i₁ → j ≠ i₂ → j ≠ i₃ → g ∣ v j)
    (hδ1 : ¬ g ∣ v i₁) (hδ2 : ¬ g ∣ v i₂) (hδ3 : ¬ g ∣ v i₃)
    (hq1 : 8 ≤ g / (Int.gcd (v i₁) g : ℤ)) (hq2 : 8 ≤ g / (Int.gcd (v i₂) g : ℤ))
    (hq3 : 8 ≤ g / (Int.gcd (v i₃) g : ℤ)) :
    ∃ t : ℝ, Lonely 14 v t := by
  have hg0 : (0 : ℤ) < g := by omega
  refine lonely14_of_three_detuned' cite v hv g hg i₁ i₂ i₃ h12 h13 h23 hdvd hδ1 hδ2 hδ3 ?_
  -- `56·Σ badCount ≤ 45·g.toNat < 56·g.toNat`
  have hG2 : 2 ≤ g.toNat := by omega
  have hb1 := coarse_badCount_bound (v i₁) g hg0 hq1
  have hb2 := coarse_badCount_bound (v i₂) g hg0 hq2
  have hb3 := coarse_badCount_bound (v i₃) g hg0 hq3
  omega

end DetunedD3
end LonelyRunner

