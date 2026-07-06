/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: opus (LRC multi-agent project, 2026-07-06-S105)
-/
import TournamentH7.LRCDivisorProtection

/-!
# The universal pinned floor: M ≥ 1/13 for any 13-coprime family (HYP-4376)

The tight-side anchor, as a clean corollary of the divisor-protection lemma (S104) at the
prime modulus `k = 13`.  A 12-runner family none of whose speeds is a multiple of `13` — in
particular EVERY residue-pinned lift `{r + 13·k_r}` (residues `{1,…,12} (mod 13)`, the shape
`residue_pinning_13` forces on tight-from-above families) — is lonely at `t = 1/13` with
margin `1/13`.  So `M ≥ 1/13` on the whole pinned locus: `1/13` is the floor the entire
gap-window analysis sits above, and it is reached exactly by the AP `{1,…,12}` (mac-mini-S12:
uniquely, BECAUSE 13 is prime — every nonzero residue is a unit).
-/

namespace LonelyRunner
namespace PinnedFloor

open DivisorProtection

/-- **The pinned floor.**  A family `v : Fin 12 → ℤ` with no speed divisible by `13` is
lonely at `t = 1/13` with margin `1/13`: every runner stays `≥ 1/13` from every integer.
Consumes `lonely_of_all_not_dvd` at `k = 13`. -/
theorem floor_of_no_dvd_13 (v : Fin 12 → ℤ) (hv : ∀ i, ¬ (13 : ℤ) ∣ v i) :
    ∀ i, ∀ m : ℤ, (1 : ℝ) / 13 ≤ |(v i : ℝ) * (1 / 13) - m| := by
  intro i m
  have hval : (v i : ℝ) * (1 / 13) - m = (v i : ℝ) / 13 - m := by ring
  rw [hval]
  have := int_far_of_not_dvd_k 13 (by norm_num) (v i) (hv i) m
  simpa using this

/-- The residue-pinned lifts are exactly the `13`-coprime families: `{r + 13·k_r}` with
`r ∈ {1,…,12}` is never a multiple of `13`, so the pinned floor applies. -/
theorem floor_of_pinned_lift (r : Fin 12 → ℤ) (k : Fin 12 → ℤ)
    (hr : ∀ i, ¬ (13 : ℤ) ∣ r i) :
    ∀ i, ∀ m : ℤ, (1 : ℝ) / 13 ≤ |((r i + 13 * k i : ℤ) : ℝ) * (1 / 13) - m| := by
  apply floor_of_no_dvd_13
  intro i
  -- 13 ∣ (r i + 13·k i) ↔ 13 ∣ r i, which is false
  intro hdvd
  apply hr i
  have : (13 : ℤ) ∣ 13 * k i := Dvd.intro (k i) rfl
  have h2 : (13 : ℤ) ∣ (r i + 13 * k i) - 13 * k i := dvd_sub hdvd this
  simpa using h2

#print axioms floor_of_no_dvd_13
#print axioms floor_of_pinned_lift

end PinnedFloor
end LonelyRunner
