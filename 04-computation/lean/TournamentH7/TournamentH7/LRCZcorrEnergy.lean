/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: kind-pasteur (LRC multi-agent project, 2026-07-11-S127)
-/
import Mathlib
import TournamentH7.LRCMcorrHyperbola

/-!
# The `ZMod q` off-diagonal energy aggregation — wiring `zcorr_percell` to the `t₂` energy (kps-S127)

`LRCMcorrHyperbola.zcorr_percell` supplies the per-cell bound `(zcorr A w − 1)·P ≤ 4K²` for the ratio
correlation `zcorr A w = #{(a,b)∈A² : a = w·b}` on `ZMod q`, drawn from death-star's combinatorial
`hyperbola_box_count`.  The consumer `MultCorrelation.offdiag_mcorr_sq_le` aggregates a per-cell bound into
the total off-diagonal energy — but it lives over an abstract `Group G`, and `zcorr` lives over the monoid
`ZMod q`.  This file supplies the **`ZMod q`-native aggregation**, closing that seam:

  `offdiag_zcorr_sq_le` : `(∀w≠1, zcorr A w ≤ M) → Σ_{w≠1} (zcorr A w)² ≤ M·(|A|² − |A|)`

for a set `A` of **units** (the identity `Σ_w zcorr A w = |A|²` needs each `b∈A` invertible, so the fiber
`{w : a = w·b}` is a singleton — exactly the character program's prime-modulus setting).  Combined with
`zcorr_percell`, the `hyperbola_box_count → zcorr_percell → energy` route is complete: the `t₂` off-diagonal
energy is bounded by the per-cell hyperbola count, kernel-pure.

`zcorr A w = mcorr A w` on units (`a·b⁻¹ = w ⟺ a = w·b`), so this is the `ZMod q` realization of
`offdiag_mcorr_sq_le`; the bridge to the abstract-group form is the units embedding
`(ZMod q)ˣ ↪ ZMod q`.

Kernel-pure: no `sorry`, no `native_decide`.  Axioms: `[propext, Classical.choice, Quot.sound]`.
-/

namespace LonelyRunner
namespace McorrHyperbola

open Finset

variable {q : ℕ} [NeZero q]

/-- **The diagonal correlation is `|A|`** — `w = 1` counts the pairs `(a,a)` (`a = 1·a`). -/
theorem zcorr_one (A : Finset (ZMod q)) : zcorr A 1 = A.card := by
  rw [zcorr]
  have heq : (A ×ˢ A).filter (fun p => p.1 = 1 * p.2) = A.image (fun a => (a, a)) := by
    ext ⟨x, y⟩
    simp only [Finset.mem_filter, Finset.mem_product, Finset.mem_image, one_mul, Prod.mk.injEq]
    constructor
    · rintro ⟨⟨hx, hy⟩, hxy⟩
      exact ⟨x, hx, rfl, hxy⟩
    · rintro ⟨a, ha, rfl, rfl⟩
      exact ⟨⟨ha, ha⟩, rfl⟩
  rw [heq, Finset.card_image_of_injective A (fun _ _ h => (Prod.ext_iff.mp h).1)]

/-- **The total correlation is `|A|²`** — on a set of units, every ordered pair `(a,b)` has exactly one
ratio `w = a·b⁻¹` with `a = w·b`. -/
theorem sum_zcorr (A : Finset (ZMod q)) (hunit : ∀ a ∈ A, IsUnit a) :
    ∑ w : ZMod q, zcorr A w = A.card ^ 2 := by
  have hfib : (A ×ˢ A).card
      = ∑ w : ZMod q, ((A ×ˢ A).filter (fun p => p.1 * p.2⁻¹ = w)).card :=
    Finset.card_eq_sum_card_fiberwise (fun _ _ => Finset.mem_univ _)
  rw [Finset.card_product, ← sq] at hfib
  rw [hfib]
  refine Finset.sum_congr rfl (fun w _ => ?_)
  rw [zcorr]
  congr 1
  refine Finset.filter_congr (fun p hp => ?_)
  obtain ⟨a, b⟩ := p
  simp only [Finset.mem_product] at hp
  have hb : IsUnit b := hunit b hp.2
  constructor
  · intro h; rw [h, mul_assoc, ZMod.mul_inv_of_unit b hb, mul_one]
  · intro h; rw [← h, mul_assoc, ZMod.inv_mul_of_unit b hb, mul_one]

/-- **The `ZMod q` off-diagonal energy bound.**  For a set `A` of units and a uniform per-cell ratio bound
`zcorr A w ≤ M` (`w ≠ 1`), the total off-diagonal `t₂` energy is `Σ_{w≠1} (zcorr A w)² ≤ M·(|A|² − |A|)`.
The `w = 1` diagonal (mass `|A|²`) is isolated; `M` is supplied per-cell by `zcorr_percell`. -/
theorem offdiag_zcorr_sq_le (A : Finset (ZMod q)) (M : ℕ) (hunit : ∀ a ∈ A, IsUnit a)
    (hM : ∀ w : ZMod q, w ≠ 1 → zcorr A w ≤ M) :
    ∑ w ∈ Finset.univ.erase 1, (zcorr A w) ^ 2 ≤ M * (A.card ^ 2 - A.card) := by
  have hbound : ∑ w ∈ Finset.univ.erase 1, (zcorr A w) ^ 2
      ≤ ∑ w ∈ Finset.univ.erase 1, M * zcorr A w := by
    apply Finset.sum_le_sum
    intro w hw
    rw [Finset.mem_erase] at hw
    rw [sq]
    exact Nat.mul_le_mul_right _ (hM w hw.1)
  calc ∑ w ∈ Finset.univ.erase 1, (zcorr A w) ^ 2
      ≤ ∑ w ∈ Finset.univ.erase 1, M * zcorr A w := hbound
    _ = M * ∑ w ∈ Finset.univ.erase 1, zcorr A w := by rw [Finset.mul_sum]
    _ = M * (A.card ^ 2 - A.card) := by
        have hadd := Finset.add_sum_erase Finset.univ (zcorr A) (Finset.mem_univ 1)
        rw [sum_zcorr A hunit, zcorr_one A] at hadd
        congr 1
        omega

/-- **The wired energy bound: per-cell hyperbola count ⟹ `t₂` off-diagonal energy.**  For a `K`-bounded set
`A` of units and a uniform hyperbola floor `P` on every off-diagonal ratio, the total off-diagonal energy is
`Σ_{w≠1} (zcorr A w)² ≤ (1 + 4K²/P)·(|A|² − |A|)` — `hyperbola_box_count → zcorr_percell → offdiag_zcorr_sq_le`
composed.  `M = 1 + 4K²/P` is the per-cell ceiling from the box count. -/
theorem zcorr_energy_of_hyperbola (A : Finset (ZMod q)) (K P : ℕ) (hP : 0 < P)
    (hunit : ∀ a ∈ A, IsUnit a) (h0 : ∀ a ∈ A, a ≠ 0) (hK : ∀ a ∈ A, HyperbolaBox.cdist a ≤ K)
    (hPmin : ∀ w : ZMod q, w ≠ 1 → ∀ h : ZMod q, h ≠ 0 →
        P ≤ HyperbolaBox.cdist h * HyperbolaBox.cdist (w * h)) :
    ∑ w ∈ Finset.univ.erase 1, (zcorr A w) ^ 2
      ≤ (1 + 4 * K * K / P) * (A.card ^ 2 - A.card) := by
  apply offdiag_zcorr_sq_le A (1 + 4 * K * K / P) hunit
  intro w hw
  -- from `(zcorr A w − 1)·P ≤ 4K²` (zcorr_percell) get `zcorr A w ≤ 1 + 4K²/P`
  have hpc := zcorr_percell A w K P h0 hK (hPmin w hw)
  rcases Nat.eq_zero_or_pos (zcorr A w) with hz0 | hz1
  · rw [hz0]; exact Nat.zero_le _
  · have hcast : ((zcorr A w : ℤ) - 1) * P = (((zcorr A w - 1) * P : ℕ) : ℤ) := by
      have h1 : 1 ≤ zcorr A w := hz1
      push_cast [Nat.cast_sub h1]; ring
    rw [hcast] at hpc
    have hnat : (zcorr A w - 1) * P ≤ 4 * K * K := by exact_mod_cast hpc
    have hdiv : zcorr A w - 1 ≤ 4 * K * K / P := (Nat.le_div_iff_mul_le hP).mpr hnat
    calc zcorr A w ≤ 4 * K * K / P + 1 := Nat.sub_le_iff_le_add.mp hdiv
      _ = 1 + 4 * K * K / P := Nat.add_comm _ _

end McorrHyperbola
end LonelyRunner
