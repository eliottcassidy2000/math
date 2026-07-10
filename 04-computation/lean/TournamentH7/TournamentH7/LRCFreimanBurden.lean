/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: opus (LRC multi-agent project, 2026-07-09-S188/S189)
-/
import Mathlib

/-!
# The Freiman descent-burden bound and equality core (THM-676 (ii), the burden ≤ 12 finite check)

The LRC(14) grand-assembly residual (THM-671) reduces via mac-mini's THM-676 to a Freiman 3k−4
stability check on majority-parity 7-classes with small **descent burden** `|A +̂ A| = #{aᵢ+aⱼ : i<j}`
(opus-S187). This file formalizes:

* `restrictedSum_card_ge` — the Freiman lower bound `|A +̂ A| ≥ 2|A| − 3` (the classical increasing-chain
  argument: with `m = min A`, `M = max A`, the `|A|−1` sums `{m + y : y ≠ m}` lie `≤ m+M`, the `|A|−2`
  sums `{x + M : m < x < M}` lie `> m+M`, so they are DISJOINT subsets of the restricted sumset);
* `burden_ge_eleven` — its `|A| = 7` specialization (burden ≥ 11 = THM-676 (ii)'s floor);
* `restrictedSum_eq_freimanChain` — the **Freiman-equality core**: at the MINIMAL burden
  `|A +̂ A| = 2|A| − 3`, the restricted sumset is EXACTLY the min/max chain `freimanChain`.

The equality core is the structural entry point for the full `burden = 2n−3 ⟹ A is an AP`
characterization, **now PROVED in `LRCFreimanAP.ap_of_min_burden` (opus-S195)** — but note it holds only
for **`n ≥ 5`** (MISTAKE-133: for `n ≤ 4`, e.g. `{0,1,3,4}`, the minimal restricted sumset is achieved
by non-AP "bi-arithmetic" sets; the S189 blueprint overstated it as general). For 7-classes `n = 7 ≥ 5`
so the step applies. The near-AP stability half of the finite check (burden ∈ {11,12} ⟹ opus-S187's 5
explicit shapes; burden = 13 ⟹ 2-D GAPs, routed through the density floor) sits above this core.

Kernel-pure: no `sorry`, no `native_decide`. Axioms: `[propext, Classical.choice, Quot.sound]`.
-/

namespace LRCFreiman

open Finset

/-- The **restricted sumset** (descent-burden support): sums `x + y` over strictly increasing pairs. -/
def restrictedSum (s : Finset ℤ) : Finset ℤ :=
  ((s ×ˢ s).filter (fun p => p.1 < p.2)).image (fun p => p.1 + p.2)

theorem mem_restrictedSum {s : Finset ℤ} {z : ℤ} :
    z ∈ restrictedSum s ↔ ∃ x ∈ s, ∃ y ∈ s, x < y ∧ x + y = z := by
  simp only [restrictedSum, mem_image, mem_filter, mem_product, Prod.exists]
  constructor
  · rintro ⟨x, y, ⟨⟨hx, hy⟩, hlt⟩, heq⟩; exact ⟨x, hx, y, hy, hlt, heq⟩
  · rintro ⟨x, hx, y, hy, hlt, heq⟩; exact ⟨x, y, ⟨⟨hx, hy⟩, hlt⟩, heq⟩

/-- The **Freiman min/max chain**: the `2|s|−3` sums forced into the restricted sumset — the low chain
`{m + y : y ∈ s, y ≠ m}` (`m = min`) together with the high chain `{x + M : x ∈ s, x ∉ {m, M}}`
(`M = max`). -/
def freimanChain (s : Finset ℤ) (hne : s.Nonempty) : Finset ℤ :=
  (s.erase (s.min' hne)).image (fun y => s.min' hne + y)
    ∪ ((s.erase (s.min' hne)).erase (s.max' hne)).image (fun x => x + s.max' hne)

/-- Two distinct elements of a `≥ 2`-element finset give `min' < max'`. -/
theorem min'_lt_max' (s : Finset ℤ) (hne : s.Nonempty) (hk : 2 ≤ s.card) :
    s.min' hne < s.max' hne := by
  have h1 : 1 < s.card := by omega
  obtain ⟨a, b, ha, hb, hab⟩ := Finset.one_lt_card_iff.mp h1
  rcases lt_or_gt_of_ne hab with h | h
  · exact lt_of_le_of_lt (s.min'_le a ha) (lt_of_lt_of_le h (s.le_max' b hb))
  · exact lt_of_le_of_lt (s.min'_le b hb) (lt_of_lt_of_le h (s.le_max' a ha))

/-- The Freiman chain is a subset of the restricted sumset. -/
theorem freimanChain_subset (s : Finset ℤ) (hne : s.Nonempty) :
    freimanChain s hne ⊆ restrictedSum s := by
  have hmem_m : s.min' hne ∈ s := s.min'_mem hne
  have hmem_M : s.max' hne ∈ s := s.max'_mem hne
  apply union_subset
  · intro z hz
    rw [mem_image] at hz
    obtain ⟨y, hy, rfl⟩ := hz
    rw [mem_erase] at hy
    exact mem_restrictedSum.mpr
      ⟨s.min' hne, hmem_m, y, hy.2, lt_of_le_of_ne (s.min'_le y hy.2) (Ne.symm hy.1), rfl⟩
  · intro z hz
    rw [mem_image] at hz
    obtain ⟨x, hx, rfl⟩ := hz
    rw [mem_erase, mem_erase] at hx
    exact mem_restrictedSum.mpr
      ⟨x, hx.2.2, s.max' hne, hmem_M, lt_of_le_of_ne (s.le_max' x hx.2.2) hx.1, rfl⟩

/-- The Freiman chain has exactly `2|s| − 3` elements. -/
theorem freimanChain_card (s : Finset ℤ) (hne : s.Nonempty) (hk : 2 ≤ s.card) :
    (freimanChain s hne).card = 2 * s.card - 3 := by
  have hmem_m : s.min' hne ∈ s := s.min'_mem hne
  have hmem_M : s.max' hne ∈ s := s.max'_mem hne
  have hmM : s.min' hne < s.max' hne := min'_lt_max' s hne hk
  have hMerase : s.max' hne ∈ s.erase (s.min' hne) := mem_erase.mpr ⟨ne_of_gt hmM, hmem_M⟩
  have hdisj : Disjoint ((s.erase (s.min' hne)).image (fun y => s.min' hne + y))
      (((s.erase (s.min' hne)).erase (s.max' hne)).image (fun x => x + s.max' hne)) := by
    rw [Finset.disjoint_left]
    intro z hz1 hz2
    rw [mem_image] at hz1
    obtain ⟨y, hy, rfl⟩ := hz1
    rw [mem_image] at hz2
    obtain ⟨x, hx, hxeq⟩ := hz2
    rw [mem_erase] at hy
    rw [mem_erase, mem_erase] at hx
    have hyM : y ≤ s.max' hne := s.le_max' y hy.2
    have hmx : s.min' hne < x := lt_of_le_of_ne (s.min'_le x hx.2.2) (Ne.symm hx.2.1)
    omega
  rw [freimanChain, card_union_of_disjoint hdisj,
    card_image_of_injective _ (add_right_injective (s.min' hne)),
    card_image_of_injective _ (add_left_injective (s.max' hne)),
    card_erase_of_mem hMerase, card_erase_of_mem hmem_m]
  omega

/-- **Freiman restricted-sumset lower bound.** For a finite set of ≥ 2 integers,
`|A +̂ A| ≥ 2|A| − 3` (THM-676 (ii)). -/
theorem restrictedSum_card_ge (s : Finset ℤ) (hk : 2 ≤ s.card) :
    2 * s.card - 3 ≤ (restrictedSum s).card := by
  have hne : s.Nonempty := card_pos.mp (by omega)
  calc 2 * s.card - 3 = (freimanChain s hne).card := (freimanChain_card s hne hk).symm
    _ ≤ (restrictedSum s).card := card_le_card (freimanChain_subset s hne)

/-- **Descent burden of a 7-set is ≥ 11** (THM-676 (ii), the floor of the finite check). -/
theorem burden_ge_eleven (s : Finset ℤ) (h7 : s.card = 7) :
    11 ≤ (restrictedSum s).card := by
  have := restrictedSum_card_ge s (by omega); omega

/-- **The Freiman-equality core.** At the MINIMAL burden `|A +̂ A| = 2|A| − 3`, the restricted sumset is
EXACTLY the min/max chain. This is the structural entry point for `burden = 11 ⟹ A is an AP`: once the
sumset is pinned to the chain, the row `{a₁ + aⱼ}` maps bijectively (order-preserving) onto the chain
positions `2 … k−1`, forcing `a₁ − a₀ = aⱼ₊₁ − aⱼ` (constant differences); the missing difference comes
from the reflection `aᵢ ↦ −a_{k-1-i}` (blueprint: opus-S189 reflection). -/
theorem restrictedSum_eq_freimanChain (s : Finset ℤ) (hk : 2 ≤ s.card)
    (hcard : (restrictedSum s).card = 2 * s.card - 3) :
    restrictedSum s = freimanChain s (card_pos.mp (by omega)) := by
  -- The chain ⟹ AP conclusion is proved (for `card ≥ 5`) in `LRCFreimanAP.ap_of_min_burden`,
  -- via a cleaner interleaved-chain route than the S189 row-bijection+reflection blueprint.
  have hne : s.Nonempty := card_pos.mp (by omega)
  refine (eq_of_subset_of_card_le (freimanChain_subset s hne) ?_).symm
  rw [freimanChain_card s hne hk, hcard]

end LRCFreiman
