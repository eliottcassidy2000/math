/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: opus (LRC multi-agent project, 2026-07-09-S188)
-/
import Mathlib

/-!
# The Freiman descent-burden lower bound (THM-675 (ii), the burden ≤ 12 finite check)

The LRC(14) grand-assembly residual (THM-671) reduces via mac-mini's THM-675 to a Freiman 3k−4
stability check on majority-parity 7-classes with small **descent burden** `|A +̂ A| = #{aᵢ+aⱼ : i<j}`
(opus-S187). This file formalizes the FLOOR of that check — the Freiman restricted-sumset lower bound

  `|A +̂ A| ≥ 2·|A| − 3`

for any finite set of integers (`|A| ≥ 2`), and specializes it to `|A| = 7` (burden ≥ 11). The proof is
the classical increasing-chain argument: with `m = min A`, `M = max A`, the `(|A|−1)` sums `{m + y : y ∈
A, y ≠ m}` all lie `≤ m + M`, and the `(|A|−2)` sums `{x + M : x ∈ A, m < x < M}` all lie `> m + M`, so
they are disjoint subsets of the restricted sumset — giving `(|A|−1) + (|A|−2) = 2|A| − 3` distinct sums.

This is the lower half of the burden check (burden ≥ 11); the equality/near-equality characterization
(burden ∈ {11,12} ⟹ near-AP; opus-S187's 5-shape finite family) is the stability half.

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

/-- **Freiman restricted-sumset lower bound.** For a finite set of ≥ 2 integers,
the restricted sumset has at least `2·|s| − 3` elements. -/
theorem restrictedSum_card_ge (s : Finset ℤ) (hk : 2 ≤ s.card) :
    2 * s.card - 3 ≤ (restrictedSum s).card := by
  have hne : s.Nonempty := card_pos.mp (by omega)
  set m := s.min' hne with hm
  set M := s.max' hne with hMdef
  have hmem_m : m ∈ s := s.min'_mem hne
  have hmem_M : M ∈ s := s.max'_mem hne
  have hcard1 : 1 < s.card := by omega
  have hmM : m < M := by
    obtain ⟨a, b, ha, hb, hab⟩ := Finset.one_lt_card_iff.mp hcard1
    rcases lt_or_gt_of_ne hab with h | h
    · exact lt_of_le_of_lt (s.min'_le a ha) (lt_of_lt_of_le h (s.le_max' b hb))
    · exact lt_of_le_of_lt (s.min'_le b hb) (lt_of_lt_of_le h (s.le_max' a ha))
  have hMerase : M ∈ s.erase m := mem_erase.mpr ⟨ne_of_gt hmM, hmem_M⟩
  -- the two chains
  set C1 := (s.erase m).image (fun y => m + y) with hC1
  set C2 := ((s.erase m).erase M).image (fun x => x + M) with hC2
  -- their cardinalities
  have e1 : C1.card = s.card - 1 := by
    rw [hC1, card_image_of_injective _ (add_right_injective m), card_erase_of_mem hmem_m]
  have e2 : C2.card = s.card - 2 := by
    rw [hC2, card_image_of_injective _ (add_left_injective M), card_erase_of_mem hMerase,
      card_erase_of_mem hmem_m]
    omega
  -- both are subsets of the restricted sumset
  have hC1sub : C1 ⊆ restrictedSum s := by
    intro z hz
    rw [hC1, mem_image] at hz
    obtain ⟨y, hy, rfl⟩ := hz
    rw [mem_erase] at hy
    exact mem_restrictedSum.mpr
      ⟨m, hmem_m, y, hy.2, lt_of_le_of_ne (s.min'_le y hy.2) (Ne.symm hy.1), rfl⟩
  have hC2sub : C2 ⊆ restrictedSum s := by
    intro z hz
    rw [hC2, mem_image] at hz
    obtain ⟨x, hx, rfl⟩ := hz
    rw [mem_erase, mem_erase] at hx
    exact mem_restrictedSum.mpr
      ⟨x, hx.2.2, M, hmem_M, lt_of_le_of_ne (s.le_max' x hx.2.2) hx.1, rfl⟩
  -- the chains are disjoint (max of C1 = m+M < min of C2)
  have hdisj : Disjoint C1 C2 := by
    rw [Finset.disjoint_left]
    intro z hz1 hz2
    rw [hC1, mem_image] at hz1
    obtain ⟨y, hy, rfl⟩ := hz1
    rw [hC2, mem_image] at hz2
    obtain ⟨x, hx, hxeq⟩ := hz2
    rw [mem_erase] at hy
    rw [mem_erase, mem_erase] at hx
    have hyM : y ≤ M := s.le_max' y hy.2
    have hmx : m < x := lt_of_le_of_ne (s.min'_le x hx.2.2) (Ne.symm hx.2.1)
    omega
  -- combine
  have hunion : (C1 ∪ C2).card = C1.card + C2.card := card_union_of_disjoint hdisj
  have hsub : C1 ∪ C2 ⊆ restrictedSum s := union_subset hC1sub hC2sub
  have hle : (C1 ∪ C2).card ≤ (restrictedSum s).card := card_le_card hsub
  omega

/-- **Descent burden of a 7-set is ≥ 11** (THM-675 (ii), the floor of the finite check). -/
theorem burden_ge_eleven (s : Finset ℤ) (h7 : s.card = 7) :
    11 ≤ (restrictedSum s).card := by
  have := restrictedSum_card_ge s (by omega)
  omega

end LRCFreiman
