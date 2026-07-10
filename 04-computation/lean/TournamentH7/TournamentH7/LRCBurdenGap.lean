/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: opus (LRC multi-agent project, 2026-07-09-S198)
-/
import Mathlib
import TournamentH7.LRCFreimanBurden

/-!
# The disjoint-block burden bound: a dominant gap forces high descent burden (LEM-016(i) tail)

LEM-016(i) (mac-mini) is the excess-1 stability `B ≤ 12 ⟹ diam ≤ 7` for 7-sets. Its `g > D/2`
(largest-gap-exceeds-half-diameter) case has a clean, exhaustion-free proof: split `A = L ⊔ R` at the
dominant gap; the three sum-blocks `L +̂ L`, `L + R` (Minkowski), `R +̂ R` occupy DISJOINT integer ranges,
so `B = |A +̂ A| ≥ |L +̂ L| + |L + R| + |R +̂ R| ≥ (2|L|−3) + (|L|+|R|−1) + (2|R|−3) = 3|A| − 7`.

For `|A| = 7` this is `B ≥ 14 > 12`, so **`B ≤ 13 ⟹` no gap exceeds `D/2`** — reducing LEM-016(i) to its
flagged all-gaps-`≤ D/2` sliver. The cross-block bound `|L + R| ≥ |L|+|R|−1` is Mathlib's
`cauchy_davenport_add_of_linearOrder_isCancelAdd`; the within-block bounds are `restrictedSum_card_ge`.

Kernel-pure: no `sorry`, no `native_decide`. Axioms: `[propext, Classical.choice, Quot.sound]`.
-/

namespace LRCFreiman

open Finset Pointwise

/-- Unconditional Freiman lower bound (the `card ≤ 1` cases are vacuous in `ℕ`). -/
theorem restrictedSum_card_ge' (s : Finset ℤ) : 2 * s.card - 3 ≤ (restrictedSum s).card := by
  rcases Nat.lt_or_ge s.card 2 with h | h
  · omega
  · exact restrictedSum_card_ge s h

/-- Restricted sumset is monotone in the set. -/
theorem restrictedSum_mono {L s : Finset ℤ} (h : L ⊆ s) : restrictedSum L ⊆ restrictedSum s := by
  intro z hz
  rw [mem_restrictedSum] at hz ⊢
  obtain ⟨x, hx, y, hy, hlt, heq⟩ := hz
  exact ⟨x, h hx, y, h hy, hlt, heq⟩

/-- Every restricted sum of `L` is `≤ 2·max L`. -/
theorem restrictedSum_le {L : Finset ℤ} (hL : L.Nonempty) {z : ℤ} (hz : z ∈ restrictedSum L) :
    z ≤ 2 * L.max' hL := by
  rw [mem_restrictedSum] at hz
  obtain ⟨x, hx, y, hy, _, rfl⟩ := hz
  have h1 := L.le_max' x hx; have h2 := L.le_max' y hy; omega

/-- Every restricted sum of `R` is `≥ 2·min R + 1` (the two smallest are distinct). -/
theorem restrictedSum_ge {R : Finset ℤ} (hR : R.Nonempty) {v : ℤ} (hv : v ∈ restrictedSum R) :
    2 * R.min' hR + 1 ≤ v := by
  rw [mem_restrictedSum] at hv
  obtain ⟨x, hx, y, hy, hlt, rfl⟩ := hv
  have h1 := R.min'_le x hx; have h2 := R.min'_le y hy; omega

/-- Every Minkowski cross sum is `≥ min L + min R`. -/
theorem add_ge {L R : Finset ℤ} (hL : L.Nonempty) (hR : R.Nonempty) {w : ℤ} (hw : w ∈ L + R) :
    L.min' hL + R.min' hR ≤ w := by
  rw [Finset.mem_add] at hw
  obtain ⟨x, hx, y, hy, rfl⟩ := hw
  have h1 := L.min'_le x hx; have h2 := R.min'_le y hy; omega

/-- Every Minkowski cross sum is `≤ max L + max R`. -/
theorem add_le {L R : Finset ℤ} (hL : L.Nonempty) (hR : R.Nonempty) {w : ℤ} (hw : w ∈ L + R) :
    w ≤ L.max' hL + R.max' hR := by
  rw [Finset.mem_add] at hw
  obtain ⟨x, hx, y, hy, rfl⟩ := hw
  have h1 := L.le_max' x hx; have h2 := R.le_max' y hy; omega

/-- The Minkowski cross set lies in the restricted sumset of the union (all `L < R`). -/
theorem add_subset_restrictedSum {L R : Finset ℤ} (hL : L.Nonempty) (hR : R.Nonempty)
    (hsep : L.max' hL < R.min' hR) : L + R ⊆ restrictedSum (L ∪ R) := by
  intro w hw
  rw [Finset.mem_add] at hw
  obtain ⟨x, hx, y, hy, rfl⟩ := hw
  rw [mem_restrictedSum]
  exact ⟨x, mem_union_left R hx, y, mem_union_right L hy,
    lt_of_le_of_lt (L.le_max' x hx) (lt_of_lt_of_le hsep (R.min'_le y hy)), rfl⟩

/-- **The three-block bound.** When `L`'s and `R`'s sum-blocks and the cross block occupy disjoint
integer ranges (the two separation conditions `hd1`, `hd2`), the restricted sumset of `L ∪ R` contains
all three disjointly, so its size is at least their total. -/
theorem three_block_card {L R : Finset ℤ} (hL : L.Nonempty) (hR : R.Nonempty)
    (hsep : L.max' hL < R.min' hR)
    (hd1 : 2 * L.max' hL < L.min' hL + R.min' hR)
    (hd2 : L.max' hL + R.max' hR ≤ 2 * R.min' hR) :
    (restrictedSum L).card + (L + R).card + (restrictedSum R).card
      ≤ (restrictedSum (L ∪ R)).card := by
  have dj_L_LR : Disjoint (restrictedSum L) (L + R) := by
    rw [Finset.disjoint_left]; intro a ha1 ha2
    have p1 := restrictedSum_le hL ha1; have p2 := add_ge hL hR ha2; omega
  have dj_LR_R : Disjoint (L + R) (restrictedSum R) := by
    rw [Finset.disjoint_left]; intro a ha1 ha2
    have p1 := add_le hL hR ha1; have p2 := restrictedSum_ge hR ha2; omega
  have dj_L_R : Disjoint (restrictedSum L) (restrictedSum R) := by
    rw [Finset.disjoint_left]; intro a ha1 ha2
    have p1 := restrictedSum_le hL ha1; have p2 := restrictedSum_ge hR ha2; omega
  have dj_union_R : Disjoint (restrictedSum L ∪ (L + R)) (restrictedSum R) :=
    Finset.disjoint_union_left.mpr ⟨dj_L_R, dj_LR_R⟩
  have hsub : (restrictedSum L ∪ (L + R)) ∪ restrictedSum R ⊆ restrictedSum (L ∪ R) := by
    refine union_subset (union_subset ?_ ?_) ?_
    · exact restrictedSum_mono subset_union_left
    · exact add_subset_restrictedSum hL hR hsep
    · exact restrictedSum_mono subset_union_right
  calc (restrictedSum L).card + (L + R).card + (restrictedSum R).card
      = (restrictedSum L ∪ (L + R)).card + (restrictedSum R).card := by
        rw [card_union_of_disjoint dj_L_LR]
    _ = ((restrictedSum L ∪ (L + R)) ∪ restrictedSum R).card := by
        rw [card_union_of_disjoint dj_union_R]
    _ ≤ (restrictedSum (L ∪ R)).card := card_le_card hsub

/-- **LEM-016(i), the dominant-gap case (7-sets), formalized.** If a 7-set `s = L ⊔ R` splits at a gap
that DOMINATES both sides' spans — `span(L) + span(R) < gap` (⟺ the gap exceeds `D/2`) — then its descent
burden is `≥ 14`. Contrapositive: **`B ≤ 13 ⟹` no gap exceeds `D/2`**, reducing LEM-016(i) to its flagged
all-gaps-`≤ D/2` sliver. Proof: `three_block_card` + the Freiman within-block bounds
(`restrictedSum_card_ge'`) + the Cauchy–Davenport cross-block bound
(`cauchy_davenport_add_of_linearOrder_isCancelAdd`); `3·7 − 7 = 14`. -/
theorem burden_ge_of_dominant_gap {s L R : Finset ℤ} (hLR : L ∪ R = s) (hdisj : Disjoint L R)
    (hL : L.Nonempty) (hR : R.Nonempty) (hsep : L.max' hL < R.min' hR)
    (hgap : (L.max' hL - L.min' hL) + (R.max' hR - R.min' hR) < R.min' hR - L.max' hL) :
    3 * s.card - 7 ≤ (restrictedSum s).card := by
  have hℓ : L.min' hL ≤ L.max' hL := L.min'_le _ (L.max'_mem hL)
  have hr : R.min' hR ≤ R.max' hR := R.min'_le _ (R.max'_mem hR)
  have hd1 : 2 * L.max' hL < L.min' hL + R.min' hR := by omega
  have hd2 : L.max' hL + R.max' hR ≤ 2 * R.min' hR := by omega
  have h3 := three_block_card hL hR hsep hd1 hd2
  rw [hLR] at h3
  have hcL := restrictedSum_card_ge' L
  have hcR := restrictedSum_card_ge' R
  have hcLR := cauchy_davenport_add_of_linearOrder_isCancelAdd hL hR
  have hcard : L.card + R.card = s.card := by
    rw [← Finset.card_union_of_disjoint hdisj, hLR]
  omega

/-- The `k = 7` specialization (LEM-016(i), `g > D/2` case): a dominant-gap 7-set has burden `≥ 14`.
(For the full 13-set the same lemma gives `≥ 3·13 − 7 = 32` — monad's THM-682 "core `B ≥ 32`" bound
on the dominant-gap branch.) -/
theorem burden_ge_fourteen_of_dominant_gap {s L R : Finset ℤ} (h7 : s.card = 7)
    (hLR : L ∪ R = s) (hdisj : Disjoint L R) (hL : L.Nonempty) (hR : R.Nonempty)
    (hsep : L.max' hL < R.min' hR)
    (hgap : (L.max' hL - L.min' hL) + (R.max' hR - R.min' hR) < R.min' hR - L.max' hL) :
    14 ≤ (restrictedSum s).card := by
  have h := burden_ge_of_dominant_gap hLR hdisj hL hR hsep hgap
  rw [h7] at h; omega

end LRCFreiman
