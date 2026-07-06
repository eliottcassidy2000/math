/-
  TournamentH7.LRCMinimalSumset — the minimal-sumset bound |S+S| ≥ 2|S|−1.

  (mac-mini-2026-07-06-S20, HYP-4482.)  The additive-combinatorics anchor of the
  Freiman frame for the density floor (opus-S112 theta-sum): safe(S,2/25) = 0
  ⟹ minimal doubling ⟹ (this + the AP characterization) the AP.  This file is
  step 2's LOWER bound: for a nonempty finite set S of integers,
  |S + S| ≥ 2|S| − 1, with the AP {1,…,n} attaining equality.

  Proof: the translates (m + S) and (S + M), m = min S, M = max S, each have
  |S| elements, both sit inside S + S, and meet only at m + M.  So their union
  has 2|S| − 1 elements inside S + S.

  Pure Finset arithmetic; no analysis.  Draft: reflection
  the-extremizer-is-stricter-than-freiman-macmini-S20.
-/
import Mathlib.Algebra.Order.Group.Finset
import Mathlib.Tactic

open Finset Pointwise

namespace LonelyRunner
namespace MinimalSumset

variable {S : Finset ℤ}

/-- **Minimal sumset bound.**  A nonempty finite set of integers `S` has
    `2·|S| − 1 ≤ |S + S|`. -/
theorem two_mul_card_sub_one_le (hS : S.Nonempty) :
    2 * S.card - 1 ≤ (S + S).card := by
  obtain ⟨m, hm, hmin⟩ := S.exists_min_image id hS
  obtain ⟨M, hM, hmax⟩ := S.exists_max_image id hS
  -- the two translates, both inside S + S
  set A := S.image (m + ·) with hA
  set B := S.image (· + M) with hB
  have hAsub : A ⊆ S + S := by
    intro x hx
    rw [hA, mem_image] at hx
    obtain ⟨s, hs, rfl⟩ := hx
    exact add_mem_add hm hs
  have hBsub : B ⊆ S + S := by
    intro x hx
    rw [hB, mem_image] at hx
    obtain ⟨s, hs, rfl⟩ := hx
    exact add_mem_add hs hM
  have hAcard : A.card = S.card :=
    card_image_of_injective _ (add_right_injective m)
  have hBcard : B.card = S.card :=
    card_image_of_injective _ (fun a b h => by simpa using h)
  -- A ∩ B = {m + M}
  have hinter : A ∩ B = {m + M} := by
    apply Subset.antisymm
    · intro x hx
      rw [mem_inter, hA, hB, mem_image, mem_image] at hx
      obtain ⟨⟨s, hs, rfl⟩, ⟨t, ht, hteq⟩⟩ := hx
      -- m + s = t + M, with m ≤ t and s ≤ M ⟹ s = M ∧ t = m
      have h1 : (id m : ℤ) ≤ id t := hmin t ht
      have h2 : (id s : ℤ) ≤ id M := hmax s hs
      simp only [id] at h1 h2
      have heq : m + s = t + M := hteq.symm
      have hsM : s = M := by omega
      rw [hsM, mem_singleton]
    · intro x hx
      rw [mem_singleton] at hx; subst hx
      rw [mem_inter]
      constructor
      · rw [hA, mem_image]; exact ⟨M, hM, rfl⟩
      · rw [hB, mem_image]; exact ⟨m, hm, by ring⟩
  -- |A ∪ B| = |A| + |B| − |A ∩ B| = 2|S| − 1, and A ∪ B ⊆ S + S
  have hunion : A ∪ B ⊆ S + S := union_subset hAsub hBsub
  have hcard_union : (A ∪ B).card = 2 * S.card - 1 := by
    have := card_union_add_card_inter A B
    rw [hinter, card_singleton, hAcard, hBcard] at this
    omega
  calc 2 * S.card - 1 = (A ∪ B).card := hcard_union.symm
    _ ≤ (S + S).card := card_le_card hunion

/- Equality holds for the AP `{1, …, n}` (sumset `{2, …, 2n}`, size `2n − 1`) —
   the density-floor extremizer is Freiman-minimal.  Verified numerically at
   n = 12 (`lrc_freiman_rigidity_macmini_S20.py`): the AP is the UNIQUE primitive
   minimal-doubling family, but `safe = 0` is STRICTER (shifted APs are
   minimal-doubling yet `safe > 0`; the extremizer needs first-term = difference,
   i.e. residue-completeness mod 13). -/

end MinimalSumset
end LonelyRunner
