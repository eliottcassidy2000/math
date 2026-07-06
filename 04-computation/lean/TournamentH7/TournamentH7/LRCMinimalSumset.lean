/-
  TournamentH7.LRCMinimalSumset — the minimal-sumset bound and its equality case.

  (mac-mini-2026-07-06-S20/S21, HYP-4482/4492.)  The additive-combinatorics
  anchor of the Freiman frame for the density floor (opus-S112 theta-sum):
  safe(S,2/25) = 0 ⟹ minimal doubling ⟹ (this + the AP characterization) the AP.

  * `two_mul_card_sub_one_le`: |S + S| ≥ 2|S| − 1 for nonempty finite `S ⊆ ℤ`.
  * `sumset_eq_translates`: at EQUALITY (|S+S| = 2|S|−1), the sumset is EXACTLY
    the two translates `(m + S) ∪ (S + M)`, m = min S, M = max S — the
    structural core of the classical "minimal doubling ⟹ arithmetic progression"
    (step 2 of the (U)-rigidity factoring).

  Proof: the translates each have |S| elements, both sit inside S + S, and meet
  only at m + M.  So their union has 2|S| − 1 elements inside S + S; equality
  forces S + S to BE that union.

  Pure Finset arithmetic; no analysis.  Draft: reflection
  the-extremizer-is-stricter-than-freiman-macmini-S20.
-/
import Mathlib.Algebra.Order.Group.Finset
import Mathlib.Tactic

open Finset Pointwise

namespace LonelyRunner
namespace MinimalSumset

variable {S : Finset ℤ}

/-- The two translates `(m + S)` and `(S + M)` (m = min, M = max) sit inside
    `S + S`, each of size `|S|`, and meet exactly at `{m + M}`. -/
private lemma translate_facts (hS : S.Nonempty) :
    let m := S.min' hS; let M := S.max' hS
    let A := S.image (m + ·); let B := S.image (· + M)
    A ⊆ S + S ∧ B ⊆ S + S ∧ A.card = S.card ∧ B.card = S.card ∧
      A ∩ B = {m + M} := by
  intro m M A B
  refine ⟨?_, ?_, ?_, ?_, ?_⟩
  · intro x hx
    rw [mem_image] at hx
    obtain ⟨s, hs, rfl⟩ := hx
    exact add_mem_add (S.min'_mem hS) hs
  · intro x hx
    rw [mem_image] at hx
    obtain ⟨s, hs, rfl⟩ := hx
    exact add_mem_add hs (S.max'_mem hS)
  · exact card_image_of_injective _ (add_right_injective m)
  · exact card_image_of_injective _ (fun a b h => by simpa using h)
  · apply Subset.antisymm
    · intro x hx
      rw [mem_inter, mem_image, mem_image] at hx
      obtain ⟨⟨s, hs, rfl⟩, ⟨t, ht, hteq⟩⟩ := hx
      have h1 : m ≤ t := S.min'_le t ht
      have h2 : s ≤ M := S.le_max' s hs
      have heq : m + s = t + M := hteq.symm
      have hsM : s = M := by omega
      rw [hsM, mem_singleton]
    · intro x hx
      rw [mem_singleton] at hx; subst hx
      rw [mem_inter]
      exact ⟨mem_image.mpr ⟨M, S.max'_mem hS, rfl⟩,
             mem_image.mpr ⟨m, S.min'_mem hS, by ring⟩⟩

/-- The union of the two translates has exactly `2|S| − 1` elements. -/
private lemma union_card (hS : S.Nonempty) :
    (S.image (S.min' hS + ·) ∪ S.image (· + S.max' hS)).card = 2 * S.card - 1 := by
  obtain ⟨_, _, hA, hB, hI⟩ := translate_facts hS
  have := card_union_add_card_inter
    (S.image (S.min' hS + ·)) (S.image (· + S.max' hS))
  rw [hI, card_singleton, hA, hB] at this
  omega

/-- **Minimal sumset bound.**  `2|S| − 1 ≤ |S + S|` for nonempty finite `S ⊆ ℤ`. -/
theorem two_mul_card_sub_one_le (hS : S.Nonempty) :
    2 * S.card - 1 ≤ (S + S).card := by
  obtain ⟨hA, hB, _, _, _⟩ := translate_facts hS
  calc 2 * S.card - 1
      = (S.image (S.min' hS + ·) ∪ S.image (· + S.max' hS)).card := (union_card hS).symm
    _ ≤ (S + S).card := card_le_card (union_subset hA hB)

/-- **Equality case (structural core of minimal doubling ⟹ AP).**  When
    `|S + S| = 2|S| − 1`, the sumset is EXACTLY the two translates
    `(min S + S) ∪ (S + max S)`. -/
theorem sumset_eq_translates (hS : S.Nonempty) (heq : (S + S).card = 2 * S.card - 1) :
    S + S = S.image (S.min' hS + ·) ∪ S.image (· + S.max' hS) := by
  obtain ⟨hA, hB, _, _, _⟩ := translate_facts hS
  refine (eq_of_subset_of_card_le (union_subset hA hB) ?_).symm
  rw [heq, union_card hS]

end MinimalSumset
end LonelyRunner
