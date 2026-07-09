/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: opus (LRC multi-agent project, 2026-07-09-S183)
-/
import Mathlib

/-!
# The interval maximizes Schur triples (LEM-014, upper bound)

The LRC density-floor resonance sum `R` (with `L(S) = (6/7)^13 + R`) is a theta-sum over the resonance
lattice `Λ = {n : Σ nᵢ vᵢ = 0}`, dominated by its minimal vectors — the height-3 **Schur triples**
`a + b = c` (opus-S182). Their count is the leading resonance coefficient, so the tight extremal
maximizes it. This file formalizes the **upper bound**

  `E₃(S) := #{(a,b) ∈ S×S : a+b ∈ S} ≤ C(k,2)`,  `k = |S|`,

for a finite set `S` of positive naturals (`0 ∉ S`). The interval `{1,…,k}` attains `C(k,2)` (e.g.
`E₃({1,…,13}) = 78`), so this bound is tight and identifies the AP as the Schur-triple maximizer.

**Proof.** Map each counted pair `(a,b)` to the 2-element set `{a, a+b} ⊆ S`. Since every `b ∈ S` is
positive, `a < a+b`, so `{a, a+b}` has two elements and the map is injective (the smaller element is `a`,
the larger is `a+b`, from which `b` is recovered). Hence `E₃(S)` is at most the number of 2-element
subsets of `S`, which is `C(k,2)`.

Kernel-pure: no `sorry`, no `native_decide`. Axioms: `[propext, Classical.choice, Quot.sound]`.
-/

namespace LRCSchurTriples

open Finset

/-- **The interval maximizes Schur triples (LEM-014, upper bound).**
For a finite set `S` of positive naturals, the number of ordered pairs `(a,b) ∈ S×S` with `a+b ∈ S`
(the Schur-triple count `E₃`) is at most `C(|S|, 2)`. Equality holds for the interval `{1,…,|S|}`. -/
theorem schurTriple_card_le (S : Finset ℕ) (hS : 0 ∉ S) :
    ((S ×ˢ S).filter (fun p => p.1 + p.2 ∈ S)).card ≤ S.card.choose 2 := by
  rw [← Finset.card_powersetCard 2 S]
  apply Finset.card_le_card_of_injOn (fun p => ({p.1, p.1 + p.2} : Finset ℕ))
  · -- the map lands in the 2-element subsets of `S`
    rintro ⟨a, b⟩ hp
    obtain ⟨hmem, habS⟩ := Finset.mem_filter.mp hp
    obtain ⟨haS, hbS⟩ := Finset.mem_product.mp hmem
    have hb : 0 < b := Nat.pos_of_ne_zero (fun h => hS (h ▸ hbS))
    show ({a, a + b} : Finset ℕ) ∈ Finset.powersetCard 2 S
    rw [Finset.mem_powersetCard]
    refine ⟨?_, ?_⟩
    · intro x hx
      rw [Finset.mem_insert, Finset.mem_singleton] at hx
      rcases hx with rfl | rfl
      · exact haS
      · exact habS
    · rw [Finset.card_pair (by omega : a ≠ a + b)]
  · -- injectivity: `{a, a+b} = {a', a'+b'}` with `b, b' > 0` forces `(a,b) = (a',b')`
    rintro ⟨a, b⟩ hp ⟨a', b'⟩ hp' heq
    dsimp only at heq
    simp only [Finset.coe_filter, Finset.mem_product, Set.mem_setOf_eq] at hp hp'
    have hb : 0 < b := Nat.pos_of_ne_zero (fun h => hS (h ▸ hp.1.2))
    have hb' : 0 < b' := Nat.pos_of_ne_zero (fun h => hS (h ▸ hp'.1.2))
    have e1 : a ∈ ({a', a' + b'} : Finset ℕ) := by rw [← heq]; simp
    have e2 : a + b ∈ ({a', a' + b'} : Finset ℕ) := by rw [← heq]; simp
    rw [Finset.mem_insert, Finset.mem_singleton] at e1 e2
    have hab : a = a' ∧ b = b' := by
      rcases e1 with h1 | h1 <;> rcases e2 with h2 | h2 <;> omega
    obtain ⟨rfl, rfl⟩ := hab
    rfl

/-- The interval `{1,…,k}` attains the bound: `E₃({1,…,k}) = C(k,2)`.
(Stated for the concrete `k = 13` LRC case; `E₃ = 78`.) -/
theorem schurTriple_interval_13 :
    (((Finset.Icc 1 13) ×ˢ (Finset.Icc 1 13)).filter (fun p => p.1 + p.2 ∈ Finset.Icc 1 13)).card
      = 78 := by decide

end LRCSchurTriples
