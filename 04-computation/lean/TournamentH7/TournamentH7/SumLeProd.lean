import Mathlib

/-
# The arithmetic kernel of the `c₃ ≤ H` SCC reduction (THM-1860)

kind-pasteur-2026-07-21-S128c135.  The Written-on-the-Wall-style candidate `c₃ ≤ H`
(3-cycles ≤ Hamiltonian paths, verified to n = 23) reduces to the strongly-connected case by
the strongly-connected-component decomposition:

  `c₃(T) = Σ_i c₃(Cᵢ)`  and  `H(T) = ∏_i H(Cᵢ)`,

so with the per-component base `c₃(Cᵢ) ≤ H(Cᵢ)` and `H(Cᵢ) ≥ 2` for every non-trivial
component, `Σ c₃(Cᵢ) ≤ Σ H(Cᵢ) ≤ ∏ H(Cᵢ) = H(T)`.  The last inequality — **sum ≤ product for a
list of naturals each ≥ 2** — is the reduction's arithmetic kernel, and is the piece the WOWII
"formalize" step asks for.  This file proves it, `sorry`-free.
-/

namespace SumLeProd

/-- **Sum ≤ product** for a list of naturals each at least `2`. -/
theorem sum_le_prod : ∀ (l : List ℕ), (∀ a ∈ l, 2 ≤ a) → l.sum ≤ l.prod
  | [], _ => by simp
  | a :: l, h => by
      have ha : 2 ≤ a := h a (by simp)
      cases l with
      | nil => simp
      | cons b t =>
          have hl : ∀ c ∈ (b :: t), 2 ≤ c := fun c hc => h c (List.mem_cons_of_mem a hc)
          have ih : (b :: t).sum ≤ (b :: t).prod := sum_le_prod (b :: t) hl
          have htp : 0 < t.prod :=
            List.prod_pos (fun c hc => by
              have := hl c (List.mem_cons_of_mem b hc); omega)
          have hb : 2 ≤ b := hl b (by simp)
          have hP : 2 ≤ (b :: t).prod := by
            simp only [List.prod_cons]; nlinarith [hb, htp]
          -- goal: (a :: b :: t).sum ≤ (a :: b :: t).prod
          simp only [List.sum_cons, List.prod_cons] at ih ⊢
          -- ih : b + t.sum ≤ b * t.prod ; hP : 2 ≤ b * t.prod
          -- abstract the product to a fresh atom P so nlinarith stays bilinear
          have key : ∀ P : ℕ, 2 ≤ P → b + t.sum ≤ P → a + (b + t.sum) ≤ a * P := by
            intro P hP2 hSP; nlinarith [ha, hP2, hSP]
          exact key (b * t.prod) hP ih

end SumLeProd
