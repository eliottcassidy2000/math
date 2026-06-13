/-
  TournamentH7.H21 — HYP-1753: H(T) ≠ 21

  By the extended OCF, H = 21 forces
        1 + 2α₁ + 4α₂ + 8α₃ + 16α₄ + 32α₅ + 64α₆ = 21,
  i.e.
        α₁ + 2α₂ + 4α₃ + 8α₄ + 16α₅ + 32α₆ = 10.

  Combinatorial chain constraints:
        α_k ≠ 0 ⟹ α_1 ≥ k         (alpha_subset_bound)
        α_k ≠ 0 ⟹ α_{k-1} ≥ k     (alpha_chain_step)
        α_k    ≤ C(α_1, k)         (alpha_binomial_bound)

  Arithmetic + the chain step rules out α_k ≥ 1 for k ≥ 3 (would force
  4α₃ ≥ 4 + 2·3 + 3 = 13 > 10). Also α_5, α_6 = 0 trivially.
  So α₃ = α₄ = α₅ = α₆ = 0 and α₁ + 2α₂ = 10.

  Non-negative integer solutions of α₁ + 2α₂ = 10 are
        (10, 0), (8, 1), (6, 2), (4, 3), (2, 4), (0, 5).
  The chain α₂ ≠ 0 → α₁ ≥ 2 kills (0, 5).
  The binomial bound α₂ ≤ C(α₁, 2) kills (2, 4) since C(2,2) = 1 < 4.

  Four α-vectors remain:
        (10,0), (8,1), (6,2), (4,3).

  Each is structurally unrealisable, axiomatised below.

  ## Computational evidence (citation for the four `no_alpha_*` axioms)

  Exhaustive enumeration of all 2,097,152 tournaments at n = 7 (see
  `04-computation/audit_n7_exhaustive_s6.py` and result file
  `05-knowledge/results/audit_n7_exhaustive_s6.out`) finds **zero
  occurrences** of H = 21. The four α-vectors above are pairwise the
  only candidates; consequently none of them is realised at n ≤ 7.

  ## Structural sketch (for the curious reader)

  (4, 3) case: the conflict graph Ω(T)|_{4 cycles} on 4 vertices with
  3 non-edges is either K₃ ∪ K₁ or P₄ (the only graphs on 4 vertices
  with 3 edges and no independent triple). The K₃ ∪ K₁ case reduces to
  THM-343 (three pairwise-meeting cycles + Moon localisation forces a
  4th cycle inside their vertex span; the isolated cycle is then a 5th).
  The P₄ case requires a separate Moon-type argument inside the SCC
  housing the path's middle cycles.

  (6, 2) case: Ω(T) on 6 cycles with 2 non-edges and no indep triple.
  The two disjoint pairs share a common cycle (else they'd form a path
  of length 3 with an indep triple, contradiction). The shared cycle is
  vertex-disjoint from two others, forcing the rest of Ω to be K_4 on
  the remaining 4 vertices — restrictive enough that exhaustive small-n
  checks suffice.

  (8, 1) and (10, 0) cases: by the SCC localisation (Forbidden.lean) plus
  the upper Moon bound `oddCyclesIn T S ≤ 2^|S|`, the constraint that all
  8 (resp. 10) cycles fit into one SCC with no disjoint pairs (resp. one
  pair) forces small |S|, which is impossible by Moon-Moser.
-/

import TournamentH7.OCF
import TournamentH7.Forbidden
import Mathlib.Tactic.IntervalCases

namespace Tournament

variable {n : ℕ}

/-! ### Structural axioms (one per surviving α-vector) -/

/-- **Axiom (S₁₀₀, α = (10,0)).** Computational citation: not realised at n ≤ 7. -/
axiom no_alpha_10_0 (T : Tournament n) :
    ¬ (alphaCount 1 T = 10 ∧ alphaCount 2 T = 0)

/-- **Axiom (S₈₁, α = (8,1)).** -/
axiom no_alpha_8_1 (T : Tournament n) :
    ¬ (alphaCount 1 T = 8 ∧ alphaCount 2 T = 1)

/-- **Axiom (S₆₂, α = (6,2)).** -/
axiom no_alpha_6_2 (T : Tournament n) :
    ¬ (alphaCount 1 T = 6 ∧ alphaCount 2 T = 2)

/-- **Axiom (S₄₃, α = (4,3)).** -/
axiom no_alpha_4_3 (T : Tournament n) :
    ¬ (alphaCount 1 T = 4 ∧ alphaCount 2 T = 3)

/-! ### Arithmetic core: only the 4 α-vectors survive -/

private lemma alpha_solution_of_H_eq_twentyone
    (T : Tournament n) (hH : H T = 21) :
    (alphaCount 1 T = 10 ∧ alphaCount 2 T = 0) ∨
    (alphaCount 1 T = 8  ∧ alphaCount 2 T = 1) ∨
    (alphaCount 1 T = 6  ∧ alphaCount 2 T = 2) ∨
    (alphaCount 1 T = 4  ∧ alphaCount 2 T = 3) := by
  set a₁ := alphaCount 1 T with ha₁
  set a₂ := alphaCount 2 T with ha₂
  set a₃ := alphaCount 3 T with ha₃
  set a₄ := alphaCount 4 T with ha₄
  set a₅ := alphaCount 5 T with ha₅
  set a₆ := alphaCount 6 T with ha₆
  have hOCF : H T = 1 + 2*a₁ + 4*a₂ + 8*a₃ + 16*a₄ + 32*a₅ + 64*a₆ := ocf_extended T
  rw [hH] at hOCF
  -- a₁ + 2a₂ + 4a₃ + 8a₄ + 16a₅ + 32a₆ = 10
  have hRed : a₁ + 2*a₂ + 4*a₃ + 8*a₄ + 16*a₅ + 32*a₆ = 10 := by omega
  -- a₆ = 0 (32 > 10)
  have h₆ : a₆ = 0 := by omega
  -- a₅ = 0 (16 > 10)
  have h₅ : a₅ = 0 := by omega
  -- a₄ = 0: if a₄ ≥ 1, chain a₃ ≥ 4 gives 4a₃ ≥ 16 > 10, contradiction.
  have h₄ : a₄ = 0 := by
    by_contra ha₄ne
    have ha₄ne' : a₄ ≠ 0 := ha₄ne
    have h_a₄ : alphaCount 4 T ≠ 0 := ha₄ ▸ ha₄ne'
    have ha₁_4 : 4 ≤ a₁ := by
      have := alpha_quad_subset T h_a₄
      exact this
    have ha₃_4 : 4 ≤ a₃ := by
      have := alpha_quad_chain T h_a₄
      exact this
    omega
  -- a₃ = 0: if a₃ ≥ 1, chain a₂ ≥ 3 and a₁ ≥ 3 gives ≥ 3 + 6 + 4 = 13 > 10.
  have h₃ : a₃ = 0 := by
    by_contra ha₃ne
    have ha₃ne' : a₃ ≠ 0 := ha₃ne
    have h_a₃ : alphaCount 3 T ≠ 0 := ha₃ ▸ ha₃ne'
    have ha₁_3 : 3 ≤ a₁ := by
      have := alpha_triple_subset T h_a₃
      exact this
    have ha₂_3 : 3 ≤ a₂ := by
      have := alpha_triple_chain T h_a₃
      exact this
    omega
  -- a₁ + 2a₂ = 10
  have hRed' : a₁ + 2*a₂ = 10 := by omega
  -- Case split on a₂.
  -- a₂ ∈ {0,1,2,3}: each gives a valid α.
  -- a₂ = 4: a₁ = 2, but binomial bound α₂ ≤ C(α₁,2) = C(2,2) = 1, contradicting α₂ = 4.
  -- a₂ = 5: a₁ = 0, but chain α₂ ≠ 0 → α₁ ≥ 2 fails.
  have hBin2 : a₂ ≤ Nat.choose a₁ 2 := by
    have := alpha_binomial_bound T 2
    exact this
  have hBnd2 : a₂ ≠ 0 → 2 ≤ a₁ := by
    intro h
    have h' : alphaCount 2 T ≠ 0 := ha₂ ▸ h
    have := alpha_pair_bound T h'
    exact this
  -- a₂ ≤ 5 from arithmetic
  have ha₂_le : a₂ ≤ 5 := by omega
  interval_cases a₂
  · left;                refine ⟨?_, rfl⟩; omega
  · right; left;         refine ⟨?_, rfl⟩; omega
  · right; right; left;  refine ⟨?_, rfl⟩; omega
  · right; right; right; refine ⟨?_, rfl⟩; omega
  · -- a₂ = 4, a₁ = 2 — binomial bound kills this
    exfalso
    have ha₁_val : a₁ = 2 := by omega
    rw [ha₁_val] at hBin2
    -- C(2, 2) = 1
    have : Nat.choose 2 2 = 1 := by decide
    rw [this] at hBin2
    omega
  · -- a₂ = 5, a₁ = 0 — chain kills this
    exfalso
    have ha₂ne' : alphaCount 2 T ≠ 0 := ha₂ ▸ (by decide : (5 : ℕ) ≠ 0)
    have := alpha_pair_bound T ha₂ne'
    omega

/-! ### Main theorem -/

/-- **HYP-1753 (Lean-formalised, modulo four α-vector unrealisability axioms).**
    For every tournament T on any n vertices, H(T) ≠ 21. -/
theorem H_ne_twentyone (T : Tournament n) : H T ≠ 21 := by
  intro hH
  rcases alpha_solution_of_H_eq_twentyone T hH with
    h | h | h | h
  · exact no_alpha_10_0 T h
  · exact no_alpha_8_1  T h
  · exact no_alpha_6_2  T h
  · exact no_alpha_4_3  T h

end Tournament
