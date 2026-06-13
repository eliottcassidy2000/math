/-
  TournamentH7.H7 — Main theorem: H(T) ≠ 7

  Formalization of THM-343 (project canon:
  `01-canon/theorems/THM-343-H7-impossible.md`).

  Proof skeleton:

  Suppose H(T) = 7. By OCF
        7 = 1 + 2α₁ + 4α₂ + 8α₃ + 16α₄
    ⟹ α₁ + 2α₂ + 4α₃ + 8α₄ = 3.

  The non-negative integer solutions of this equation are
    (3, 0, 0, 0)  and  (1, 1, 0, 0)
  but (1, 1, 0, 0) is excluded by the subset bound `α₂ ≠ 0 → α₁ ≥ 2`.
  So α₁ = 3 and α_k = 0 for k ≥ 2.

  Hence T has exactly three odd directed cycles, pairwise vertex-meeting.
  By A4 they all lie in a single SCC `S` of T with `|S| = s ≥ 3`.

  Case-split on s:
    • s = 3: A5a gives oddCyclesIn = 1 ≠ 3 — contradiction.
    • s = 4: A5b gives oddCyclesIn = 2 ≠ 3 — contradiction.
    • s = 5: A3 gives oddCyclesIn ≥ 4 > 3 — contradiction.
    • s ≥ 6: A2 gives oddCyclesIn ≥ s − 2 ≥ 4 > 3 — contradiction.

  Q.E.D.
-/

import TournamentH7.OCF
import Mathlib.Tactic.IntervalCases
import Mathlib.Algebra.Ring.Parity

namespace Tournament

variable {n : ℕ}

/-! ### Arithmetic core: 7 = 1 + 2α₁ + 4α₂ + 8α₃ + 16α₄ + subset bound ⟹ α₁ = 3 -/

private lemma alpha_solution_of_H_eq_seven
    (a₁ a₂ a₃ a₄ : ℕ)
    (hSum : 1 + 2 * a₁ + 4 * a₂ + 8 * a₃ + 16 * a₄ = 7)
    (hBnd : a₂ ≠ 0 → 2 ≤ a₁) :
    a₁ = 3 ∧ a₂ = 0 ∧ a₃ = 0 ∧ a₄ = 0 := by
  -- The reduced sum.
  have hRed : a₁ + 2 * a₂ + 4 * a₃ + 8 * a₄ = 3 := by omega
  -- From hRed, a₃ = a₄ = 0 (since 4·1 = 4 > 3 and 8·1 = 8 > 3).
  have h₄ : a₄ = 0 := by omega
  have h₃ : a₃ = 0 := by omega
  -- So a₁ + 2 a₂ = 3.
  have hRed' : a₁ + 2 * a₂ = 3 := by omega
  -- Two non-neg solutions: (3,0) and (1,1). The bound rules out (1,1).
  by_cases ha₂ : a₂ = 0
  · refine ⟨?_, ha₂, h₃, h₄⟩
    omega
  · -- a₂ ≠ 0, so a₁ ≥ 2. But hRed' forces a₂ = 1 and a₁ = 1, contradicting a₁ ≥ 2.
    have ha₁_ge2 : 2 ≤ a₁ := hBnd ha₂
    have : a₂ ≥ 1 := Nat.one_le_iff_ne_zero.mpr ha₂
    exfalso; omega

/-! ### Main theorem -/

/-- **THM-343 (formalised).**  For every tournament T on any n vertices,
    the Hamiltonian path count H(T) is never equal to 7. -/
theorem H_ne_seven (T : Tournament n) : H T ≠ 7 := by
  intro hH
  -- Step 1: OCF.
  have hOCF : 1 + 2 * alphaCount 1 T + 4 * alphaCount 2 T
              + 8 * alphaCount 3 T + 16 * alphaCount 4 T = 7 := by
    have := ocf T; omega
  -- Step 2: subset bound for k = 2 (specialised).
  have hBnd : alphaCount 2 T ≠ 0 → 2 ≤ alphaCount 1 T := by
    intro h
    have h12 : (1 : ℕ) ≤ 2 := by decide
    exact alpha_subset_bound T 2 h12 h
  -- Step 3: solve for α₁, α₂, α₃, α₄.
  obtain ⟨hα₁, hα₂, hα₃, hα₄⟩ :=
    alpha_solution_of_H_eq_seven _ _ _ _ hOCF hBnd
  -- Step 4: localise the three cycles into a single SCC S.
  obtain ⟨S, hSCC, hS3, hOdd⟩ := omegaTriangleLocalises T hα₁ hα₂
  -- Step 5: case-split on s = |S|.
  rcases Nat.lt_or_ge S.card 5 with hsmall | hlarge
  · -- Cases s = 3 or s = 4.
    have hcases : S.card = 3 ∨ S.card = 4 := by omega
    rcases hcases with h3 | h4
    · -- s = 3
      have h := oddCyclesIn_size3 T S hSCC h3
      omega
    · -- s = 4
      have h := oddCyclesIn_size4 T S hSCC h4
      omega
  · -- Cases s ≥ 5.
    rcases Nat.even_or_odd S.card with hev | hodd
    · -- s even ⟹ s ≥ 6 (since s ≥ 5 even ⟹ s ≥ 6).
      have hs6 : 6 ≤ S.card := by
        rcases Nat.lt_or_ge S.card 6 with hlt | hge
        · exfalso
          -- s ≥ 5, s < 6, so s = 5; but 5 is not even.
          have h5 : S.card = 5 := by omega
          rw [h5] at hev
          exact (by decide : ¬ Even 5) hev
        · exact hge
      have hM := moonMoser T S hSCC (by omega)
      omega
    · -- s odd ∧ s ≥ 5: Moon–Camion gives oddCyclesIn ≥ s − 1 ≥ 4.
      have hMC := moonCamion_oddSize T S hSCC hodd hlarge
      omega

end Tournament
