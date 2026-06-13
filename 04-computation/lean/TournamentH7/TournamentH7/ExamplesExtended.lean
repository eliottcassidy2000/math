/-
  TournamentH7.ExamplesExtended — Worked examples using the extended
  framework (forbidden values, N_min, iso classes, etc.)

  All examples use only project-axiomatic facts; many are FULLY proved.
-/

import TournamentH7.Examples
import TournamentH7.ForbiddenHCounting
import TournamentH7.HSpectrumExtended
import TournamentH7.IsomorphismClasses
import TournamentH7.PaleyAxiomatic

namespace Tournament

/-! ### Example 1: forbidden-trio + parity at n = 7 -/

/-- Every tournament on 7 vertices has H ∈ {1, 3, 5, 9, 11, 13, 15, 17, 19, …, 189}
    (odd, not in {7, 21, 63}). -/
example (T : Tournament 7) :
    Odd (H T) ∧ H T ≤ 189 ∧ H T ≠ 7 ∧ H T ≠ 21 ∧ H T ≠ 63 := by
  obtain ⟨_, h_le, h_odd, h_ne7, h_ne21, h_ne63⟩ := H_spectrum_n7 T
  exact ⟨h_odd, h_le, h_ne7, h_ne21, h_ne63⟩

/-! ### Example 2: N_min(k) corollaries -/

/-- A tournament with H = 5 has only 2 odd cycles (α_1 = 2) and no
    independent pairs. -/
example (T : Tournament 7) (hH : H T = 5) :
    alphaCount 1 T = 2 ∧ alphaCount 2 T = 0 := by
  obtain ⟨h1, h2, _, _⟩ := alpha_solution_H5 T hH
  exact ⟨h1, h2⟩

/-- A tournament with H = 3 is the trivial transitive — has exactly 1 odd cycle. -/
example (T : Tournament 7) (hH : H T = 3) : alphaCount 1 T = 1 := by
  obtain ⟨h1, _, _, _⟩ := alpha_solution_H3 T hH
  exact h1

/-! ### Example 3: H ≥ 27 needed for α_3 ≥ 1 -/

/-- A tournament with H ≤ 26 has α_3 = 0 (no vertex-disjoint odd triple). -/
example (T : Tournament 7) (hH : H T ≤ 26) : alphaCount 3 T = 0 :=
  H_lt_27_no_alpha3 T (by omega)

/-! ### Example 4: iso-class structure -/

/-- A000568(7) = 456. -/
example : numIsoClasses 7 = 456 := numIsoClasses_7

/-- numNS at n = 7 is 448 (= 456 − 8). -/
example : numNS 7 = 448 := by simp [numNS]

/-! ### Example 5: Paley(7) -/

/-- Paley(7) is regular (on 7 vertices with score (3, 3, …, 3)). -/
example : IsRegular paley_7.T := paley_7.isRegular

/-- Paley(7) has H = 189. -/
example : H paley_7.T = 189 := H_paley_7

/-- Paley(7) is the H-maximiser at n = 7. -/
example (T : Tournament 7) : H T ≤ H paley_7.T := by
  rw [H_paley_7]
  exact paley_7_maximises_H T

/-! ### Example 6: Forbidden trio implies α_3 = 0 -/

/-- If H is forbidden (7 or 21), then α_3 = 0. -/
example (T : Tournament n) (h : H T = 7 ∨ H T = 21) : alphaCount 3 T = 0 :=
  forbidden_pair_alpha3_zero T h

/-! ### Example 7: H = 21 is forbidden -/

/-- H = 21 has no realization (HYP-1753). -/
example (T : Tournament n) : H T ≠ 21 := H_ne_twentyone T

end Tournament
