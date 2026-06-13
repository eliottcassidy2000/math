/-
  TournamentH7.HSpectrumSmallN — H-spectrum at n = 3, 4, 5, 6

  Project-canon maximum H values (OEIS A038375):
    n=1: 1
    n=2: 1
    n=3: 3
    n=4: 5
    n=5: 15
    n=6: 45
    n=7: 189 (proved via Paley(7))

  We axiomatize these via `maxH_n` axioms (computationally verified by
  exhaustive enumeration), then derive H-spectrum constraints.
-/

import TournamentH7.HSpectrumClean

namespace Tournament

/-! ### Axiomatic maximum H values -/

/-- **Axiom (A038375 verified by exhaustive enumeration).** -/
axiom maxH_3 : ∀ T : Tournament 3, H T ≤ 3
axiom maxH_4 : ∀ T : Tournament 4, H T ≤ 5
axiom maxH_5 : ∀ T : Tournament 5, H T ≤ 15
axiom maxH_6 : ∀ T : Tournament 6, H T ≤ 45

/-! ### H-spectrum at n = 3 -/

/-- H-spectrum at n=3: H ∈ {1, 3}. -/
theorem H_spectrum_n3_full (T : Tournament 3) :
    1 ≤ H T ∧ H T ≤ 3 ∧ Odd (H T) ∧ H T ≠ 7 ∧ H T ≠ 21 := by
  refine ⟨?_, ?_, ?_, ?_, ?_⟩
  · exact H_ge_one T
  · exact maxH_3 T
  · exact H_odd_from_ocf T
  · exact H_ne_seven T
  · exact H_ne_twentyone T

/-- H at n=3 is 1 or 3. -/
theorem H_n3_eq_one_or_three (T : Tournament 3) : H T = 1 ∨ H T = 3 := by
  have ⟨h1, h2, hodd, _, _⟩ := H_spectrum_n3_full T
  -- H is odd, ≥ 1, ≤ 3.  Odd values in [1, 3]: 1, 3.
  rcases hodd with ⟨k, hk⟩
  omega

/-! ### H-spectrum at n = 4 -/

/-- H at n=4 is 1, 3, or 5. -/
theorem H_n4_eq_135 (T : Tournament 4) : H T = 1 ∨ H T = 3 ∨ H T = 5 := by
  have h1 := H_ge_one T
  have h2 := maxH_4 T
  have hodd := H_odd_from_ocf T
  -- H ∈ {1, 3, 5}
  rcases hodd with ⟨k, hk⟩
  have h_ne7 := H_ne_seven T
  omega

/-! ### H-spectrum at n = 5 -/

/-- H at n=5 is 1, 3, 5, 9, 11, 13, or 15 (skipping forbidden 7). -/
theorem H_n5_in_spectrum (T : Tournament 5) :
    H T = 1 ∨ H T = 3 ∨ H T = 5 ∨ H T = 9 ∨ H T = 11 ∨ H T = 13 ∨ H T = 15 := by
  have h1 := H_ge_one T
  have h2 := maxH_5 T
  have hodd := H_odd_from_ocf T
  have h_ne7 := H_ne_seven T
  rcases hodd with ⟨k, hk⟩
  omega

/-! ### H-spectrum at n = 6 -/

/-- H at n=6 is odd in [1, 45] and ∉ {7, 21}. -/
theorem H_n6_constraints (T : Tournament 6) :
    1 ≤ H T ∧ H T ≤ 45 ∧ Odd (H T) ∧ H T ≠ 7 ∧ H T ≠ 21 :=
  ⟨H_ge_one T, maxH_6 T, H_odd_from_ocf T, H_ne_seven T, H_ne_twentyone T⟩

end Tournament
