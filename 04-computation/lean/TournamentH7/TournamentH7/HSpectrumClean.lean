/-
  TournamentH7.HSpectrumClean — H-spectrum theorems with minimal axioms

  This module collects H-spectrum results using ONLY OCF + structural
  axioms.  No `redei_existence` or `redei_parity` dependency.

  All theorems here have axiom dependencies in the family
    {propext, Classical.choice, Quot.sound,
     ocf, alpha_*, no_alpha_*, moonMoser, moonCamion_oddSize,
     omegaTriangleLocalises, oddCyclesIn_*, paley_7_maximises_H, ...}
  but NOT {redei_existence, redei_parity}.
-/

import TournamentH7.RedeiFromOCF
import TournamentH7.H7
import TournamentH7.H21
import TournamentH7.H63
import TournamentH7.PaleyAxiomatic

namespace Tournament

variable {n : ℕ}

/-! ### Cleaner H-spectrum at n = 7 -/

/-- The H-spectrum at n = 7: {1, 3, 5, 9, 11, 13, ..., 189} ∖ {7, 21, 63}.

    Uses OCF-derived Rédei facts (no redei_existence/parity axiom). -/
theorem H_spectrum_n7_clean (T : Tournament 7) :
    1 ≤ H T ∧ H T ≤ 189 ∧ Odd (H T) ∧ H T ≠ 7 ∧ H T ≠ 21 ∧ H T ≠ 63 := by
  refine ⟨?_, ?_, ?_, ?_, ?_, ?_⟩
  · exact H_ge_one T
  · exact paley_7_maximises_H T
  · exact H_odd_from_ocf T
  · exact H_ne_seven T
  · exact H_ne_twentyone T
  · exact H_ne_sixtythree_le_seven (by omega) T

/-! ### Universal H-spectrum constraints (any n) -/

/-- For every tournament T on any n vertices:
    H(T) ≥ 1 ∧ H(T) is odd ∧ H(T) ≠ 7 ∧ H(T) ≠ 21. -/
theorem H_spectrum_universal (T : Tournament n) :
    1 ≤ H T ∧ Odd (H T) ∧ H T ≠ 7 ∧ H T ≠ 21 := by
  exact ⟨H_ge_one T, H_odd_from_ocf T, H_ne_seven T, H_ne_twentyone T⟩

/-- H ∉ {7, 21} universally; H ∉ {63} for n ≤ 7. -/
theorem H_universal_constraints (T : Tournament n) :
    (∀ m, Even m → H T ≠ m) ∧ H T ≠ 7 ∧ H T ≠ 21 := by
  refine ⟨?_, H_ne_seven T, H_ne_twentyone T⟩
  intro m hm h
  have hodd := H_odd_from_ocf T
  rw [h] at hodd
  exact (Nat.not_odd_iff_even.mpr hm) hodd

/-! ### H = 1 ⟺ "trivial transitive class" -/

/-- If H(T) = 1, then α_k = 0 for all k ≥ 1.

    Proof: H = 1 = 1 + 2α_1 + 4α_2 + … ⟹ 2α_1 + … = 0 ⟹ all α_k = 0. -/
theorem alpha_solution_H1 (T : Tournament n) (hH : H T = 1) :
    alphaCount 1 T = 0 ∧ alphaCount 2 T = 0
       ∧ alphaCount 3 T = 0 ∧ alphaCount 4 T = 0 := by
  have hocf := ocf T
  refine ⟨?_, ?_, ?_, ?_⟩ <;> omega

/-! ### Concrete even-H non-realizability -/

theorem H_ne_zero_clean (T : Tournament n) : H T ≠ 0 := H_ne_zero T

theorem H_ne_two_clean (T : Tournament n) : H T ≠ 2 :=
  H_ne_two_from_ocf T

theorem H_ne_four (T : Tournament n) : H T ≠ 4 := by
  intro h
  have hodd := H_odd_from_ocf T
  rw [h] at hodd
  exact (by decide : ¬ Odd 4) hodd

theorem H_ne_six (T : Tournament n) : H T ≠ 6 := by
  intro h
  have hodd := H_odd_from_ocf T
  rw [h] at hodd
  exact (by decide : ¬ Odd 6) hodd

theorem H_ne_eight (T : Tournament n) : H T ≠ 8 := by
  intro h
  have hodd := H_odd_from_ocf T
  rw [h] at hodd
  exact (by decide : ¬ Odd 8) hodd

theorem H_ne_ten (T : Tournament n) : H T ≠ 10 := by
  intro h
  have hodd := H_odd_from_ocf T
  rw [h] at hodd
  exact (by decide : ¬ Odd 10) hodd

/-! ### Combined forbidden set -/

/-- For every tournament T: H(T) ∉ {0, 2, 4, 6, 7, 8, 10, 21}. -/
theorem H_in_forbidden_set (T : Tournament n) :
    H T ≠ 0 ∧ H T ≠ 2 ∧ H T ≠ 4 ∧ H T ≠ 6 ∧ H T ≠ 7 ∧ H T ≠ 8 ∧
    H T ≠ 10 ∧ H T ≠ 21 :=
  ⟨H_ne_zero_clean T, H_ne_two_clean T, H_ne_four T, H_ne_six T,
   H_ne_seven T, H_ne_eight T, H_ne_ten T, H_ne_twentyone T⟩

end Tournament
