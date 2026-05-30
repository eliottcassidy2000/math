/-
  TournamentH7.HSpectrumExtended — Extended H-spectrum theorems

  Combines the universal forbidden-H theorems with the Paley(7) maximum
  to give sharper H-spectrum constraints at n = 7.
-/

import TournamentH7.HSpectrum
import TournamentH7.PaleyAxiomatic
import TournamentH7.Redei
import TournamentH7.RedeiFromOCF
import TournamentH7.ForbiddenHCounting

namespace Tournament

/-! ### H-spectrum at n = 7 -/

/-- **Theorem.** At n = 7, H(T) ∈ {1, 3, 5, …, 189} ∖ {7, 21, 63}.

    The lower bound (≥ 1) is Rédei's theorem.
    The parity (odd) is Rédei's theorem.
    The upper bound (≤ 189) is Paley(7) maximality.
    The exclusions (7, 21, 63) are THM-343, HYP-1753, and the finite n=7 verification. -/
theorem H_spectrum_n7 (T : Tournament 7) :
    1 ≤ H T ∧ H T ≤ 189 ∧ Odd (H T) ∧ H T ≠ 7 ∧ H T ≠ 21 ∧ H T ≠ 63 := by
  refine ⟨?_, ?_, ?_, ?_, ?_, ?_⟩
  · exact H_ge_one T
  · exact paley_7_maximises_H T
  · exact H_odd_from_ocf T
  · exact H_ne_seven T
  · exact H_ne_twentyone T
  · exact H_ne_sixtythree_le_seven (by omega) T

/-! ### H-spectrum at small n -/

/-- For n = 3: H(T) ∈ {1, 3}. -/
theorem H_spectrum_n3 (T : Tournament 3) :
    1 ≤ H T ∧ Odd (H T) ∧ H T ≠ 7 ∧ H T ≠ 21 := by
  refine ⟨?_, ?_, ?_, ?_⟩
  · exact H_ge_one T
  · exact H_odd_from_ocf T
  · exact H_ne_seven T
  · exact H_ne_twentyone T

/-! ### Connection to ForbiddenHCounting phase transitions -/

/-- For H(T) < 27, no triple of disjoint odd cycles.  So H = 3, 5, 9, 11,
    13, 15, 17, 19, 23, 25 only allow α-tuples with α_3 = 0. -/
theorem small_H_implies_alpha3_zero {n : ℕ} (T : Tournament n) (hH : H T < 27) :
    alphaCount 3 T = 0 :=
  H_lt_27_no_alpha3 T hH

/-- For H(T) < 81, no quadruple of disjoint odd cycles. -/
theorem H_lt_81_implies_alpha4_zero {n : ℕ} (T : Tournament n) (hH : H T < 81) :
    alphaCount 4 T = 0 :=
  H_lt_81_no_alpha4 T hH

end Tournament
