/- LRCSevenModuli.lean -- opus-2026-07-17-S360 (HYP-7430).
   THE EXISTENCE STEP, REDUCED TO SEVEN MODULI.

   The denominator sieve (`sieve_one_div`) gives a lonely time `1/q` whenever
   some `q ≤ 14` divides no speed.  Naively that is twelve checks,
   `q = 2 … 14`.  It reduces to SEVEN: every `q ∈ [2,14]` divides some
   element of `{8,9,10,11,12,13,14}` (2∣8, 3∣9, 4∣8, 5∣10, 6∣12, 7∣14, and
   8…14 divide themselves), and a multiple of `q'` is a multiple of `q`.  So
   if ANY `q ∈ [2,14]` fails to divide a speed, one of the seven fails too.

   SEVEN IS EXACTLY RIGHT, not six: `{8,…,14}` is the unique minimum
   covering set, and the six-number window `{9,…,14}` misses precisely
   `q = 8` (no multiple of 8 lies in `9…14`).  Verified by exhaustive search
   over subsets of `{2,…,14}`.

   Kernel-pure target: no sorry, no native_decide. -/
import Mathlib
import TournamentH7.LonelyRunner

namespace LonelyRunner

/-- Every modulus in `[2,14]` divides one of the seven `8 … 14`. -/
theorem seven_covers (q : ℕ) (hq2 : 2 ≤ q) (hq14 : q ≤ 14) :
    ∃ q' ∈ ({8, 9, 10, 11, 12, 13, 14} : Finset ℕ), q ∣ q' := by
  interval_cases q <;> simp

/-- **THE EXISTENCE STEP, SEVEN CHECKS.**  If one of the seven moduli
`8 … 14` divides no speed, the family has a lonely time. -/
theorem lonely_of_seven_moduli {ι : Type*} (v : ι → ℤ) (q : ℕ)
    (hq : q ∈ ({8, 9, 10, 11, 12, 13, 14} : Finset ℕ))
    (hdiv : ∀ i, ¬ ((q : ℤ) ∣ v i)) :
    ∃ t : ℝ, Lonely 14 v t := by
  refine ⟨(1 : ℝ) / q, sieve_one_div 14 q v ?_ ?_ hdiv⟩
  · fin_cases hq <;> norm_num
  · fin_cases hq <;> norm_num

/-- **THE REDUCTION.**  If all seven moduli divide some speed, then so does
every `q ∈ [2,14]` — the twelve sieve conditions collapse to seven. -/
theorem all_moduli_of_seven {ι : Type*} (v : ι → ℤ)
    (h : ∀ q' ∈ ({8, 9, 10, 11, 12, 13, 14} : Finset ℕ), ∃ i, (q' : ℤ) ∣ v i)
    (q : ℕ) (hq2 : 2 ≤ q) (hq14 : q ≤ 14) :
    ∃ i, (q : ℤ) ∣ v i := by
  obtain ⟨q', hq'mem, hqq'⟩ := seven_covers q hq2 hq14
  obtain ⟨i, hi⟩ := h q' hq'mem
  exact ⟨i, dvd_trans (Int.natCast_dvd_natCast.mpr hqq') hi⟩

/-- Contrapositive: a counterexample must be blocked at all seven moduli. -/
theorem counterexample_needs_seven {ι : Type*} (v : ι → ℤ)
    (hcex : ∀ t : ℝ, ¬ Lonely 14 v t) :
    ∀ q ∈ ({8, 9, 10, 11, 12, 13, 14} : Finset ℕ), ∃ i, (q : ℤ) ∣ v i := by
  intro q hq
  by_contra h
  push_neg at h
  exact absurd (lonely_of_seven_moduli v q hq h) (by push_neg; exact hcex)

/-! ## Axiom audit -/
#print axioms seven_covers
#print axioms lonely_of_seven_moduli
#print axioms all_moduli_of_seven
#print axioms counterexample_needs_seven

end LonelyRunner
