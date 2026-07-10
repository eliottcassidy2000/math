/-
  TournamentH7.LRCCommonResidue — THE COMMON-RESIDUE DISPATCH (THM-682(a)).
  monad-explorer-2026-07-09-S12 (HYP-5827).

  If all thirteen speeds share a NONTRIVIAL residue class mod d (some d ≥ 2, some a
  with d ∤ a, and d ∣ v i − a for every i), then the family is lonely — with clearance
  ≥ 1/4, far above 1/14.

  Mechanism: write g = gcd(a, d), a = g·a', d = g·d'.  Since d ∤ a we get d' ≥ 2, and
  gcd(a', d') = 1.  Choose c with a'·c ≡ ⌊d'/2⌋ (mod d') (Bezout).  At t = c/d every
  runner's phase is  v i · t = a'c/d' + (integer):  ALL THIRTEEN PHASES COINCIDE at
  ⌊d'/2⌋/d' ∈ [1/4, 1/2] — the quarter window (LRCDetunedDispatch.quarter_window).

  This is STRONGER than the paper form of THM-682(a): no covering or primitivity
  hypothesis is needed (those only sharpen the clearance from 1/4 to 8/17).

  Affine sibling of THM-668 (`lonely14_of_detuned`): there d divides all but one VALUE;
  here d divides all DIFFERENCES v i − a.

  Kernel-pure: no native_decide, no sorry, no axioms beyond propext/choice/Quot.sound.
-/
import Mathlib
import TournamentH7.LonelyRunner
import TournamentH7.LRCDetunedDispatch

namespace LonelyRunner
namespace CommonResidue

open LonelyRunner

/-- **The common-residue dispatch** (THM-682(a), Lean form).  A nontrivial common
residue class mod `d ≥ 2` forces loneliness at `t = c/d` with clearance ≥ 1/4:
all phases coincide in the quarter window. -/
theorem lonely_of_common_residue (v : Fin 13 → ℤ) (d a : ℤ)
    (hd : 2 ≤ d) (hna : ¬ d ∣ a) (hres : ∀ i, d ∣ (v i - a)) :
    ∃ t : ℝ, Lonely 14 v t := by
  -- the gcd split a = g·a', d = g·d'
  set g : ℤ := (Int.gcd a d : ℤ) with hgdef
  have hga : g ∣ a := Int.gcd_dvd_left a d
  have hgd : g ∣ d := Int.gcd_dvd_right a d
  have hdpos : (0 : ℤ) < d := by omega
  have hgnat : 0 < Int.gcd a d := by
    apply Nat.pos_of_ne_zero
    intro h0
    have : d = 0 := (Int.gcd_eq_zero_iff.mp h0).2
    omega
  have hgpos : (0 : ℤ) < g := by rw [hgdef]; exact_mod_cast hgnat
  set a' : ℤ := a / g with ha'def
  set d' : ℤ := d / g with hd'def
  have hA : g * a' = a := Int.mul_ediv_cancel' hga
  have hD : g * d' = d := Int.mul_ediv_cancel' hgd
  have hd'pos : (0 : ℤ) < d' := by
    rcases lt_trichotomy d' 0 with h | h | h
    · nlinarith
    · rw [h, mul_zero] at hD; omega
    · exact h
  -- d' = 1 would give d = g ∣ a, contradiction
  have hd'2 : (2 : ℤ) ≤ d' := by
    rcases (by omega : d' = 1 ∨ 2 ≤ d') with h1 | h2
    · exfalso
      apply hna
      rw [h1, mul_one] at hD
      rw [← hD]
      exact hga
    · exact h2
  -- coprimality of the reduced pair
  have hcop : Int.gcd a' d' = 1 := by
    have := Int.gcd_div_gcd_div_gcd (i := a) (j := d) hgnat
    simpa [ha'def, hd'def, hgdef] using this
  obtain ⟨u, w, huw⟩ := Int.isCoprime_iff_gcd_eq_one.mpr hcop
  -- the half-point target and the Bezout multiplier
  set m₀ : ℤ := d' / 2 with hm0def
  have hm0 : 2 * m₀ ≤ d' ∧ d' ≤ 4 * m₀ := by
    have := Int.mul_ediv_add_emod d' 2
    have hmod : 0 ≤ d' % 2 ∧ d' % 2 < 2 := ⟨Int.emod_nonneg d' (by norm_num),
      Int.emod_lt_of_pos d' (by norm_num)⟩
    constructor <;> omega
  set c : ℤ := u * m₀ with hcdef
  set W : ℤ := w * m₀ with hWdef
  -- the exact Bezout identity: a'·c = m₀ − W·d'
  have key : a' * c = m₀ - W * d' := by
    have : (u * a' + w * d') * m₀ = m₀ := by rw [huw, one_mul]
    calc a' * c = (u * a') * m₀ := by rw [hcdef]; ring
    _ = m₀ - (w * d') * m₀ := by linarith [this]
    _ = m₀ - W * d' := by rw [hWdef]; ring
  -- the witness time
  refine ⟨(c : ℝ) / (d : ℝ), ?_⟩
  intro i m
  obtain ⟨k, hk⟩ := hres i
  have hvi : v i = a + d * k := by linarith [hk]
  -- real-cast nonvanishing
  have hdR : ((d : ℝ)) ≠ 0 := by
    have : (0 : ℝ) < (d : ℝ) := by exact_mod_cast hdpos
    linarith
  have hd'R : ((d' : ℝ)) ≠ 0 := by
    have : (0 : ℝ) < (d' : ℝ) := by exact_mod_cast hd'pos
    linarith
  have hgR : ((g : ℝ)) ≠ 0 := by
    have : (0 : ℝ) < (g : ℝ) := by exact_mod_cast hgpos
    linarith
  -- the phase identity: v i · t − m = m₀/d' − (W − k·c + m)
  have hmain : (v i : ℝ) * ((c : ℝ) / (d : ℝ)) - (m : ℝ)
      = (m₀ : ℝ) / (d' : ℝ) - ((W - k * c + m : ℤ) : ℝ) := by
    have hAR : ((g : ℝ)) * ((a' : ℝ)) = (a : ℝ) := by exact_mod_cast hA
    have hDR : ((g : ℝ)) * ((d' : ℝ)) = (d : ℝ) := by exact_mod_cast hD
    have hKR : ((a' : ℝ)) * ((c : ℝ)) = (m₀ : ℝ) - (W : ℝ) * (d' : ℝ) := by
      exact_mod_cast key
    have hVR : ((v i : ℤ) : ℝ) = (a : ℝ) + (d : ℝ) * (k : ℝ) := by exact_mod_cast hvi
    push_cast
    rw [hVR]
    field_simp
    linear_combination (-(c : ℝ) * (d' : ℝ)) * hAR + ((c : ℝ) * (a' : ℝ)) * hDR
      + (d : ℝ) * hKR
  rw [hmain]
  calc (1 : ℝ) / 14 ≤ 1 / 4 := by norm_num
  _ ≤ |(m₀ : ℝ) / (d' : ℝ) - ((W - k * c + m : ℤ) : ℝ)| :=
      DetunedDispatch.quarter_window d' m₀ hd'pos hm0.2 hm0.1 _

end CommonResidue
end LonelyRunner

#print axioms LonelyRunner.CommonResidue.lonely_of_common_residue
