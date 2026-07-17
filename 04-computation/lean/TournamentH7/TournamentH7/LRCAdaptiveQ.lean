/-
  TournamentH7.LRCAdaptiveQ — THE ADAPTIVE-Q PIGEONHOLE
  (death-star-2026-07-17-S45, HYP-7211; open item 3's named attack after the S44
  refutation-with-replacement).

  Recon truth: over 400 ladder-7 families, the worst number of ALL-COPRIME deep-
  carrying moduli in the whole q-window is ONE (the raw worst case of 161 carriers
  was a composite-`v_bot` artifact — every even q resonates with 24).  The
  formal frame, with the component count as an honest named parameter:

  * `farey_separation` — two distinct rationals `p₁/q₁ ≠ p₂/q₂` in a closed
    interval of width `w` force `1 ≤ w·q₁·q₂` (cross-multiplication; the
    denominator-bound Farey spacing in exact integers).
  * `carrier_unique_per_component` — on SUPER-LADDER strata (`D² < 7·v_top` for
    the window bound `D`), a deep component (an interval of width `≤ 1/(7·v_top)`
    hosting the pinned t's) admits at most ONE carrier modulus `q ≤ D`: two
    carriers in one component would violate Farey separation.
  * `adaptive_q_pigeonhole` — THE SKELETON: if every deep-carrying coprime
    modulus lands in one of `K` components (K = the tower-component count; the
    two-choice ladder bound gives `K ≤ 64·v_bot`, the recon says `K ≤ 2` — the
    gap is the named refinement), and the window holds more than `K` coprime
    moduli, then SOME coprime modulus in the window is deep-free — and THM-950's
    census criterion there needs only `liveCount > 0`.

  The live floor at the chosen modulus remains the program's nucleus — every
  route meets it; this module reduces "∃ good q" to it on super-ladder strata.

  Kernel-pure: no `sorry`, no `native_decide`, no new axioms.
-/
import Mathlib
import TournamentH7.LRCDeepCount

namespace LonelyRunner
namespace LRC14Concrete

open Finset

/-- **Farey separation**: distinct rationals in a width-`w` interval force
`1 ≤ w·q₁·q₂` — in exact integers, `|p₁·q₂ − p₂·q₁| ≥ 1` and the mean-value
bound. -/
theorem farey_separation (p₁ q₁ p₂ q₂ : ℕ) (w : ℚ)
    (hq₁ : 0 < q₁) (hq₂ : 0 < q₂)
    (hne : (p₁ : ℚ) / q₁ ≠ (p₂ : ℚ) / q₂)
    (hw : |(p₁ : ℚ) / q₁ - (p₂ : ℚ) / q₂| ≤ w) :
    1 ≤ w * q₁ * q₂ := by
  have hq₁Q : (0 : ℚ) < (q₁ : ℚ) := by exact_mod_cast hq₁
  have hq₂Q : (0 : ℚ) < (q₂ : ℚ) := by exact_mod_cast hq₂
  -- the cross-difference is a nonzero integer
  have hcross : (p₁ : ℤ) * q₂ ≠ (p₂ : ℤ) * q₁ := by
    intro hEq
    apply hne
    have hEqQ : (p₁ : ℚ) * q₂ = (p₂ : ℚ) * q₁ := by exact_mod_cast hEq
    rw [div_eq_div_iff (ne_of_gt hq₁Q) (ne_of_gt hq₂Q)]
    exact hEqQ
  have habs1 : (1 : ℤ) ≤ |(p₁ : ℤ) * q₂ - (p₂ : ℤ) * q₁| := by
    have : (p₁ : ℤ) * q₂ - (p₂ : ℤ) * q₁ ≠ 0 := sub_ne_zero.mpr hcross
    exact Int.one_le_abs this
  have habs1Q : (1 : ℚ) ≤ |(p₁ : ℚ) * q₂ - (p₂ : ℚ) * q₁| := by
    have h := habs1
    have hcast : ((|(p₁ : ℤ) * q₂ - (p₂ : ℤ) * q₁| : ℤ) : ℚ)
        = |(p₁ : ℚ) * q₂ - (p₂ : ℚ) * q₁| := by
      rw [Int.cast_abs]
      push_cast
      ring_nf
    exact_mod_cast h
  -- |p₁/q₁ − p₂/q₂| = |p₁q₂ − p₂q₁|/(q₁q₂)
  have hdiff : |(p₁ : ℚ) / q₁ - (p₂ : ℚ) / q₂|
      = |(p₁ : ℚ) * q₂ - (p₂ : ℚ) * q₁| / (q₁ * q₂) := by
    rw [div_sub_div _ _ (ne_of_gt hq₁Q) (ne_of_gt hq₂Q), abs_div,
      abs_of_pos (mul_pos hq₁Q hq₂Q)]
    have hnum : (p₁ : ℚ) * q₂ - q₁ * p₂ = (p₁ : ℚ) * q₂ - p₂ * q₁ := by ring
    rw [hnum]
  rw [hdiff] at hw
  have hpos : (0 : ℚ) < (q₁ : ℚ) * q₂ := by positivity
  have h1 : |(p₁ : ℚ) * q₂ - (p₂ : ℚ) * q₁| ≤ w * (q₁ * q₂) := by
    rw [div_le_iff₀ hpos] at hw
    linarith
  calc (1 : ℚ) ≤ |(p₁ : ℚ) * q₂ - (p₂ : ℚ) * q₁| := habs1Q
    _ ≤ w * (q₁ * q₂) := h1
    _ = w * q₁ * q₂ := by ring

/-- **One carrier per component on super-ladders**: if two moduli `q₁ ≠ q₂ ≤ D`
carry deep multipliers whose rational instants lie in a common interval of width
`1/(7·V)` with `D² < 7·V`, they coincide — the component pins a unique carrier. -/
theorem carrier_unique_per_component (p₁ q₁ p₂ q₂ D V : ℕ)
    (hq₁ : 0 < q₁) (hq₂ : 0 < q₂) (hD₁ : q₁ ≤ D) (hD₂ : q₂ ≤ D)
    (hV : 0 < V) (hsuper : (D : ℚ) * D < 7 * V)
    (hwidth : |(p₁ : ℚ) / q₁ - (p₂ : ℚ) / q₂| ≤ 1 / (7 * V))
    (hne : (p₁ : ℚ) / q₁ ≠ (p₂ : ℚ) / q₂) : False := by
  have h := farey_separation p₁ q₁ p₂ q₂ (1 / (7 * V)) hq₁ hq₂ hne hwidth
  have hVQ : (0 : ℚ) < (V : ℚ) := by exact_mod_cast hV
  have hq₁D : (q₁ : ℚ) ≤ (D : ℚ) := by exact_mod_cast hD₁
  have hq₂D : (q₂ : ℚ) ≤ (D : ℚ) := by exact_mod_cast hD₂
  have hq₁0 : (0 : ℚ) < (q₁ : ℚ) := by exact_mod_cast hq₁
  have hq₂0 : (0 : ℚ) < (q₂ : ℚ) := by exact_mod_cast hq₂
  have hbound : (1 : ℚ) / (7 * V) * q₁ * q₂ ≤ (D : ℚ) * D / (7 * V) := by
    have h1 : (q₁ : ℚ) * q₂ ≤ (D : ℚ) * D := by nlinarith
    rw [div_mul_eq_mul_div, div_mul_eq_mul_div, one_mul]
    apply div_le_div_of_nonneg_right _ (by positivity)
    · exact h1
  have hlt : (D : ℚ) * D / (7 * V) < 1 := by
    rw [div_lt_one (by positivity)]
    exact hsuper
  linarith [h, hbound, hlt]

/-- **THE ADAPTIVE-Q PIGEONHOLE** (the skeleton, with the component count as an
honest named parameter): if the deep-carrying moduli of a window inject into `K`
components and the window holds more than `K` moduli, some modulus is deep-free.
At that modulus THM-950's census criterion needs only `liveCount > 0`. -/
theorem adaptive_q_pigeonhole (v : Fin 13 → ℤ) (window : Finset ℕ) (K : ℕ)
    (comp : ℕ → ℕ)
    (hcarriers : ∀ q ∈ window,
      (0 < (((Finset.Ioo 0 q).filter fun p => 6 ≤ bandCount v q p).card)) →
      comp q < K)
    (hcomp_inj : ∀ q₁ ∈ window, ∀ q₂ ∈ window,
      (0 < (((Finset.Ioo 0 q₁).filter fun p => 6 ≤ bandCount v q₁ p).card)) →
      (0 < (((Finset.Ioo 0 q₂).filter fun p => 6 ≤ bandCount v q₂ p).card)) →
      comp q₁ = comp q₂ → q₁ = q₂)
    (hsize : K < window.card) :
    ∃ q ∈ window,
      (((Finset.Ioo 0 q).filter fun p => 6 ≤ bandCount v q p).card) = 0 := by
  by_contra hcon
  push Not at hcon
  -- every window modulus is a carrier; comp injects them into Fin K — too small
  have hall : ∀ q ∈ window,
      0 < (((Finset.Ioo 0 q).filter fun p => 6 ≤ bandCount v q p).card) := by
    intro q hq
    have := hcon q hq
    omega
  have hmap : ∀ q ∈ window, comp q ∈ Finset.range K := by
    intro q hq
    rw [Finset.mem_range]
    exact hcarriers q hq (hall q hq)
  have hinj : ∀ q₁ ∈ window, ∀ q₂ ∈ window, comp q₁ = comp q₂ → q₁ = q₂ := by
    intro q₁ h₁ q₂ h₂ hEq
    exact hcomp_inj q₁ h₁ q₂ h₂ (hall q₁ h₁) (hall q₂ h₂) hEq
  have hle := Finset.card_le_card_of_injOn comp hmap hinj
  rw [Finset.card_range] at hle
  omega

/-! ## Axiom audit -/
#print axioms farey_separation
#print axioms carrier_unique_per_component
#print axioms adaptive_q_pigeonhole

end LRC14Concrete
end LonelyRunner
