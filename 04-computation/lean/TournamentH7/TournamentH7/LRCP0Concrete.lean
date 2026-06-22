/-
  TournamentH7.LRCP0Concrete -- the CONCRETE cover atom `p0(E) = meas(coverSet E)`
  and the concrete wide bound `hp0cap` (kind-pasteur-2026-06-22-S31, HYP-2839).

  The skeleton (`LRCFourteenSkeleton`) carries `p0 : List ℤ → ℝ` opaquely.  This
  file gives it a concrete Lebesgue-measure definition

      p0 E := (slowμ (coverSet E)).toReal              -- slowμ = volume on one period

  so that `hp0cap` (`p0 E ≤ cap`) becomes a genuine statement about the measurable
  cover event `coverSet E` (the all-6-inner-sectors-hit set = `measS7`).  The
  elementary cores from `LRCCoverBound` then apply directly:
  `p0` is a probability in `[0,1]`, vanishes for `< 6` distinct speeds, is monotone
  in the speed list, and the wide bound reduces to the resonance bound
  `p0 ≤ p0_decorr` (the analytic residual).  Sorry-free.

  To wire the skeleton, set `LRC14.p0 E := DenseCovers.p0 E` (a one-line change);
  every `p0`-node then refers to this concrete atom.
-/

import TournamentH7.LRCCoverBound

namespace LonelyRunner
namespace DenseCovers

open MeasureTheory

/-- **The concrete cover atom** `p0(E) = meas{x : the phases hit all 6 inner
sectors}`, a real number in `[0,1]`. -/
noncomputable def p0 (E : List ℤ) : ℝ := (slowμ (coverSet E)).toReal

theorem p0_nonneg (E : List ℤ) : 0 ≤ p0 E := ENNReal.toReal_nonneg

/-- `p0` is a probability: `p0 E ≤ 1`. -/
theorem p0_le_one (E : List ℤ) : p0 E ≤ 1 := by
  have h : slowμ (coverSet E) ≤ 1 := by
    calc slowμ (coverSet E) ≤ slowμ Set.univ := measure_mono (Set.subset_univ _)
      _ = 1 := measure_univ
  calc p0 E = (slowμ (coverSet E)).toReal := rfl
    _ ≤ (1 : ENNReal).toReal := ENNReal.toReal_mono ENNReal.one_ne_top h
    _ = 1 := ENNReal.toReal_one

/-- **`p0` vanishes for `< 6` distinct speeds** (the small-k pigeonhole): six
disjoint inner sectors cannot be covered by fewer than six speeds. -/
theorem p0_eq_zero_of_card_lt_six {E : List ℤ} (h : E.toFinset.card < 6) : p0 E = 0 := by
  rw [p0, slowμ_coverSet_eq_zero_of_card_lt_six h, ENNReal.toReal_zero]

/-- **`p0` is monotone** in the speed list. -/
theorem p0_mono {E E' : List ℤ} (h : ∀ e ∈ E, e ∈ E') : p0 E ≤ p0 E' :=
  ENNReal.toReal_mono (measure_ne_top _ _) (slowμ_coverSet_mono h)

/-- **hp0cap for `< 6` speeds is unconditional:** `p0 E = 0 ≤ cap` for any `cap ≥ 0`. -/
theorem p0_le_cap_of_card_lt_six {E : List ℤ} {cap : ℝ} (hcap : 0 ≤ cap)
    (h : E.toFinset.card < 6) : p0 E ≤ cap := by
  rw [p0_eq_zero_of_card_lt_six h]; exact hcap

/-- **The concrete wide bound `hp0cap`**, parameterized by the cap.  For `< 6`
speeds it is unconditional; for `≥ 6` it is the resonance residual
`p0 ≤ p0_decorr ≤ Q < cap`.  This packages hp0cap as a concrete statement about
`p0 = meas(coverSet)`, ready to discharge the skeleton's `WideBound`. -/
def WideBoundConcrete (cap : ℝ) (E : List ℤ) : Prop := p0 E ≤ cap

/-- hp0cap holds once the resonance/finite/rational chain is supplied for the
`≥ 6` case (the `< 6` case is the pigeonhole, free). -/
theorem wideBoundConcrete_of_decorrelation {cap : ℝ} (E : List ℤ) (hcap : 0 ≤ cap)
    (hbig : 6 ≤ E.toFinset.card → ∃ p0decorr Q : ℝ,
      p0 E ≤ p0decorr ∧ p0decorr ≤ Q ∧ Q ≤ cap) :
    WideBoundConcrete cap E := by
  unfold WideBoundConcrete
  rcases lt_or_ge E.toFinset.card 6 with hsmall | hge
  · exact p0_le_cap_of_card_lt_six hcap hsmall
  · obtain ⟨p0decorr, Q, hr, hf, hm⟩ := hbig hge
    linarith

/-! ## Axiom audit -/

#print axioms p0
#print axioms p0_le_one
#print axioms p0_eq_zero_of_card_lt_six
#print axioms wideBoundConcrete_of_decorrelation

end DenseCovers
end LonelyRunner
