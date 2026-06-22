/-
  TournamentH7.LRCApex7Floor -- the apex-7 obstruction fragment (Idea 3, the owner's
  inspiration; mac-mini S34 HYP-2876 proved it analytically), kind-pasteur S31d.

  THE cleanest fragment of why `14 = 2·7` composite makes LRC(14) hard: in the
  rational-witness (residue) route a covering set is certified at denominator `D`
  by a residue `a` with `‖s·a/D‖ ≥ 1/14` for every speed `s`.  At `D = 14`, any
  runner whose speed is a multiple of `14` sits EXACTLY on the observer for every
  `a`:

      14 ∣ s  ⟹  ‖s · (a/14)‖ = 0   (the runner is never lonely),

  because `s·a/14 = (s/14)·a ∈ ℤ`.  So `D = 14` never certifies a covering set
  (which forces such a runner) — the minimal witness denominator is `≥ 15`.  This
  is precisely the apex prime `7` (`= 14/2`) acting: the observer absorbs the
  composite scale.  Sorry-free, elementary `Int.fract`.
-/

import TournamentH7.LRCMreachConcrete

namespace LonelyRunner
namespace Apex7

open LRC14Concrete

/-- **The apex-7 collapse.**  A multiple of `14` lands exactly on the observer at
denominator `14`: `‖s · (a/14)‖ = 0` for every residue `a`. -/
theorem nearInt_eq_zero_of_dvd (s : ℤ) (hs : (14 : ℤ) ∣ s) (a : ℤ) :
    nearInt ((s : ℝ) * ((a : ℝ) / 14)) = 0 := by
  obtain ⟨m, hm⟩ := hs
  have heq : (s : ℝ) * ((a : ℝ) / 14) = ((m * a : ℤ) : ℝ) := by
    rw [hm]; push_cast; ring
  rw [heq, nearInt, Int.fract_intCast]
  norm_num

/-- The multiple-of-14 runner is strictly inside the danger radius `1/14`: it can
never witness loneliness at `D = 14`. -/
theorem on_observer_lt_of_dvd (s : ℤ) (hs : (14 : ℤ) ∣ s) (a : ℤ) :
    nearInt ((s : ℝ) * ((a : ℝ) / 14)) < 1 / 14 := by
  rw [nearInt_eq_zero_of_dvd s hs a]; norm_num

/-- **D = 14 never certifies.**  If the speed set `S` contains a multiple of `14`,
then for EVERY residue `a` some runner is within `1/14` of the observer at time
`a/14` — so no `a/14` is `14`-lonely.  Hence the minimal certifying denominator is
`≥ 15` (the apex-7 floor). -/
theorem D14_never_certifies (S : Finset ℤ) (s : ℤ) (hsS : s ∈ S) (hs : (14 : ℤ) ∣ s)
    (a : ℤ) : ∃ t ∈ S, nearInt ((t : ℝ) * ((a : ℝ) / 14)) < 1 / 14 :=
  ⟨s, hsS, on_observer_lt_of_dvd s hs a⟩

/-! ## Axiom audit -/

#print axioms nearInt_eq_zero_of_dvd
#print axioms D14_never_certifies

end Apex7
end LonelyRunner
