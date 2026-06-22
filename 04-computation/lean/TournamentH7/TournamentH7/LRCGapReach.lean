/-
  TournamentH7.LRCGapReach -- the geometric REACH CORE of THM-527 Part A
  (kind-pasteur-2026-06-22-S31): the elementary "1/7 gap ⟹ 1/14 margin" fact.

  This is the heart of "criterion C": why a `> 1/7` gap in the cluster teeth
  produces an actual `≥ 1/14` loneliness margin.  In the slow-fast picture the
  parked runners create "danger teeth" of radius `1/14` centred at the cluster
  phases `c = frac(e·x)`; a fast-phase window free of all teeth has its midpoint
  `≥ 1/14` from every tooth modulo 1.

  Phrasing the gap as "no integer-translate of any tooth lands in the open
  interval `(a, a+g)`" makes this a one-line interval fact -- no cyclic-distance
  or wrap-around case analysis.  Sorry-free.  (The remaining Part-A content is the
  `Vmax`-ruler embedding of `φ` into a real time and the equidistribution
  `ρ_K → ρ*`; this lemma is the geometric core they both feed.)
-/

import Mathlib.Tactic

namespace LonelyRunner
namespace GapReach

/-- **The reach core.**  If the open interval `(a, a+g)` contains no integer
translate of any tooth `c ∈ C`, then the midpoint `a + g/2` is at distance `≥ g/2`
from every tooth modulo `1`: `∀ c ∈ C, ∀ n : ℤ, g/2 ≤ |(a + g/2) - (c + n)|`.
(Equivalently `‖(a+g/2) - c‖_{ℝ/ℤ} ≥ g/2` for every tooth.) -/
theorem margin_ge_of_free_interval (C : Finset ℝ) (a g : ℝ)
    (hfree : ∀ c ∈ C, ∀ n : ℤ, (c + (n : ℝ)) ∉ Set.Ioo a (a + g)) :
    ∀ c ∈ C, ∀ n : ℤ, g / 2 ≤ |(a + g / 2) - (c + (n : ℝ))| := by
  intro c hc n
  by_contra h
  push_neg at h
  rw [abs_lt] at h
  obtain ⟨h1, h2⟩ := h
  exact hfree c hc n ⟨by linarith, by linarith⟩

/-- **"1/7 ⟹ 1/14".**  A free gap wider than `1/7` puts the midpoint strictly more
than `1/14` from every tooth modulo `1` — the loneliness margin certifying
`M(S) ≥ 1/14` once the fast phase `φ = a + g/2` is realized at a `Vmax`-ruler time.
This is the exact reach inequality of THM-527 Part A's criterion C. -/
theorem margin_gt_one_div_14_of_gap (C : Finset ℝ) (a g : ℝ) (hg : 1 / 7 < g)
    (hfree : ∀ c ∈ C, ∀ n : ℤ, (c + (n : ℝ)) ∉ Set.Ioo a (a + g)) :
    ∀ c ∈ C, ∀ n : ℤ, 1 / 14 < |(a + g / 2) - (c + (n : ℝ))| := by
  intro c hc n
  have hmid := margin_ge_of_free_interval C a g hfree c hc n
  linarith

/-- The midpoint of a `> 1/7` free gap exists and witnesses a uniform `> 1/14`
loneliness margin: a single real `φ` (`= a + g/2`) works for all teeth. -/
theorem exists_lonely_phase_of_gap (C : Finset ℝ) (a g : ℝ) (hg : 1 / 7 < g)
    (hfree : ∀ c ∈ C, ∀ n : ℤ, (c + (n : ℝ)) ∉ Set.Ioo a (a + g)) :
    ∃ φ : ℝ, ∀ c ∈ C, ∀ n : ℤ, 1 / 14 < |φ - (c + (n : ℝ))| :=
  ⟨a + g / 2, margin_gt_one_div_14_of_gap C a g hg hfree⟩

/-! ## Axiom audit -/

#print axioms margin_ge_of_free_interval
#print axioms margin_gt_one_div_14_of_gap
#print axioms exists_lonely_phase_of_gap

end GapReach
end LonelyRunner
