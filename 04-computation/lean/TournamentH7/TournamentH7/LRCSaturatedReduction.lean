/-
  TournamentH7.LRCSaturatedReduction — the SIEVE reduction of LRC(14) to the saturated core.
  (kind-pasteur-2026-07-07-S56, HYP-4737; formalizes leg 1 of the S55 decomposition.)

  opus-S131's sieve dichotomy, as a clean Lean node.  For any `q ∈ {2,…,14}`, if no speed is
  divisible by `q`, then `t = 1/q` is `14`-lonely (`‖v/q‖ ≥ 1/q ≥ 1/14`).  Contrapositive
  (`counterexample_needs_all_divisors`): a counterexample to LRC(14) must be **saturated** — a
  multiple of *every* `q ∈ {2,…,14}`.  Hence LRC(14) reduces to: every saturated 13-family is
  lonely.  This is the non-analytic half of the decomposition
  `LRC(14) = sieve[GREEN] ⊕ coarse-reduction[GREEN] ⊕ decorrelation[open, saturated core]`.

  Kernel-pure (propext, Classical.choice, Quot.sound); no sorry, no native_decide.
-/
import Mathlib
import TournamentH7.LonelyRunner
import TournamentH7.LRCFourteenSkeleton

namespace LonelyRunner
namespace LRC14

/-- A 13-speed family is **saturated** if some speed is divisible by every `q ∈ {2,…,14}`.
By the sieve, an LRC(14) counterexample must be saturated. -/
def Saturated (v : Fin 13 → ℤ) : Prop :=
  ∀ q : ℕ, 2 ≤ q → q ≤ 14 → ∃ i, (q : ℤ) ∣ v i

/-- **THE SIEVE REDUCTION.**  If every *saturated* nonzero 13-family is lonely, then LRC(14)
holds — because a counterexample is forced to be saturated by
`counterexample_needs_all_divisors`.  So the non-saturated families are discharged (sieve,
GREEN), and the whole content of LRC(14) is the **saturated core**. -/
theorem lrc14_of_saturated_lonely
    (hsat : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → Saturated v → ∃ t : ℝ, Lonely 14 v t) :
    LRC14Statement := by
  intro v hv
  by_contra h
  push_neg at h
  have hsatv : Saturated v := fun q hq2 hq14 =>
    counterexample_needs_all_divisors 14 v h q hq2 hq14
  obtain ⟨t, ht⟩ := hsat v hv hsatv
  exact h t ht

/-- The converse is immediate, so the reduction is an equivalence: LRC(14) holds **iff** every
saturated nonzero 13-family is lonely. -/
theorem lrc14_iff_saturated_lonely :
    LRC14Statement ↔
      (∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → Saturated v → ∃ t : ℝ, Lonely 14 v t) := by
  constructor
  · intro hlrc v hv _; exact hlrc v hv
  · exact lrc14_of_saturated_lonely

#print axioms lrc14_of_saturated_lonely
#print axioms lrc14_iff_saturated_lonely

end LRC14
end LonelyRunner
