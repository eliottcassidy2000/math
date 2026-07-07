/-
  TournamentH7.LRCMod25Floor — THE mod-25 2/25-FLOOR CERTIFICATE
  (kind-pasteur-2026-07-06-S41, HYP-4567).

  Pins the lower-bound half of the LRC(14) crux.  After opus-S120's architecture
  audit, (G) [first gap (1/13, 2/25) empty at N=12] reduces to the FREIMAN step:
  a gap member is an (N−2)-AP + exactly 2 defects; the ≤2-defect world is swept
  empty, so the residual is "rule out ≥3 defects".  Since LRC(13) is proved (owner
  directive), every 12-speed family has M ≥ 1/13, so the Freiman step is EQUIVALENT
  to the lower bound

        longest-AP-subset(v) ≤ 9  (≥ 3 defects)   ⟹   M(v) ≥ 2/25.

  kps-S41 reduces the SHARP core of that lower bound to a mod-25 covering fact:
  a ≥3-defect family with no multiple of 25 admits a rotation `c ∈ (ℤ/25)*` taking
  every speed off the forbidden band `{0, ±1}` mod 25 (verified: ~all such families
  clear this way; the exceptions contain a multiple of 25 and clear at small
  denominators with M ≫ 2/25).  This file gives the Lean certificate that such a
  rotation forces the 2/25 floor — a direct instance of the `rational_point_margin`
  atom at `s = 25`, `μ = 2`, so the loose conclusion `∃ t, ∀ i m, 2/25 ≤ |vᵢ t − m|`
  is machine-checkable from a decidable integer hypothesis.

  Kernel-pure; no `sorry`, no `native_decide`.
-/
import Mathlib
import TournamentH7.LRCHarmonicGate

namespace LonelyRunner
namespace Mod25Floor

open HarmonicGate

/-- **THE mod-25 2/25-FLOOR CERTIFICATE (the Freiman-residual lower bound).**
If some `c` rotates every speed off the forbidden band `{0, ±1}` mod 25 — i.e.
`2 ≤ (v i · c) % 25 ≤ 23` for every runner — then at `t = c/25` every runner is at
distance `≥ 2/25` from every integer, so the family is LOOSE (`M ≥ 2/25`, above the
gap).  Direct instance of `rational_point_margin` (`s = 25`, `μ = 2`). -/
theorem mod25_covering_floor {ι : Type*} (v : ι → ℤ) (c : ℤ)
    (hcov : ∀ i, 2 ≤ (v i * c) % 25 ∧ (v i * c) % 25 ≤ 23) :
    ∀ i, ∀ m : ℤ, (2 : ℝ) / 25 ≤ |(v i : ℝ) * ((c : ℝ) / 25) - m| :=
  rational_point_margin v c 25 2 (by norm_num) hcov

/-- **Existential loose form.**  A mod-25 rotation `c` off `{0, ±1}` gives a common
`2/25`-margin point `t = c/25` — exactly the loose branch of the dichotomy. -/
theorem loose_of_mod25_covering {ι : Type*} (v : ι → ℤ) (c : ℤ)
    (hcov : ∀ i, 2 ≤ (v i * c) % 25 ∧ (v i * c) % 25 ≤ 23) :
    ∃ t : ℝ, ∀ i, ∀ m : ℤ, (2 : ℝ) / 25 ≤ |(v i : ℝ) * t - m| :=
  ⟨(c : ℝ) / 25, mod25_covering_floor v c hcov⟩

#print axioms mod25_covering_floor
#print axioms loose_of_mod25_covering

end Mod25Floor
end LonelyRunner
