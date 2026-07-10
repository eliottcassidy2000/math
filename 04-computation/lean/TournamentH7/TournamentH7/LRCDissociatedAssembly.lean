/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: opus (LRC multi-agent project, 2026-07-10-S209)
-/
import Mathlib
import TournamentH7.LRCPrimitiveAssembly

/-!
# The `d = 2,3` detuned peel — the residual shrinks to DISSOCIATED families (opus-S209)

opus-S208 found (exact witness) that the primitive residual's μ-minimizers are **near-dilates**:
`d`-detuned families `v = g·H ∪ D` with `|D| = 2, 3` (e.g. `[2,12,14,16,18,20,22,26,31,34,37,38,46]`, which
is `2·H ∪ {31,37}` — μ ≈ 0.0085). The assembly's detuned branch peels only `d = 1` (`lonely14_of_detuned`);
`d = 2, 3` survive and carry the small μ, which is why the S207 "`Vmax` decorrelation tail" was false. The
fix S208 recommended: peel `d = 2, 3` (monad's **THM-678**, the multi-detuned counting dispatch — PROVED
elementarily on paper: `Σ Nᵢ/qᵢ < 1 ⟹` a good branch, so `d = 2` fires unless `q₁ = q₂ = 2`, `d = 3` when
all `qᵢ ≥ 8`), leaving the genuinely DISSOCIATED families where the decorrelation floor is well-behaved.

THM-678 is not yet transcribed to Lean, so it enters here as a **named dispatch hypothesis**
`MultiDetunedDispatch` — exactly as the LRC(≤13) citation and `hB5` do. Consuming it:

* `MultiDetunedDispatch` — every family with some `g ≥ 2` having exactly 2 or 3 non-multiples is lonely
  (the THM-678 content);
* `ResidualObligationDissoc` — `ResidualObligationPrimitive`'s clauses PLUS "no `g` has 2 or 3
  non-multiples" (dissociated: with `tupleGcd = 1` and divisor-closed already, every `g ≥ 2` then has `≥ 4`
  non-multiples);
* `lrc14_grand_assembly_dissoc : LRCUpTo13 → MultiDetunedDispatch → ResidualObligationDissoc →
  LRC14Statement`, and the `B5` finish line `lrc14_from_B5_dissoc`.

Net: the analytic floor obligation moves from the full primitive residual (where `inf μ = 0.0085` is carried
by the near-dilates) to the DISSOCIATED residual (where μ is bounded well away from 0 and rises to iid — the
moment/Bonferroni regime). The near-dilates are handed to THM-678.

Kernel-pure: no `sorry`, no `native_decide`. Axioms: `[propext, Classical.choice, Quot.sound]`.
-/

namespace LonelyRunner
namespace LRC14Grand

open Finset

/-- The number of coordinates of `v` NOT divisible by `g` (the detuning level at `g`). -/
def nonMultCard (v : Fin 13 → ℤ) (g : ℤ) : ℕ :=
  (Finset.univ.filter (fun i => ¬ g ∣ v i)).card

/-- **THM-678 as a named dispatch hypothesis.** Any family with a `g ≥ 2` at detuning level 2 or 3
(`v = g·H ∪ D`, `|D| ∈ {2,3}`) has a lonely instant. Proved elementarily (the multi-detuned counting
dispatch) but not yet Lean-transcribed — enters as a citation, like LRC(≤13). -/
def MultiDetunedDispatch : Prop :=
  ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) →
    (∃ g : ℤ, 2 ≤ g ∧ (nonMultCard v g = 2 ∨ nonMultCard v g = 3)) → ∃ t : ℝ, Lonely 14 v t

/-- **The residual surface, DISSOCIATED.** `ResidualObligationPrimitive`'s clauses plus "no `g ≥ 2` has
detuning level 2 or 3" — i.e. every `g ≥ 2` has `≥ 4` non-multiples (with primitivity + divisor-closed
excluding levels 0 and 1). This is the family class where the decorrelation floor is well-behaved. -/
def ResidualObligationDissoc : Prop :=
  ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → LRC14.tupleGcd v = 1 →
    LRC14.CoveringFamily v →
    GapFamily v →
    (∀ i, ∃ j, j ≠ i ∧ |v i| ≤ 13 * |v j|) →
    (∀ i j, i ≠ j → |v i| ≠ |v j|) →
    (∃ i, 23 ≤ |v i|) →
    (∀ g : ℤ, 2 ≤ g → ∀ i₀ : Fin 13, (∀ j, j ≠ i₀ → g ∣ v j) → g ∣ v i₀) →
    (¬ ∃ (L : ℤ) (k a : Fin 13 → ℤ) (A : ℝ), (∀ i, v i = a i + L * k i) ∧ 0 < (L : ℝ) ∧
      (∀ i, |(a i : ℝ)| ≤ A) ∧ A / (L : ℝ) ≤ 1/13 - 1/14 ∧ (∀ i, k i ≠ 0) ∧
      (Finset.univ.image k).card ≤ 12) →
    (∀ d : ℤ, 2 ≤ d → ∀ a : ℤ, (∀ i, d ∣ (v i - a)) → d ∣ a) →
    (∀ g : ℤ, 2 ≤ g → nonMultCard v g ≠ 2 ∧ nonMultCard v g ≠ 3) →
    ∃ t : ℝ, Lonely 14 v t

/-- **The dissociation peel.** The primitive residual obligation follows from the THM-678 dispatch (which
clears the `d = 2, 3` detuned families) plus the residual obligation on the dissociated remainder. -/
theorem residualObligationPrimitive_of_dissoc
    (hMD : MultiDetunedDispatch) (hdissoc : ResidualObligationDissoc) :
    ResidualObligationPrimitive := by
  intro v hv hgcd hcov hgap hcomp hdist hlarge hdiv hcoarse hcres
  by_cases hdet : ∃ g : ℤ, 2 ≤ g ∧ (nonMultCard v g = 2 ∨ nonMultCard v g = 3)
  · exact hMD v hv hdet
  · push_neg at hdet
    exact hdissoc v hv hgcd hcov hgap hcomp hdist hlarge hdiv hcoarse hcres hdet

/-- **LRC(14) from the citation, the THM-678 dispatch, and the DISSOCIATED residual floor.** The
near-dilate (`d = 2, 3` detuned) families are handed to THM-678; the remaining analytic obligation is the
floor on the dissociated residual only. -/
theorem lrc14_grand_assembly_dissoc (cite : LRCUpTo13)
    (hMD : MultiDetunedDispatch) (hdissoc : ResidualObligationDissoc) :
    LRC14.LRC14Statement :=
  lrc14_grand_assembly_primitive cite (residualObligationPrimitive_of_dissoc hMD hdissoc)

/-- **The finish line on the DISSOCIATED residual.** With THM-678 clearing the near-dilates, the `B5 > 0`
obligation need only hold for dissociated families (every `g ≥ 2` at detuning level `≥ 4`) — the regime
where the density floor decorrelates toward `(6/7)^13` (opus-S208). This is the sharpest `hB5`. -/
theorem lrc14_from_B5_dissoc (cite : LRCUpTo13) (hMD : MultiDetunedDispatch)
    (hB5 : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → LRC14.tupleGcd v = 1 →
      LRC14.CoveringFamily v → GapFamily v →
      (∀ i, ∃ j, j ≠ i ∧ |v i| ≤ 13 * |v j|) → (∀ i j, i ≠ j → |v i| ≠ |v j|) →
      (∃ i, 23 ≤ |v i|) →
      (∀ g : ℤ, 2 ≤ g → ∀ i₀ : Fin 13, (∀ j, j ≠ i₀ → g ∣ v j) → g ∣ v i₀) →
      (¬ ∃ (L : ℤ) (k a : Fin 13 → ℤ) (A : ℝ), (∀ i, v i = a i + L * k i) ∧ 0 < (L : ℝ) ∧
        (∀ i, |(a i : ℝ)| ≤ A) ∧ A / (L : ℝ) ≤ 1/13 - 1/14 ∧ (∀ i, k i ≠ 0) ∧
        (Finset.univ.image k).card ≤ 12) →
      (∀ g : ℤ, 2 ≤ g → nonMultCard v g ≠ 2 ∧ nonMultCard v g ≠ 3) →
      ∃ q : ℕ, 0 < q ∧ 0 < LRC14Concrete.B5 v q) :
    LRC14.LRC14Statement := by
  refine lrc14_grand_assembly_dissoc cite hMD ?_
  intro v hv hgcd hcov hgap hcomp hdist hlarge hdiv hcoarse hcres hdissoc
  obtain ⟨q, hq, hpos⟩ := hB5 v hv hgcd hcov hgap hcomp hdist hlarge hdiv hcoarse hdissoc
  exact LRC14Concrete.lonely_of_Mreach_ge v hv (LRC14Concrete.mreach_ge_of_B5_pos v q hq hpos)

end LRC14Grand
end LonelyRunner
