/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: opus (LRC multi-agent project, 2026-07-10-S206)
-/
import Mathlib
import TournamentH7.LRC14GrandAssembly
import TournamentH7.LRCPrimitivePeel
import TournamentH7.LRCResidualFromLedger

/-!
# The grand assembly with the PRIMITIVITY PEEL (opus-S206)

`lrc14_grand_assembly` derives `LRC14Statement` from the LRC(≤13) citation plus `ResidualObligation`.
opus-S205 showed that obligation's domain **admits non-primitive dilates**: the machine-checked witness
`v = 2·[1,2,3,4,5,6,7,8,9,11,12,13,20]` satisfies every residual clause with `gcd = 2` and `μ = 1/980`,
while its core has `Vmax = 20 ≤ 22` (already window-censused). Since `α ↦ c·α` is measure-preserving,
`μ(c·w) = μ(w)` exactly, and cores range over the window census (which contains near-APs with `μ → 0`).
Hence `inf μ = 0` over the residual as stated — **no uniform measure floor exists**, and every proof
attempt routing through one is unreachable.

This file removes the dilates. Because `c·w` is lonely at `t/c` whenever `w` is lonely at `t`
(`lonely_scale`, packaged as `PrimitivePeel.lrc14_of_primitive`), LRC(14) reduces to PRIMITIVE families,
and the residual may be restated with `tupleGcd v = 1`:

* `ResidualObligationPrimitive` — `ResidualObligation`'s clauses **plus** `tupleGcd v = 1`;
* `lrc14_grand_assembly_primitive : LRCUpTo13 → ResidualObligationPrimitive → LRC14Statement`;
* `residualPrimitive_of_residual : ResidualObligation → ResidualObligationPrimitive` (so the primitive
  surface is **strictly weaker**: it subsumes the old one, never the reverse).

After the peel the residual is dilate-free and the measure-floor target is well-posed (adversarial descent
over the full primitive residual predicate gives `min μ ≈ 0.0094`, vs iid `(6/7)^13 ≈ 0.1348`).

The branch cascade is identical to `lrc14_grand_assembly`'s (non-covering ∣ window-≤22 ∣ ratio-≤13 ∣
dominant peel ∣ repeated ∣ detuned ∣ coarse ∣ common-residue); only the residual leg differs, now carrying
primitivity.

Kernel-pure: no `sorry`, no `native_decide`. Axioms: `[propext, Classical.choice, Quot.sound]`.
-/

namespace LonelyRunner
namespace LRC14Grand

open LonelyRunner

/-- **The residual surface, dilate-free.** `ResidualObligation`'s clauses plus primitivity
(`tupleGcd v = 1`). Strictly weaker than `ResidualObligation` — the dilates are gone. -/
def ResidualObligationPrimitive : Prop :=
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
    ∃ t : ℝ, Lonely 14 v t

/-- The primitive surface is weaker: the old obligation implies it (just forget `tupleGcd v = 1`). -/
theorem residualPrimitive_of_residual (h : ResidualObligation) : ResidualObligationPrimitive :=
  fun v hv _ hcov hgap hcomp hdist hlarge hdiv hcoarse hcres =>
    h v hv hcov hgap hcomp hdist hlarge hdiv hcoarse hcres

/-- **THE GRAND ASSEMBLY WITH THE PRIMITIVITY PEEL.** `LRC14Statement` from the LRC(≤13) citation and
the DILATE-FREE residual obligation. -/
theorem lrc14_grand_assembly_primitive (cite : LRCUpTo13)
    (hresP : ResidualObligationPrimitive) : LRC14.LRC14Statement := by
  apply LRC14.PrimitivePeel.lrc14_of_primitive
  intro v hv hgcd
  -- (1) non-covering: the denominator sieve at the missing modulus
  by_cases hcov : LRC14.CoveringFamily v
  swap
  · rw [LRC14.CoveringFamily] at hcov
    push_neg at hcov
    obtain ⟨q, h2, h14, hdiv⟩ := hcov
    exact ⟨(1 : ℝ) / q, sieve_one_div 14 q v h14 (by omega) hdiv⟩
  -- (4) bounded window ≤ 22 — kernel-pure six-witness pigeonhole (LEM-024)
  by_cases hwin : ∀ i, |v i| ≤ 22
  · exact WindowData.hwindowW_closed 22 cite Window22Census.hdistinct22_kernel v hv hwin
  -- (2) all-comparable: ratio ≤ 13 throughout
  by_cases hgap : GapFamily v
  swap
  · rw [GapFamily, not_not] at hgap
    obtain ⟨imax, -, hmax⟩ :=
      Finset.exists_max_image (Finset.univ) (fun i => |v i|) ⟨0, Finset.mem_univ 0⟩
    obtain ⟨jmin, -, hmin⟩ :=
      Finset.exists_min_image (Finset.univ) (fun i => |v i|) ⟨0, Finset.mem_univ 0⟩
    exact ⟨_, LRC14.spread13_lonely v (|v jmin|) (|v imax|)
      (abs_pos.mpr (hv jmin))
      (fun i => hmin i (Finset.mem_univ i))
      (fun i => hmax i (Finset.mem_univ i))
      (hgap imax jmin)⟩
  -- (3) dominant runner
  by_cases hdom : ∃ i, ∀ j, j ≠ i → 13 * |v j| < |v i|
  · exact LRC14.hdom_discharged cite v hv hcov hdom
  -- (5) repeated absolute speed
  by_cases hrep : ∃ i j, i ≠ j ∧ |v i| = |v j|
  · obtain ⟨i, j, hij, heq⟩ := hrep
    obtain ⟨t, ht⟩ := lonely14_of_repeat cite (fun i => |v i|)
      (fun i => by simpa using hv i) hij heq
    exact ⟨t, lonely_of_abs v t ht⟩
  -- (6) detuned harmonic
  by_cases hdet : ∃ g : ℤ, 2 ≤ g ∧ ∃ i₀ : Fin 13, (∀ j, j ≠ i₀ → g ∣ v j) ∧ ¬ g ∣ v i₀
  · obtain ⟨g, hg2, i₀, hH, hnd⟩ := hdet
    exact DetunedDispatch.lonely14_of_detuned cite v hv g hg2 i₀ hH hnd
  -- (7) multi-scale: the coarse reduction
  by_cases hms : ∃ (L : ℤ) (k a : Fin 13 → ℤ) (A : ℝ), (∀ i, v i = a i + L * k i) ∧
      0 < (L : ℝ) ∧ (∀ i, |(a i : ℝ)| ≤ A) ∧ A / (L : ℝ) ≤ 1/13 - 1/14 ∧
      (∀ i, k i ≠ 0) ∧ (Finset.univ.image k).card ≤ 12
  · obtain ⟨L, k, a, A, hdecomp, hL, ha, hbudget, hkne, hcard⟩ := hms
    exact CoarseReduction.lonely14_of_coarse_le12 cite v k a L A hdecomp hL ha hbudget hkne hcard
  -- (8) nontrivial common residue
  by_cases hcr : ∃ d : ℤ, 2 ≤ d ∧ ∃ a : ℤ, ¬ d ∣ a ∧ ∀ i, d ∣ (v i - a)
  · obtain ⟨d, hd2, a, hna, hcres⟩ := hcr
    exact CommonResidue.lonely_of_common_residue v d a hd2 hna hcres
  -- the residual class — now carrying `tupleGcd v = 1`
  · push_neg at hdom hrep hwin hdet
    obtain ⟨iw, hiw⟩ := hwin
    refine hresP v hv hgcd hcov hgap ?_ ?_ ⟨iw, by omega⟩ ?_ hms ?_
    · intro i
      obtain ⟨j, hji, hj⟩ := hdom i
      exact ⟨j, hji, hj⟩
    · intro i j hij
      exact hrep i j hij
    · intro g hg2 i₀ hH
      exact hdet g hg2 i₀ hH
    · intro d hd2 a hall
      by_contra hnda
      exact hcr ⟨d, hd2, a, hnda, hall⟩

/-- **Subsumption.** The primitive assembly recovers the original one: anyone holding the OLD
`ResidualObligation` still gets `LRC14Statement`. So switching to the primitive surface loses nothing
and strictly shrinks what must be proved (the dilates are gone). -/
theorem lrc14_grand_assembly_of_residual (cite : LRCUpTo13) (h : ResidualObligation) :
    LRC14.LRC14Statement :=
  lrc14_grand_assembly_primitive cite (residualPrimitive_of_residual h)

/-! ## The finish line, with primitivity: `hB5` restated on the dilate-free residual -/

/-- The `B5 > 0` certifier discharges the DILATE-FREE residual. Identical to
`residualObligation_of_B5` (kps-S123) except the supplier now also receives `tupleGcd v = 1`. -/
theorem residualObligationPrimitive_of_B5
    (hB5 : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → LRC14.tupleGcd v = 1 →
      LRC14.CoveringFamily v → GapFamily v →
      (∀ i, ∃ j, j ≠ i ∧ |v i| ≤ 13 * |v j|) → (∀ i j, i ≠ j → |v i| ≠ |v j|) →
      (∃ i, 23 ≤ |v i|) →
      (∀ g : ℤ, 2 ≤ g → ∀ i₀ : Fin 13, (∀ j, j ≠ i₀ → g ∣ v j) → g ∣ v i₀) →
      (¬ ∃ (L : ℤ) (k a : Fin 13 → ℤ) (A : ℝ), (∀ i, v i = a i + L * k i) ∧ 0 < (L : ℝ) ∧
        (∀ i, |(a i : ℝ)| ≤ A) ∧ A / (L : ℝ) ≤ 1/13 - 1/14 ∧ (∀ i, k i ≠ 0) ∧
        (Finset.univ.image k).card ≤ 12) →
      ∃ q : ℕ, 0 < q ∧ 0 < LRC14Concrete.B5 v q) :
    ResidualObligationPrimitive := by
  intro v hv hgcd hcov hgap hcomp hdist hlarge hdiv hcoarse _hcres
  obtain ⟨q, hq, hpos⟩ := hB5 v hv hgcd hcov hgap hcomp hdist hlarge hdiv hcoarse
  exact LRC14Concrete.lonely_of_Mreach_ge v hv (LRC14Concrete.mreach_ge_of_B5_pos v q hq hpos)

/-- **The finish line with the primitivity peel.** `LRC14Statement` from the LRC(≤13) citation and the
`B5 > 0` obligation on the DILATE-FREE residual. Strictly weaker than `lrc14_from_B5`: the supplier now
also gets `tupleGcd v = 1`, so the non-primitive dilates (which forced `inf μ = 0`, opus-S205) need not be
handled at all. -/
theorem lrc14_from_B5_primitive (cite : LRCUpTo13)
    (hB5 : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → LRC14.tupleGcd v = 1 →
      LRC14.CoveringFamily v → GapFamily v →
      (∀ i, ∃ j, j ≠ i ∧ |v i| ≤ 13 * |v j|) → (∀ i j, i ≠ j → |v i| ≠ |v j|) →
      (∃ i, 23 ≤ |v i|) →
      (∀ g : ℤ, 2 ≤ g → ∀ i₀ : Fin 13, (∀ j, j ≠ i₀ → g ∣ v j) → g ∣ v i₀) →
      (¬ ∃ (L : ℤ) (k a : Fin 13 → ℤ) (A : ℝ), (∀ i, v i = a i + L * k i) ∧ 0 < (L : ℝ) ∧
        (∀ i, |(a i : ℝ)| ≤ A) ∧ A / (L : ℝ) ≤ 1/13 - 1/14 ∧ (∀ i, k i ≠ 0) ∧
        (Finset.univ.image k).card ≤ 12) →
      ∃ q : ℕ, 0 < q ∧ 0 < LRC14Concrete.B5 v q) :
    LRC14.LRC14Statement :=
  lrc14_grand_assembly_primitive cite (residualObligationPrimitive_of_B5 hB5)

end LRC14Grand
end LonelyRunner
